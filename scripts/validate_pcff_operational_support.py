"""J5 checked retained-input inspection; incomplete operational gates exit nonzero."""

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from pcff_j5_reference import charge_raw, raw_inventory, read_source, resolve_raw

from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.workflows.bundle import system_from
from island.workflows.storage import decode, encode


def write(path, data):
    with path.open("x") as stream:
        json.dump(encode(data), stream, indent=2, allow_nan=False)
        stream.write("\n")


def operational_gate(cases, required):
    names = [c.get("case") for c in cases]
    return (
        bool(required)
        and Counter(names) == Counter(required)
        and all(
            not c.get("error")
            and all(
                c.get(k) is True
                for k in (
                    "typed",
                    "charges_complete",
                    "model_complete",
                    "source_verified",
                    "independent_selection_verified",
                    "whole_system_verified",
                    "relocated_bundle_verified",
                )
            )
            for c in cases
        )
    )


def inspect_case(c, declaration, source, raw_rows, directory):
    p = Path(c["retained_input"]["path"])
    assert hashlib.sha256(p.read_bytes()).hexdigest() == c["retained_input"]["sha256"]
    system = system_from(decode(json.loads(p.read_text())))
    typing = type_pcff_atoms(system, source, profile=declaration["profile"])
    assert Counter(typing.payload["assignments"].values()) == c["expected_counts"]
    result = inspect_pcff_operational_support(
        system,
        typing,
        resolution_policy=declaration["resolution_policy"],
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    write(directory / "inspection.json", result)
    queries = result["interaction_queries"]
    assert Counter((q["family"], tuple(q["sites"])) for q in queries) == Counter(
        raw_inventory(system)
    )
    oracle = []
    for q in queries:
        if q["resolution"]["status"] == "not_applicable":
            continue
        r = resolve_raw(raw_rows, q["family"], q["supplied_types"])
        a = q["resolution"]
        same = a["status"] == r["status"]
        if a["status"] == "assigned" and r["status"] == "assigned":
            same = (
                same
                and len(a["normalized_values"]) == len(r["values"])
                and all(
                    math.isclose(
                        x,
                        y,
                        rel_tol=declaration["tolerances"]["rtol"],
                        abs_tol=declaration["tolerances"]["energy_atol_kj_mol"],
                    )
                    for x, y in zip(a["normalized_values"], r["values"], strict=True)
                )
            )
            same = (
                same
                and sorted({v["record_id"] for v in a["selected"]}) == r["selected_ids"]
            )
        oracle.append({"request_id": q["request_id"], "agrees": same, "reference": r})
    charges = charge_raw(raw_rows, system, typing.payload["assignments"])
    assert charges["complete"] == result["charges_complete"]
    native = result["native_charge_record"]["native_charge_record"]
    max_charge_difference = max(
        abs(float(v) - native["partial_charges"][i])
        for i, v in charges["partial_charges"].items()
    )
    assert max_charge_difference <= declaration["tolerances"]["native_charge_e"]
    charges["maximum_native_difference_e"] = max_charge_difference
    write(
        directory / "raw-reference.json", {"interactions": oracle, "charges": charges}
    )
    old = json.loads((p.parent / "outcome.json").read_text())
    identities = {
        k: result[k]
        for k in (
            "typing_identity",
            "charge_identity",
            "assignment_identity",
            "model_identity",
        )
    }
    preserved = all(
        v == old[k] for k, v in identities.items() if v is not None and k in old
    )
    assert preserved
    return {
        "case": c["name"],
        "typed": typing.complete,
        "charges_complete": result["charges_complete"],
        "model_complete": result["model_complete"],
        "source_verified": True,
        "independent_selection_verified": all(q["agrees"] for q in oracle),
        "independent_requests": len(oracle),
        "inventory_verified": True,
        "historical_identities_preserved": preserved,
        **identities,
        "query_status_counts": dict(
            Counter(q["resolution"]["status"] for q in queries)
        ),
        "unresolved": [
            {
                "request_id": q["request_id"],
                "types": q["supplied_types"],
                "resolution": q["resolution"],
                "dependencies": q["equilibrium_dependencies"],
            }
            for q in queries
            if q["resolution"]["status"] not in ("assigned", "not_applicable")
            or not q["dependencies_complete"]
        ],
        "charge_components": charges["components"],
        "missing_charge_bonds": charges["missing"],
        "whole_system_verified": False,
        "relocated_bundle_verified": False,
        "blocked_gate": "Incomplete native model: no whole-system evaluation or bundle publication",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--declaration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-full-source", action="store_true")
    args = parser.parse_args()
    declaration = json.loads(args.declaration.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / "declaration.json", declaration)
    report = {
        "schema": "island_j5_operational_receipt_v1",
        "cases": [],
        "operational_gate": False,
        "full_source_gate": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        raw = args.source.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == declaration["source"]["sha256"]
        source = load_pcff_source(args.source)
        assert source.identity == declaration["source"]
        rows = read_source(raw)
        for c in declaration["cases"]:
            directory = args.output / c["name"]
            directory.mkdir()
            try:
                result = inspect_case(c, declaration, source, rows, directory)
            except Exception as error:  # noqa: BLE001 -- retain diagnostic failures
                result = {
                    "case": c["name"],
                    "error": f"{type(error).__name__}: {error}",
                }
            report["cases"].append(result)
            print(
                c["name"],
                result.get("query_status_counts", result.get("error")),
                flush=True,
            )
        report["operational_gate"] = operational_gate(
            report["cases"], [c["name"] for c in declaration["cases"]]
        )
    except Exception as error:  # noqa: BLE001 -- retain diagnostic failures
        report["error"] = f"{type(error).__name__}: {error}"
    write(args.output / "report.json", report)
    return (
        0
        if report["operational_gate"]
        and (not args.require_full_source or report["full_source_gate"])
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
