"""Checked retained PSMILES graphs through the opt-in J12 native resolver.

The independently authored Decimal reader checks source selections and charges.
Missing rows remain failures; this CLI cannot authorize a diagnostic model.
"""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from pcff_j5_reference import charge_raw, read_source, resolve_raw

from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY
from island.workflows import storage
from island.workflows.bundle import system_from


def inspect(source_path, declaration, output):
    d = storage.read_json(declaration)
    source = load_pcff_source(source_path)
    policy = d["resolution_policy"]
    if source.identity["sha256"] != d["source"]["sha256"] or policy not in (
        MSI_POLICY,
        COMPATIBILITY_POLICY,
    ):
        raise ValueError("Declaration/source/resolver mismatch")
    output.mkdir(parents=True, exist_ok=False)
    storage.publish(output / "declaration.json", storage.json_bytes(d))
    raw = read_source(source.raw)
    cases = []
    for case in [d["case"], *d["blocked_cases"]]:
        path = Path(case["system_path"])
        if storage.checksum(path.read_bytes()) != case["system_sha256"]:
            raise ValueError("Retained system checksum mismatch: " + case["name"])
        system = system_from(storage.decode(storage.read_json(path)))
        before = system.to_dict()
        typing = type_pcff_atoms(system, source, profile=d["typing_profile"])
        result = {
            "case": case["name"],
            "typing_complete": typing.complete,
            "typing_identity": typing.identity,
            "system_sha256": case["system_sha256"],
        }
        if typing.complete:
            payload = inspect_pcff_operational_support(
                system,
                typing,
                resolution_policy=policy,
                special_pairs=special_pair_policy(
                    lj=d["special_pairs"]["lj"], coulomb=d["special_pairs"]["coulomb"]
                ),
            )
            charge = charge_raw(raw, system, typing.assignments)
            native = payload["native_charge_record"]["native_charge_record"]
            if charge["complete"] != payload["charges_complete"]:
                raise ValueError("Independent component-charge acceptance mismatch")
            # Independently selected raw rows, oriented endpoint contributions and
            # Decimal sums; the native partial vector is checked even on failure.
            np.testing.assert_allclose(
                [
                    float(charge["partial_charges"][i])
                    for i in sorted(typing.assignments)
                ],
                [native["partial_charges"][i] for i in sorted(typing.assignments)],
                atol=d["tolerances"]["charge_atol_e"],
                rtol=0,
            )
            checks, blockers = [], []
            for q in payload["interaction_queries"]:
                n = q["resolution"]
                if n["status"] == "not_applicable":
                    continue
                r = resolve_raw(raw, q["family"], q["supplied_types"], guarded_msi=True)
                if r["status"] != n["status"]:
                    raise ValueError(
                        "Independent selection status mismatch: " + q["request_id"]
                    )
                if r["status"] == "assigned":
                    if r["selected_ids"] != sorted(
                        {x["record_id"] for x in n["selected"]}
                    ):
                        raise ValueError(
                            "Independent source row mismatch: " + q["request_id"]
                        )
                    np.testing.assert_allclose(
                        r["values"],
                        n["normalized_values"],
                        atol=d["tolerances"]["coefficient_atol"],
                        rtol=0,
                    )
                else:
                    blockers.append({"request": q, "independent_search": r})
                checks.append(
                    {
                        "id": q["request_id"],
                        "status": r["status"],
                        "selected_ids": r.get("selected_ids", []),
                    }
                )
            storage.publish(
                output / (case["name"] + "-inspection.json"),
                storage.json_bytes(storage.encode(payload)),
            )
            result.update(
                native_charges_complete=payload["charges_complete"],
                model_complete=payload["model_complete"],
                model_identity=payload["model_identity"],
                independent_charge=charge,
                independent_selection_checks=checks,
                blockers=blockers,
                model_diagnostics=payload["model_diagnostics"],
            )
        if (
            system.to_dict() != before
            or storage.checksum(path.read_bytes()) != case["system_sha256"]
        ):
            raise ValueError("Retained input mutation")
        cases.append(result)
    receipt = {
        "schema": "island_j13_converter_receipt_v1"
        if policy == COMPATIBILITY_POLICY
        else "island_j12_resolver_receipt_v1",
        "source": source.identity,
        "resolution_policy": policy,
        "cases": cases,
        "counts": {
            k: sum(bool(c.get(k)) for c in cases)
            for k in ("typing_complete", "native_charges_complete", "model_complete")
        },
        "denominator": len(cases),
        "source_missing_by_family": dict(
            Counter(
                b["request"]["family"]
                for c in cases
                for b in c.get("blockers", [])
                if b["independent_search"]["status"] == "missing"
            )
        ),
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    storage.publish(output / "report.json", storage.json_bytes(storage.encode(receipt)))
    print(receipt["counts"])
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j12_declaration.json"),
    )
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    r = inspect(a.source, a.declaration, a.output)
    return 1 if a.require_full_source and not r["full_source_complete"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
