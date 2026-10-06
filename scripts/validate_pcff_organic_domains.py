"""Execute J4 declared source checks; retain incomplete native models explicitly."""

import argparse
import json
import time
from collections import Counter
from pathlib import Path

from validate_pcff_expanded import independent_charges

from island.chemistry import from_smiles
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    assign_pcff_source_types,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.domain_coverage import pcff_domain_coverage
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.forcefields.pcff.organic_domains import PROFILE_NAME
from island.workflows.bundle import system_data
from island.workflows.storage import checksum, encode, json_bytes, publish


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    p.add_argument(
        "--controls-only",
        action="store_true",
        help="Record declared near misses; no parameterization or numerical acceptance",
    )
    a = p.parse_args()
    root = a.output
    root.mkdir(parents=True, exist_ok=False)
    declaration = json.loads(
        Path("docs/evidence/phase_4j4_declaration.json").read_text()
    )
    src = load_pcff_source(a.source)
    assert src.identity == declaration["source"]
    publish(root / "declaration.json", json_bytes(declaration))
    report = {
        "source": src.identity,
        "profile": PROFILE_NAME,
        "cases": [],
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if a.controls_only:
        for smiles in declaration["negative_controls"]:
            row = {"smiles": smiles}
            try:
                system = from_smiles(smiles, random_seed=declaration["seed"])
                typing = type_pcff_atoms(system, src, profile=PROFILE_NAME)
                row.update(
                    typed=typing.complete,
                    assignments=typing.payload["assignments"],
                    diagnostics=typing.payload["diagnostics"],
                )
                if typing.complete:
                    charges = assign_automatic_pcff_charges(
                        system, typing, resolution_policy=DOMAIN_POLICY
                    )
                    row.update(
                        charges_complete=charges.complete,
                        charge_diagnostics=charges.payload["native_charge_record"][
                            "diagnostics"
                        ],
                    )
            except Exception as error:  # noqa: BLE001 -- preserve unsupported construction
                row["error"] = type(error).__name__ + ": " + str(error)
            report["cases"].append(row)
        publish(root / "controls.json", json_bytes(report))
        return 1 if a.require_full_source else 0
    for case in declaration["cases"]:
        out = {
            **case,
            "typed": False,
            "charges_complete": False,
            "model_complete": False,
        }
        work = root / case["name"]
        work.mkdir()
        start = time.monotonic()
        try:
            s = from_smiles(case["smiles"], random_seed=declaration["seed"])
            publish(work / "system.json", json_bytes(encode(system_data(s))))
            t = type_pcff_atoms(s, src, profile=PROFILE_NAME)
            publish(work / "typing.json", t.json_text.encode())
            out.update(
                typed=t.complete,
                typing_identity=t.identity,
                typing_diagnostics=t.payload["diagnostics"],
                actual_counts=dict(Counter(t.assignments.values())),
            )
            assert out["actual_counts"] == case["expected_counts"]
            explicit = assign_pcff_source_types(
                s,
                src,
                t.assignments,
                profile=PROFILE_NAME,
                provenance="J4 validated automatic labels, not independent automatic typing",
            )
            assert explicit.assignments == t.assignments
            q = assign_automatic_pcff_charges(s, t, resolution_policy=DOMAIN_POLICY)
            publish(work / "charges.json", q.json_text.encode())
            out.update(
                charges_complete=q.complete,
                charge_identity=q.identity,
                charge_diagnostics=q.payload["native_charge_record"]["diagnostics"],
            )
            if not q.complete:
                continue
            try:
                ref, rows = independent_charges(src.raw, s, t.assignments)
                error = max(abs(ref[i] - q.charges[i]) for i in ref)
                assert error <= 1e-12
                out["independent_charge"] = {"max_error": error, "rows": rows}
            except ValueError as err:
                out["independent_charge_unavailable"] = str(err)
            assignment = assign_pcff_parameters(
                s, t, q, resolution_policy=DOMAIN_POLICY
            )
            publish(work / "assignment.json", assignment.json_text.encode())
            model = define_pcff_model(
                assignment,
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
            publish(work / "model.json", model.json_text.encode())
            out.update(
                model_complete=model.payload["model_definition_complete"],
                model_identity=model.identity,
                diagnostics=model.payload["diagnostics"],
                coverage=assignment.payload["coverage"],
            )
        except Exception as err:  # noqa: BLE001 -- diagnostic matrix preserves failures
            out["error"] = type(err).__name__ + ": " + str(err)
        finally:
            out["wall_seconds"] = time.monotonic() - start
            publish(work / "outcome.json", json_bytes(out))
            report["cases"].append(out)
            print(
                case["name"],
                out["typed"],
                out["charges_complete"],
                out["model_complete"],
                out.get("error", ""),
                flush=True,
            )
    report["source_unchanged"] = (
        checksum(a.source.read_bytes()) == declaration["source"]["sha256"]
    )
    report["counts"] = {
        key: sum(c[key] for c in report["cases"])
        for key in ("typed", "charges_complete", "model_complete")
    }
    report["denominator"] = len(report["cases"])
    publish(root / "report.json", json_bytes(report))
    publish(
        root / "coverage.json",
        json_bytes(pcff_domain_coverage(src, profile=PROFILE_NAME)),
    )
    # Typing evidence alone is never a numerical/domain or full-source gate.
    return (
        1
        if a.require_full_source
        else int(
            not report["source_unchanged"]
            or any(not c["typed"] or "error" in c for c in report["cases"])
        )
    )


if __name__ == "__main__":
    raise SystemExit(main())
