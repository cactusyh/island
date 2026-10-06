"""Declared final-graph domains: strict eligibility separated from source diagnostics."""

import argparse
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields import PCFFOptions
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_source,
    prepare_pcff_polymer,
    select_pcff_profile,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.forcefields.pcff.registry import FUSED_BENZENOID
from island.workflows import storage
from island.workflows.bundle import system_data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--declaration", default="docs/evidence/phase_4j10_declaration.json")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    d = storage.read_json(a.declaration)
    a.output.mkdir(parents=True, exist_ok=False)
    storage.publish(a.output / "declaration.json", storage.json_bytes(d))
    source = load_pcff_source(a.source)
    results = []
    for case in d["matrix"] + d["charged_cases"]:
        out = a.output / case["name"]
        out.mkdir()
        result = {
            "case": case["name"],
            "typing_complete": False,
            "native_charge_complete": False,
            "raw_parameters_complete": False,
            "strict_prepared": False,
            "independent_numerical_verified": False,
        }
        try:
            s = (
                build_linear_polymer(case["psmiles"], dp=case["dp"], random_seed=2026)
                if "psmiles" in case
                else from_smiles(case["smiles"], random_seed=2026)
            )
            before = s.to_dict()
            storage.publish(
                out / "system.json", storage.json_bytes(storage.encode(system_data(s)))
            )
            t = type_pcff_atoms(s, source, profile="island_pcff_source_graph_v4")
            storage.publish(out / "typing.json", t.json_text.encode())
            result["typing_complete"] = t.complete
            result["typing_diagnostics"] = t.payload["diagnostics"]
            if t.complete:
                q = assign_automatic_pcff_charges(s, t, resolution_policy=DOMAIN_POLICY)
                storage.publish(out / "charge.json", q.json_text.encode())
                result["native_charge_complete"] = q.complete
                result["native_charge_observation"] = q.payload["native_charge_record"]
                if q.complete:
                    assignment = assign_pcff_parameters(
                        s, t, q, resolution_policy=DOMAIN_POLICY
                    )
                    storage.publish(
                        out / "assignment.json", assignment.json_text.encode()
                    )
                    result["raw_parameters_complete"] = assignment.payload[
                        "parameter_coverage_complete"
                    ]
                    result["interaction_blockers"] = [
                        r
                        for r in assignment.payload["assignments"]
                        if r["status"] not in ("assigned", "not_applicable")
                    ]
                    model = define_pcff_model(
                        assignment,
                        special_pairs=special_pair_policy(
                            lj=(0, 0, 1), coulomb=(0, 0, 1)
                        ),
                    )
                    storage.publish(out / "model.json", model.json_text.encode())
                    result["model_definition_complete"] = model.payload[
                        "model_definition_complete"
                    ]
                    result["model_diagnostics"] = model.payload["diagnostics"]
            # Strict preparation independently requires the operational domain
            # and raw source-backed terms; diagnostic completeness cannot authorize it.
            try:
                options = PCFFOptions(
                    a.source,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v1",
                    source_profile=select_pcff_profile(
                        FUSED_BENZENOID, sha256=source.identity["sha256"]
                    ),
                )
                strict = prepare_pcff_polymer(s, options)
                storage.publish(
                    out / "strict-final-graph.json", strict.json_text.encode()
                )
                result.update(
                    strict_prepared=True,
                    strict_model=strict.prepared.native_result.identity,
                    strict_prepared_identity=strict.prepared.identity,
                )
            except Exception as exc:  # noqa: BLE001 -- diagnostic ledger never produces an executable substitute
                result["strict_blocker"] = type(exc).__name__ + ": " + str(exc)
            if s.to_dict() != before:
                raise ValueError("Input mutation")
            result["nonmutation"] = True
        except Exception as exc:  # noqa: BLE001 -- retain every required case and exact failure
            result["execution_error"] = type(exc).__name__ + ": " + str(exc)
        storage.publish(out / "outcome.json", storage.json_bytes(result))
        results.append(result)
    report = {
        "schema": "island_j10_parity_diagnosis_v1",
        "cases": results,
        "source": source.identity,
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    storage.publish(a.output / "report.json", storage.json_bytes(report))
    print(
        [
            (
                r["case"],
                r["typing_complete"],
                r["native_charge_complete"],
                r["raw_parameters_complete"],
                r["strict_prepared"],
            )
            for r in results
        ]
    )
    return (
        1
        if a.require_full_source or not all(r["strict_prepared"] for r in results)
        else 0
    )


if __name__ == "__main__":
    raise SystemExit(main())
