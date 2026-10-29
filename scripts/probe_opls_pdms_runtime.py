"""Probe real OPLS preparation of PDMS; no parameter fallback or PDMS MD."""

import argparse
import json
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from math import isfinite, sqrt
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.exceptions import OPLSAssignmentError, ValidationError
from island.forcefields.external_pdms import LIMITS, digest, pdms_final_graph
from island.forcefields.oplsaa import load_oplsaa_source
from island.forcefields.oplsaa.adapter import check_foyer_installation
from island.forcefields.preparation import (
    ForceFieldRequest,
    OPLSOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.graph import final_graph
from island.minimization import minimize_geometry


def fixtures():
    for name, smiles in (
        ("ethane_control", "CC"),
        ("butane_control", "CCCC"),
        ("benzene_control", "c1ccccc1"),
        ("disiloxane_SiH3", "[SiH3]O[SiH3]"),
        ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C"),
    ):
        yield name, smiles, from_smiles(smiles, random_seed=2026)
    psmiles = "[*:1]O[Si](C)(C)[*:2]"
    for dp in (1, 2, 4):
        yield (
            f"TMS_capped_PDMS_DP{dp}",
            psmiles,
            build_linear_polymer(
                psmiles,
                dp=dp,
                random_seed=2026,
                add_hydrogens=True,
                head_end_group="[*][Si](C)(C)C",
                tail_end_group="[*]C",
            ),
        )


def probe(xml):
    source = load_oplsaa_source(xml)
    environment = check_foyer_installation()[3]
    xml_types, _ = source.entries
    silicon_types = [n for n, r in xml_types.items() if r.get("element") == "Si"]
    rows = []
    for name, smiles, system in fixtures():
        before = deepcopy(system.to_dict())
        silicon = [i for i, a in system.topology.sites.items() if a.element == "Si"]
        native_graph, generic_error = None, None
        if silicon:
            # The existing separately versioned Si-aware chemical validator is
            # an audit check only. It is not substituted into the native API.
            validated_graph = pdms_final_graph(system)
            graph_id, graph_schema = digest(validated_graph), validated_graph["schema"]
            try:
                native_graph = final_graph(system)
            except ValidationError as error:
                generic_error = {
                    "exception_type": type(error).__name__,
                    "message": str(error),
                    "probe_classification": "typing_blocked",
                    "reason": "historical_generic_final_graph_does_not_support_Si",
                }
        else:
            native_graph = final_graph(system)
            graph_id = native_graph.identity
            graph_schema = native_graph.payload["schema"]
        row = {
            "fixture": name,
            "input": smiles,
            "validated_graph_identity": graph_id,
            "validated_graph_schema": graph_schema,
            "final_graph_valid": True,
            "generic_final_graph_error": generic_error,
            "elements": dict(Counter(a.element for a in system.topology.sites.values())),
            "atom_count": len(system.topology.sites),
            "bond_count": len(system.topology.bonds),
            "silicon_site_ids": silicon,
            "prepared": False,
            "evaluated": False,
            "minimization_run": False,
            "MD_run": False,
        }
        try:
            prepared = prepare_forcefield(
                system,
                ForceFieldRequest("oplsaa", OPLSOptions(xml), final_graph=native_graph),
            )
        except OPLSAssignmentError as error:
            # This is the probe's classification of the native exception, not a
            # fabricated status field on the original native exception.
            if (
                not silicon
                or str(error) != "Unsupported or inconsistent element/atomic number"
            ):
                raise
            row.update(
                probe_classification="typing_blocked",
                blocking_stage="native_opls_element_validation",
                exception_type=type(error).__name__,
                exception_message=str(error),
                downstream="not_attempted_no_prepared_forcefield",
            )
        else:
            row["prepared"] = True
            row["prepared_identity"] = prepared.identity
            if silicon:
                # This probe cannot promote newly accepted Si to a qualified
                # domain. Independent source/numerical gates remain mandatory.
                row["probe_classification"] = "unexpected_native_Si_preparation"
                row["downstream"] = "not_attempted_pending_independent_qualification"
            else:
                evaluator = create_evaluator(system, prepared)
                energy = evaluator.evaluate()
                row.update(
                    probe_classification="control_evaluated",
                    evaluated=True,
                    energy_kJ_per_mol=energy.potential_energy,
                    energy_components=dict(energy.energy_components),
                    forces_finite=all(
                        isfinite(v) for xyz in energy.forces.values() for v in xyz
                    ),
                    fmax_kJ_per_mol_A=max(
                        sqrt(sum(v * v for v in xyz)) for xyz in energy.forces.values()
                    ),
                    backend=energy.backend_name,
                    backend_version=energy.backend_version,
                    platform=energy.platform,
                )
                minimum = minimize_geometry(system, evaluator)
                row["minimization_run"] = True
                row["minimization"] = {
                    "converged": minimum.converged,
                    "termination_reason": minimum.termination_reason,
                    "iterations": minimum.iterations,
                    "evaluations": minimum.evaluations,
                    "final_evaluation_verified": minimum.final_evaluation_verified,
                    "final_energy_kJ_per_mol": minimum.final_evaluation.potential_energy,
                }
        if system.to_dict() != before:
            raise ValueError("Probe mutated the input system")
        row["input_unchanged"] = True
        rows.append(row)
    report = {
        "schema": "island_opls_pdms_runtime_probe_v1",
        "native_source_identity": source.identity,
        "native_silicon_type_labels": silicon_types,
        "environment": environment,
        "cases": rows,
        "native_PDMS_executable": False,
        "PDMS_MD_runs": 0,
        "new_execution_profiles": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
        "scope": "Actual public preparation calls; native non-Si controls evaluated/minimized. No external model substitution or PDMS MD.",
    }
    controls = [r for r in rows if not r["silicon_site_ids"]]
    pdms = [r for r in rows if r["silicon_site_ids"]]
    report["boundary_probe_passed"] = (
        not silicon_types
        and len(pdms) == 5
        and all(r["probe_classification"] == "typing_blocked" for r in pdms)
        and all(
            r["probe_classification"] == "control_evaluated"
            and r["forces_finite"]
            and r["minimization"]["converged"]
            and r["minimization"]["final_evaluation_verified"]
            for r in controls
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-pdms-executable", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output file already exists")
    try:
        report = probe(args.xml)
    except Exception as error:  # noqa: BLE001 -- report failures, never substitute data
        print(
            json.dumps(
                {
                    "status": "probe_failed",
                    "error_type": type(error).__name__,
                    "message": str(error),
                    "acceptance_gate": 1,
                }
            )
        )
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(
        json.dumps(
            {
                "boundary_probe_passed": report["boundary_probe_passed"],
                "native_PDMS_executable": report["native_PDMS_executable"],
                "acceptance_gate": report["acceptance_gate"],
                "report_sha256": sha256(args.output.read_bytes()).hexdigest(),
            }
        )
    )
    return 1 if args.require_pdms_executable or not report["boundary_probe_passed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
