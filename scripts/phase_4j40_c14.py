"""Reproduce backend-scoped decisions from unchanged, pinned C13 artifacts."""

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

from island.analysis.backend_readiness import SCOPES, assess_periodic_backend_readiness
from island.analysis.periodic_preflight import PeriodicPreflightPolicy, preflight_box
from scripts.phase_4j40_c11_numeric import require, sha, signed, write
from scripts.phase_4j40_c13 import validate_summary

C13_RECEIPT = "4d394a73291a535f9c03f6770376780f060edb49be6e132ec0198638efae06e2"
CELLS = ("gaff2-pe", "opls-peo", "opls-resp-aromatic", "opls-network")


def load_source(directory):
    require(sha(directory / "receipt.json") == C13_RECEIPT, "Changed C13 receipt")
    receipt = json.loads((directory / "receipt.json").read_text())
    require(
        sha(directory / "manifest.json") == receipt["manifest_sha256"],
        "Changed C13 manifest",
    )
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, digest in manifest["files"].items():
        require(sha(directory / name) == digest, "Changed C13 source: " + name)
    summary = json.loads((directory / "summary.json").read_text())
    validate_summary(summary)
    require(receipt["summary_identity"] == summary["identity"], "Changed C13 summary")
    require(
        receipt["native_evaluations"] == 300
        and len(summary["rows"]) == 108
        and sum(r["existing_energy_force_bounds_passed"] for r in summary["rows"])
        == 13,
        "Changed C13 outcomes",
    )
    return receipt, summary


def build_report(directory):
    _receipt, summary = load_source(directory)
    inventory = json.loads((directory / "artifact-inventory.json").read_text())
    require(
        sum(name.endswith("/point.json") for name in inventory) == 300,
        "Incomplete native inventory",
    )
    backends = {}
    with tarfile.open(directory / "artifacts.tar.xz") as archive:
        for cell in CELLS:
            name = f"study/source-originals/{cell}/backend.json"
            raw = archive.extractfile(name).read()
            require(
                hashlib.sha256(raw).hexdigest() == inventory[name],
                "Changed backend source",
            )
            backends[cell] = (json.loads(raw), inventory[name])
        additional = json.load(archive.extractfile("study/additional-matrix.json"))
    with tarfile.open(directory / "validation.tar.xz") as archive:
        engines = json.load(archive.extractfile("validation/backend-inventory.json"))
        engines["lammps_version_and_build"] = (
            archive.extractfile("validation/lammps-build.txt").read().decode()
        )
    engines = {
        **engines,
        "openmm_version": additional[0]["openmm"]["version"],
        "openmm_platform": additional[0]["openmm"]["platform"],
        "comparison_semantics": "Native LAMMPS and explicitly labelled analytical special-pair composite; not a production LAMMPS repair",
    }
    cells = []
    for cell, (backend, source_hash) in backends.items():
        bindings = {
            "final_graph": backend["final_graph_identity"],
            "parameters": backend["assignment_identity"],
            "charges": backend["charge_identity"],
            "backend": signed(engines)["identity"],
            "protocol": signed(
                {
                    "c13_policy": summary["policy"],
                    "cell": cell,
                    "source_hash": source_hash,
                }
            )["identity"],
        }
        preflight = preflight_box(
            initial_box_nm=[v / 10 for v in backend["config"]["box_lengths"]],
            cutoff_nm=backend["config"]["nonbonded_cutoff"] / 10,
            policy=PeriodicPreflightPolicy(safety_margin_nm=0.05, neighbor_skin_nm=0.2),
        )
        gates = {}

        def gate(
            key,
            status,
            reason,
            *,
            gates=gates,
            bindings=bindings,
            source_hash=source_hash,
        ):
            gates[key] = {
                "name": key,
                "scope": SCOPES[key],
                "bindings": bindings,
                "status": status,
                "reason": reason,
                "evidence_sha256": [source_hash, C13_RECEIPT],
            }

        gate(
            "openmm.exact_graph_parameters",
            "passed",
            "Exact retained graph/assignment/source identities match pinned C13 bytes; no parameter changes",
        )
        gate(
            "openmm.validated_charge_origin",
            "not_measured",
            "Retained charge metadata and provenance preserved; C14 does not independently requalify raw charge-generation evidence",
        )
        gate(
            "openmm.cutoff_box_preflight",
            "blocked" if preflight["status"] == "blocked" else "not_measured",
            preflight["diagnostic"]
            + "; no independently established minimum NPT box bound",
        )
        gate(
            "openmm.internal_energy_force",
            "not_measured",
            "Native single-point software execution exists; complete actual-model internal force/energy qualification not established by cross-engine execution",
        )
        gate(
            "openmm.pressure_estimator_controls",
            "unavailable",
            "pressure_unavailable: full periodic reciprocal virial not independently measured; finite-volume agreement alone is insufficient",
        )
        gate(
            "openmm.checkpoint_continuation",
            "not_measured",
            "No new C13/C14 qualified periodic dynamics/checkpoint continuation; historical software tests do not attest this exact long protocol",
        )
        gate(
            "openmm.sufficient_sampling",
            "blocked",
            "Single-point diagnostics provide no independent trajectory samples, burn-in, ESS or replica convergence",
        )
        gate(
            "lammps_export.exact_frame",
            "passed",
            "C13 native data decks exported and executed for retained single frames only; production exporter remains conservative",
        )
        gate(
            "lammps_kspace.native_execution",
            "passed",
            "Real pinned KSPACE executable run in C13; execution is not physical equivalence",
        )
        gate(
            "cross_engine.actual_polymer_systems",
            "passed",
            "Actual retained polymer frames are present in the unchanged comparison matrix",
        )
        gate(
            "cross_engine.numerical_policy",
            "blocked",
            "C13 numerical-equivalence policy remains candidate_not_qualified; no widened bounds",
        )
        gate(
            "cross_engine.component_energy",
            "blocked",
            "C13 component/total residuals fail unchanged bounds across actual polymer matrix",
        )
        gate(
            "cross_engine.per_atom_force",
            "blocked",
            "Actual polymer force closure across the retained matrix is unresolved; per-atom residuals preserved",
        )
        gate(
            "cross_engine.virial",
            "unavailable",
            "Full independent OpenMM periodic reciprocal virial tensor unavailable",
        )
        gate(
            "cross_engine.pressure",
            "unavailable",
            "Three-way actual-polymer pressure qualification not established",
        )
        decision = assess_periodic_backend_readiness(bindings=bindings, checks=gates)
        cells.append(
            {
                "cell": cell,
                "backend_source_sha256": source_hash,
                "retained_charge_metadata": {
                    k: backend["graph_charges"].get(k)
                    for k in (
                        "charge_origin",
                        "charge_method",
                        "provenance",
                        "evidence",
                    )
                },
                "retained_config_provenance": backend["config"]["provenance"],
                "preflight": preflight,
                "decision": decision,
            }
        )
    return signed(
        {
            "schema": "island.j40-c14-comparison-report.v1",
            "c13_receipt_sha256": C13_RECEIPT,
            "c13_summary_identity": summary["identity"],
            "source_hashes": {
                name: sha(directory / name)
                for name in (
                    "receipt.json",
                    "manifest.json",
                    "summary.json",
                    "artifacts.tar.xz",
                    "validation.tar.xz",
                    "artifact-inventory.json",
                )
            },
            "native_evaluations": 300,
            "comparison_count": 108,
            "passing_comparisons": 13,
            "systems_tested": sorted({r["case"] for r in additional}),
            "engines": engines,
            "comparisons": summary["rows"],
            "inherited_policy": summary["policy"],
            "new_qualified_equivalence_policy": None,
            "cells": cells,
            "known_limitations": summary["convergence_limitations"],
            "new_md_execution": "not_executed",
            "workflow_scope": "software_tests_only",
            "cross_engine_equivalence_status": "ewald_unvalidated",
            "pressure_validation_status": "pressure_unavailable",
            "scientific_equilibrium_status": "not_established",
            "scientific_npt_acceptance": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    write(args.output, build_report(args.source))
    return 1  # Scientific acceptance remains blocked; report construction succeeded.


if __name__ == "__main__":
    raise SystemExit(main())
