"""Bounded actual AmberTools acceptance with synthetic NON-scientific charges.

Declared matrix: PE50 GAFF (302), PE100 GAFF2 (602), PEO20 GAFF2 (142).
Fixed seeds=2026; alternating sorted-site +/-0.01 e, neutral; charge tol=1e-4.
Reference acceptance: energies/components atol=2e-5 kJ/mol, forces atol=2e-5
kJ/(mol*A), rtol=2e-6. No changes after observing results.
Workflow smoke: PEO20 only, fmax=0.1, 500 iterations/1000 evaluations,
300 K, friction=5/ps, 0.1 fs, 4 steps in two 2-step segments, 4 calls/segment.
"""

import argparse
import os
import resource
import time
from pathlib import Path

import numpy as np

from island.analysis import analyze_workflow, export_analysis
from island.builders import build_linear_polymer
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import (
    AmberToolsOptions,
    AmberToolsParameterizationEngine,
    import_amber_prmtop,
)
from island.minimization import MinimizationOptions
from island.workflows import WorkflowConfig, start_prepared_workflow, storage
from island.workflows.bundle import system_data

try:
    from scripts.validate_singlepoint_references import reference
except ModuleNotFoundError:
    from validate_singlepoint_references import reference

CASES = (
    ("pe50", "[*]CC[*]", 50, "gaff", 302),
    ("pe100", "[*]CC[*]", 100, "gaff2", 602),
    ("peo20", "[*]CCO[*]", 20, "gaff2", 142),
)
SOURCE = "Synthetic alternating sorted-site +/-0.01 e; software acceptance only; no scientific charge model"


def acceptance_gates(rows):
    """Separate software acceptance gates; missing cases never pass."""
    indexed = {row["case"]: row for row in rows}
    parameterization = len(indexed) == len(CASES) == len(rows) and all(
        indexed.get(case[0], {}).get("status")
        == "parameterization_and_numerical_acceptance_passed"
        for case in CASES
    )
    target = indexed.get("peo20", {})
    workflow = target.get("workflow", {})
    return {
        "parameterization_numerical": parameterization,
        "workflow_dynamics_analysis": bool(
            parameterization
            and workflow.get("status") == "completed"
            and workflow.get("accepted_step") == 4
            and target.get("analyzed_samples") == 5
            and target.get("analysis_verified") is True
        ),
    }


def run(root, amberhome):
    os.environ["PATH"] = (
        str(amberhome / "bin") + os.pathsep + os.environ.get("PATH", "")
    )
    root.mkdir(parents=True, exist_ok=False)
    storage.publish(
        root / "declared-settings.json",
        storage.json_bytes(
            {
                "matrix": CASES,
                "charge_source": SOURCE,
                "max_atoms": 1000,
                "charge_tolerance": 1e-4,
                "timeout_seconds_per_tool": 600,
                "rtol": 2e-6,
                "atol": 2e-5,
                "workflow": __doc__,
            }
        ),
    )
    rows = []
    for name, psmiles, dp, family, expected in CASES:
        directory = root / name
        directory.mkdir()
        row = {
            "case": name,
            "force_field": family,
            "expected_atoms": expected,
            "charge_source": SOURCE,
            "status": "failed",
            "timings_seconds": {},
        }
        t = time.perf_counter()
        try:
            print(f"{name}: building", flush=True)
            system = build_linear_polymer(
                psmiles,
                dp=dp,
                coordinate_method="local_templates",
                template_seed=2026,
                assembly_seed=2026,
            )
            row["timings_seconds"]["construction"] = time.perf_counter() - t
            row["actual_atoms"] = system.number_of_sites
            ids = sorted(system.topology.sites)
            charges = {s: (0.01 if i % 2 == 0 else -0.01) for i, s in enumerate(ids)}
            assert len(ids) == expected and len(ids) % 2 == 0
            storage.publish(
                directory / "input-system.json",
                storage.json_bytes(storage.encode(system_data(system))),
            )
            storage.publish(
                directory / "charges.json",
                storage.json_bytes({"source": SOURCE, "charges": charges}),
            )
            before = system.to_dict()
            t = time.perf_counter()
            print(f"{name}: parameterizing {len(ids)} sites", flush=True)
            prepared = AmberToolsParameterizationEngine().parameterize(
                system,
                AmberToolsOptions(
                    family,
                    "provided",
                    charges,
                    max_atoms=1000,
                    charge_source=SOURCE,
                    amberhome=amberhome,
                    work_root=directory,
                    retain_success_artifacts=True,
                ),
            )
            row["timings_seconds"]["parameterization_total"] = time.perf_counter() - t
            assert before == system.to_dict()
            r = prepared.record
            artifacts = Path(r["artifact_dir"])
            row["record_signature"] = prepared.record_signature
            row["import_signature"] = prepared.imported_result.result_signature
            row["timings_seconds"].update(
                {s["stage"]: s["elapsed_seconds"] for s in r["stages"]}
            )
            row["charge_max_absolute_error_e"] = max(
                abs(charges[s] - a.charge)
                for s, a in prepared.imported_result.charge_result.assignments.items()
            )
            assert row["charge_max_absolute_error_e"] <= 1e-5
            assert (
                r["charge_outcome"]["qm_run"] is False
                and not (artifacts / "sqm.out").exists()
            )
            storage.publish(
                directory / "preparation.json",
                storage.json_bytes(
                    storage.encode(
                        {
                            "record": dict(r),
                            "record_signature": prepared.record_signature,
                            "import_source": prepared.imported_result.source,
                            "import_provenance": dict(
                                prepared.imported_result.provenance
                            ),
                        }
                    )
                ),
            )
            t = time.perf_counter()
            import_amber_prmtop(
                system,
                artifacts / "result.prmtop",
                dict(prepared.imported_result.mapping),
                source="Independent timing reimport",
                force_field=family,
                charge_method="provided",
                charge_tolerance=1e-4,
            )
            row["timings_seconds"]["independent_reimport_including_validation"] = (
                time.perf_counter() - t
            )
            t = time.perf_counter()
            prepared.validate_integrity(system)
            row["timings_seconds"]["preparation_integrity_recheck"] = (
                time.perf_counter() - t
            )
            t = time.perf_counter()
            evaluator = OpenMMSinglePointEvaluator(
                system, prepared.imported_result, platform="Reference"
            )
            row["timings_seconds"]["openmm_binding"] = time.perf_counter() - t
            t = time.perf_counter()
            result = evaluator.evaluate()
            row["timings_seconds"]["fresh_singlepoint"] = time.perf_counter() - t
            coordinates = {s: system.coordinates.get(s) for s in ids}
            t = time.perf_counter()
            energy, forces, components = reference(
                artifacts / "result.prmtop",
                prepared.imported_result.mapping,
                coordinates,
            )
            row["timings_seconds"]["independent_prmtop_reference"] = (
                time.perf_counter() - t
            )
            observed = {
                "bond": result.energy_components["bond"],
                "angle": result.energy_components["angle"],
                "torsion": result.energy_components["proper_torsion"]
                + result.energy_components["periodic_improper"],
                "nonbonded": result.energy_components["nonbonded"],
            }
            comparisons = {}
            for label, a, b in [
                ("energy", result.potential_energy, energy),
                ("forces", [result.forces[s] for s in ids], [forces[s] for s in ids]),
                *[(key, observed[key], value) for key, value in components.items()],
            ]:
                a, b = np.array(a), np.array(b)
                comparisons[label] = {
                    "max_absolute_error": float(np.max(abs(a - b))),
                    "max_relative_error_scale_floor_1e-12": float(
                        np.max(abs(a - b) / np.maximum(abs(b), 1e-12))
                    ),
                }
                np.testing.assert_allclose(a, b, atol=2e-5, rtol=2e-6)
            row["comparisons"] = comparisons
            row["energy_kj_mol"] = result.potential_energy
            row["max_force_kj_mol_angstrom"] = max(
                float(np.linalg.norm(f)) for f in result.forces.values()
            )
            row["status"] = "parameterization_and_numerical_acceptance_passed"
            row["artifact_bytes"] = sum(
                p.stat().st_size for p in artifacts.iterdir() if p.is_file()
            )
            row["tool_package"] = r["tool_package"]
            row["force_field_data"] = r["force_field_data"]
            if name == "peo20":
                print("peo20: bounded workflow smoke", flush=True)
                t = time.perf_counter()
                cfg = WorkflowConfig(
                    psmiles,
                    dp,
                    str(root / "workflow"),
                    family,
                    "provided",
                    provided_charges=charges,
                    charge_source=SOURCE,
                    max_atoms=1000,
                    minimization=MinimizationOptions(
                        force_tolerance=0.1, max_iterations=500, max_evaluations=1000
                    ),
                    total_steps=4,
                    segment_steps=2,
                    recording_interval=1,
                    max_evaluations_per_segment=4,
                    max_frames_per_segment=3,
                )
                manifest = start_prepared_workflow(
                    cfg,
                    system,
                    prepared,
                    artifacts,
                    evidence="synthetic_software_test",
                    segments=2,
                )
                row["timings_seconds"]["workflow"] = time.perf_counter() - t
                row["workflow"] = {
                    "status": manifest["status"],
                    "accepted_step": manifest["accepted_step"],
                    "stages": manifest["stages"],
                }
                if manifest["segments"]:
                    report = analyze_workflow(root / "workflow")
                    export_analysis(report, root / "analysis")
                    row["analyzed_samples"] = report.payload["sample_count"]
                    row["analysis_verified"] = True
        except Exception as error:  # noqa: BLE001 -- retain actual acceptance failures
            row["failure"] = f"{type(error).__name__}: {error}"
            print(row["failure"], flush=True)
        row["process_high_water_rss_kib_cumulative"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss
        row["child_high_water_rss_kib_cumulative"] = resource.getrusage(
            resource.RUSAGE_CHILDREN
        ).ru_maxrss
        rows.append(row)
        storage.publish(
            root / "report.json",
            storage.json_bytes(
                {
                    "cases": rows,
                    "gates": acceptance_gates(rows),
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                }
            ),
            replace=True,
        )
        print(f"{name}: {row['status']}", flush=True)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--amberhome", required=True, type=Path)
    p.add_argument(
        "--require-workflow",
        action="store_true",
        help="Require completed dynamics and validated analysis as well as all numerical cases",
    )
    args = p.parse_args(argv)
    rows = run(args.output, args.amberhome)
    gates = acceptance_gates(rows)
    requested = (
        "workflow_dynamics_analysis"
        if args.require_workflow
        else "parameterization_numerical"
    )
    storage.publish(
        args.output / "report.json",
        storage.json_bytes(
            {
                "cases": rows,
                "gates": gates,
                "requested_gate": requested,
                "requested_gate_passed": gates[requested],
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
        replace=True,
    )
    return 0 if gates[requested] else 1


if __name__ == "__main__":
    raise SystemExit(main())
