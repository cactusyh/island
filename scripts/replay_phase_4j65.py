"""Retain the bounded campaign and replay native intake/fitting after relocation."""

import argparse
import shutil
import tempfile
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, pin, require, signed, write
from island.pcff_qm_pilot.bundle import pack, unpack
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j63 import sample
from scripts.analyze_phase_4j65 import (
    controls,
    independent_value,
    primary_gates,
    trial,
    validation,
)

REPORTS = (
    "targeted-controls.json",
    "primary-numerical-gates.json",
    "fitting-dataset.json",
    "fitting-numerical-checks.json",
    "trial-fit.json",
    "block-uncertainty.json",
    "untouched-validation.json",
)


def retain(repo, root, work):
    with tempfile.TemporaryDirectory(prefix="j65-native-collect-") as tmp:
        target = Path(tmp)
        for name in ("requests", "outputs", "processes", "resources"):
            if (work / name).exists():
                shutil.copytree(work / name, target / name)
        for name in (
            "scientific-plan.json",
            "jobs.json",
            "ledger.json",
            "attempts.json",
        ):
            shutil.copyfile(work / name, target / name)
        write(root / "native-file-inventory.json", inventory(target))
        pack(target, root / "native-outputs.tar.xz")
    plan = read(root / "plan/scientific-plan.json")
    attempts = read(work / "attempts.json")["calls"]
    native = read(work / "ledger.json")["calls"]
    require(len(attempts) <= 32, "J65_attempt_cap")
    refs = {}
    checks = {}
    resources = {}
    for key, row in native.items():
        if row["status"] == "completed":
            refs[key] = sample(work, key)
            require(
                refs[key]["graph_identity"] == plan["graph_identity"],
                "J65_chemical_scope",
            )
            checks[key] = independent_value(work, key, refs[key])
            file = work / "resources" / (key + "-usage.json")
            if file.exists():
                resources[key] = read(file)
    write(
        root / "normalized-native-controls-and-observations.json",
        signed({"samples": refs, "protocols_not_pooled": True}),
    )
    write(root / "independent-native-checks.json", signed({"samples": checks}))
    write(
        root / "compute-accounting.json",
        signed(
            {
                "counted_attempts": len(attempts),
                "cap": 32,
                "successful_native_gradients": len(refs),
                "attempts": attempts,
                "native_gradient_starts": len(native),
                "resources": resources,
                "elapsed_seconds": {
                    k: r.get("elapsed_seconds") for k, r in native.items()
                },
                "unexecuted_requests": [
                    k for k in read(work / "jobs.json")["jobs"] if k not in native
                ],
                "new_SiH_QM": 0,
                "PCFF_MD": False,
                "new_registered_profiles": [],
                "reused_J64_native_executions": 7,
                "reused_records_have_new_execution_identity": False,
            }
        ),
    )
    fit = read(root / "trial-fit.json")
    tested = (root / "untouched-validation.json").exists()
    passed = (
        tested
        and read(root / "untouched-validation.json")["method_specific_target_passed"]
    )
    write(
        root / "fitting-decision.json",
        signed(
            {
                "scientific_plan_identity": plan["identity"],
                "native_calculation_success": all(
                    r["status"] == "completed" for r in attempts.values()
                ),
                "internal_numerical_consistency": fit["numerical_eligible"],
                "fit_executed": fit["fit_executed"],
                "trial_physical_coefficients": 58 if fit["fit_executed"] else 0,
                "energy_offset_is_not_physical_coefficient": True,
                "training_performance_passed": fit.get(
                    "training_performance_passed", False
                ),
                "uncertainty_qualified": fit.get("uncertainty_qualified", False),
                "untouched_validation_executed": tested,
                "method_specific_target_reproduction": "passed"
                if passed
                else "not_accepted",
                "accepted_method_specific_research_candidates": [fit["identity"]]
                if passed
                else [],
                "force_qualified_transferable_candidates": [],
                "validation_activation_reason": "Frozen training/numerical/identifiability/uncertainty gates passed"
                if fit["validation_authorized"]
                else "Training candidate did not pass every frozen gate; see trial-fit and block-uncertainty, reservation preserved.",
                "basis_sensitivity": "Historical failed/unestablished qualifications retained; QZVP controls separate from denseTZVP target fitting.",
                "chemical_scope": "Exact neutral-singlet explicit-H27atom HMDS, graph d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd, declared PBE0/def2TZVP dense-grid target only. Other sites are nearly symmetry-equivalent; no independent conformational basin, SiH or cappedDP2 transfer.",
                "parameter_completeness": "incomplete",
                "backend_correctness": "bounded evaluator/importer controls only",
                "physical_validation": "not_established",
                "cross_engine_equivalence": "not_claimed",
                "production_validated": False,
                "simulation_readiness": "not_established",
                "ordinary_API": 'forcefield="pcff"',
                "automatic_assignment": False,
                "original_author_artifact_status": "not_yet_supplied",
                "confirmation_authority": "28 user-attested author confirmations preserved",
                "source_rows_charges_or_equilibria_changed": False,
                "capped_DP2_dependency": "Separate exact graph/source/charge inventory and data for interior sio(osi,osi,c3,c3), O-Si-O and Si-O-Si-O; HMDS cannot supply these environments.",
            }
        ),
    )


def replay(repo, root):
    with tempfile.TemporaryDirectory(prefix="j65-fresh-relocation-") as tmp:
        temp = Path(tmp)
        work = temp / "native"
        output = temp / "evidence"
        unpack(root / "native-outputs.tar.xz", work)
        require(
            read(work / "scientific-plan.json")
            == read(root / "plan/scientific-plan.json")
            and read(work / "jobs.json") == read(root / "plan/jobs.json"),
            "J65_native_frozen_request_contract_mutation",
        )
        requests = read(root / "plan/requests.json")
        actual_requests = {
            str(p.relative_to(work / "requests"))[:-5]: read(p)
            for p in (work / "requests").rglob("*.json")
        }
        require(actual_requests == requests, "J65_native_request_inventory_mutation")
        require(
            set(read(work / "ledger.json")["calls"]) <= set(requests),
            "J65_unplanned_native_request",
        )
        require(
            inventory(work) == read(root / "native-file-inventory.json"),
            "J65_native_inventory_mutation",
        )
        shutil.copytree(root / "plan", output / "plan")
        raw = read(root / "normalized-native-controls-and-observations.json")["samples"]
        checks = read(root / "independent-native-checks.json")["samples"]
        for key, expected in raw.items():
            value = sample(work, key)
            require(value == expected, "J65_native_sample_mutation")
            require(
                independent_value(work, key, value) == checks[key],
                "J65_independent_units_mutation",
            )
        attempts = read(work / "attempts.json")["calls"]
        require(len(attempts) <= 32, "J65_attempt_cap")
        for file in (work / "processes").rglob("result.json"):
            record = read(file)
            require(
                all(
                    pin(file.parent / name) == sha
                    for name, sha in record["files"].items()
                ),
                "J65_process_capture_mutation",
            )
        controls(repo, output, work)
        primary_gates(repo, output, work)
        if (root / "trial-fit.json").exists():
            trial(repo, output, work)
        if (root / "untouched-validation.json").exists():
            validation(repo, output, work)
        for name in REPORTS:
            require(
                (root / name).exists() == (output / name).exists(),
                "J65_replay_report_coverage",
                file=name,
            )
            if (root / name).exists():
                require(
                    read(root / name) == read(output / name),
                    "J65_numerical_or_fit_replay_changed",
                    file=name,
                )
        if (root / "actual-fitting-system.npz").exists():
            with (
                np.load(root / "actual-fitting-system.npz", allow_pickle=False) as a,
                np.load(output / "actual-fitting-system.npz", allow_pickle=False) as b,
            ):
                require(
                    a.files == b.files
                    and all(np.array_equal(a[n], b[n]) for n in a.files),
                    "J65_actual_Jacobian_replay",
                )
        return {
            "relocated_native_directory": True,
            "QM_invoked": False,
            "native_gradients": len(raw),
            "counted_attempts": len(attempts),
            "exact_replay": True,
            "fit_executed": read(root / "trial-fit.json")["fit_executed"]
            if (root / "trial-fit.json").exists()
            else False,
        }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=["retain", "replay"])
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j65"))
    p.add_argument("--work", type=Path)
    a = p.parse_args()
    repo, root = a.repo.resolve(), a.evidence.resolve()
    if a.operation == "retain":
        retain(repo, root, a.work.resolve())
    else:
        print(replay(repo, root))
