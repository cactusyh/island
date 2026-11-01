"""Retain and independently replay J63 native HMDS evidence without invoking QM."""

import argparse
import json
import lzma
import shutil
import tempfile
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, require, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.bundle import pack, unpack
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from scripts.analyze_phase_4j63 import analyze, independent_native, read, sample
from scripts.prepare_phase_4j63 import optimizer_outcome, stationary


def retain(work, evidence):
    ledger = read(work / "ledger.json")
    with tempfile.TemporaryDirectory(prefix="j63-retain-") as tmp:
        out = Path(tmp)
        for name in ("requests", "outputs", "processes", "optimization"):
            if (work / name).exists():
                shutil.copytree(work / name, out / name)
        for name in ("ledger.json", "scientific-plan.json"):
            shutil.copyfile(work / name, out / name)
        write(evidence / "native-file-inventory.json", inventory(out))
        pack(out, evidence / "native-outputs.tar.xz")
    samples = {
        k: sample(work, k)
        for k, v in ledger["calls"].items()
        if v["status"] == "completed"
    }
    independent = {k: independent_native(work, k, s) for k, s in samples.items()}
    write(
        evidence / "normalized-reference.json",
        signed(
            {
                "schema": "island.j63-native-HMDS-reference.v1",
                "samples": samples,
                "units": {"energy": "kcal/mol", "forces": "kcal/mol/angstrom"},
                "conversion": CONVERSION,
                "origin": "executed native Psi4 gradients; raw native outputs archived separately",
                "unsuccessful_attempts": {
                    k: v
                    for k, v in ledger["calls"].items()
                    if v["status"] != "completed"
                },
                "protocols_kept_separate": True,
            }
        ),
    )
    write(
        evidence / "independent-native-checks.json",
        signed(
            {
                "samples": independent,
                "samples_reconstructed": len(independent),
                "independent_path": "Parse @DF-RKS text energy and all81 explicitly indexed ASCII gradient entries; independent constants/sign conversion, no production importer or npy reader.",
                "maximum_force_difference": max(
                    r["maximum_force_difference_kcal_mol_A"]
                    for r in independent.values()
                ),
            }
        ),
    )


def reconstruct(repo, evidence):
    with tempfile.TemporaryDirectory(prefix="j63-relocated-") as tmp:
        work, output = Path(tmp) / "native", Path(tmp) / "evidence"
        unpack(evidence / "native-outputs.tar.xz", work)
        require(
            inventory(work) == read(evidence / "native-file-inventory.json"),
            "J63_native_inventory_mutation",
        )
        shutil.copytree(evidence / "plan", output / "plan")
        expected = read(evidence / "normalized-reference.json")
        actual = {k: sample(work, k) for k in expected["samples"]}
        require(actual == expected["samples"], "J63_native_reconstruction_changed")
        independent = {k: independent_native(work, k, v) for k, v in actual.items()}
        require(
            independent == read(evidence / "independent-native-checks.json")["samples"],
            "J63_independent_native_check_changed",
        )
        center = stationary(repo, output, work)
        require(
            center == read(evidence / "stationary-center.json"),
            "J63_stationarity_reconstruction_changed",
        )
        write(output / "normalized-reference.json", expected)
        if (evidence / "optimizer-model-diagnosis.json").exists():
            from scripts.trajectory_phase_4j63 import diagnose

            diagnosis = diagnose(repo, output)
            saved = read(evidence / "optimizer-model-diagnosis.json")
            require(
                diagnosis["samples"] == saved["samples"]
                and diagnosis["actual_energy_force_Jacobian"]["rank"]
                == saved["actual_energy_force_Jacobian"]["rank"]
                and np.isclose(
                    diagnosis["optimistic_unconstrained_lower_bound"][
                        "weighted_squared_residual"
                    ],
                    saved["optimistic_unconstrained_lower_bound"][
                        "weighted_squared_residual"
                    ],
                    rtol=1e-8,
                ),
                "J63_actual_trajectory_model_diagnosis_changed",
            )
        paths = None
        if (evidence / "deformation-reference-results.json").exists():
            paths = analyze(
                load_model(repo),
                read(evidence / "plan/scientific-plan.json"),
                center["center"],
                actual,
            )
            require(
                paths == read(evidence / "deformation-reference-results.json"),
                "J63_path_reconstruction_changed",
            )
            write(output / "deformation-reference-results.json", paths)
            write(output / "normalized-reference.json", expected)
            from scripts.curvature_phase_4j63 import audit as curvature_audit
            from scripts.design_phase_4j63 import design

            proposal, compatibility = design(repo, output)
            saved_design = read(evidence / "joint-design.json")
            require(
                proposal["current"]["rank"] == saved_design["current"]["rank"]
                and np.isclose(
                    proposal["current"]["scaled_condition"],
                    saved_design["current"]["scaled_condition"],
                    rtol=1e-8,
                )
                and [v["id"] for v in proposal["selected_additional_configurations"]]
                == [
                    v["id"] for v in saved_design["selected_additional_configurations"]
                ],
                "J63_design_reconstruction_changed",
            )
            saved_compatibility = read(evidence / "source-baseline-compatibility.json")
            require(
                np.isclose(
                    compatibility["optimistic_joint_column_space_lower_bound"][
                        "weighted_squared_lower_bound"
                    ],
                    saved_compatibility["optimistic_joint_column_space_lower_bound"][
                        "weighted_squared_lower_bound"
                    ],
                    rtol=1e-8,
                ),
                "J63_baseline_lower_bound_changed",
            )
            curvature = curvature_audit(repo, output)
            saved_curvature = read(evidence / "source-baseline-path-curvature.json")
            for key, value in curvature["paths"].items():
                for h, row in value["steps"].items():
                    require(
                        np.isclose(
                            row["QM_minus_model_curvature"],
                            saved_curvature["paths"][key]["steps"][h][
                                "QM_minus_model_curvature"
                            ],
                            atol=1e-8,
                            rtol=1e-8,
                        ),
                        "J63_unknown_null_curvature_changed",
                    )
        ledger = read(work / "ledger.json")["calls"]
        require(
            set(actual)
            == {k for k, row in ledger.items() if row["status"] == "completed"},
            "J63_completed_native_sample_coverage",
        )
        plan = read(evidence / "plan/scientific-plan.json")
        executed_sha = read(evidence / "plan/execution-pins.json")["files"][
            "src/island/pcff_qm_pilot/hmds_engine.py"
        ]
        coordinator_sha = sha256(
            lzma.decompress(
                (evidence / "plan/executed-phase-coordinator.py.xz").read_bytes()
            )
        ).hexdigest()
        for command_file in (work / "processes").rglob("command.json"):
            require(
                read(command_file)["coordinator_sha256"] == coordinator_sha,
                "J63_executed_coordinator_changed",
            )
        for key in actual:
            raw = read(work / "outputs" / key / "result.json")
            basis = plan["basis_composition"][raw["basis"].lower()]
            require(
                raw["basis_functions"] == basis["nbf"]
                and raw["basis_shells"] == basis["nshell"]
                and raw["adapter_sha256"] == executed_sha
                and raw["adapter"] == "island.psi4-hmds-j63.v1",
                "J63_executed_adapter_or_contracted_basis_changed",
            )
        require(
            len(ledger) <= 80
            and sum(k.startswith("optimization/") for k in ledger) <= 40,
            "J63_budget_changed",
        )
        # All optimizer input geometries are checked against their independent native records.
        opt, _ = optimizer_outcome(work)
        if opt["status"] == "converged":
            hist = read(work / "optimization/history.json")
            require(
                len(hist["energy"]) == len(opt["calls"]),
                "J63_optimizer_history_coverage",
            )
            require(
                np.allclose(
                    hist["energy"],
                    [
                        read(work / "outputs" / k / "result.json")["energy_hartree"]
                        for k in opt["calls"]
                    ],
                    atol=1e-12,
                    rtol=0,
                ),
                "J63_optimizer_history_energy",
            )
        return {
            "native_samples_reconstructed": len(actual),
            "attempted_gradients": len(ledger),
            "stationary_candidate": center["stationary_candidate"],
            "optimizer_status": opt["status"],
            "bounded_path_precision_passed": None
            if paths is None
            else paths["bounded_precision_passed"],
            "QM_invoked": False,
            "fit_invoked": False,
            "relocated": True,
        }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("retain", "replay"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    p.add_argument("--work", type=Path)
    a = p.parse_args()
    if a.operation == "retain":
        retain(a.work, a.evidence)
    else:
        print(json.dumps(reconstruct(a.repo, a.evidence), sort_keys=True))
