"""Retain/reconstruct real J62 molecular gradients and nested model trials."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, require, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.angular_fit import Problem, metrics
from island.pcff_qm_pilot.bundle import pack, unpack
from scripts.analyze_phase_4j62 import controls, hmds_report, sample
from scripts.audit_phase_4j62 import build


def read(path):
    return json.loads(Path(path).read_text())


def retain(work, evidence):
    work, evidence = Path(work), Path(evidence)
    ledger = read(work / "ledger.json")
    with tempfile.TemporaryDirectory(prefix="j62-retain-") as tmp:
        out = Path(tmp)
        for name in ("requests", "outputs", "processes", "hmds-optimization"):
            if (work / name).exists():
                shutil.copytree(work / name, out / name)
        for name in ("ledger.json", "scientific-plan.json"):
            shutil.copyfile(work / name, out / name)
        write(evidence / "native-file-inventory.json", inventory(out))
        pack(out, evidence / "native-outputs.tar.xz")
    write(
        evidence / "normalized-reference.json",
        signed(
            {
                "schema": "island.j62-native-reference.v1",
                "samples": {
                    k: sample(work, k)
                    for k, v in ledger["calls"].items()
                    if v["status"] == "completed"
                },
                "units": {"energy": "kcal/mol", "forces": "kcal/mol/angstrom"},
                "conversion": CONVERSION,
                "origin": "executed native Psi4 gradients; no synthetic physical data",
                "unsuccessful_attempts": {
                    k: v
                    for k, v in ledger["calls"].items()
                    if v["status"] != "completed"
                },
            }
        ),
    )


def reconstruct(repo, evidence):
    repo, evidence = Path(repo), Path(evidence)
    with tempfile.TemporaryDirectory(prefix="j62-relocated-") as tmp:
        work, output = Path(tmp) / "native", Path(tmp) / "evidence"
        unpack(evidence / "native-outputs.tar.xz", work)
        require(
            inventory(work) == read(evidence / "native-file-inventory.json"),
            "J62_native_inventory_mutation",
        )
        shutil.copytree(evidence / "plan", output / "plan")
        for name in ("HMDS-center.json", "nonstretch-projection.json"):
            shutil.copyfile(evidence / name, output / name)
        normalized = read(evidence / "normalized-reference.json")
        reproduced = {k: sample(work, k) for k in normalized["samples"]}
        require(
            reproduced == normalized["samples"], "J62_native_reconstruction_changed"
        )
        control = controls(repo, output, work)
        require(
            control == read(evidence / "angular-torsional-control.json"),
            "J62_control_reconstruction_changed",
        )
        hmds_report(output, work)
        require(
            read(output / "HMDS-reference-results.json")
            == read(evidence / "HMDS-reference-results.json"),
            "J62_HMDS_reconstruction_changed",
        )
        base, initial = build(repo)
        audit = read(evidence / "training-weighting-tradeoff.json")
        residual, _ = base.system(initial)
        require(
            np.isclose(residual @ residual, audit["J61_soft_optimum_reconstructed_SS"]),
            "J62_J61_soft_optimum_changed",
        )
        feasible = []
        for row in audit["trials"]:
            if "parameters" not in row:
                continue
            values = np.array([row["parameters"][k] for k in base.names])
            r, _ = base.system(values)
            f = base.stationary.evaluate(values)["forces"]
            h, _ = base.mm_curvature(values)
            stationarity = float(np.sqrt(np.mean(f * f)))
            curvature = float(
                np.linalg.norm(h - base.curvature) / np.linalg.norm(base.curvature)
            )
            require(
                np.isclose(r @ r, row["unmodified_objective_SS"])
                and np.isclose(stationarity, row["stationary_force_RMS"])
                and np.isclose(curvature, row["curvature_relative_error"]),
                "J62_tradeoff_reconstruction_changed",
            )
            witness = stationarity <= 0.05 and curvature <= 0.05
            require(
                witness == row["stationarity_and_curvature_feasible"],
                "J62_false_feasibility_witness",
            )
            feasible.append(witness)
        require(
            any(feasible) == audit["stationarity_and_curvature_feasible_witness_found"],
            "J62_feasibility_claim_changed",
        )
        fine = control["steps"]["0.002"]
        p = Problem(
            base,
            control["Q_combined"],
            fine["combined_symmetric_projected"],
            fine["Cartesian_response"],
        )
        trials = {}
        for name in ("M0", "M1", "M2", "M3"):
            t = read(evidence / "nested" / f"{name}.json")
            v = np.array([t["parameters"][k] for k in p.names])
            require(
                np.all(v >= p.lower) and np.all(v <= p.upper_bounds),
                "J62_parameter_bound_mutation",
            )
            residual, jac = p.system(v)
            require(
                np.isclose(
                    residual @ residual, t["weighted_SS"], rtol=1e-10, atol=1e-8
                ),
                "J62_trial_objective_changed",
            )
            training = metrics(
                p.configurations, p.training, v, t["training_reference_offset"]
            )
            require(
                np.isclose(training["energy_RMSE"], t["training"]["energy_RMSE"])
                and np.isclose(training["force_RMS"], t["training"]["force_RMS"]),
                "J62_trial_training_changed",
            )
            n = t["parameter_count"]
            j = jac[:, :n]
            sv = np.linalg.svd(j / np.linalg.norm(j, axis=0), compute_uv=False)
            rank = int(sum(sv > sv[0] * 1e-10))
            condition = float(sv[0] / sv[-1])
            require(
                rank == t["rank"]
                and np.isclose(condition, t["column_scaled_condition"]),
                "J62_trial_rank_changed",
            )
            h, _, response = p.curvature(v)
            require(
                np.allclose(h, t["projected_curvature"], rtol=1e-10, atol=1e-8)
                and np.allclose(
                    response[:, 8:],
                    t["nonstretch_Cartesian_response"],
                    rtol=1e-10,
                    atol=1e-8,
                ),
                "J62_model_response_changed",
            )
            stationary = p.stationary.evaluate(v)["forces"]
            error = np.sqrt(
                np.mean(
                    (response[:, 8:] - np.array(fine["Cartesian_response"])) ** 2,
                    axis=0,
                )
            )
            require(
                np.isclose(np.sqrt(np.mean(stationary**2)), t["stationary_force_RMS"])
                and np.allclose(
                    error,
                    t["nonstretch_response_column_RMS_errors"],
                    rtol=1e-10,
                    atol=1e-8,
                )
                and np.isclose(
                    np.linalg.norm(h - p.target) / np.linalg.norm(p.target),
                    t["projected_curvature_relative_error"],
                ),
                "J62_diagnostic_metric_changed",
            )
            norms = np.linalg.norm(j, axis=0)
            _, singular, vt = np.linalg.svd(j / norms, full_matrices=False)
            inverse = np.where(
                singular > singular[0] * 1e-10, 1 / np.maximum(singular, 1e-150) ** 2, 0
            )
            covariance = (
                (vt.T * inverse)
                @ vt
                / np.outer(norms, norms)
                * (residual @ residual)
                / max(len(residual) - n, 1)
            )
            sigma = np.sqrt(np.maximum(np.diag(covariance), 0))
            correlation = np.divide(
                covariance,
                np.outer(sigma, sigma),
                out=np.zeros_like(covariance),
                where=np.outer(sigma, sigma) > 0,
            )
            require(
                np.allclose(
                    sigma, t["residual_covariance_sigma_proxy"], rtol=1e-8, atol=1e-8
                )
                and np.allclose(
                    correlation,
                    t["residual_covariance_correlation_proxy"],
                    rtol=1e-8,
                    atol=1e-8,
                ),
                "J62_uncertainty_proxy_changed",
            )
            trials[name] = {
                "parameters": n,
                "rank": rank,
                "weighted_SS": t["weighted_SS"],
                "training_qualifies": t["training_qualifies"],
            }
        return {
            "native_samples_reconstructed": len(reproduced),
            "control_identity": control["identity"],
            "HMDS_geometry_identity": read(evidence / "HMDS-center.json")[
                "geometry_identity"
            ],
            "trials": trials,
            "QM_invoked": False,
            "refitting_invoked": False,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("retain", "replay"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j62")
    )
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.operation == "retain":
        retain(args.work, args.evidence)
    else:
        print(json.dumps(reconstruct(args.repo, args.evidence), sort_keys=True))
