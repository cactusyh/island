"""Retain/reconstruct the actual J61 native campaign without invoking QM."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, inventory, require, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.bundle import pack, unpack
from island.pcff_qm_pilot.campaign import ENGINE_BOHR_ANGSTROM
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.stretch_fit import Problem, metrics
from scripts.analyze_phase_4j61 import curvature_report
from scripts.phase_4j61 import imported, read


def bound_protocol(sample):
    """Resolve redundant SCF metadata from explicit engine settings, never defaults.

    J60 places scf_type in settings; J61 also repeats it at protocol top level.
    An inconsistent repetition or absent engine setting is rejected.
    """
    protocol = dict(sample["protocol"])
    settings = protocol["settings"]
    require("scf_type" in settings, "J61_unbound_reference_SCF_type")
    require(
        protocol.get("scf_type", settings["scf_type"]) == settings["scf_type"],
        "J61_conflicting_reference_SCF_type",
    )
    protocol["scf_type"] = settings["scf_type"]
    return protocol


def retain(work, evidence):
    work, evidence = Path(work), Path(evidence)
    with tempfile.TemporaryDirectory(prefix="j61-retain-") as temp:
        root = Path(temp)
        for name in ("requests", "outputs", "optimization", "logs", "attempts"):
            if (work / name).exists():
                shutil.copytree(work / name, root / name)
        for name in ("ledger.json", "scientific-plan.json"):
            shutil.copyfile(work / name, root / name)
        write(evidence / "native-file-inventory.json", inventory(root))
        pack(root, evidence / "native-outputs.tar.xz")
    keys = list(read(work / "ledger.json")["calls"])
    write(
        evidence / "normalized-reference.json",
        signed(
            {
                "schema": "island.j61-native-reference.v1",
                "samples": {key: imported(work, key) for key in keys},
                "units": {"energy": "kcal/mol", "forces": "kcal/mol/angstrom"},
                "conversion": CONVERSION,
                "engine_coordinate_bohr_angstrom": ENGINE_BOHR_ANGSTROM,
                "charge": 0,
                "multiplicity": 1,
                "origin": "executed native Psi4 calculations; no synthetic QM",
            }
        ),
    )


def reconstruct(repo, evidence):
    """Re-read native matrices, reconstruct curvature, then evaluate retained fit.

    Relocated evidence paths are supported. Bootstrap and QM are never rerun.
    """
    repo, evidence = Path(repo), Path(evidence)
    with tempfile.TemporaryDirectory(prefix="j61-relocated-") as temp:
        temp = Path(temp)
        work, output = temp / "native", temp / "evidence"
        unpack(evidence / "native-outputs.tar.xz", work)
        require(
            inventory(work) == read(evidence / "native-file-inventory.json"),
            "J61_native_inventory_mutation",
        )
        shutil.copytree(evidence / "plan", output / "plan")
        for name in ("stationary-center.json", "stretching-subspace.json"):
            shutil.copyfile(evidence / name, output / name)
        reference, samples = curvature_report(work, output)
        require(
            reference == read(evidence / "projected-Hessian.json"),
            "J61_curvature_reconstruction_changed",
        )
        normalized = read(evidence / "normalized-reference.json")
        reconstructed = {key: imported(work, key) for key in normalized["samples"]}
        require(
            reconstructed == normalized["samples"], "J61_native_reconstruction_changed"
        )
        summary = {
            "native_samples_reconstructed": len(reconstructed),
            "projected_Hessian_identity": reference["identity"],
            "stationary_geometry_identity": read(evidence / "stationary-center.json")[
                "geometry_identity"
            ],
            "fit_reconstructed": False,
            "QM_invoked": False,
        }
        if (evidence / "baseline-versus-joint-fit.json").exists():
            old = repo / "docs/evidence/phase_4j60"
            design = read(old / "joint-target-design.json")
            prior = read(old / "stages/tight-d-joint-fit/normalized-reference.json")
            training = [
                {**prior["samples"]["cc5/" + s["id"]], "key": "J60/" + s["id"]}
                for s in design["selected_training"]
            ]
            training.extend(samples[f"fresh/training-{i}"] for i in range(4))
            validation = [samples[f"fresh/validation-{i}"] for i in range(4)]
            model = load_model(repo)
            center = read(evidence / "stationary-center.json")["center"]
            require(
                all(
                    bound_protocol(s) == bound_protocol(center)
                    for s in [*training, *validation]
                ),
                "J61_undeclared_mixed_reference_protocol",
            )
            q = read(evidence / "stretching-subspace.json")["Q"]
            target = reference["steps"]["0.002"]["projected_symmetric_kcal_mol_A2"]
            problem = Problem(model, training, center, q, target)
            fit = read(evidence / "baseline-versus-joint-fit.json")
            for label in (
                "baseline_adjustable_joint",
                "coupling_only_refit_same_constraints",
            ):
                candidate = fit[label]
                p = np.array([candidate["parameters"][k] for k in problem.names])
                require(
                    np.all(p >= problem.lower) and np.all(p <= problem.upper_bounds),
                    "J61_retained_parameter_outside_bounds",
                )
                residual, jacobian = problem.system(p)
                h, _ = problem.mm_curvature(p)
                require(
                    np.isclose(
                        residual @ residual,
                        candidate["weighted_residual_SS"],
                        rtol=1e-9,
                        atol=1e-7,
                    ),
                    "J61_joint_objective_changed",
                )
                require(
                    np.allclose(h, candidate["projected_curvature"], atol=1e-8),
                    "J61_joint_curvature_changed",
                )
                actual_stationary = problem.stationary.evaluate(p)
                require(
                    np.isclose(
                        np.sqrt(np.mean(actual_stationary["forces"] ** 2)),
                        candidate["stationary_force_RMS"],
                        atol=1e-9,
                    )
                    and np.isclose(
                        np.linalg.norm(h - problem.curvature)
                        / np.linalg.norm(problem.curvature),
                        candidate["projected_curvature_relative_Frobenius"],
                        atol=1e-9,
                    ),
                    "J61_stationary_curvature_error_changed",
                )
                for split, records in (
                    ("training", training),
                    ("fresh_validation", validation),
                ):
                    actual = metrics(
                        model, records, p, candidate["training_reference_offset"]
                    )
                    expected = candidate[split]
                    require(
                        np.isclose(
                            actual["energy_RMSE"], expected["energy_RMSE"], atol=1e-8
                        )
                        and np.isclose(
                            actual["force_RMS"], expected["force_RMS"], atol=1e-8
                        ),
                        "J61_joint_errors_changed",
                    )
                if label == "baseline_adjustable_joint":
                    require(
                        np.allclose(
                            jacobian,
                            read(evidence / "joint-parameter-jacobian.json")["matrix"],
                            atol=1e-8,
                        ),
                        "J61_parameter_Jacobian_changed",
                    )
            summary.update(
                fit_reconstructed=True,
                parameter_checksum=digest(
                    fit["baseline_adjustable_joint"]["parameters"]
                ),
                accepted_coefficients=len(fit["accepted_coefficients"]),
            )
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("retain", "replay"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j61")
    )
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.operation == "retain":
        require(args.work is not None, "work_required")
        retain(args.work, args.evidence)
    else:
        print(json.dumps(reconstruct(args.repo, args.evidence), sort_keys=True))


if __name__ == "__main__":
    main()
