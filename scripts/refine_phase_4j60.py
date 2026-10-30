"""Versioned tight-d reference refinement and conditional identifiable pilot fit.

Reuses the existing campaign/native importer/solver; preserves all primary
requests and failures. The tight-d plan predates all refinement execution.
"""

import argparse
import json
import shutil
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.hamiltonian_analysis import convergence
from island.pcff_qm_pilot.hamiltonian_fit import (
    bootstrap_training,
    discrepancy_diagnostics,
    joint_fit,
)
from scripts.finalize_phase_4j60 import run_stage, stage_plan
from scripts.residuals_phase_4j60 import residuals

MODULE = "island.pcff_qm_pilot.hamiltonian_refinement_engine"


def summarize(repo, evidence):
    """Interpret the executed fit; never alter a parameter or rerun QM."""
    from importlib.metadata import version

    repo, evidence = Path(repo), Path(evidence)
    path = evidence / "experimental-joint-fit.json"
    require(path.exists(), "executed_fit_required")
    fit = json.loads(path.read_text())
    require(
        fit["held_out_thresholds_passed"] is False,
        "J60_rejection_interpretation_changed",
    )
    data_path = evidence / "stages/tight-d-joint-fit/normalized-reference.json"
    data = json.loads(data_path.read_text())
    design = json.loads((evidence / "joint-target-design.json").read_text())
    model = load_model(repo)
    training = [
        {**data["samples"]["cc5/" + s["id"]], **s, "conformer_identity": s["id"]}
        for s in design["selected_training"]
    ]
    diagnostics = signed(discrepancy_diagnostics(model, training, fit))
    require(
        [p["probe"] for p in diagnostics["probes"][:2]]
        == ["quartic_bond:osi,sio:delta^1", "quartic_bond:hsi,sio:delta^1"],
        "J60_baseline_probe_order_changed",
    )
    write(evidence / "model-discrepancy-diagnostics.json", diagnostics)
    bootstrap = json.loads((evidence / "training-bootstrap.json").read_text())
    specs = json.loads((evidence / "parameter-specifications.json").read_text())
    sigma = dict(zip(bootstrap["columns"], bootstrap["sigma"], strict=True))
    weak = [k for k, v in fit["coefficients"].items() if sigma[k] > 0.5 * abs(v)]
    provenance = signed(
        {
            "schema": "island.j60-rejected-research-fit-provenance.v1",
            "origin": "fitted_external experimental research-extension candidates; not original/native source rows",
            "accepted_parameters": [],
            "validated_scope": "none; the two independent control conformers reject this fit",
            "numerical_scope": "exact nine-atom neutral singlet Si-H disiloxane graph only; no methyl-PDMS transfer",
            "graph_identity": model["graph_identity"],
            "source_model_identity": model["identity"],
            "dataset_identity": data["identity"],
            "dataset_sha256": pin(data_path),
            "design_identity": design["identity"],
            "training_jobs": {
                s["id"]: data["samples"]["cc5/" + s["id"]]["job_identity"]
                for s in design["selected_training"]
            },
            "untouched_validation_jobs": {
                s["id"]: data["samples"]["cc5/" + s["id"]]["job_identity"]
                for s in design["held_out"]
            },
            "coefficient_checksum": digest(fit["coefficients"]),
            "fit_record_sha256": pin(path),
            "parameter_specifications_sha256": pin(
                evidence / "parameter-specifications.json"
            ),
            "optimizer": {
                **specs["optimizer"],
                "scipy_version": version("scipy"),
                "numpy_version": version("numpy"),
            },
            "parameters": [
                {
                    "name": p["name"],
                    "value": fit["coefficients"][p["name"]],
                    "units": p["units"],
                    "bootstrap_sigma": sigma[p["name"]],
                    "accepted": False,
                }
                for p in specs["parameters"]
            ],
            "weak_relative_uncertainty_parameters": weak,
            "statistical_uncertainty_scope": "Residual-only Cartesian weighted covariance is a linearized proxy, not a formal independent-force confidence interval. Whole-geometry bootstrap excludes source, method, basis and transfer-systematic errors.",
            "model_discrepancy_identity": diagnostics["identity"],
            "source_rows_changed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(evidence / "fitted-candidate-provenance.json", provenance)
    previous = evidence / "initial-compute-fitting-decision.json"
    if not previous.exists():
        shutil.copyfile(evidence / "fitting-decision.json", previous)
    final = json.loads(previous.read_text())
    final.pop("identity")
    final["accepted_coefficients"] = []
    final["model_discrepancy_identity"] = diagnostics["identity"]
    final["fitting_provenance_identity"] = provenance["identity"]
    final["smallest_next_requirement"] = (
        f"The complete20-column fit leaves training force RMS {fit['training']['force_RMS']:.9g} and untouched RMS {fit['untouched_validation']['force_RMS']:.9g} kcal/mol/A with no active bounds. Training-only source-baseline probes project {diagnostics['probes'][0]['fraction_of_training_weighted_residual_SS_projected']:.3%} (Si-O linear stretch) and {diagnostics['probes'][1]['fraction_of_training_weighted_residual_SS_projected']:.3%} (Si-H linear stretch) individually; these overlapping projections are not additive and create no new coefficients. Independently constrain the Si-O and Si-H equilibrium-force baselines with a PBE0/cc-pV(5+d)Z stationary-center geometry and a projected eight-bond stretching Hessian (2 Si-O,6 Si-H), or matching independently justified numerical baseline rows, before a new versioned joint research fit. Any baseline refit must propagate its equilibrium values through all Class-II cross terms and reserve new conformers because the original hold-outs have now been examined. Repeating the same coupling-only fit or adding its grid points cannot remove the demonstrated training residual. No extra author engine-setting confirmation is needed for the declared research Hamiltonian."
    )
    final = signed(final)
    write(evidence / "fitting-decision.json", final)
    return final


def new_stage(parent, geoms, protocols, points, budget, authority, role):
    raw = stage_plan(parent, geoms, protocols, points, budget, authority)
    return signed(
        {
            **{k: v for k, v in raw.items() if k != "identity"},
            "adapter_module": MODULE,
            "stage_role": role,
        }
    )


def decision(repo, evidence, reference, data, fit=None, stages=None):
    primary = json.loads((evidence / "primary-fitting-decision.json").read_text())
    gates = {
        r["effect"]: r["predeclared_gate_passed"]
        for r in reference["comparisons"]
        if r["predeclared_gate_passed"] is not None
    }
    basis = next(r for r in reference["comparisons"] if r["effect"] == "cardinal_QZ_5Z")
    updated = {
        **{k: v for k, v in primary.items() if k != "identity"},
        "selected_reference": "PBE0/cc-pV(5+d)Z versioned tight-d pilot protocol",
        "reference_gates": gates,
        "basis_gate_numbers": basis,
        "primary_reference_gates_preserved": primary["reference_gates"],
        "primary_reference_file": "primary-reference-convergence.json",
        "stages_executed": stages or [],
        "fit_completed": fit is not None,
        "fit_validated": False,
        "fit_readiness": "experimental_fit_executed"
        if fit
        else "blocked_reference_accuracy",
        "fit_completion": "completed" if fit else "not_run",
        "fit_validation": "not_run"
        if not fit
        else "failed"
        if not fit["held_out_thresholds_passed"]
        else "not_established",
        "QM_reference_adequacy": "passed_declared_pilot_criteria"
        if fit
        else "not_established",
        "fitted_coefficients": list(fit["coefficients"]) if fit else [],
    }
    if fit:
        updated["decision"] = (
            "experimental_fit_completed"
            if fit["held_out_thresholds_passed"]
            else "experimental_fit_rejected_by_untouched_validation"
        )
        updated["smallest_next_requirement"] = (
            "Use the recorded untouched-conformer residuals, coefficient correlations and reference sensitivities to test this exact Si-H research model; conformer-group uncertainty and protocol validation remain required. No Si-C transfer is authorized."
        )
    else:
        updated["decision"] = "no_fit"
        bad = {k: v for k, v in gates.items() if not v}
        if "cardinal_QZ_5Z" in bad:
            updated["smallest_next_requirement"] = (
                f"Tight-d QZ/5Z still fails unchanged criteria: max relative deformation energy difference {basis['max_relative_energy_difference']:.9g} kcal/mol; gradient RMS {basis['max_gradient_RMS_difference']:.9g} kcal/mol/A. The exact worst bond/angle corners and mixed-derivative differences in refinement-convergence.json specify a bounded cc-pV(6+d)Z or independent equivalent derivative reference requirement; no missing author convention remains for the declared research Hamiltonian."
            )
        elif bad:
            updated["smallest_next_requirement"] = (
                "Resolve the exact retained corrected-basis numerical contrast: "
                + ",".join(bad)
                + ". Preserve the existing thresholds and failed outputs."
            )
        else:
            updated["smallest_next_requirement"] = (
                "The separate corrected-basis held-out reference-accuracy report identifies the two conformers requiring converged orbital-basis derivatives before fitting; no physical coefficient is inferred."
            )
    return signed(updated)


def basis_effects(primary, corrected):
    records = []
    for protocol in ("ccq", "cc5"):
        a0, b0 = (
            primary["samples"][protocol + "/center"],
            corrected["samples"][protocol + "/center"],
        )
        for key, b in corrected["samples"].items():
            if not key.startswith(protocol + "/"):
                continue
            a = primary["samples"][key]
            df = np.array(b["forces"]) - np.array(a["forces"])
            rel = df - (np.array(b0["forces"]) - np.array(a0["forces"]))
            records.append(
                {
                    "geometry": key,
                    "basis_a": a["protocol"]["basis"],
                    "basis_b": b["protocol"]["basis"],
                    "relative_energy_difference": (b["energy"] - b0["energy"])
                    - (a["energy"] - a0["energy"]),
                    "gradient_RMS_difference": float(np.sqrt(np.mean(df * df))),
                    "deformation_gradient_RMS_difference": float(
                        np.sqrt(np.mean(rel * rel))
                    ),
                }
            )
    return signed(
        {
            "schema": "island.j60-tight-d-sensitivity.v1",
            "records": records,
            "interpretation": "Fixed-cardinality tight-d sensitivity, separately from the corrected QZ/5Z cardinal test; no CBS claim.",
        }
    )


def held_out_accuracy(fit_data, basis_data, design, criteria):
    """Reference adequacy on reserved conformers, with no fitted values involved."""
    errors = []
    for s in design["held_out"]:
        a, b = (fit_data["samples"][p + "/" + s["id"]] for p in ("ccq", "cc5"))
        ca, cb = (basis_data["samples"][p + "/center"] for p in ("ccq", "cc5"))
        errors.append(
            {
                "id": s["id"],
                "relative_energy_difference": (b["energy"] - cb["energy"])
                - (a["energy"] - ca["energy"]),
                "gradient_RMS_difference": float(
                    np.sqrt(
                        np.mean((np.array(a["forces"]) - np.array(b["forces"])) ** 2)
                    )
                ),
            }
        )
    return {
        "records": errors,
        "predeclared_basis_gate_passed": bool(errors)
        and all(
            abs(e["relative_energy_difference"]) <= criteria["relative_energy_kcal_mol"]
            and e["gradient_RMS_difference"] <= criteria["gradient_RMS_kcal_mol_A"]
            for e in errors
        ),
    }


def reference_gates_pass(reference):
    gates = {r["effect"]: r for r in reference["comparisons"]}
    required = ("cardinal_QZ_5Z", "numerical_grid", "SCF_threshold")
    return all(
        k in gates and gates[k]["predeclared_gate_passed"] is True for k in required
    ) and all(
        set(gates[k]["mixed_derivative_differences"]) == {"angle-angle", "bond-bond"}
        for k in ("cardinal_QZ_5Z", "numerical_grid")
    )


def replay_fit(repo, evidence, stages):
    path = evidence / "experimental-joint-fit.json"
    if not path.exists():
        require(
            not json.loads((evidence / "fitting-decision.json").read_text())[
                "fit_completed"
            ],
            "fit_claim_without_coefficients",
        )
        return
    require("tight-d-joint-fit" in stages, "fit_native_data_missing")
    fd = stages["tight-d-joint-fit"]
    design = json.loads((evidence / "joint-target-design.json").read_text())
    model = load_model(repo)
    require("tight-d-numerics" in stages, "fit_numerical_controls_missing")
    corrected = signed(
        {
            **{k: v for k, v in stages["tight-d-basis"].items() if k != "identity"},
            "samples": {
                **stages["tight-d-basis"]["samples"],
                **stages["tight-d-numerics"]["samples"],
            },
        }
    )
    reference = convergence(repo, evidence, corrected)
    require(reference_gates_pass(reference), "reconstructed_reference_gate_failed")
    plan = json.loads((evidence / "plan/reference-plan.json").read_text())
    accuracy = held_out_accuracy(
        fd, stages["tight-d-basis"], design, plan["accuracy_criteria"]["basis"]
    )
    require(
        accuracy
        == json.loads(
            (evidence / "tight-d-held-out-reference-accuracy.json").read_text()
        ),
        "held_out_reference_reconstruction_mutation",
    )
    require(accuracy["predeclared_basis_gate_passed"], "held_out_reference_inadequate")

    def sample(s):
        return {
            **fd["samples"]["cc5/" + s["id"]],
            "split": s["split"],
            "graph_identity": model["graph_identity"],
            "atom_order": model["atom_order"],
            "conformer_identity": s["id"],
        }

    training = [sample(s) for s in design["selected_training"]]
    validation = [sample(s) for s in design["held_out"]]
    criteria = json.loads((evidence / "plan/reference-plan.json").read_text())[
        "fit_criteria"
    ]
    result = joint_fit(
        model, training, validation, reference_gate_passed=True, criteria=criteria
    )
    bootstrap = bootstrap_training(model, training)
    require(
        signed(bootstrap)
        == json.loads((evidence / "training-bootstrap.json").read_text()),
        "bootstrap_reconstruction_mutation",
    )
    values = np.array([result["coefficients"][k] for k in bootstrap["columns"]])
    result["bootstrap_uncertainty_gate_passed"] = bool(
        bootstrap["valid_replicates"] >= 90
        and np.all(np.array(bootstrap["sigma"]) <= 0.5 * np.abs(values))
        and np.all(np.abs(values) < 1000 - 1e-6)
    )
    require(
        signed(result) == json.loads(path.read_text()),
        "joint_fit_reconstruction_mutation",
    )
    diagnostics = signed(discrepancy_diagnostics(model, training, result))
    require(
        diagnostics
        == json.loads((evidence / "model-discrepancy-diagnostics.json").read_text()),
        "model_discrepancy_reconstruction_mutation",
    )
    provenance = json.loads((evidence / "fitted-candidate-provenance.json").read_text())
    require(
        provenance["identity"]
        == digest({k: v for k, v in provenance.items() if k != "identity"})
        and provenance["dataset_identity"] == fd["identity"]
        and provenance["dataset_sha256"]
        == pin(evidence / "stages/tight-d-joint-fit/normalized-reference.json")
        and provenance["coefficient_checksum"] == digest(result["coefficients"])
        and provenance["fit_record_sha256"] == pin(path)
        and provenance["accepted_parameters"] == [],
        "fitted_provenance_mutation",
    )
    for p in provenance["parameters"]:
        require(
            p["value"] == result["coefficients"][p["name"]] and p["accepted"] is False,
            "coefficient_promotion_or_mutation",
        )


def refine(repo, evidence, work, python):
    repo, evidence, work = Path(repo), Path(evidence), Path(work)
    plan = json.loads((evidence / "plan/tight-d-refinement-plan.json").read_text())
    parent = json.loads((evidence / "plan/reference-plan.json").read_text())
    preserve = evidence / "primary-fitting-decision.json"
    if not preserve.exists():
        shutil.copyfile(evidence / "fitting-decision.json", preserve)
    primary = json.loads(preserve.read_text())
    require(
        primary["reference_gates"]["cardinal_QZ_5Z"] is False,
        "refinement_rationale_changed",
    )
    root = evidence / "stages/tight-d-basis"
    sp = new_stage(
        parent,
        plan["geometries"],
        plan["protocols"],
        plan["protocol_points"],
        plan["budget"],
        pin(evidence / "plan/tight-d-refinement-plan.json"),
        "corrected_cardinal_reference",
    )
    data = run_stage(root, work / "basis", python, sp)
    write(root / "real-residual-subtraction.json", residuals(repo, root))
    primary_data = json.loads((evidence / "normalized-reference.json").read_text())
    write(evidence / "tight-d-basis-effects.json", basis_effects(primary_data, data))
    reference = convergence(repo, evidence, data)
    write(evidence / "tight-d-cardinal-convergence.json", reference)
    stages = [{"stage": "tight-d-basis", "jobs": 18, "identity": data["identity"]}]
    combined = data
    fit = None
    if reference["cardinal_basis_gate_passed"]:
        labels = list(plan["geometries"])
        three = ["center", "angle-angle-0000", "bond-bond-0024"]
        protocol = plan["protocols"]["cc5"]
        numerical = {
            "grid": {**protocol, "grid": [99, 590]},
            "tight": {**protocol, "thresholds": [1e-12, 1e-10]},
        }
        nplan = new_stage(
            parent,
            plan["geometries"],
            numerical,
            {"grid": labels, "tight": three},
            {**parent["budget"], "maximum_gradient_jobs": 12},
            pin(evidence / "plan/tight-d-refinement-plan.json"),
            "corrected_basis_numerical_controls",
        )
        nd = run_stage(
            evidence / "stages/tight-d-numerics", work / "numerics", python, nplan
        )
        combined = signed(
            {
                **{k: v for k, v in data.items() if k != "identity"},
                "samples": {**data["samples"], **nd["samples"]},
            }
        )
        reference = convergence(repo, evidence, combined)
        stages.append(
            {"stage": "tight-d-numerics", "jobs": 12, "identity": nd["identity"]}
        )
        passed = reference_gates_pass(reference)
        if passed:
            d = json.loads((evidence / "joint-target-design.json").read_text())
            training = d["selected_training"]
            validation = d["held_out"]
            allpoints = training + validation
            fplan = new_stage(
                parent,
                {s["id"]: s["coordinates"] for s in allpoints},
                plan["protocols"],
                {
                    "cc5": [s["id"] for s in allpoints],
                    "ccq": [s["id"] for s in validation],
                },
                {**parent["budget"], "maximum_gradient_jobs": 22},
                pin(evidence / "plan/tight-d-refinement-plan.json"),
                "versioned_joint_fit_reference_v2",
            )
            write(evidence / "plan/conditional-fit-data-plan-v2.json", fplan)
            fd = run_stage(
                evidence / "stages/tight-d-joint-fit", work / "fit", python, fplan
            )
            stages.append(
                {"stage": "tight-d-joint-fit", "jobs": 22, "identity": fd["identity"]}
            )
            accuracy = held_out_accuracy(
                fd, data, d, parent["accuracy_criteria"]["basis"]
            )
            write(
                evidence / "tight-d-held-out-reference-accuracy.json",
                accuracy,
            )
            if accuracy["predeclared_basis_gate_passed"]:
                model = load_model(repo)

                def sample(s):
                    return {
                        **fd["samples"]["cc5/" + s["id"]],
                        "split": s["split"],
                        "graph_identity": model["graph_identity"],
                        "atom_order": model["atom_order"],
                        "conformer_identity": s["id"],
                    }

                fit = joint_fit(
                    model,
                    [sample(s) for s in training],
                    [sample(s) for s in validation],
                    reference_gate_passed=True,
                    criteria=parent["fit_criteria"],
                )
                bootstrap = bootstrap_training(model, [sample(s) for s in training])
                write(evidence / "training-bootstrap.json", signed(bootstrap))
                values = np.array(
                    [fit["coefficients"][k] for k in bootstrap["columns"]]
                )
                fit["bootstrap_uncertainty_gate_passed"] = bool(
                    bootstrap["valid_replicates"] >= 90
                    and np.all(np.array(bootstrap["sigma"]) <= 0.5 * np.abs(values))
                    and np.all(np.abs(values) < 1000 - 1e-6)
                )
                write(evidence / "experimental-joint-fit.json", signed(fit))
    write(evidence / "tight-d-combined-normalized-reference.json", combined)
    write(evidence / "refinement-convergence.json", reference)
    write(evidence / "reference-convergence.json", reference)
    final = decision(repo, evidence, reference, combined, fit, stages)
    write(evidence / "fitting-decision.json", final)
    return summarize(repo, evidence) if fit else final


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j60"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path, required=True)
    p.add_argument("--summarize-only", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            summarize(a.repo, a.evidence)
            if a.summarize_only
            else refine(a.repo, a.evidence, a.work, a.python),
            sort_keys=True,
        )
    )
