"""Finish actual J60 numerical work; conditional stages reuse the same campaign.

Call after the primary bounded campaign returns. This never changes its plan,
requests, failed gates or native results. No MD or profile registration occurs.
"""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, pin, require, signed, write
from island.pcff_qm_pilot.bundle import pack
from island.pcff_qm_pilot.hamiltonian import IDENTITY, load_model
from island.pcff_qm_pilot.hamiltonian_analysis import convergence, normalize
from island.pcff_qm_pilot.hamiltonian_campaign import execute
from island.pcff_qm_pilot.hamiltonian_fit import joint_fit
from scripts.independent_phase_4j59 import check as independent_check


def retain(work, evidence):
    with tempfile.TemporaryDirectory(prefix="j60-retain-") as temp:
        temp = Path(temp)
        for name in ("requests", "outputs", "logs", "attempts"):
            if (work / name).exists():
                shutil.copytree(work / name, temp / name)
        if (work / "execution-summary.json").exists():
            shutil.copyfile(
                work / "execution-summary.json", temp / "execution-summary.json"
            )
        write(evidence / "native-file-inventory.json", inventory(temp))
        pack(temp, evidence / "native-outputs.tar.xz")


def stage_plan(parent, geometries, protocols, points, budget, authority):
    return signed(
        {
            **{
                k: v
                for k, v in parent.items()
                if k
                not in (
                    "identity",
                    "geometries",
                    "protocols",
                    "protocol_points",
                    "budget",
                )
            },
            "geometries": geometries,
            "protocols": protocols,
            "protocol_points": points,
            "budget": budget,
            "additive_stage_authority": authority,
        }
    )


def run_stage(evidence, work, python, plan):
    path = evidence / "plan/reference-plan.json"
    if path.exists():
        require(json.loads(path.read_text()) == plan, "stage_resume_mutation")
    else:
        write(path, plan)
    rows = execute(evidence, work, python)
    require(
        all(r["status"] in ("completed", "verified_resume_skip") for r in rows),
        "stage_calculations_failed_preserved",
    )
    data = normalize(evidence, work)
    write(evidence / "normalized-reference.json", data)
    write(
        evidence / "independent-native-checks.json",
        independent_check(work, evidence / "normalized-reference.json"),
    )
    retain(work, evidence)
    return data


def finish(repo, evidence, work, python):
    repo, evidence, work = Path(repo), Path(evidence), Path(work)
    data = normalize(evidence, work)
    write(evidence / "normalized-reference.json", data)
    write(
        evidence / "independent-native-checks.json",
        independent_check(work, evidence / "normalized-reference.json"),
    )
    result = convergence(repo, evidence, data)
    write(evidence / "primary-reference-convergence.json", result)
    primary = json.loads((evidence / "plan/reference-plan.json").read_text())
    gates = {
        r["effect"]: r["predeclared_gate_passed"]
        for r in result["comparisons"]
        if r["predeclared_gate_passed"] is not None
    }
    required = ("cardinal_QZ_5Z", "numerical_grid", "SCF_threshold")
    combined = data
    fit = None
    stage_info = []
    if all(gates[k] for k in required):
        amendment = json.loads(
            (evidence / "plan/grid-derivative-amendment.json").read_text()
        )
        labels = amendment["additional_jobs"]
        stage = evidence / "stages/grid-derivatives"
        gridplan = stage_plan(
            primary,
            {k: primary["geometries"][k] for k in labels},
            {"grid": amendment["protocol"]},
            {"grid": labels},
            {**amendment["budget"], "maximum_gradient_jobs": 6},
            pin(evidence / "plan/grid-derivative-amendment.json"),
        )
        gd = run_stage(
            stage, work.parent / "phase4j60-grid-derivatives", python, gridplan
        )
        combined = signed(
            {
                **{k: v for k, v in data.items() if k != "identity"},
                "samples": {**data["samples"], **gd["samples"]},
            }
        )
        result = convergence(repo, evidence, combined)
        stage_info.append(
            {"stage": "grid-derivatives", "jobs": 6, "identity": gd["identity"]}
        )
        write(evidence / "combined-normalized-reference.json", combined)
        gates = {
            r["effect"]: r["predeclared_gate_passed"]
            for r in result["comparisons"]
            if r["predeclared_gate_passed"] is not None
        }
    if all(gates[k] for k in required):
        conditional = json.loads(
            (evidence / "plan/conditional-fit-data-plan.json").read_text()
        )
        training = conditional["training"]
        validation = conditional["validation"]
        samples = training + validation
        geoms = {s["id"]: s["coordinates"] for s in samples}
        fp = stage_plan(
            primary,
            geoms,
            {
                "cc5": conditional["protocol"],
                "ccq": conditional["validation_basis_control"],
            },
            {"cc5": [s["id"] for s in samples], "ccq": [s["id"] for s in validation]},
            {**conditional["budget"], "maximum_gradient_jobs": 22},
            pin(evidence / "plan/conditional-fit-data-plan.json"),
        )
        fd = run_stage(
            evidence / "stages/joint-fit-data",
            work.parent / "phase4j60-fit-data",
            python,
            fp,
        )
        stage_info.append(
            {"stage": "joint-fit-data", "jobs": 22, "identity": fd["identity"]}
        )
        model = load_model(repo)

        def ref(s):
            row = fd["samples"]["cc5/" + s["id"]]
            return {
                **row,
                "split": s["split"],
                "graph_identity": model["graph_identity"],
                "atom_order": model["atom_order"],
                "conformer_identity": s["id"],
            }

        # Basis accuracy on held-out geometries is assessed before any coefficient
        # inference; reference accuracy checks do not select/tune fitting weights.
        differences = []
        for s in validation:
            a, b = fd["samples"]["ccq/" + s["id"]], fd["samples"]["cc5/" + s["id"]]
            ca, cb = data["samples"]["ccq/center"], data["samples"]["cc5/center"]
            differences.append(
                {
                    "id": s["id"],
                    "relative_energy_difference": (b["energy"] - cb["energy"])
                    - (a["energy"] - ca["energy"]),
                    "gradient_RMS_difference": float(
                        np.sqrt(
                            np.mean(
                                (np.array(a["forces"]) - np.array(b["forces"])) ** 2
                            )
                        )
                    ),
                }
            )
        held_ref = all(
            abs(s["relative_energy_difference"]) <= 0.05
            and s["gradient_RMS_difference"] <= 0.5
            for s in differences
        )
        write(
            evidence / "held-out-reference-accuracy.json",
            {"records": differences, "predeclared_basis_gate_passed": held_ref},
        )
        if held_ref:
            fit = joint_fit(
                model,
                [ref(s) for s in training],
                [ref(s) for s in validation],
                reference_gate_passed=True,
                criteria=primary["fit_criteria"],
            )
            # These are experimental numbers even if held-out errors pass. The
            # required uncertainty/protocol checks remain explicit, never promoted.
            write(evidence / "experimental-joint-fit.json", signed(fit))
    write(evidence / "reference-convergence.json", result)
    from scripts.residuals_phase_4j60 import residuals

    write(evidence / "real-residual-subtraction.json", residuals(repo, evidence))
    retain(work, evidence)
    design = json.loads((evidence / "joint-target-design.json").read_text())
    basis = next(r for r in result["comparisons"] if r["effect"] == "cardinal_QZ_5Z")
    if fit:
        decision = (
            "experimental_fit_completed"
            if fit["held_out_thresholds_passed"]
            else "experimental_fit_rejected_by_untouched_validation"
        )
        next_requirement = "Conformer-group uncertainty/bootstrap and protocol sensitivity at the exact training and untouched conformers; no methyl transfer. If held-out errors fail, inspect recorded residual structure before a new versioned model/data plan."
    else:
        decision = "no_fit"
        failed = {k: v for k, v in gates.items() if not v}
        if not gates["cardinal_QZ_5Z"]:
            worst = max(
                basis["records"], key=lambda r: abs(r["relative_energy_difference"])
            )
            next_requirement = f"At the existing exact revised-center {worst['geometry']} geometry, bind a demonstrably adequate additional orbital basis/augmentation control against the retained cc-pV5Z reference. The QZ/5Z maximum relative-energy discrepancy is {basis['max_relative_energy_difference']:.9g} kcal/mol and maximum gradient RMS discrepancy is {basis['max_gradient_RMS_difference']:.9g} kcal/mol/A. Use the signed mixed-derivative errors to select the smallest follow-up; no author coefficient or pair convention is required to evaluate this research Hamiltonian."
        elif failed:
            next_requirement = (
                "Resolve the exact failing numerical protocol contrast(s): "
                + ",".join(failed)
                + " using retained native gradients and unchanged thresholds before the reserved fitting stage."
            )
        else:
            next_requirement = "The held-out reference-accuracy file identifies the exact distinct-conformer basis discrepancies; obtain converged values for those two geometries before joint inference."
    final = signed(
        {
            "schema": "island.j60-fitting-decision.v1",
            "decision": decision,
            "hamiltonian": IDENTITY,
            "functional_accounting_complete": True,
            "independent_unknown_coefficients": 20,
            "complete_design_rank": design["selected_complete_design"]["rank"],
            "complete_design_condition": design["selected_complete_design"][
                "condition"
            ],
            "reference_gates": gates,
            "basis_gate_numbers": basis,
            "stages_executed": stage_info,
            "fitted_coefficients": list(fit["coefficients"]) if fit else [],
            "fit_completed": fit is not None,
            "fit_validated": False,
            "fit_readiness": "experimental_fit_executed"
            if fit
            else "blocked_reference_accuracy",
            "fit_completion": "completed" if fit else "not_run",
            "fit_validation": "not_established" if fit else "not_run",
            "QM_reference_adequacy": "passed_declared_pilot_criteria"
            if fit
            else "not_established",
            "smallest_next_requirement": next_requirement,
            "scope": "9-atom Si-H disiloxane research extension only; no methyl-PDMS coefficient transfer",
            "parameter_completeness": "incomplete",
            "typing_status": "complete",
            "charge_status": "complete",
            "QM_execution": "passed",
            "QM_data_integrity": "passed",
            "lammps_export_status": "not_run",
            "lammps_runtime_status": "not_run",
            "synthetic_lammps_static_kernels": "passed",
            "openmm_build_status": "not_run",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(evidence / "fitting-decision.json", final)
    return final


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j60"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(finish(a.repo, a.evidence, a.work, a.python), sort_keys=True))
