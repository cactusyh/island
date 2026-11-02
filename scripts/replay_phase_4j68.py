"""Offline physical-array, complete-MM and bound-free certificate replay; no solves."""

import argparse
import tempfile
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, pin, require
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_augmentation import identity, solver_state_diagnostics
from island.pcff_qm_pilot.hmds_constrained import dual_lower_bound, physical_problem
from island.pcff_qm_pilot.hmds_local_baseline import (
    derivative_check,
    extend,
    fixed_families,
    span_sensitivity,
)
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.hmds_span import (
    certificate,
    check,
    classification,
    projection,
    system,
)
from scripts.analyze_phase_4j68 import (
    FROZEN,
    HYPOTHESES,
    complete_check,
    prediction_report,
    reconstruct,
)
from scripts.replay_phase_4j66 import compare_system
from scripts.validate_phase_4j67_receipt import semantics as j67_semantics

BASE = "873e1cc7d5d20f0d6f7c7e5345c3841e0892011c"
BRANCH = "codex/phase-4j68-hmds-model-span-baseline-diagnosis"


def parent_semantics(receipt, manifest, decision):
    j67_semantics(receipt, manifest, training_feasible=False)
    expected = {
        "MD_steps": 0,
        "accepted_candidates": [],
        "accepted_coefficients": 0,
        "author_confirmations": "28 user-attested",
        "automatic_assignment": False,
        "candidate_qualified": False,
        "chemical_transfer": "not_established",
        "constrained_trial_executed": False,
        "explicit_trial_vector_complete": False,
        "new_QM": 4,
        "original_author_artifact_status": "not_yet_supplied",
        "parameter_completeness": "incomplete",
        "production_validated": False,
        "scientific_validation": "not_established",
        "simulation_readiness": "not_established",
        "training_feasible": False,
        "uncertainty_qualified": False,
        "validation_QM": 0,
        "validation_activated": False,
    }
    identity(decision)
    require(set(decision) == set(expected) | {"identity"}, "J68_parent_decision_schema")
    for k, v in expected.items():
        require(
            type(decision[k]) is type(v) and decision[k] == v,
            "J68_parent_scientific_semantics",
            key=k,
        )
    nested = receipt["numerical_reconstruction"]
    required = {
        "native_attempts": 4,
        "new_gradients": 4,
        "training_configurations": 26,
        "protected_configurations": 22,
        "bootstrap_draws": 64,
        "eligible": 1,
        "training_feasible": False,
        "candidate_qualified": False,
        "validation_untouched": True,
        "offline": True,
        "MD_steps": 0,
    }
    require(
        nested == required
        and all(type(nested[k]) is type(v) for k, v in required.items()),
        "J68_parent_nested_semantics",
    )
    return True


def close(actual, previous, label):
    if isinstance(actual, dict):
        require(set(actual) == set(previous), label)
        for k in actual:
            close(actual[k], previous[k], label + "." + str(k))
    elif isinstance(actual, list):
        require(len(actual) == len(previous), label)
        for a, b in zip(actual, previous, strict=True):
            close(a, b, label)
    elif isinstance(actual, (float, np.floating)):
        require(
            np.isfinite(actual)
            and abs(actual - previous) <= 1e-8 * max(1, abs(actual)),
            label,
        )
    else:
        require(type(actual) is type(previous) and actual == previous, label)


def replay(repo, root, deep=False):
    parent = repo / "docs/evidence/phase_4j67"
    parent_semantics(
        *(read(parent / n) for n in ("receipt.json", "manifest.json", "decision.json"))
    )
    if deep:
        from scripts.replay_phase_4j67 import replay as parent_replay

        require(
            parent_replay(repo, parent)["validation_untouched"], "J68_parent_replay"
        )
    require(
        pin(root / "plan/scientific-plan.json") == FROZEN
        and pin(root / "plan/baseline-hypotheses.json") == HYPOTHESES,
        "J68_frozen_contract",
    )
    plan, model, rows, data, lo, hi = reconstruct(repo, root)
    p, protected = physical_problem(data), physical_problem(data, np.arange(22))
    expected = read(root / "model-span-projection.json")
    actual = projection(p)
    actual["checks"] = check(p, protected, actual["values"], lo, hi)
    actual["predictions"] = prediction_report(
        data, actual["values"], [r["id"] for r in rows]
    )
    actual["complete_evaluator_check"] = complete_check(
        model, rows, data, np.array(actual["values"])
    )
    close(
        actual,
        {k: v for k, v in expected.items() if k != "identity"},
        "J68_projection_replay",
    )
    attempts = read(root / "solver-attempts.json")["attempts"]
    require(
        len(attempts) <= plan["solver_call_cap"]
        and len({r["name"] for r in attempts}) == len(attempts),
        "J68_solver_budget",
    )
    outcomes = {}
    comparison = read(root / "executed-comparison.json")
    violations = read(root / "original-bound-violations.json")
    specs = read(parent / "plan/scientific-plan.json")["parameter_specifications"]
    original_values = None
    with tempfile.TemporaryDirectory(prefix="J68-offline-") as temp:
        directory = root / "solves"
        if not directory.exists():
            directory = Path(temp) / "solves"
            unpack(root / "solver-states.tar.xz", directory)
            require(
                inventory(directory) == read(root / "solver-file-inventory.json"),
                "J68_solver_archive",
            )
        for a in attempts:
            name = a["name"]
            file = directory / (name + ".json")
            if not file.exists():
                require(a["status"] != "completed", "J68_missing_completed_state")
                continue
            record = read(file)
            identity(record)
            raw = read(directory / (name + "-raw.json"))
            identity(raw)
            require(
                record["solution"] == {k: v for k, v in raw.items() if k != "identity"},
                "J68_raw_state_overwritten",
            )
            variant = next(
                (
                    v
                    for v in ("bridge_angle", "SiO_stretch")
                    if name.startswith(v + "-")
                ),
                None,
            )
            d = extend(model, rows, data, variant) if variant else data
            pp, pr = physical_problem(d), physical_problem(d, np.arange(22))
            s = system(pp, lo, hi, pr, a["bounded"])
            compare_system(directory / (name + "-canonical.npz"), s)
            checks = check(
                pp, pr, record["solution"]["values"], lo, hi, record["solution"]["rho"]
            )
            close(checks, record["checks"], "J68_original_unit_replay")
            state = solver_state_diagnostics(s, record["solution"])
            close(state, record["raw_state_diagnostics"], "J68_raw_closure_replay")
            if a["bounded"]:
                worst = max(
                    (pp, pr),
                    key=lambda prob: check(
                        prob, None, record["solution"]["values"], lo, hi
                    )["required_rho"],
                )
                cert = dual_lower_bound(s, record["solution"], worst, lo, hi)
                cert["certified"] = True
            else:
                cert = certificate(s, record["solution"], pp, pr, lo, hi)
            close(cert, record["certificate"], "J68_checked_certificate_replay")
            require(
                classification(checks, cert) == record["classification"],
                "J68_classification",
            )
            close(
                prediction_report(
                    d, record["solution"]["values"], [r["id"] for r in rows]
                ),
                record["predictions"],
                "J68_prediction_replay",
            )
            if a["backend"] == "clarabel":
                close(
                    complete_check(
                        model, rows, d, np.array(record["solution"]["values"]), variant
                    ),
                    record["complete_evaluator_check"],
                    "J68_complete_MM_replay",
                )
            if a["status"] == "failed":
                recovery = read(directory / (name + "-postprocessing-recovery.json"))
                require(
                    not recovery["solver_reexecuted"]
                    and recovery["original_failure_log_sha256"]
                    == pin(directory / (name + ".log")),
                    "J68_failed_attempt_lost",
                )
            summary = comparison["models"][name]
            for k, v in {
                "rho": record["solution"]["rho"],
                "lower": cert["lower_bound"],
                "upper": cert["primal_upper_bound"],
                "status": record["solution"]["status"],
                "classification": record["classification"],
                "parameter_values": record["solution"]["values"],
                "parameter_names": d["names"] + ["energy_offset"],
                "performance_limits_satisfied": checks["performance_limits_satisfied"],
                "original_bounds_satisfied": checks["original_bounds_satisfied"],
            }.items():
                close(v, summary[k], "J68_comparison_" + k)
            for subset, c in checks["subsets"].items():
                for metric in summary["subsets"][subset]:
                    close(
                        c[metric],
                        summary["subsets"][subset][metric],
                        "J68_comparison_metrics",
                    )
            actual_violations = [
                {
                    **v,
                    "name": specs[v["index"]]["name"],
                    "units": specs[v["index"]]["units"],
                }
                for v in checks["outside_original_bounds"]
            ]
            close(
                actual_violations,
                violations[name],
                "J68_original_bound_violation_report",
            )
            if name == "bound-free-Clarabel":
                original_values = np.array(record["solution"]["values"])
            outcomes[name] = record["classification"]
    localization = read(root / "fixed-family-localization.json")
    for row, old in zip(rows, localization["rows"], strict=True):
        fixed = fixed_families(model, row["coordinates"])
        close(fixed["families"], old["families"], "J68_fixed_family_replay")
        i = rows.index(row)
        adjustable = {}
        for j, name in enumerate(data["names"]):
            family = name.split(":")[0] if not name.startswith("delta_K2") else name
            adjustable.setdefault(family, np.zeros((27, 3)))
            adjustable[family] += data["force"][i, :, :, j] * original_values[j]
        close(
            {k: v.tolist() for k, v in adjustable.items()},
            old["adjustable_family_forces"],
            "J68_adjustable_family_replay",
        )
    sensitivity = read(root / "baseline-pre-fit-sensitivities.json")
    for variant in ("bridge_angle", "SiO_stretch"):
        amended = physical_problem(extend(model, rows, data, variant))
        close(
            span_sensitivity(p, amended),
            sensitivity[variant]["span"],
            "J68_independent_added_span",
        )
        for row, old in zip(
            rows, sensitivity[variant]["derivative_checks"], strict=True
        ):
            close(
                derivative_check(model, row["coordinates"], variant),
                {k: v for k, v in old.items() if k != "id"},
                "J68_derivative_replay",
            )
    from scripts.summarize_phase_4j68 import remaining_localization

    amended = comparison["models"]["SiO_stretch-clarabel"]
    remaining = remaining_localization(
        model,
        rows,
        data,
        amended["parameter_values"],
        amended,
        localization["rows"],
        expected,
    )
    close(
        remaining,
        {
            k: v
            for k, v in read(root / "remaining-residual-localization.json").items()
            if k != "identity"
        },
        "J68_remaining_atom_localization",
    )
    mapping = read(root / "dataset-and-producers.json")
    require(
        mapping["graph_identity"] == model["graph_identity"]
        and mapping["atom_order"] == list(range(1, 28))
        and mapping["configuration_order"] == [r["id"] for r in rows]
        and mapping["parameter_specifications"] == specs
        and mapping["energy_only_intercepts"] == 1,
        "J68_dataset_semantics",
    )
    require(
        mapping["units"]
        == {
            "energy": "kcal/mol",
            "force": "kcal/mol/angstrom",
            "native_gradient": "hartree/bohr",
            "coordinates": "angstrom",
        }
        and mapping["charge"] == 0
        and mapping["multiplicity"] == 1
        and mapping["new_execution_records_for_historical_QM"] == 0,
        "J68_dataset_units_scope",
    )
    require(
        mapping["observations"]
        == [
            {
                "id": r["id"],
                "producer": r["producer"],
                "native_energy_hartree_decimal": r["native_energy_hartree_decimal"],
                "native_center_energy_hartree_decimal": r[
                    "native_center_energy_hartree_decimal"
                ],
                "parent_record_index": i,
            }
            for i, r in enumerate(rows)
        ],
        "J68_original_producer_report",
    )
    decision = read(root / "decision.json")
    require(
        decision["new_QM"]
        == decision["MD_steps"]
        == decision["validation_QM"]
        == decision["accepted_coefficients"]
        == 0
        and not decision["validation_activated"],
        "J68_no_native_or_validation",
    )
    return {
        "observations": 26,
        "protected": 22,
        "physical_columns": 58,
        "intercept_columns": 1,
        "research_solver_calls": len(attempts),
        "failed_processes_retained": sum(a["status"] != "completed" for a in attempts),
        "classifications": outcomes,
        "new_QM": 0,
        "MD_steps": 0,
        "validation_untouched": True,
        "accepted_coefficients": 0,
        "offline": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j68")
    )
    parser.add_argument("--deep-parent", action="store_true")
    a = parser.parse_args()
    print(replay(a.repo.resolve(), a.evidence.resolve(), a.deep_parent))
