"""Offline J66 physical-array, canonical-problem and complete-MM replay."""

import argparse
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse

from island.pcff_completion import inventory, require
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_constrained import canonical, physical_problem, verify
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_fit import prediction_uncertainty
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j63 import sample as native_sample
from scripts.analyze_phase_4j65 import assemble, independent_value
from scripts.analyze_phase_4j66 import full_evaluator_check, reconstruct


def compare_system(file, expected):
    with np.load(file, allow_pickle=False) as saved:
        for name in ("A", "P"):
            actual = sparse.csc_matrix(
                (
                    saved[name + "_data"],
                    saved[name + "_indices"],
                    saved[name + "_indptr"],
                ),
                shape=tuple(saved[name + "_shape"]),
            )
            require(
                actual.shape == expected[name].shape
                and np.max(abs((actual - expected[name]).toarray()), initial=0) < 1e-11,
                "J66_canonical_matrix_mutation",
            )
        for name in ("b", "q", "scale"):
            require(
                np.allclose(saved[name], expected[name], atol=1e-11, rtol=0),
                "J66_scaling_rhs_mutation",
            )


def verify_record(record, problem, lo, hi):
    result = record["solution"]
    old = record["original_unit_verification"]
    actual = verify(problem, result["values"], lo, hi, result["rho"])
    for key in (
        "energy_RMSE",
        "force_RMS",
        "maximum_absolute_energy_error",
        "maximum_local_force_error",
        "weighted_objective",
    ):
        require(abs(actual[key] - old[key]) < 1e-9, "J66_original_unit_replay")
    require(
        actual["original_limits_satisfied"] == old["original_limits_satisfied"]
        and actual["active_bound_indices"] == old["active_bound_indices"],
        "J66_physical_gate_mutation",
    )
    return actual


def replay(repo, root, native=True):
    if native:
        parent = repo / "docs/evidence/phase_4j65"
        with tempfile.TemporaryDirectory(prefix="j66-native-intake-") as temp:
            work = Path(temp) / "native"
            unpack(parent / "native-outputs.tar.xz", work)
            require(
                inventory(work) == read(parent / "native-file-inventory.json"),
                "J66_native_inventory_mutation",
            )
            raw = read(parent / "normalized-native-controls-and-observations.json")[
                "samples"
            ]
            for key, expected in raw.items():
                actual = native_sample(work, key)
                require(actual == expected, "J66_native_sample_mutation")
                independent_value(work, key, actual)
            reconstructed, missing = assemble(
                repo, work, read(parent / "plan/scientific-plan.json")
            )
            require(
                not missing
                and reconstructed
                == read(parent / "fitting-dataset.json")["observations"],
                "J66_native_training_reconstruction",
            )
    with tempfile.TemporaryDirectory(prefix="j66-offline-relocated-") as tmp:
        target = Path(tmp)
        (target / "plan").mkdir()
        (target / "plan/scientific-plan.json").write_bytes(
            (root / "plan/scientific-plan.json").read_bytes()
        )
        plan, model, observations, data, old = reconstruct(repo, target)
        with (
            np.load(root / "physical-arrays.npz", allow_pickle=False) as a,
            np.load(target / "physical-arrays.npz", allow_pickle=False) as b,
        ):
            require(
                a.files == b.files and all(np.array_equal(a[k], b[k]) for k in a.files),
                "J66_physical_array_or_unit_mutation",
            )
        lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
        hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
        p = physical_problem(data)
        for name, options in (
            ("nominal-rho", {}),
            ("nominal-rho-SCS", {}),
            ("unbounded-rho", {"bounded": False}),
            ("force-only-rho", {"force_only": True}),
            ("constrained-epigraph", {"kind": "least_squares"}),
            ("interior-rho", {"interior": 1e-6}),
        ):
            compare_system(
                root / (name + "-canonical.npz"), canonical(p, lo, hi, **options)
            )
            verify_record(read(root / (name + ".json")), p, lo, hi)
        trial = read(root / "constrained-epigraph.json")
        direct = full_evaluator_check(
            model, observations, data, trial["solution"]["values"]
        )
        require(
            direct["maximum_local_force_error"] <= 2
            and direct["force_RMS"] <= 0.5
            and direct["energy_RMSE"] <= 0.1,
            "J66_complete_evaluator_gate",
        )
        summary = read(root / "block-uncertainty.json")
        if (root / "bootstrap-solves.tar.xz").exists():
            directory = target / "bootstrap"
            unpack(root / "bootstrap-solves.tar.xz", directory)
            require(
                inventory(directory) == read(root / "bootstrap-file-inventory.json"),
                "J66_bootstrap_archive_mutation",
            )
        else:
            directory = root / "bootstrap"
        draws = []
        for original in old["attempts"]:
            number = original["replicate"]
            name = f"replicate-{number:02d}"
            row = read(directory / (name + ".json"))
            require(
                row["sampled_blocks"] == original["sampled_blocks"]
                and row["replicate"] == number,
                "J66_bootstrap_draw_mutation",
            )
            indices = [
                i
                for b in original["sampled_blocks"]
                for i, r in enumerate(observations)
                if r["block"] == b
            ]
            require(indices == row["indices"], "J66_bootstrap_mapping_mutation")
            sample = physical_problem(data, indices)
            compare_system(
                directory / (name + "-rho-canonical.npz"), canonical(sample, lo, hi)
            )
            feasible = verify_record(row["feasibility"], sample, lo, hi)
            design = metrics(sample["w"])
            failures = []
            if design["rank"] != 59:
                failures.append("rank")
            if design["scaled_condition"] > 10000:
                failures.append("conditioning")
            if "fit" in row:
                require(
                    feasible["original_limits_satisfied"], "J66_fit_activation_mutation"
                )
                compare_system(
                    directory / (name + "-fit-canonical.npz"),
                    canonical(sample, lo, hi, kind="least_squares"),
                )
                fit_check = verify_record(row["fit"], sample, lo, hi)
                if row["fit"]["solution"]["status"] != "Solved":
                    failures.append("solver_status")
                if not fit_check["original_limits_satisfied"]:
                    failures.append("performance_or_bounds")
                if fit_check["active_bound_indices"]:
                    failures.append("active_bounds")
                require(
                    row["values"] == row["fit"]["solution"]["values"],
                    "J66_selected_vector_mutation",
                )
            else:
                require(
                    not feasible["original_limits_satisfied"],
                    "J66_feasible_draw_dropped",
                )
                failures.append("no_verified_feasible_witness")
            require(
                failures == row["failure_categories"]
                and row["eligible"] == (not failures),
                "J66_draw_qualification_mutation",
            )
            draws.append(row)
        require(
            len(draws) == 64
            and sum(r["eligible"] for r in draws) == summary["eligible_replicates"]
            and summary["replacement_resamples"] == 0,
            "J66_bootstrap_failure_retention",
        )
        survivors = [r["values"] for r in draws if r["eligible"]]
        if len(survivors) >= 2:
            spread = prediction_uncertainty(
                data,
                survivors,
                trial["solution"]["values"],
                plan["configuration_order"],
            )
            require(
                all(
                    abs(spread[k] - summary["prediction_spread"][k]) < 1e-9
                    for k in ("energy_RMS95", "force_RMS95")
                ),
                "J66_prediction_uncertainty_mutation",
            )
        require(
            summary["eligible_fraction"] == len(survivors) / 64
            and not summary["qualified"],
            "J66_uncertainty_promoted",
        )
        decision = read(root / "decision.json")
        require(
            decision["new_QM"] == decision["validation_QM"] == decision["MD_steps"] == 0
            and not decision["automatic_assignment"]
            and not decision["production_validated"],
            "J66_reserved_validation_or_scope_violation",
        )
        return {
            "native_reconstructed": native,
            "configurations": 22,
            "atoms": 27,
            "physical_columns": 58,
            "bootstrap_draws": 64,
            "eligible": summary["eligible_replicates"],
            "offline": True,
            "new_QM": 0,
            "MD_steps": 0,
            "training_feasible": True,
            "validation_untouched": True,
        }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j66"))
    p.add_argument("--skip-native", action="store_true")
    a = p.parse_args()
    print(replay(a.repo.resolve(), a.evidence.resolve(), not a.skip_native))
