"""Retain raw/normalized evidence and reconstruct it in a fresh relocated tree."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

from island.pcff_completion import inventory, require, signed, write
from island.pcff_qm_pilot.bundle import pack, unpack
from island.pcff_qm_pilot.hmds_precision import decimal_energy, read
from scripts.analyze_phase_4j63 import independent_native, sample
from scripts.analyze_phase_4j64 import analyze

REPORTS = (
    "dense-grid-decision.json",
    "QZVP-center-decision.json",
    "torsion-step-control.json",
    "rigid-rotation-control.json",
    "baseline-calibration.json",
    "updated-joint-design.json",
)


def samples(work):
    ledger = read(work / "ledger.json")["calls"]
    return {
        key: sample(work, key)
        for key, r in ledger.items()
        if r["status"] == "completed"
    }


def independently_check(work, actual):
    checks = {}
    for key, value in actual.items():
        checks[key] = {
            **independent_native(work, key, value),
            "native_total_energy_hartree_decimal": str(
                decimal_energy(work / "outputs" / key)[0]
            ),
        }
    return checks


def retain(work, root):
    with tempfile.TemporaryDirectory(prefix="j64-collect-") as tmp:
        directory = Path(tmp)
        for name in ("requests", "outputs", "processes", "resources"):
            if (work / name).exists():
                shutil.copytree(work / name, directory / name)
        for name in (
            "scientific-plan.json",
            "jobs.json",
            "ledger.json",
            "attempts.json",
        ):
            shutil.copyfile(work / name, directory / name)
        write(root / "native-file-inventory.json", inventory(directory))
        pack(directory, root / "native-outputs.tar.xz")
    actual = samples(work)
    write(
        root / "normalized-reference.json",
        signed(
            {
                "samples": actual,
                "origin": "executed real Psi4; native evidence separate; protocols never pooled into a fitting dataset",
                "protocols_kept_separate": True,
            }
        ),
    )
    write(
        root / "independent-native-checks.json",
        signed(
            {
                "samples": independently_check(work, actual),
                "reconstruction": "Independent ASCII gradient parsing, decimal native total energy, independent unit/sign conversion",
            }
        ),
    )
    attempted = read(work / "attempts.json")["calls"]
    write(
        root / "compute-accounting.json",
        signed(
            {
                "counted_process_attempts": len(attempted),
                "native_gradient_starts": len(read(work / "ledger.json")["calls"]),
                "successful_gradients": len(actual),
                "failed_process_attempts": {
                    k: r for k, r in attempted.items() if r["status"] != "completed"
                },
                "maximum": 24,
                "new_SiH_QM": 0,
                "joint_design_QM": 0,
                "validation_QM": 0,
                "full_or_numerical_Hessian": False,
                "native_resources": {
                    k: read(work / "resources" / (k + "-usage.json")) for k in actual
                },
                "gradient_elapsed_seconds": {
                    k: read(work / "outputs" / k / "result.json")["elapsed_seconds"]
                    for k in actual
                },
            }
        ),
    )


def replay(repo, root):
    with tempfile.TemporaryDirectory(prefix="j64-relocated-") as tmp:
        work, output = Path(tmp) / "native", Path(tmp) / "reconstructed"
        unpack(root / "native-outputs.tar.xz", work)
        require(
            inventory(work) == read(root / "native-file-inventory.json"),
            "J64_native_inventory_mutation",
        )
        actual = samples(work)
        require(
            actual == read(root / "normalized-reference.json")["samples"],
            "J64_native_reconstruction_changed",
        )
        require(
            independently_check(work, actual)
            == read(root / "independent-native-checks.json")["samples"],
            "J64_independent_reconstruction_changed",
        )
        shutil.copytree(root / "plan", output / "plan")
        for process in (work / "processes").rglob("result.json"):
            value = read(process)
            for name, sha in value["files"].items():
                from island.pcff_completion import pin

                require(
                    pin(process.parent / name) == sha, "J64_process_evidence_mutation"
                )
        attempts = read(work / "attempts.json")["calls"]
        require(len(attempts) <= 24, "J64_attempt_budget_mutation")
        for key in actual:
            require(key in read(work / "jobs.json")["jobs"], "J64_unplanned_native_job")
            require(
                not any(label in key for label in ("joint", "validation", "SiH")),
                "J64_forbidden_campaign",
            )
        status = analyze(repo, output, work)
        from scripts.independent_phase_4j64 import check

        independent = check(repo, output, work)
        require(
            independent == read(root / "independent-path-and-baseline-checks.json"),
            "J64_independent_path_or_scalar_replay_changed",
        )
        for name in REPORTS:
            require(
                (root / name).exists() == (output / name).exists(),
                "J64_report_coverage_changed",
                file=name,
            )
            if (root / name).exists():
                require(
                    read(root / name) == read(output / name),
                    "J64_numerical_replay_changed",
                    file=name,
                )
        return {
            "relocated": True,
            "QM_invoked": False,
            "native_gradients": len(actual),
            "counted_attempts": len(attempts),
            "objective_status": status,
            "fresh_work_directory": True,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["retain", "replay"])
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j64")
    )
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.operation == "retain":
        retain(args.work.resolve(), args.evidence.resolve())
    else:
        print(
            json.dumps(
                replay(args.repo.resolve(), args.evidence.resolve()), sort_keys=True
            )
        )
