"""Offline source/physical-array replay; no converter execution or solver import."""

import argparse
from pathlib import Path

from island.pcff_completion import inventory, pin, require
from island.pcff_qm_pilot.hmds_precision import read
from scripts.audit_phase_4j69 import ROOT, converter_text, numerical
from scripts.replay_phase_4j68 import close
from scripts.validate_phase_4j68_receipt import decision_semantics, semantics


def replay(repo, root, deep=False):
    parent = repo / "docs/evidence/phase_4j68"
    semantics(read(parent / "receipt.json"), read(parent / "manifest.json"))
    decision_semantics(read(parent / "decision.json"))
    if deep:
        from scripts.replay_phase_4j68 import replay as parent_replay

        require(
            parent_replay(repo, parent, True)["validation_untouched"],
            "J69_parent_replay",
        )
    actual = numerical(repo, root)
    for name, result in actual.items():
        close(result, read(root / name), "J69_offline." + name)
    converter = read(root / "converter/attempt.json")
    require(
        converter["input_files"] == inventory(root / "converter/input")
        and converter["executable_sha256"] == pin(root / "converter/msi2lmp"),
        "J69_converter_execution_inventory",
    )
    from island.forcefields.pcff.interop import load_msi2lmp_car_mdf
    from island.pcff_qm_pilot.hmds_hamiltonian import load_model
    from island.workflows import storage

    system = load_msi2lmp_car_mdf(root / "converter/input")
    payload = storage.decode(storage.read_json(root / "converter/input/system.json"))
    model = load_model(repo)
    original = read(repo / "docs/evidence/phase_4j67/fitting-dataset.json")[
        "observations"
    ][0]
    require(
        payload["types"] == {int(k): v for k, v in model["typing"].items()}
        and payload["charges"] == {int(k): v for k, v in model["charges"].items()}
        and [system.coordinates.get(i).tolist() for i in range(1, 28)]
        == original["coordinates"],
        "J69_converter_physical_inputs",
    )
    require(
        converter["returncode"] == 17
        and converter["ignore_used"] is False
        and converter["new_QM"] == converter["MD_steps"] == 0,
        "J69_converter_completeness",
    )
    require(
        "Unable to find endbonddihedral data for hc c3 sio c3"
        in converter_text(root, "stdout.log"),
        "J69_converter_missing_row",
    )
    return {
        "offline": True,
        "observations": 26,
        "protected": 22,
        "physical_columns": 58,
        "energy_only_intercepts": 1,
        "source_available_subset_closed": True,
        "source_complete": False,
        "new_QM": 0,
        "MD_steps": 0,
        "validation_untouched": True,
        "hypotheses_executed": 0,
        "accepted_coefficients": 0,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=ROOT)
    p.add_argument("--deep-parent", action="store_true")
    a = p.parse_args()
    print(replay(a.repo.resolve(), a.evidence.resolve(), a.deep_parent))
