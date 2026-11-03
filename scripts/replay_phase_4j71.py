"""Offline J71 replay; original PDFs, converter and solver packages are unnecessary."""

import argparse
from pathlib import Path

from island.pcff_completion import require
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.literature_intake import semantic_declarations
from scripts.extract_phase_4j71 import ROOT, human_table, reports


def replay(repo, root, deep_parent=False):
    expected = reports(repo, root)
    for name, value in expected.items():
        require(read(root / name) == value, "J71_semantic_reconstruction", path=name)
    require(
        (root / "source-evidence-table.md").read_text() == human_table(expected),
        "J71_human_table_mutation",
    )
    semantic_declarations(read(root / "decision.json"))
    if deep_parent:
        from scripts.validate_phase_4j70_receipt import validate

        validate(repo, repo / "docs/evidence/phase_4j70")
    result = expected["decision.json"]
    return {
        "offline": True,
        "extracted_candidate_rows": result["extracted_candidate_rows"],
        "nonblank_numeric_cells": result["extracted_nonblank_numeric_cells"],
        "J70_instances": 378,
        "resolved_instances": 0,
        "decision": result["decision"],
        "parameter_completeness": "incomplete",
        "converter_executed": False,
        "static_complete_closure": False,
        "new_QM": 0,
        "QM_fits": 0,
        "MD_steps": 0,
        "accepted_coefficients": 0,
        "validation_untouched": True,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--deep-parent", action="store_true")
    args = p.parse_args()
    print(replay(args.repo.resolve(), args.root.resolve(), args.deep_parent))
