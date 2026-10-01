"""Opt-in read-only existing PE bundle comparisons, never new MD."""

import os

import pytest


@pytest.mark.archived_amber_reference
def test_retained_pe_analysis(tmp_path):
    root = os.environ.get("ISLAND_ANALYSIS_REFERENCE_ROOT")
    if not root:
        pytest.skip(
            "Set ISLAND_ANALYSIS_REFERENCE_ROOT to retained Phase 4E7.1 bundles"
        )
    from pathlib import Path

    from scripts.validate_trajectory_analysis import validate

    for name in ("live-pe-relocated", "archive-pe-relocated"):
        result = validate(Path(root) / name, tmp_path / name)
        assert result["samples"] == 7
        assert result["steps"] == [0, 10, 20, 30, 40, 50, 60]
