"""Read real validated synthetic saved records; no mocked scientific success."""

import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from test_workflow import config, start
from test_workflow import prepared as prepared_fixture

from island.analysis import AnalysisOptions, analyze_workflow, export_analysis
from island.exceptions import AnalysisError
from island.workflows import resume_workflow, start_prepared_workflow, storage

prepared = prepared_fixture


def hashes(root):
    return {
        str(p.relative_to(root)): storage.checksum(p.read_bytes())
        for p in root.rglob("*")
        if p.is_file()
    }


def test_saved_analysis_export_ownership(tmp_path, prepared):
    root = tmp_path / "run"
    start(root, prepared)
    before = hashes(root)
    report = analyze_workflow(root)
    payload = report.payload
    assert payload["source"]["status"] == "paused"
    assert [r["step"] for r in payload["frames"]] == [0, 1, 2]
    assert payload["frames"][0]["geometry"]["end_to_end_distance"] is None
    # Independent pair-distance identity, not a second centroid implementation.
    from island.workflows import read_workflow_frames

    for frame, row in zip(read_workflow_frames(root), payload["frames"], strict=True):
        xyz = np.array([frame.coordinates[s] for s in sorted(frame.coordinates)])
        mass = np.array(row["geometry"]["selected_masses"])
        rg2 = sum(
            mass[i] * mass[j] * np.dot(xyz[i] - xyz[j], xyz[i] - xyz[j])
            for i in range(len(xyz))
            for j in range(len(xyz))
        ) / (2 * mass.sum() ** 2)
        assert row["geometry"]["radius_of_gyration"] ** 2 == pytest.approx(
            rg2, abs=1e-12, rel=1e-12
        )
    output = export_analysis(report, tmp_path / "analysis")
    data = json.loads((output / "report.json").read_text())
    assert data == payload
    rows = list(csv.DictReader((output / "frames.csv").open()))
    assert rows[0]["Re_angstrom"] == "" and rows[0]["weighting"] == "mass"
    assert (
        float(rows[0]["Rg_angstrom"])
        == data["frames"][0]["geometry"]["radius_of_gyration"]
    )
    payload["frames"].clear()
    assert report.payload["sample_count"] == 3
    with pytest.raises(AnalysisError):
        export_analysis(report, output)
    with pytest.raises(AnalysisError):
        export_analysis(report, root / "analysis")
    assert hashes(root) == before


def test_unequal_spacing_windows(tmp_path, prepared):
    root = tmp_path / "run"
    cfg = replace(
        config(root), total_steps=5, recording_interval=2, max_frames_per_segment=2
    )
    start_prepared_workflow(
        cfg, *prepared, evidence="synthetic_software_test", segments=3
    )
    report = analyze_workflow(root)
    assert [r["step"] for r in report.payload["frames"]] == [0, 2, 4, 5]
    selected = analyze_workflow(
        root, AnalysisOptions(start_step=1, end_time_ps=0.0005, stride=2)
    )
    assert [r["step"] for r in selected.payload["frames"]] == [2, 5]
    with pytest.raises(AnalysisError, match="No retained"):
        analyze_workflow(root, AnalysisOptions(start_step=100))


@pytest.mark.parametrize("offset", [5e-9, -5e-9])
def test_tolerated_boundary_uses_earlier_frame(tmp_path, prepared, monkeypatch, offset):
    from island.evaluation import OpenMMEvaluationSession

    root = tmp_path / "run"
    start(root, prepared)
    before = analyze_workflow(root).payload["frames"][-1]
    original = OpenMMEvaluationSession.evaluate_fresh
    calls = 0

    def shifted(self, *a, **kw):
        nonlocal calls
        calls += 1
        r = original(self, *a, **kw)
        if calls == 1:
            components = dict(r.energy_components)
            components[next(iter(components))] += offset
            return replace(
                r,
                potential_energy=r.potential_energy + offset,
                energy_components=components,
            )
        return r

    with monkeypatch.context() as patch:
        patch.setattr(OpenMMEvaluationSession, "evaluate_fresh", shifted)
        resume_workflow(root)
    p = analyze_workflow(root).payload
    assert len(p["frames"]) == 5
    assert p["frames"][2] == before


def test_coherent_snapshot_during_publication(tmp_path, prepared, monkeypatch):
    root = tmp_path / "run"
    start(root, prepared)
    manifest_bytes = (root / "manifest.json").read_bytes()
    original = Path.read_bytes
    published = False

    def read(path):
        nonlocal published
        data = original(path)
        if path == root / "manifest.json" and not published:
            published = True
            resume_workflow(root)
        return data

    monkeypatch.setattr(Path, "read_bytes", read)
    p = analyze_workflow(root).payload
    assert published
    assert p["source"]["accepted_step"] == 2
    assert p["source"]["manifest_file_sha256"] == storage.checksum(manifest_bytes)
    assert [r["step"] for r in p["frames"]] == [0, 1, 2]
    assert json.loads(original(root / "manifest.json"))["payload"]["accepted_step"] == 4


def test_corrupt_snapshot_and_failed_export(tmp_path, prepared, monkeypatch):
    root = tmp_path / "run"
    m = start(root, prepared)
    report = analyze_workflow(root)
    before = hashes(root)

    def fail(*a, **kw):
        raise OSError("injected publication failure")

    monkeypatch.setattr(storage, "publish", fail)
    with pytest.raises(AnalysisError, match="injected"):
        export_analysis(report, tmp_path / "failed-output")
    assert not (tmp_path / "failed-output").exists()
    assert hashes(root) == before
    path = root / m["segments"][0]
    path.write_bytes(path.read_bytes() + b"corrupt")
    with pytest.raises(AnalysisError):
        analyze_workflow(root)


def test_failed_history_and_no_scientific_execution(tmp_path, prepared, monkeypatch):
    from island.evaluation import OpenMMSinglePointEvaluator
    from island.workflows import chain

    root = tmp_path / "run"
    m = start(root, prepared)
    m["status"] = "stage_failed"
    chain._manifest(root, m)

    def forbidden(*a, **kw):
        raise AssertionError("Analysis invoked a scientific engine")

    for name in (
        "inspect_chain",
        "minimize_geometry",
        "initialize_velocities",
        "run_dynamics_segment",
        "resume_dynamics",
    ):
        monkeypatch.setattr(chain, name, forbidden)
    monkeypatch.setattr(
        chain.AmberToolsParameterizationEngine, "parameterize", forbidden
    )
    monkeypatch.setattr(OpenMMSinglePointEvaluator, "__init__", forbidden)
    report = analyze_workflow(root)
    assert report.payload["source"]["status"] == "stage_failed"
    assert report.payload["sample_count"] == 3


def test_cli_subprocess(tmp_path, prepared):
    import subprocess
    import sys

    root = tmp_path / "run"
    start(root, prepared)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.analysis",
            str(root),
            "--output",
            str(tmp_path / "export"),
            "--selection",
            "explicit",
            "--site-ids",
            "125",
            "101",
            "--endpoints",
            "101",
            "125",
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "3 retained samples" in completed.stdout
    p = json.loads((tmp_path / "export" / "report.json").read_text())
    assert p["frames"][0]["geometry"]["selected_ids"] == [101, 125]


def test_analysis_rejects_consistently_checksummed_trajectory_mismatch(
    tmp_path, prepared
):
    from test_workflow_consistency import alternative, publish_alternative

    root = tmp_path / "run"
    m = start(root, prepared)
    candidate = alternative(root, m, "seed")
    publish_alternative(root, m, candidate)
    with pytest.raises(AnalysisError, match="initialization|Initialization|velocit"):
        analyze_workflow(root)


def test_missing_bundle_has_no_invented_frames(tmp_path, prepared):
    from island.workflows.chain import _manifest, _new_manifest

    root = tmp_path / "empty"
    root.mkdir()
    m = _new_manifest(config(root), "synthetic_software_test")
    _manifest(root, m)
    with pytest.raises(AnalysisError, match="No accepted retained frames"):
        analyze_workflow(root)
