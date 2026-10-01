"""Public report/input boundary regressions from independent review."""

import json
from copy import deepcopy

import numpy as np
import pytest
from test_workflow import prepared as prepared_fixture
from test_workflow import start

from island.analysis import (
    AnalysisOptions,
    AnalysisReport,
    analyze_workflow,
    export_analysis,
    geometry_metrics,
)
from island.exceptions import AnalysisError

prepared = prepared_fixture


@pytest.fixture
def report(tmp_path, prepared):
    start(tmp_path / "run", prepared)
    return analyze_workflow(tmp_path / "run")


@pytest.mark.parametrize("change", ["count", "rg", "readiness", "json", "missing"])
def test_reconstructed_invalid_report_before_publication(tmp_path, report, change):
    p = report.payload
    if change == "count":
        p["sample_count"] = 999
    elif change == "rg":
        p["frames"][0]["geometry"]["radius_of_gyration"] = -10
    elif change == "readiness":
        p["production_validated"] = True
    raw = (
        "not json"
        if change == "json"
        else "{}"
        if change == "missing"
        else json.dumps(p)
    )
    output = tmp_path / "export"
    with pytest.raises(AnalysisError):
        export_analysis(AnalysisReport(raw), output)
    assert not output.exists()


@pytest.mark.parametrize("value", [True, np.bool_(False)])
def test_mixed_boolean_coordinates(value):
    with pytest.raises(AnalysisError):
        geometry_metrics({1: (value, 0.0, 0.0), 2: (0.0, 0.0, 0.0)}, {1: 1.0, 2: 1.0})


def test_large_steps_and_unrepresentable_time():
    assert AnalysisOptions(start_step=10**400).start_step == 10**400
    with pytest.raises(AnalysisError):
        AnalysisOptions(start_time_ps=10**400)


def test_valid_reconstruction_owned(tmp_path, report):
    p = report.payload
    original = deepcopy(p)
    reconstructed = AnalysisReport(json.dumps(p))
    reconstructed.validate_integrity()
    p["frames"].clear()
    assert reconstructed.payload == original
    output = export_analysis(reconstructed, tmp_path / "export")
    assert json.loads((output / "report.json").read_text()) == original


def test_strict_json_and_structural_fields(tmp_path, report):
    from island.workflows import storage

    p = report.payload
    malformed = [
        "not json",
        "{}",
        "[]",
        "null",
        '{"a":1,"a":2}',
        report._json.replace('"sample_count": 3', '"sample_count": NaN'),
        report._json.replace('"sample_count": 3', '"sample_count": 1e999'),
    ]
    for raw in malformed:
        with pytest.raises(AnalysisError):
            export_analysis(AnalysisReport(raw), tmp_path / "never")
        assert not (tmp_path / "never").exists()
    # Duplicates within a valid report are rejected too.
    raw = json.dumps(p).replace(
        '"sample_count": 3', '"sample_count": 3, "sample_count": 3'
    )
    with pytest.raises(AnalysisError, match="Duplicate"):
        AnalysisReport(raw).validate_integrity()
    assert (
        storage.checksum((tmp_path / "run" / "manifest.json").read_bytes())
        == p["source"]["manifest_file_sha256"]
    )


def test_geometry_and_frame_consistency(tmp_path, report):
    mutations = [
        lambda p: p["frames"].clear(),
        lambda p: p["frames"][1].update(step=0),
        lambda p: p["frames"][1].update(time_ps=0),
        lambda p: p["options"].update(start_step=2),
        lambda p: p["frames"][0].update(source_segment="unknown.json"),
        lambda p: p["source"]["segments"][0].update(checksum="invalid"),
        lambda p: p["frames"][0].update(total_energy_kj_mol=999),
        lambda p: p["frames"][0].update(temperature_kelvin=-1),
        lambda p: p["source"].update(production_validated=True),
        lambda p: p.update(simulation_readiness="ready"),
        lambda p: p["frames"][0]["geometry"].update(selected_ids=[101, 101]),
        lambda p: p["frames"][0]["geometry"].update(selected_masses=[1]),
        lambda p: p["frames"][0]["geometry"].update(normalized_weights=[1, 1, 1, 1]),
        lambda p: p["frames"][0]["geometry"].update(coordinate_unit="nm"),
        lambda p: p["frames"][0]["geometry"].update(center=[0, 0]),
        lambda p: p["frames"][0]["geometry"].update(eigenvalues=[3, 2, 1]),
        lambda p: p["frames"][0]["geometry"].update(kappa_squared=2),
        lambda p: p["frames"][0]["geometry"].update(radius_of_gyration=0.123),
        lambda p: p["frames"][0]["geometry"]["gyration_tensor"][0].__setitem__(1, 999),
        lambda p: p["frames"][0]["geometry"].update(end_to_end_distance=1),
        lambda p: p["frames"][0]["geometry"].update(re_over_rg=1),
        lambda p: p["frames"][0]["geometry"].update(diagnostics=[]),
    ]
    for change in mutations:
        p = report.payload
        change(p)
        with pytest.raises(AnalysisError):
            export_analysis(AnalysisReport(json.dumps(p)), tmp_path / "never")
        assert not (tmp_path / "never").exists()


def test_valid_degenerate_reports_and_missing_source(tmp_path, prepared):
    import shutil

    from island.analysis import GeometryOptions

    start(tmp_path / "run", prepared)
    for endpoints in (None, (101, 101), (101, 109)):
        report = analyze_workflow(
            tmp_path / "run",
            AnalysisOptions(
                geometry=GeometryOptions("explicit", (101,), endpoints=endpoints)
            ),
        )
        g = report.payload["frames"][0]["geometry"]
        assert (
            g["radius_of_gyration"] == 0
            and g["kappa_squared"] is None
            and g["re_over_rg"] is None
        )
        report.validate_integrity()
        export_analysis(AnalysisReport(report._json), tmp_path / str(endpoints))
    shutil.rmtree(tmp_path / "run")
    # Source lifetime is independent of report export lifetime.
    export_analysis(report, tmp_path / "without-source")


def test_existing_output_path_and_cleanup_failures(tmp_path, report, monkeypatch):
    import island.analysis.export as export_module
    from island.workflows import storage

    output = tmp_path / "existing"
    output.mkdir()
    (output / "sentinel").write_bytes(b"unchanged")
    with pytest.raises(AnalysisError):
        export_analysis(report, output)
    assert (output / "sentinel").read_bytes() == b"unchanged"
    with pytest.raises(AnalysisError):
        export_analysis(report, object())

    def fail_publish(*a, **kw):
        raise OSError("original publication failure")

    def fail_cleanup(*a, **kw):
        raise OSError("cleanup failure")

    monkeypatch.setattr(storage, "publish", fail_publish)
    monkeypatch.setattr(export_module.shutil, "rmtree", fail_cleanup)
    with pytest.raises(AnalysisError, match="original publication failure") as caught:
        export_analysis(report, tmp_path / "failed")
    assert isinstance(caught.value.__cause__, OSError)
    assert "cleanup failure" in caught.value.__cause__.__notes__[0]
    assert not (tmp_path / "failed" / "COMPLETE.json").exists()


def test_numeric_arrays_still_supported():
    for value in (
        (1.0, 0.0, 0.0),
        np.array([1, 0, 0]),
        np.array([True, False, False], dtype=float),
    ):
        result = geometry_metrics({1: value, 2: [0.0, 0.0, 0.0]}, {1: 1, 2: 1})
        assert result.radius_of_gyration == 0.5
