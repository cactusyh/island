"""Independent geometrical identities, numerical limits and dependency isolation."""

import json
import subprocess
import sys

import numpy as np
import pytest

from island.analysis import GeometryOptions, geometry_metrics
from island.exceptions import AnalysisError


def test_two_unequal_masses_and_uniform():
    coordinates = {90: (4, 0, 0), 7: (0, 0, 0)}
    masses = {7: 1.0, 90: 3.0}
    result = geometry_metrics(
        coordinates, masses, options=GeometryOptions(endpoints=(7, 90))
    )
    assert result.selected_ids == (7, 90)
    assert result.radius_of_gyration**2 == pytest.approx(1 * 3 / (1 + 3) ** 2 * 4**2)
    assert result.center == (3, 0, 0)
    assert result.kappa_squared == 1
    assert result.end_to_end_distance == 4
    assert result.re_over_rg == pytest.approx(4 / np.sqrt(3))
    uniform = geometry_metrics(
        coordinates, masses, options=GeometryOptions(weighting="uniform")
    )
    assert uniform.radius_of_gyration == 2
    assert result == geometry_metrics(
        dict(reversed(list(coordinates.items()))),
        masses,
        options=GeometryOptions(endpoints=(7, 90)),
    )
    coordinates[90] = (99, 0, 0)
    data = result.to_dict()
    data["center"] = [99] * 3
    assert result.center == (3, 0, 0)


def test_rotation_translation_and_isotropic():
    xyz = np.vstack([np.eye(3), -np.eye(3)])
    ids = (30, 2, 71, 19, 6, 100)
    masses = dict.fromkeys(ids, 1.0)
    before = geometry_metrics(dict(zip(ids, xyz)), masses)
    assert before.radius_of_gyration == pytest.approx(1)
    assert before.eigenvalues == pytest.approx([1 / 3] * 3)
    assert before.kappa_squared == pytest.approx(0, abs=1e-15)
    xyz[0] *= 2
    a = geometry_metrics(dict(zip(ids, xyz)), masses)
    q, _ = np.linalg.qr(np.array([[1.0, 2, 3], [3, 2, 4], [5, 2, 1]]))
    b = geometry_metrics(dict(zip(ids, xyz @ q.T + [7, -9, 4])), masses)
    np.testing.assert_allclose(
        b.gyration_tensor, q @ a.gyration_tensor @ q.T, atol=1e-14
    )
    np.testing.assert_allclose(a.eigenvalues, b.eigenvalues, atol=1e-14)
    assert b.radius_of_gyration == pytest.approx(a.radius_of_gyration)
    assert b.kappa_squared == pytest.approx(a.kappa_squared)


def test_selections_zero_spread_and_endpoints_independent():
    xyz = {8: (0, 0, 0), 22: (0, 0, 0), 91: (4, 0, 0)}
    masses = {8: 12.0, 22: 16.0, 91: 1.0}
    numbers = {8: 6, 22: 8, 91: 1}
    heavy = geometry_metrics(
        xyz,
        masses,
        atomic_numbers=numbers,
        options=GeometryOptions(selection="heavy", endpoints=(8, 91)),
    )
    assert heavy.selected_ids == (8, 22)
    assert heavy.radius_of_gyration == 0
    assert heavy.kappa_squared is None and heavy.re_over_rg is None
    assert heavy.end_to_end_distance == 4
    single = geometry_metrics(
        xyz, masses, options=GeometryOptions("explicit", (91,), endpoints=(91, 91))
    )
    assert single.end_to_end_distance == 0 and single.kappa_squared is None
    assert "same_site_endpoints: Re=0" in single.diagnostics
    json.dumps(single.to_dict(), allow_nan=False)


@pytest.mark.parametrize(
    "change",
    [
        {"coordinates": {8: (1, 2)}},
        {"coordinates": {8: (float("nan"), 0, 0)}},
        {"coordinates": {8: (float("inf"), 0, 0)}},
        {"coordinates": {9: (0, 0, 0)}},
        {"masses": {8: 0}},
        {"masses": {8: -1}},
        {"masses": {8: float("inf")}},
        {"coordinate_unit": "nm"},
        {"mass_unit": "kg"},
        {"options": GeometryOptions("explicit", (999,))},
        {"options": GeometryOptions("heavy")},
        {"options": GeometryOptions(endpoints=(8, 9))},
    ],
)
def test_invalid_data(change):
    args = {"coordinates": {8: (0, 0, 0)}, "masses": {8: 1}}
    args.update(change)
    with pytest.raises(AnalysisError):
        geometry_metrics(**args)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"site_ids": (1,)},
        {"selection": "explicit", "site_ids": ()},
        {"selection": "explicit", "site_ids": (1, 1)},
        {"endpoints": (True, 1)},
        {"weighting": "volume"},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(AnalysisError):
        GeometryOptions(**kwargs)


def test_overflow_and_eigenvalue_policy(monkeypatch):
    with pytest.raises(AnalysisError):
        geometry_metrics({1: (1e308, 0, 0), 2: (-1e308, 0, 0)}, {1: 1, 2: 1})
    for small, accepted in [(-1e-16, True), (-0.01, False)]:
        monkeypatch.setattr(
            np.linalg, "eigvalsh", lambda x, small=small: np.array([small, 0.0, 1.0])
        )
        if accepted:
            result = geometry_metrics({1: (0, 0, 0), 2: (1, 0, 0)}, {1: 1, 2: 1})
            assert result.eigenvalues[0] == 0
            assert any("tiny_negative" in s for s in result.diagnostics)
        else:
            with pytest.raises(AnalysisError, match="Substantial negative"):
                geometry_metrics({1: (0, 0, 0), 2: (1, 0, 0)}, {1: 1, 2: 1})


def test_numpy_only_import():
    code = """
import sys
class Block:
    def find_spec(self, name, *args):
        if name.split('.')[0] in {'rdkit', 'openmm', 'parmed', 'scipy'}:
            raise AssertionError('optional dependency imported: '+name)
sys.meta_path.insert(0, Block())
from island.analysis import geometry_metrics
assert geometry_metrics({7:(0,0,0)}, {7:1}).radius_of_gyration == 0
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True)


def test_provenance_endpoints_and_missing_diagnostic():
    pytest.importorskip("rdkit")
    from island.analysis.workflow import resolve_endpoints
    from island.builders import build_linear_polymer

    system = build_linear_polymer("[*]CC[*]", dp=3, coordinate_method="local_templates")
    ids, convention, reason = resolve_endpoints(system)
    assert ids == (
        system.metadata["polymer"]["head_site_id"],
        system.metadata["polymer"]["tail_site_id"],
    )
    assert convention == "builder_head_tail" and reason is None
    system.metadata["polymer"]["head_site_id"] = next(
        s for s, a in system.topology.sites.items() if a.atomic_number == 1
    )
    ids, convention, reason = resolve_endpoints(system)
    assert ids is None and "provenance" in reason
    system.metadata.pop("polymer")
    assert resolve_endpoints(system)[0] is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"stride": 0},
        {"stride": True},
        {"start_step": 0.5},
        {"start_time_ps": float("inf")},
        {"start_time_ps": -1},
        {"start_step": 3, "end_step": 2},
    ],
)
def test_analysis_windows_invalid(kwargs):
    from island.analysis import AnalysisOptions

    with pytest.raises(AnalysisError):
        AnalysisOptions(**kwargs)


def test_empty_heavy_and_bool_coordinates():
    with pytest.raises(AnalysisError):
        geometry_metrics(
            {1: (0, 0, 0)},
            {1: 1},
            atomic_numbers={1: 1},
            options=GeometryOptions(selection="heavy"),
        )
    with pytest.raises(AnalysisError):
        geometry_metrics({1: (True, False, False)}, {1: 1})
