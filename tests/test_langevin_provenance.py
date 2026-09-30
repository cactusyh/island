"""Dynamics state changes retire prior coordinate producers without rewriting history."""

from copy import deepcopy

import pytest

pytest.importorskip("scipy")
pytest.importorskip("rdkit")
from test_coordinate_provenance import Target
from test_dynamics import Harmonic
from test_langevin import options

from island.builders import build_linear_polymer
from island.conformations import ConformationResult
from island.dynamics import run_langevin
from island.minimization import minimize_geometry


def test_minimization_dynamics_and_later_producers():
    system = build_linear_polymer(
        "[*:1]CC[*:2]",
        dp=3,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    system.metadata["ambertools_preparation"] = {"historical_signature": "unchanged"}
    minimized = minimize_geometry(system, Target(system)).to_system(system)
    before = deepcopy(minimized.to_dict())
    velocities = {site: (0.1, -0.1, 0.2) for site in system.topology.sites}
    result = run_langevin(minimized, Harmonic(0), velocities, options(2))
    moved = result.to_system(minimized)
    moved.validate()
    assert minimized.to_dict() == before
    assert (
        moved.metadata["coordinate_source"]
        == moved.metadata["polymer"]["coordinates"]
        == "langevin_dynamics"
    )
    assert (
        moved.metadata["polymer"]["coordinate_generation"]["coordinate_source"]
        == "langevin_dynamics"
    )
    assert "local_minimization" not in moved.metadata
    assert (
        moved.metadata["coordinate_history"][-1]["records"]["local_minimization"]
        == minimized.metadata["local_minimization"]
    )
    repeated = run_langevin(moved, Harmonic(0), velocities, options(2)).to_system(moved)
    assert (
        repeated.metadata["coordinate_history"][-1]["records"]["dynamics"]
        == moved.metadata["dynamics"]
    )
    assert repeated.metadata["coordinate_source"] == "langevin_dynamics"
    next_minimum = minimize_geometry(moved, Target(moved)).to_system(moved)
    assert "dynamics" not in next_minimum.metadata
    assert (
        next_minimum.metadata["coordinate_history"][-1]["records"]["dynamics"]
        == moved.metadata["dynamics"]
    )
    conformation = ConformationResult(
        moved.coordinates.copy(), "test", 1, True, 1, 0, 0
    )
    conformed = conformation.apply_to(moved)
    assert "dynamics" not in conformed.metadata
    assert (
        conformed.metadata["coordinate_history"][-1]["records"]["dynamics"]
        == moved.metadata["dynamics"]
    )
    for output in (moved, next_minimum, conformed):
        output.validate()
        assert (
            output.metadata["ambertools_preparation"]
            == system.metadata["ambertools_preparation"]
        )
        assert output.topology.sites == system.topology.sites
