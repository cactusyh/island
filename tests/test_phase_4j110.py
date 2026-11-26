"""Adaptive diagnostics must preserve honest nonconvergence and resumable evidence."""

import json
import os

import numpy as np
import pytest

from island.exceptions import AnalysisError
from island.workflows.recipe import run_recipe
from island.workflows.recipe_config import PolymerRecipe
from island.workflows.recipe_equilibration import RecipeEquilibration, assess, structure


def recipe():
    return {
        "architecture": "homopolymer",
        "polymer": {"repeat_units": {"A": "*CC*"}, "dp": 3},
        "forcefield": {"name": "gaff2", "charge": "gasteiger"},
        "chains": 2,
        "density_g_cm3": 0.3,
        "md": {
            "nvt_steps": 20,
            "npt_steps": 20,
            "production_steps": 20,
            "sample_every": 10,
            "segment_steps": 20,
        },
        "equilibration": {"max_segments": 2, "steps_per_segment": 40},
    }


def stationary():
    rng = np.random.default_rng(173)
    return [
        {
            "temperature": 300 + rng.normal(),
            "potential": -1000 + rng.normal(),
            "density": 0.95 + rng.normal() * 0.001,
            "pressure": 1 + rng.normal(),
            "mean_bond_length": 1.4 + rng.normal() * 0.001,
            "mean_rg2": 10 + rng.normal() * 0.01,
            "minimum_bond_length": 1,
            "maximum_bond_length": 1.6,
            "winding_rank": 0,
            "components": 2,
        }
        for _ in range(1024)
    ]


def test_statistics_reject_drift_short_and_target_mismatch():
    rows = stationary()
    cfg = RecipeEquilibration()
    assert assess(rows, cfg, 300, 1)["passed"]
    assert not assess(rows[:20], cfg, 300, 1)["passed"]
    assert not assess(rows, cfg, 450, 1)["passed"]
    for i, row in enumerate(rows):
        row["density"] += i / len(rows)
    report = assess(rows, cfg, 300, 1)
    assert not report["passed"] and any("density" in s for s in report["failures"])


def test_structural_checks_and_periodic_unwrapping():
    xyz = np.array([[9.5, 1, 1], [0.5, 1, 1], [1.5, 1, 1]])
    first = structure(xyz, np.array([10] * 3), np.ones(3), [(0, 1), (1, 2)])
    second = structure(
        xyz + [10, -10, 20], np.array([10] * 3), np.ones(3), [(0, 1), (1, 2)]
    )
    assert first == second and first["mean_rg2"] == pytest.approx(2 / 3)
    rows = stationary()
    rows[-1]["winding_rank"] = 1
    assert (
        "periodic_winding_changed"
        in assess(rows, RecipeEquilibration(), 300, 1, network=True)["failures"]
    )
    rows[-1]["maximum_bond_length"] = 5
    assert "bond_geometry" in assess(rows, RecipeEquilibration(), 300, 1)["failures"]


@pytest.mark.parametrize(
    "config",
    [
        {"max_segments": 0},
        {"steps_per_segment": True},
        {"policy": {"minimum_frames": 1}},
        {"consecutive_passes": 11},
        {"pressure_tolerance_atm": float("nan")},
    ],
)
def test_invalid_budgets_rejected(config):
    with pytest.raises((ValueError, AnalysisError)):
        RecipeEquilibration(**config)


def test_real_budget_exhaustion_resume_and_tamper(tmp_path):
    raw = recipe()
    kwargs = {
        "lammps": os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        "packmol": os.environ["ISLAND_PACKMOL"],
    }
    root = tmp_path / "run"
    paused = run_recipe(raw, root, stop_after="equilibration", **kwargs)
    assert paused["status"] == "paused"
    before = (root / "equilibration-0001/final.restart").read_bytes()
    result = run_recipe(raw, root, **kwargs)
    assert result["status"] == "budget_exhausted" and result["equilibrated"] is False
    assert result["equilibration"]["stop_reason"] == "budget_exhausted"
    assert result["equilibration"]["segments"] == 2
    assert not (root / "production").exists()
    assert (root / result["checkpoint"]).exists()
    assert before == (root / "equilibration-0001/final.restart").read_bytes()
    assert run_recipe(raw, root, **kwargs) == result
    raw["equilibration"]["max_segments"] = 3
    with pytest.raises(ValueError, match="Recipe differs"):
        run_recipe(raw, root, **kwargs)
    raw["equilibration"]["max_segments"] = 2
    path = root / "equilibration-0001/observations.json"
    values = json.loads(path.read_text())
    values[0]["temperature"] = 300
    path.write_text(json.dumps(values))
    with pytest.raises(ValueError, match="Checkpoint content"):
        run_recipe(raw, root, **kwargs)


def test_explicit_exploratory_production_preserves_nonconvergence(tmp_path):
    raw = recipe()
    raw["equilibration"]["production_on_exhaustion"] = True
    result = run_recipe(
        raw,
        tmp_path / "run",
        lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        packmol=os.environ["ISLAND_PACKMOL"],
    )
    assert result["status"] == "completed" and result["equilibrated"] is False
    assert result["checkpoint"] == "production/final.restart"
    assert result["equilibration"]["stop_reason"] == "budget_exhausted"


def test_sample_alignment_required():
    raw = recipe()
    raw["equilibration"]["steps_per_segment"] = 33
    with pytest.raises(ValueError, match="align"):
        PolymerRecipe.from_dict(raw)


def test_controller_waits_for_sample_floor_and_consecutive_passes(
    tmp_path, monkeypatch
):
    from island.workflows import recipe_equilibration as controller
    from island.workflows import recipe_md
    from island.workflows.recipe import checkpoint, write_json

    raw = recipe()
    raw["equilibration"].update(max_segments=5, steps_per_segment=5120)
    cfg = PolymerRecipe.from_dict(raw)
    root = tmp_path / "controller"
    root.mkdir()
    with checkpoint(root, "npt") as (work, exists):
        assert not exists
        (work / "final.restart").write_text("unit-test restart placeholder")
    calls = []

    def segment(source, work, md, steps, ensemble, seed, lammps, **kwargs):
        calls.append(steps)
        (work / "final.restart").write_text("unit-test restart placeholder")
        write_json(
            work / "result.json", {"steps": steps, "time_fs": steps * md.timestep_fs}
        )

    monkeypatch.setattr(recipe_md, "run_segment", segment)
    monkeypatch.setattr(controller, "observations", lambda *args: stationary()[:512])
    report, source, outputs, paused = controller.run(
        root, cfg, root / "npt", "unused-unit-test"
    )
    assert report["converged"] and report["stop_reason"] == "diagnostics_passed"
    assert len(calls) == len(outputs) == 3 and source.name == "equilibration-0003"
    assert report["consecutive_passes"] == 2 and not paused
    again = controller.run(root, cfg, root / "npt", "unused-unit-test")
    assert again[0] == report and len(calls) == 3
