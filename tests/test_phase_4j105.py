"""Coverage reporting must retain failed cases and exercise final polymer graphs."""

import runpy
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.topology import molecule
from island.workflows.recipe_config import PolymerRecipe

benchmark = runpy.run_path(
    str(
        Path(__file__).resolve().parents[1]
        / "scripts/benchmark_polymer_architectures.py"
    )
)


def test_panel_constructs_connected_block_and_branch_graphs():
    from rdkit import Chem

    cases = benchmark["cases"]()
    assert len(cases) == len({case["name"] for case in cases}) == 54
    graphs = {}
    for case in cases:
        recipe = PolymerRecipe.from_dict(case["recipe"])
        key = (case["chemistry"], recipe.architecture)
        if key in graphs:
            continue
        graph = normalize(recipe.polymer_spec())
        assert not graph.blockers, (key, graph.blockers)
        mol = molecule(graph)
        assert len(Chem.GetMolFrags(mol)) == 1
        assert all(a.GetAtomicNum() > 0 for a in mol.GetAtoms())
        graphs[key] = mol
    for chemistry in benchmark["REPEATS"]:
        # Replacing one C-H with ethyl adds C2H4 to the final explicit-H graph.
        assert graphs[chemistry, "branched"].GetNumAtoms() == (
            graphs[chemistry, "homopolymer"].GetNumAtoms() + 6
        )


def test_report_does_not_drop_early_failures_or_promote_assignment_to_md():
    records = [
        {
            "name": "ok",
            "chemistry": "PE",
            "architecture": "block",
            "field": {"name": "opls"},
            "graph": True,
            "packed": True,
            "assignment": True,
            "single_point": True,
            "md": True,
        },
        {
            "name": "failed-md",
            "chemistry": "PA6",
            "architecture": "block",
            "field": {"name": "opls"},
            "graph": True,
            "packed": True,
            "assignment": True,
            "single_point": True,
            "md": False,
        },
        {
            "name": "failed-graph",
            "chemistry": "PA6",
            "architecture": "branched",
            "field": {"name": "pcff"},
            "graph": False,
            "packed": False,
            "assignment": False,
            "single_point": False,
            "md": False,
        },
    ]
    summary = benchmark["summarize"](records, 54)
    assert not summary["complete"]
    assert summary["equilibration"] == "unmeasured"
    assert summary["groups"]["all"] == {
        "attempted": 3,
        "graph": 2,
        "packed": 2,
        "assignment": 2,
        "single_point": 2,
        "md": 1,
    }
    assert summary["groups"]["field:opls"]["attempted"] == 2
    assert summary["groups"]["chemistry:PA6"]["md"] == 0


def test_incomplete_or_nonfinite_single_point_is_not_counted(tmp_path):
    path = tmp_path / "start.txt"
    for contents in ("", "nan", "inf", "truncated output"):
        path.write_text(contents)
        assert not benchmark["finite_single_point"]([path])
    path.write_text("-12.34\n")
    assert benchmark["finite_single_point"]([path])


@pytest.mark.parametrize("failure_count, geometric", [(1, True), (4, True), (1, False)])
def test_packing_retries_conformers_at_fixed_density_and_retains_failures(
    tmp_path, monkeypatch, failure_count, geometric
):
    import json

    import numpy as np
    from rdkit import Chem

    from island.epoxy import dynamics
    from island.workflows import recipe as workflow

    case = benchmark["cases"]()[0]["recipe"]
    recipe = PolymerRecipe.from_dict(case)
    seen = []

    def build(recipe, out, *, random_coordinates, conformer_seed):
        seen.append((random_coordinates, conformer_seed))
        (out / "chain.json").write_text("{}")
        return Chem.AddHs(Chem.MolFromSmiles("CC"))

    def pack(monomers, counts, density, seed, executable, out):
        assert density == recipe.density_g_cm3
        assert seed == recipe.seed and counts == [recipe.chains]
        for name in ("monomer0.xyz", "pack.in", "packed.xyz"):
            (out / name).write_text("retained input")
        if len(seen) <= failure_count:
            (out / "pack.log").write_text(
                "ERROR: Packmol was unable to put the molecules"
                if geometric
                else "Invalid input keyword"
            )
            raise RuntimeError("geometry did not fit")
        (out / "pack.log").write_text("Success!")
        return np.zeros((8 * recipe.chains, 3)), np.ones(3) * 20

    monkeypatch.setattr(workflow, "_build_chain", build)
    monkeypatch.setattr(dynamics, "pack", pack)
    if failure_count == 4 or not geometric:
        with pytest.raises(RuntimeError, match="geometry did not fit"):
            workflow._pack_chains(recipe, tmp_path, "packmol")
    else:
        workflow._pack_chains(recipe, tmp_path, "packmol")
        assert (tmp_path / "coordinates.npz").is_file()
    attempts = json.loads((tmp_path / "packing-attempts.json").read_text())
    assert len(attempts) == (min(failure_count + 1, 4) if geometric else 1)
    assert seen[0] == (False, recipe.seed)
    if geometric:
        assert seen[1] == (True, recipe.seed)
    assert all((tmp_path / a["directory"] / "pack.log").is_file() for a in attempts)
    assert sum(a["status"] == "failed" for a in attempts) == failure_count


def test_random_coordinate_retry_keeps_graph_and_stereochemistry(tmp_path):
    from rdkit import Chem

    from island.workflows.recipe import _build_chain

    case = benchmark["cases"]()[0]["recipe"]
    case["polymer"] = {"repeat_units": {"A": "*C[C@H](F)O*"}, "dp": 3}
    recipe = PolymerRecipe.from_dict(case)
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    original = _build_chain(recipe, first)
    retried = _build_chain(
        recipe, second, random_coordinates=True, conformer_seed=104902
    )
    assert Chem.MolToSmiles(original) == Chem.MolToSmiles(retried)
    assert (first / "chain.json").read_bytes() == (second / "chain.json").read_bytes()
