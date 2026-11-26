"""Independent chain conformers, mixtures and reproducible length/sequence draws."""

import json
import os
from copy import deepcopy

import numpy as np
import pytest

from island.workflows.recipe import run_recipe
from island.workflows.recipe_chains import expand
from island.workflows.recipe_config import PolymerRecipe


def mixture():
    return {
        "architecture": "mixture",
        "forcefield": {"name": "gaff2", "charge": "gasteiger"},
        "density_g_cm3": 0.2,
        "seed": 173,
        "components": [
            {
                "name": "PE",
                "polymer": {"repeat_units": {"A": "*CC*"}},
                "lengths": [3, 5],
            },
            {
                "name": "block",
                "architecture": "block",
                "polymer": {
                    "repeat_units": {"A": "*CC*", "B": "*CCO*"},
                    "blocks": [["A", 2], ["B", 2]],
                },
            },
            {
                "name": "branch",
                "architecture": "branched",
                "polymer": {
                    "repeat_units": {"A": "*CC*"},
                    "dp": 4,
                    "branches": [{"anchor_site": 4, "psmiles": "*CC"}],
                },
            },
        ],
        "md": {
            "nvt_steps": 20,
            "npt_steps": 20,
            "production_steps": 20,
            "sample_every": 10,
            "segment_steps": 20,
        },
    }


def test_distributions_are_repeatable_and_detached():
    raw = mixture()
    raw["components"] = [
        {
            "name": "random",
            "architecture": "copolymer",
            "chains": 8,
            "polymer": {"repeat_units": {"A": "*CC*", "B": "*CCO*"}},
            "length_distribution": {"values": [5, 10], "probabilities": [0.4, 0.6]},
            "sequence_probabilities": {"A": 0.6, "B": 0.4},
        }
    ]
    saved = deepcopy(raw)
    recipe = PolymerRecipe.from_dict(raw)
    first, second = expand(recipe), expand(recipe)
    assert first == second and raw == saved
    assert len({r.seed for _, r in first}) == 8
    assert {r.polymer_spec().degree_of_polymerization for _, r in first} == {5, 10}
    raw["seed"] = 719
    assert expand(PolymerRecipe.from_dict(raw)) != first


@pytest.mark.parametrize(
    "change",
    [
        {"lengths": [3]},
        {"length_distribution": {"values": [4], "probabilities": [0.8]}},
        {"architecture": "crosslinked"},
        {"lengths": [True, 3]},
    ],
)
def test_ambiguous_components_rejected(change):
    raw = mixture()
    raw["components"][0].update(chains=2, **change)
    with pytest.raises(ValueError):
        PolymerRecipe.from_dict(raw)


def test_real_mixture_preserves_chain_and_component_map(tmp_path):
    raw = mixture()
    kwargs = {
        "lammps": os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        "packmol": os.environ["ISLAND_PACKMOL"],
    }
    result = run_recipe(raw, tmp_path / "mix", **kwargs)
    assert result["status"] == "completed"
    root = tmp_path / "mix/packing"
    chains = json.loads((root / "chains.json").read_text())
    sites = json.loads((root / "site-map.json").read_text())
    assert [r["dp"] for r in chains] == [3, 5, 4, 4]
    assert len({r["conformer_seed"] for r in chains}) == 4
    assert len(sites) == sum(r["atoms"] for r in chains)
    assert {r["component"] for r in sites} == {"PE", "block", "branch"}
    for chain in chains:
        assert [
            s["chain_site_id"] for s in sites if s["chain"] == chain["chain"]
        ] == list(range(1, chain["atoms"] + 1))
    assert run_recipe(raw, tmp_path / "mix", **kwargs) == result


def test_equal_chemistry_has_independent_conformations(tmp_path):
    from rdkit import Chem

    raw = mixture()
    raw.pop("components")
    raw.update(
        architecture="homopolymer",
        polymer={"repeat_units": {"A": "*CC*"}, "dp": 10},
        chains=2,
    )
    run_recipe(
        raw,
        tmp_path / "pe",
        lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        packmol=os.environ["ISLAND_PACKMOL"],
        stop_after="packing",
    )
    mols = [
        Chem.MolFromXYZFile(str(tmp_path / f"pe/packing/attempt-000/monomer{i}.xyz"))
        for i in (0, 1)
    ]
    coordinates = [m.GetConformer().GetPositions() for m in mols]
    distances = [np.linalg.norm(x[:, None] - x[None, :], axis=2) for x in coordinates]
    assert not np.allclose(*distances, atol=1e-3)
