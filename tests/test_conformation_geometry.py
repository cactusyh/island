"""Physical-coordinate regressions independent of stored stereochemical tags."""

from itertools import combinations

import numpy as np
import pytest
from rdkit import Chem

from island import SimulationBox
from island.builders import build_linear_polymer
from island.chemistry import to_rdkit
from island.conformations import (
    FixedDistanceStericPolicy,
    UniformTorsionSampler,
    VanDerWaalsStericPolicy,
    generate_polymer_conformation,
)
from island.conformations.sterics import build_excluded_pairs
from island.exceptions import UnsupportedConformationError


def geometry_cip(system):
    converted = to_rdkit(system)
    # to_rdkit restores metadata tags: explicitly discard them before measuring.
    Chem.RemoveStereochemistry(converted.mol)
    Chem.AssignStereochemistryFrom3D(converted.mol, replaceExistingTags=True)
    return {
        site_id: atom.GetProp("_CIPCode")
        for site_id, index in converted.site_id_to_rdkit_index.items()
        if (atom := converted.mol.GetAtomWithIdx(index)).HasProp("_CIPCode")
    }


def local_geometry(system):
    lengths = {
        bond.key: np.linalg.norm(
            system.coordinates.get(bond.site1) - system.coordinates.get(bond.site2)
        )
        for bond in system.topology.bonds.values()
    }
    cosines = {}
    for center in system.topology.sites:
        for left, right in combinations(sorted(system.topology.neighbors(center)), 2):
            a = system.coordinates.get(left) - system.coordinates.get(center)
            b = system.coordinates.get(right) - system.coordinates.get(center)
            cosines[left, center, right] = np.dot(a, b) / (
                np.linalg.norm(a) * np.linalg.norm(b)
            )
    return lengths, cosines


@pytest.mark.parametrize("seed", [31, 2026])
@pytest.mark.parametrize(
    "psmiles,tacticity",
    [
        ("[*:1]N[C@H](F)C[*:2]", "syndiotactic"),
        # Attachment-centered chirality was inverted by the previous algorithm.
        ("[*:1]C[C@H](F)[*:2]", None),
        ("[*:1][C@H](F)C[*:2]", None),
    ],
)
def test_actual_3d_chirality_and_all_bond_angles_survive(psmiles, tacticity, seed):
    system = build_linear_polymer(psmiles, dp=6, tacticity=tacticity)
    before = geometry_cip(system)
    assert len(before) >= 5
    if tacticity:
        assert list(before.values()) == ["R", "S", "R", "S", "R", "S"]
    generated = generate_polymer_conformation(system, seed=seed).apply_to(system)
    assert geometry_cip(generated) == before
    for original, actual in zip(
        local_geometry(system), local_geometry(generated), strict=True
    ):
        assert actual == pytest.approx(original, abs=1e-10)


def test_two_dimensional_source_is_rejected_without_mutation():
    system = build_linear_polymer(
        "[*:1]N[C@H](F)C[*:2]", dp=6, tacticity="syndiotactic", generate_3d=False
    )
    before = system.to_dict()
    assert geometry_cip(system) == {}
    with pytest.raises(UnsupportedConformationError, match="3D source"):
        generate_polymer_conformation(system)
    assert system.to_dict() == before


def test_flattened_source_is_detected_even_with_3d_metadata():
    system = build_linear_polymer("[*]CC[*]", dp=3)
    for site_id in system.topology.sites:
        point = system.coordinates.get(site_id)
        point[2] = 0.0
        system.coordinates.set(site_id, point)
    with pytest.raises(UnsupportedConformationError, match="Degenerate tetrahedral"):
        generate_polymer_conformation(system)


def test_planar_aromatic_repeat_is_not_rejected_as_tetrahedral():
    system = build_linear_polymer("[*]c1ccccc1[*]", dp=3)
    assert generate_polymer_conformation(system).success


def test_rigid_motion_of_source_does_not_change_the_sampled_internal_geometry():
    system = build_linear_polymer("[*]CC[*]", dp=5)
    moved = system.copy()
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    for site_id in moved.topology.sites:
        moved.coordinates.set(
            site_id, rotation @ system.coordinates.get(site_id) + [8, 2, -3]
        )
    a = generate_polymer_conformation(system, seed=31)
    b = generate_polymer_conformation(moved, seed=31)
    assert a.attempts == b.attempts
    assert a.rejected_trials == b.rejected_trials
    for site_id in system.topology.sites:
        assert np.allclose(
            b.coordinates.get(site_id), rotation @ a.coordinates.get(site_id)
        )


def test_exclusions_do_not_depend_on_cached_derived_interactions():
    system = build_linear_polymer("[*]CC[*]", dp=4)
    expected = build_excluded_pairs(system)
    expected14 = build_excluded_pairs(system, exclude_one_four=True)
    system.topology.angles.clear()
    system.topology.dihedrals.clear()
    assert build_excluded_pairs(system) == expected
    assert build_excluded_pairs(system, exclude_one_four=True) == expected14
    result = generate_polymer_conformation(system, exclude_one_four=True)
    assert result.metadata["one_four_pairs_checked"] is False
    assert result.metadata["excluded_pairs"] == "1-2, 1-3 and 1-4"


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, 0])
@pytest.mark.parametrize("policy", [FixedDistanceStericPolicy, VanDerWaalsStericPolicy])
def test_steric_cutoffs_must_be_finite_and_positive(policy, value):
    with pytest.raises(ValueError):
        policy(value)


def test_invalid_sampler_output_is_rejected():
    class NonFiniteSampler(UniformTorsionSampler):
        def sample(self, rng, *, repeat_index, trial):
            return float("nan")

    system = build_linear_polymer("[*]CC[*]", dp=3)
    with pytest.raises(ValueError, match="finite angle"):
        generate_polymer_conformation(system, torsion_sampler=NonFiniteSampler())


def test_periodic_source_is_rejected():
    system = build_linear_polymer("[*]CC[*]", dp=3)
    system.box = SimulationBox(30, 30, 30)
    with pytest.raises(UnsupportedConformationError, match="Periodic"):
        generate_polymer_conformation(system)


def test_bond_angle_override_is_not_silently_ignored():
    system = build_linear_polymer("[*]CC[*]", dp=3)
    with pytest.raises(UnsupportedConformationError, match="Bond-angle resampling"):
        generate_polymer_conformation(system, bond_angle_degrees=109.5)
