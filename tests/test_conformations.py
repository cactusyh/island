from copy import deepcopy

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from island import AtomSite, BeadSite, Coordinates, MolecularSystem, Topology
from island.builders import build_linear_polymer
from island.chemistry import to_rdkit
from island.conformations import (
    ConformationGenerator,
    ConformationResult,
    FixedDistanceStericPolicy,
    RetryPolicy,
    SelfAvoidingRandomWalkGenerator,
    UniformTorsionSampler,
    VanDerWaalsStericPolicy,
    generate_polymer_conformation,
)
from island.exceptions import (
    ConformationGenerationError,
    UnsupportedConformationError,
)


def coordinate_array(system: MolecularSystem, coordinates: Coordinates) -> np.ndarray:
    return np.array([coordinates.get(site_id) for site_id in system.topology.sites])


def chemical_snapshot(system: MolecularSystem) -> dict[str, object]:
    snapshot = deepcopy(system.to_dict())
    snapshot.pop("coordinates")
    polymer = snapshot.get("metadata", {}).get("polymer")
    if isinstance(polymer, dict):
        polymer.pop("coordinates", None)
        polymer.pop("coordinate_generation", None)
    return snapshot


def controlled_cip_sequence(system: MolecularSystem) -> list[str]:
    converted = to_rdkit(system)
    Chem.AssignStereochemistry(converted.mol, cleanIt=True, force=True)
    sequence: list[tuple[int, str]] = []
    for site_id, rdkit_index in converted.site_id_to_rdkit_index.items():
        site = system.topology.get_site(site_id)
        if site.metadata.get("controllable_stereocenter"):
            sequence.append(
                (
                    site.metadata["repeat_unit_index"],
                    converted.mol.GetAtomWithIdx(rdkit_index).GetProp("_CIPCode"),
                )
            )
    return [state for _, state in sorted(sequence)]


def test_general_conformation_api_is_structured_and_non_destructive() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=12, generate_3d=True)
    before_chemistry = chemical_snapshot(system)
    before_coordinates = coordinate_array(system, system.coordinates)

    result = generate_polymer_conformation(system, method="random_walk", seed=2026)

    assert isinstance(result, ConformationResult)
    assert result.success is True
    assert result.method == "self_avoiding_random_walk"
    assert result.seed == 2026
    assert chemical_snapshot(system) == before_chemistry
    assert np.allclose(coordinate_array(system, system.coordinates), before_coordinates)

    generated = result.apply_to(system)
    assert generated is not system
    assert chemical_snapshot(generated) == before_chemistry
    assert np.allclose(
        coordinate_array(generated, generated.coordinates),
        coordinate_array(system, result.coordinates),
    )
    assert not np.allclose(
        coordinate_array(generated, generated.coordinates), before_coordinates
    )


def test_random_walk_is_deterministic_and_seed_sensitive() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=15, generate_3d=True)
    first = generate_polymer_conformation(system, seed=12)
    second = generate_polymer_conformation(system, seed=12)
    different = generate_polymer_conformation(system, seed=13)
    assert np.allclose(
        coordinate_array(system, first.coordinates),
        coordinate_array(system, second.coordinates),
    )
    assert first.metadata == second.metadata
    assert first.attempts == second.attempts
    assert not np.allclose(
        coordinate_array(system, first.coordinates),
        coordinate_array(system, different.coordinates),
    )


def test_polyethylene_has_no_nonexcluded_clash_below_fixed_threshold() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=20, generate_3d=True)
    threshold = 1.0
    result = generate_polymer_conformation(
        system,
        seed=2026,
        steric_policy=FixedDistanceStericPolicy(threshold),
    )
    assert result.minimum_nonbonded_distance is not None
    assert result.minimum_nonbonded_distance >= threshold
    assert result.attempts >= 19
    assert result.rejected_trials == result.attempts - 19
    assert result.metadata["excluded_pairs"] == "1-2 and 1-3"
    assert result.metadata["one_four_pairs_checked"] is True


def test_heteroatom_polymer_generation() -> None:
    system = build_linear_polymer("[*]CO[*]", dp=10, generate_3d=True)
    result = generate_polymer_conformation(system, seed=7)
    generated = result.apply_to(system)
    assert generated.number_of_sites == system.number_of_sites
    assert len(generated.topology.connected_components()) == 1
    assert result.minimum_nonbonded_distance is not None


def test_tacticity_and_stereochemical_provenance_survive_coordinate_generation() -> (
    None
):
    system = build_linear_polymer(
        "[*:1]N[C@H](F)C[*:2]",
        dp=6,
        tacticity="syndiotactic",
        generate_3d=True,
    )
    expected = ["R", "S", "R", "S", "R", "S"]
    before_chemistry = chemical_snapshot(system)
    result = generate_polymer_conformation(system, seed=31)
    generated = result.apply_to(system)
    assert chemical_snapshot(generated) == before_chemistry
    assert controlled_cip_sequence(generated) == expected
    assert generated.metadata["polymer"]["stereochemical_sequence"] == expected
    for site_id, original_site in system.topology.sites.items():
        assert generated.topology.get_site(site_id).metadata == original_site.metadata


def test_ring_containing_repeat_units_remain_rigid() -> None:
    system = build_linear_polymer("[*]c1ccccc1[*]", dp=4, generate_3d=True)
    result = generate_polymer_conformation(system, seed=99)
    ring_ids = [
        site.id
        for site in system.topology.sites.values()
        if site.metadata["repeat_unit_index"] == 2 and site.element == "C"
    ]
    original = np.array([system.coordinates.get(site_id) for site_id in ring_ids])
    generated = np.array([result.coordinates.get(site_id) for site_id in ring_ids])
    original_distances = np.linalg.norm(original[:, None] - original[None, :], axis=2)
    generated_distances = np.linalg.norm(
        generated[:, None] - generated[None, :], axis=2
    )
    assert np.allclose(original_distances, generated_distances)
    assert result.metadata["repeat_units_are_rigid"] is True
    assert result.metadata["ring_intersection_check"] is False


def test_retry_failure_contains_diagnostics() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=True)
    generator = SelfAvoidingRandomWalkGenerator(
        seed=1,
        steric_policy=FixedDistanceStericPolicy(100.0),
        retry_policy=RetryPolicy(
            attempts_per_unit=2, rollback_units=1, max_rollbacks=1
        ),
    )
    with pytest.raises(ConformationGenerationError, match="repeat 1") as error:
        generator.generate(system)
    assert error.value.repeat_index == 1
    assert error.value.attempts == 4
    assert error.value.rejected_trials == 4
    assert error.value.rollback_count == 1


def test_successful_generation_can_report_rollback() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=20, generate_3d=True)
    result = generate_polymer_conformation(
        system,
        seed=2026,
        steric_policy=FixedDistanceStericPolicy(1.2),
        retry_policy=RetryPolicy(
            attempts_per_unit=5, rollback_units=3, max_rollbacks=20
        ),
    )
    assert result.rollback_count >= 1
    assert result.rejected_trials > 0
    assert result.minimum_nonbonded_distance >= 1.2


def test_pe_dp50_scalability_smoke() -> None:
    # Supply genuine 3D geometry explicitly: the default whole-chain ETKDG
    # distance-matrix initialization can fail at DP50. This tests the walker,
    # not a local-template or end-to-end long-chain building workflow.
    system = build_linear_polymer("[*]CC[*]", dp=50, generate_3d=False)
    converted = to_rdkit(system)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 2026
    parameters.useRandomCoords = True
    assert AllChem.EmbedMolecule(converted.mol, parameters) == 0
    conformer = converted.mol.GetConformer()
    for site_id, index in converted.site_id_to_rdkit_index.items():
        point = conformer.GetAtomPosition(index)
        system.coordinates.set(site_id, [point.x, point.y, point.z])
    system.metadata["polymer"]["coordinates"] = (
        "test_etkdg_random_coordinate_initialization"
    )
    result = generate_polymer_conformation(system, seed=2026)
    assert result.success
    assert len(result.coordinates) == system.number_of_sites
    assert result.minimum_nonbonded_distance >= 0.8


def test_unsupported_inputs_fail_clearly() -> None:
    topology = Topology()
    topology.add_site(BeadSite(id=1, name="B", mass=10.0, bead_type="B"))
    coarse = MolecularSystem(
        topology, Coordinates({1: [0, 0, 0]}), representation="coarse_grained"
    )
    with pytest.raises(UnsupportedConformationError, match="atomistic"):
        generate_polymer_conformation(coarse)

    atom_topology = Topology()
    atom_topology.add_site(
        AtomSite(id=1, name="C", mass=12.0, element="C", atomic_number=6)
    )
    molecule = MolecularSystem(atom_topology, Coordinates({1: [0, 0, 0]}))
    with pytest.raises(UnsupportedConformationError, match="polymer metadata"):
        generate_polymer_conformation(molecule)
    with pytest.raises(UnsupportedConformationError, match="Unsupported"):
        generate_polymer_conformation(molecule, method="unknown")


def test_sampler_policy_and_interface_types_are_public() -> None:
    assert issubclass(SelfAvoidingRandomWalkGenerator, ConformationGenerator)
    assert isinstance(UniformTorsionSampler(), UniformTorsionSampler)
    system = build_linear_polymer("[*]CO[*]", dp=3, generate_3d=True)
    vdw = VanDerWaalsStericPolicy(scale=0.4)
    assert vdw.minimum_distance(system, 1, 2) > 0
