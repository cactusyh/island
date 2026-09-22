from copy import deepcopy

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from island.builders import build_linear_polymer
from island.chemistry import to_rdkit
from island.conformations import (
    ConformationResult,
    LocalTemplateConformationGenerator,
    generate_polymer_conformation,
)
from island.conformations.local_validation import (
    authoritative_cip_assignments,
    validate_coordinate_stereochemistry,
)
from island.core import Coordinates
from island.exceptions import UnsupportedConformationError

ATTACHMENT_CHIRAL_REPEATS = (
    "[*:1]C[C@H](F)[*:2]",
    "[*:1][C@H](F)C[*:2]",
)


def coordinate_cip_by_site(system):
    converted = to_rdkit(system)
    Chem.RemoveStereochemistry(converted.mol)
    Chem.AssignStereochemistryFrom3D(converted.mol, replaceExistingTags=True)
    return {
        site_id: converted.mol.GetAtomWithIdx(index).GetProp("_CIPCode")
        for site_id, index in converted.site_id_to_rdkit_index.items()
        if converted.mol.GetAtomWithIdx(index).HasProp("_CIPCode")
    }


def reverse_topology_insertion_order(system):
    reordered = system.copy()
    reordered.topology.sites = dict(reversed(list(reordered.topology.sites.items())))
    reordered.topology.bonds = dict(reversed(list(reordered.topology.bonds.items())))
    reordered.topology.angles = dict(reversed(list(reordered.topology.angles.items())))
    reordered.topology.dihedrals = dict(
        reversed(list(reordered.topology.dihedrals.items()))
    )
    return reordered


@pytest.mark.parametrize("psmiles", ATTACHMENT_CHIRAL_REPEATS)
@pytest.mark.parametrize("template_seed,assembly_seed", [(1, 2), (31, 44), (2026, 7)])
@pytest.mark.parametrize("reverse_order", [False, True])
def test_unspecified_tacticity_preserves_all_authoritative_cip_sites(
    psmiles, template_seed, assembly_seed, reverse_order
):
    source = build_linear_polymer(psmiles, dp=6, generate_3d=False)
    expected = authoritative_cip_assignments(source)
    assert len(expected) == 5
    if reverse_order:
        source = reverse_topology_insertion_order(source)
    result = LocalTemplateConformationGenerator(
        template_seed=template_seed, assembly_seed=assembly_seed
    ).generate(source)
    generated = result.apply_to(source)
    assert coordinate_cip_by_site(generated) == expected
    validation = result.metadata["stereochemistry_validation"]
    assert validation["status"] == "validated"
    assert validation["validated_site_ids"] == sorted(expected)


@pytest.mark.parametrize(
    ("tacticity", "expected"),
    [
        ("isotactic", ["R"] * 6),
        ("syndiotactic", ["R", "S", "R", "S", "R", "S"]),
    ],
)
def test_tacticity_validation_remains_coordinate_derived(tacticity, expected):
    system = build_linear_polymer(
        "[*:1]N[C@H](F)C[*:2]",
        dp=6,
        tacticity=tacticity,
        coordinate_method="local_templates",
        template_seed=12,
        assembly_seed=21,
    )
    assert system.metadata["polymer"]["stereochemical_sequence"] == expected
    assert list(coordinate_cip_by_site(system).values()) == expected
    assert (
        system.metadata["polymer"]["coordinate_generation"][
            "stereochemistry_validation"
        ]["status"]
        == "validated"
    )


def test_stereochemistry_mismatch_diagnostics_identify_site_repeat_and_states():
    system = build_linear_polymer(
        ATTACHMENT_CHIRAL_REPEATS[0],
        dp=6,
        coordinate_method="local_templates",
    )
    expected = authoritative_cip_assignments(system)
    reflected = {
        site_id: system.coordinates.get(site_id) * np.array([-1.0, 1.0, 1.0])
        for site_id in system.topology.sites
    }
    with pytest.raises(
        UnsupportedConformationError,
        match=r"site_id=\d+, repeat_index=\d+, expected=[RS], observed=[RS]",
    ):
        validate_coordinate_stereochemistry(system, reflected, expected)


def test_apply_to_updates_coordinate_provenance_and_enables_random_walk():
    source = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    original = deepcopy(source.to_dict())
    result = LocalTemplateConformationGenerator(
        template_seed=17, assembly_seed=18
    ).generate(source)
    generated = result.apply_to(source)
    assert source.to_dict() == original
    assert generated.metadata["polymer"]["coordinates"] == (
        "local_templates_etkdg_incremental"
    )
    provenance = generated.metadata["polymer"]["coordinate_generation"]
    assert provenance["method"] == "local_templates_self_avoiding_random_walk"
    assert provenance["template_seed"] == 17
    assert provenance["assembly_seed"] == 18
    assert provenance["attempts"] == result.attempts
    assert generate_polymer_conformation(generated, seed=2026).success


def test_apply_to_in_place_validates_before_mutating():
    system = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    before = deepcopy(system.to_dict())
    invalid = ConformationResult(
        coordinates=Coordinates({next(iter(system.topology.sites)): [0, 0, 0]}),
        method="invalid",
        seed=1,
        success=True,
        attempts=0,
        rejected_trials=0,
        rollback_count=0,
    )
    with pytest.raises(Exception, match="missing coordinates"):
        invalid.apply_to(system, copy=False)
    assert system.to_dict() == before


def test_direct_generator_rejects_missing_required_hydrogens_before_embedding(
    monkeypatch,
):
    source = build_linear_polymer(
        "[*]CC[*]", dp=3, generate_3d=False, add_hydrogens=False
    )
    calls = 0

    def instrumented(*args, **kwargs):
        nonlocal calls
        calls += 1
        return 0

    monkeypatch.setattr(AllChem, "EmbedMolecule", instrumented)
    with pytest.raises(UnsupportedConformationError, match="explicit sites"):
        LocalTemplateConformationGenerator().generate(source)
    assert calls == 0


def test_genuinely_hydrogen_free_graph_is_allowed():
    source = build_linear_polymer(
        "[*:1][Si](F)(F)[Si](F)(F)[*:2]",
        dp=1,
        generate_3d=False,
        add_hydrogens=False,
    )
    assert all(site.element != "H" for site in source.topology.sites.values())
    result = LocalTemplateConformationGenerator().generate(source)
    assert result.success
    validation = result.metadata["stereochemistry_validation"]
    assert validation["status"] == "not_applicable"
    assert result.metadata["stereochemistry_validated_from_3d"] is False
