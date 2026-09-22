"""Sequence-driven construction of finite, linear polymers from PSMILES."""

import math
import random
from collections.abc import Mapping, Sequence
from numbers import Real
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

from island.builders._tacticity import prepare_stereochemical_sequence
from island.builders.sequence import PolymerSequence
from island.chemistry._stereo_assignment import assign_cip_sequence, verify_cip_sequence
from island.chemistry.psmiles import RepeatUnit, parse_psmiles
from island.chemistry.rdkit_adapter import from_rdkit
from island.chemistry.stereochemistry import (
    StereochemicalSequence,
    require_controllable_stereocenter,
)
from island.core import MolecularSystem
from island.exceptions import (
    EmbeddingError,
    InvalidRepeatUnitError,
    PolymerBuildError,
)

_CHAIN_ID_PROPERTY = "_island_chain_id"
_REPEAT_INDEX_PROPERTY = "_island_repeat_unit_index"
_REPEAT_TYPE_PROPERTY = "_island_repeat_unit_type"
_SOURCE_INDEX_PROPERTY = "_island_source_repeat_atom_index"
_GENERATED_HYDROGEN_PROPERTY = "_island_generated_hydrogen"
_STEREO_STATE_PROPERTY = "_island_repeat_unit_stereochemical_state"
_CONTROLLED_CENTER_PROPERTY = "_island_controllable_stereocenter"


def build_polymer_from_sequence(
    repeat_units: Mapping[str, str | RepeatUnit],
    sequence: Sequence[str] | PolymerSequence,
    *,
    polymer_type: str = "sequence_defined",
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
    tacticity: str | None = None,
    stereo_seed: int = 2026,
    atactic_fraction: float = 0.5,
) -> MolecularSystem:
    """Build a finite linear polymer from an explicit repeat-unit sequence."""
    _validate_build_options(random_seed, chain_id, polymer_type)
    if isinstance(sequence, str):
        raise PolymerBuildError(
            "sequence must be a sequence of repeat-unit labels, not a bare string"
        )
    polymer_sequence = (
        sequence
        if isinstance(sequence, PolymerSequence)
        else PolymerSequence(tuple(sequence))
    )
    library = _parse_repeat_unit_library(repeat_units)
    missing = set(polymer_sequence.repeat_unit_types) - set(library)
    if missing:
        raise PolymerBuildError(
            f"Sequence references undefined repeat-unit types: {sorted(missing)}"
        )
    stereochemical_sequence = prepare_stereochemical_sequence(
        library,
        polymer_sequence,
        tacticity=tacticity,
        stereo_seed=stereo_seed,
        atactic_fraction=atactic_fraction,
    )
    ordered_units = [library[identity] for identity in polymer_sequence.identities]
    polymer, unit_atom_maps, stereo_atom_indices = _assemble_heavy_atom_graph(
        ordered_units,
        polymer_sequence.identities,
        chain_id,
        stereochemical_sequence=stereochemical_sequence,
    )
    try:
        Chem.SanitizeMol(polymer)
    except Exception as error:
        raise PolymerBuildError(
            "The assembled head-to-tail polymer graph is chemically invalid"
        ) from error
    if stereochemical_sequence is not None:
        assign_cip_sequence(
            polymer, stereo_atom_indices, stereochemical_sequence.states
        )

    head_atom_index = unit_atom_maps[0][ordered_units[0].head.neighbor_atom_index]
    tail_atom_index = unit_atom_maps[-1][ordered_units[-1].tail.neighbor_atom_index]
    if add_hydrogens:
        polymer = Chem.AddHs(polymer)
        _annotate_generated_hydrogens(polymer)
    if stereochemical_sequence is not None:
        verify_cip_sequence(
            polymer, stereo_atom_indices, stereochemical_sequence.states
        )

    coordinate_kind = _generate_coordinates(
        polymer, generate_3d=generate_3d, random_seed=random_seed
    )
    conversion = from_rdkit(polymer)
    system = conversion.system
    _transfer_atom_provenance(system, polymer, conversion.rdkit_index_to_site_id)
    system.metadata["polymer"] = {
        "architecture": "linear",
        "polymer_type": polymer_type,
        "sequence": list(polymer_sequence.identities),
        "repeat_unit_definitions": {
            identity: repeat.original_psmiles for identity, repeat in library.items()
        },
        "degree_of_polymerization": polymer_sequence.number_of_repeat_units,
        "number_of_repeat_units": polymer_sequence.number_of_repeat_units,
        "composition_counts": polymer_sequence.composition_counts,
        "number_of_inter_repeat_unit_bonds": (
            polymer_sequence.number_of_repeat_units - 1
        ),
        "chain_id": chain_id,
        "head_site_id": conversion.rdkit_index_to_site_id[head_atom_index],
        "tail_site_id": conversion.rdkit_index_to_site_id[tail_atom_index],
        "coordinates": coordinate_kind,
    }
    if polymer_sequence.generation_metadata:
        system.metadata["polymer"]["sequence_generation"] = dict(
            polymer_sequence.generation_metadata
        )
    if stereochemical_sequence is not None:
        system.metadata["polymer"].update(
            {
                "tacticity": stereochemical_sequence.tacticity,
                "stereo_seed": stereochemical_sequence.stereo_seed,
                "atactic_fraction": stereochemical_sequence.inverted_fraction,
                "stereochemical_sequence": list(stereochemical_sequence.states),
                "stereochemical_reference_state": (
                    stereochemical_sequence.reference_state
                ),
                "stereochemistry_convention": "final_graph_absolute_cip",
            }
        )
    return system


def build_linear_polymer(
    psmiles: str,
    *,
    dp: int,
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
    tacticity: str | None = None,
    stereo_seed: int = 2026,
    atactic_fraction: float = 0.5,
) -> MolecularSystem:
    """Build a finite linear homopolymer with ``dp`` total repeat units."""
    _validate_dp(dp)
    system = build_polymer_from_sequence(
        {"A": psmiles},
        PolymerSequence(
            ("A",) * dp,
            {"strategy": "homopolymer", "degree_of_polymerization": dp},
        ),
        polymer_type="homopolymer",
        generate_3d=generate_3d,
        random_seed=random_seed,
        add_hydrogens=add_hydrogens,
        chain_id=chain_id,
        tacticity=tacticity,
        stereo_seed=stereo_seed,
        atactic_fraction=atactic_fraction,
    )
    system.metadata["polymer"]["source_psmiles"] = psmiles
    return system


def build_alternating_copolymer(
    psmiles_a: str,
    psmiles_b: str,
    *,
    dp: int,
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
) -> MolecularSystem:
    """Build A-B-A-B-... with ``dp`` total repeat units, beginning with A."""
    _validate_dp(dp)
    identities = tuple("A" if index % 2 == 0 else "B" for index in range(dp))
    return build_polymer_from_sequence(
        {"A": psmiles_a, "B": psmiles_b},
        PolymerSequence(
            identities,
            {"strategy": "alternating", "starts_with": "A"},
        ),
        polymer_type="alternating_copolymer",
        generate_3d=generate_3d,
        random_seed=random_seed,
        add_hydrogens=add_hydrogens,
        chain_id=chain_id,
    )


def build_block_copolymer(
    blocks: Sequence[tuple[str, int]],
    *,
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
) -> MolecularSystem:
    """Build a linear polymer from any number of ordered PSMILES/count blocks."""
    if not blocks:
        raise PolymerBuildError("At least one polymer block is required")
    definitions: dict[str, str] = {}
    labels_by_psmiles: dict[str, str] = {}
    identities: list[str] = []
    block_metadata: list[dict[str, Any]] = []
    for block_index, block in enumerate(blocks):
        if not isinstance(block, (tuple, list)) or len(block) != 2:
            raise PolymerBuildError("Each block must be a (psmiles, repeat_count) pair")
        psmiles, count = block
        _validate_dp(count, name=f"Block {block_index} repeat count")
        if not isinstance(psmiles, str) or not psmiles:
            raise PolymerBuildError("Each block PSMILES must be a non-empty string")
        if psmiles not in labels_by_psmiles:
            label = _default_repeat_label(len(labels_by_psmiles))
            labels_by_psmiles[psmiles] = label
            definitions[label] = psmiles
        label = labels_by_psmiles[psmiles]
        identities.extend([label] * count)
        block_metadata.append({"repeat_unit_type": label, "count": count})
    return build_polymer_from_sequence(
        definitions,
        PolymerSequence(
            tuple(identities), {"strategy": "block", "blocks": block_metadata}
        ),
        polymer_type="block_copolymer",
        generate_3d=generate_3d,
        random_seed=random_seed,
        add_hydrogens=add_hydrogens,
        chain_id=chain_id,
    )


def build_random_copolymer(
    repeat_units: Mapping[str, str | RepeatUnit],
    *,
    dp: int,
    fractions: Mapping[str, float],
    sequence_seed: int = 2026,
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
) -> MolecularSystem:
    """Build a random linear copolymer using largest-remainder allocation.

    Fractions must be non-negative and sum to one. Integer counts use the
    largest-remainder method, with repeat-unit insertion order breaking ties.
    Only a local ``random.Random(sequence_seed)`` instance shuffles the sequence.
    """
    _validate_dp(dp)
    if not isinstance(sequence_seed, int) or isinstance(sequence_seed, bool):
        raise PolymerBuildError("sequence_seed must be an integer")
    counts = _allocate_composition_counts(repeat_units, fractions, dp)
    identities = [
        identity for identity in repeat_units for _ in range(counts[identity])
    ]
    random.Random(sequence_seed).shuffle(identities)
    return build_polymer_from_sequence(
        repeat_units,
        PolymerSequence(
            tuple(identities),
            {
                "strategy": "random",
                "sequence_seed": sequence_seed,
                "requested_fractions": dict(fractions),
                "allocated_counts": counts,
                "allocation_method": "largest_remainder",
            },
        ),
        polymer_type="random_copolymer",
        generate_3d=generate_3d,
        random_seed=random_seed,
        add_hydrogens=add_hydrogens,
        chain_id=chain_id,
    )


def _parse_repeat_unit_library(
    definitions: Mapping[str, str | RepeatUnit],
) -> dict[str, RepeatUnit]:
    if not definitions:
        raise PolymerBuildError("At least one repeat-unit definition is required")
    library: dict[str, RepeatUnit] = {}
    for identity, definition in definitions.items():
        if not isinstance(identity, str) or not identity:
            raise PolymerBuildError("Repeat-unit identities must be non-empty strings")
        if isinstance(definition, RepeatUnit):
            repeat = definition
        elif isinstance(definition, str):
            repeat = parse_psmiles(definition)
        else:
            raise PolymerBuildError(
                f"Repeat-unit definition {identity!r} must be PSMILES or RepeatUnit"
            )
        _validate_attachment_chemistry(repeat)
        library[identity] = repeat
    return library


def _validate_dp(dp: int, *, name: str = "Degree of polymerization") -> None:
    if not isinstance(dp, int) or isinstance(dp, bool) or dp < 1:
        raise PolymerBuildError(f"{name} must be an integer >= 1")


def _validate_build_options(random_seed: int, chain_id: str, polymer_type: str) -> None:
    if not isinstance(random_seed, int) or isinstance(random_seed, bool):
        raise PolymerBuildError("random_seed must be an integer")
    if not isinstance(chain_id, str) or not chain_id:
        raise PolymerBuildError("chain_id must be a non-empty string")
    if not isinstance(polymer_type, str) or not polymer_type:
        raise PolymerBuildError("polymer_type must be a non-empty string")


def _validate_attachment_chemistry(repeat_unit: RepeatUnit) -> None:
    for point in repeat_unit.attachment_points:
        if point.bond_type != Chem.BondType.SINGLE or point.bond_order != 1.0:
            raise InvalidRepeatUnitError(
                "Phase 3.5 supports single-bond polymer attachment points only; "
                f"the {point.role} attachment has order {point.bond_order}"
            )


def _assemble_heavy_atom_graph(
    repeat_units: Sequence[RepeatUnit],
    repeat_unit_types: Sequence[str],
    chain_id: str,
    *,
    stereochemical_sequence: StereochemicalSequence | None,
) -> tuple[Chem.Mol, list[dict[int, int]], list[int]]:
    editable = Chem.RWMol()
    unit_atom_maps: list[dict[int, int]] = []
    stereo_atom_indices: list[int] = []
    for repeat_index, (repeat, repeat_type) in enumerate(
        zip(repeat_units, repeat_unit_types, strict=True)
    ):
        source = repeat.mol
        atom_map: dict[int, int] = {}
        stereo_state = (
            stereochemical_sequence.states[repeat_index]
            if stereochemical_sequence is not None
            else None
        )
        stereo_source_index = (
            require_controllable_stereocenter(repeat).atom_index
            if stereochemical_sequence is not None
            else None
        )
        for atom in source.GetAtoms():
            if atom.GetAtomicNum() == 0:
                continue
            source_index = atom.GetIdx()
            copied_atom = Chem.Atom(atom)
            copied_atom.SetAtomMapNum(0)
            copied_atom.SetProp(_CHAIN_ID_PROPERTY, chain_id)
            copied_atom.SetIntProp(_REPEAT_INDEX_PROPERTY, repeat_index)
            copied_atom.SetProp(_REPEAT_TYPE_PROPERTY, repeat_type)
            copied_atom.SetIntProp(_SOURCE_INDEX_PROPERTY, source_index)
            if stereo_state is not None:
                copied_atom.SetProp(_STEREO_STATE_PROPERTY, stereo_state)
            atom_map[source_index] = editable.AddAtom(copied_atom)
            if source_index == stereo_source_index:
                copied_atom_index = atom_map[source_index]
                editable.GetAtomWithIdx(copied_atom_index).SetBoolProp(
                    _CONTROLLED_CENTER_PROPERTY, True
                )
                stereo_atom_indices.append(copied_atom_index)
        for bond in source.GetBonds():
            begin = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if begin not in atom_map or end not in atom_map:
                continue
            editable.AddBond(atom_map[begin], atom_map[end], bond.GetBondType())
            copied_bond = editable.GetBondBetweenAtoms(atom_map[begin], atom_map[end])
            copied_bond.SetIsAromatic(bond.GetIsAromatic())
        unit_atom_maps.append(atom_map)

    for repeat_index in range(len(repeat_units) - 1):
        tail = unit_atom_maps[repeat_index][
            repeat_units[repeat_index].tail.neighbor_atom_index
        ]
        head = unit_atom_maps[repeat_index + 1][
            repeat_units[repeat_index + 1].head.neighbor_atom_index
        ]
        if editable.GetBondBetweenAtoms(tail, head) is not None:
            raise PolymerBuildError("Polymerization would create a duplicate bond")
        editable.AddBond(tail, head, Chem.BondType.SINGLE)
    return editable.GetMol(), unit_atom_maps, stereo_atom_indices


def _annotate_generated_hydrogens(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1 or atom.HasProp(_SOURCE_INDEX_PROPERTY):
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 1:
            raise PolymerBuildError(
                "A generated hydrogen must have one bonded neighbor"
            )
        parent = neighbors[0]
        atom.SetProp(_CHAIN_ID_PROPERTY, parent.GetProp(_CHAIN_ID_PROPERTY))
        atom.SetIntProp(
            _REPEAT_INDEX_PROPERTY, parent.GetIntProp(_REPEAT_INDEX_PROPERTY)
        )
        atom.SetProp(_REPEAT_TYPE_PROPERTY, parent.GetProp(_REPEAT_TYPE_PROPERTY))
        if parent.HasProp(_STEREO_STATE_PROPERTY):
            atom.SetProp(_STEREO_STATE_PROPERTY, parent.GetProp(_STEREO_STATE_PROPERTY))
        atom.SetBoolProp(_GENERATED_HYDROGEN_PROPERTY, True)


def _generate_coordinates(mol: Chem.Mol, *, generate_3d: bool, random_seed: int) -> str:
    if generate_3d:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = random_seed
        try:
            embedding_status = AllChem.EmbedMolecule(mol, parameters)
        except Exception as error:
            raise EmbeddingError("RDKit failed to embed the linear polymer") from error
        if embedding_status != 0:
            raise EmbeddingError("RDKit failed to embed the linear polymer")
        return "rdkit_etkdg_v3_initial_conformation"
    rdDepictor.Compute2DCoords(mol)
    return "rdkit_2d"


def _transfer_atom_provenance(
    system: MolecularSystem,
    mol: Chem.Mol,
    rdkit_index_to_site_id: dict[int, int],
) -> None:
    for atom in mol.GetAtoms():
        metadata: dict[str, Any] = {
            "chain_id": atom.GetProp(_CHAIN_ID_PROPERTY),
            "repeat_unit_index": atom.GetIntProp(_REPEAT_INDEX_PROPERTY),
            "repeat_unit_type": atom.GetProp(_REPEAT_TYPE_PROPERTY),
        }
        if atom.HasProp(_SOURCE_INDEX_PROPERTY):
            metadata["source_repeat_atom_index"] = atom.GetIntProp(
                _SOURCE_INDEX_PROPERTY
            )
        if atom.HasProp(_GENERATED_HYDROGEN_PROPERTY):
            metadata["generated_hydrogen"] = atom.GetBoolProp(
                _GENERATED_HYDROGEN_PROPERTY
            )
        if atom.HasProp(_STEREO_STATE_PROPERTY):
            metadata["repeat_unit_stereochemical_state"] = atom.GetProp(
                _STEREO_STATE_PROPERTY
            )
        if atom.HasProp(_CONTROLLED_CENTER_PROPERTY):
            metadata["controllable_stereocenter"] = atom.GetBoolProp(
                _CONTROLLED_CENTER_PROPERTY
            )
        site_id = rdkit_index_to_site_id[atom.GetIdx()]
        system.topology.get_site(site_id).metadata.update(metadata)


def _allocate_composition_counts(
    repeat_units: Mapping[str, str | RepeatUnit],
    fractions: Mapping[str, float],
    dp: int,
) -> dict[str, int]:
    if not repeat_units:
        raise PolymerBuildError("At least one repeat-unit definition is required")
    if set(fractions) != set(repeat_units):
        raise PolymerBuildError(
            "Fraction labels must exactly match repeat-unit definition labels"
        )
    values: dict[str, float] = {}
    for identity in repeat_units:
        value = fractions[identity]
        if not isinstance(value, Real) or isinstance(value, bool) or value < 0:
            raise PolymerBuildError("Copolymer fractions must be non-negative numbers")
        values[identity] = float(value)
    if not math.isclose(sum(values.values()), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise PolymerBuildError("Copolymer fractions must sum to 1.0")

    raw_counts = {identity: values[identity] * dp for identity in repeat_units}
    counts = {identity: math.floor(raw_counts[identity]) for identity in repeat_units}
    remaining = dp - sum(counts.values())
    ranked = sorted(
        repeat_units,
        key=lambda identity: raw_counts[identity] - counts[identity],
        reverse=True,
    )
    for identity in ranked[:remaining]:
        counts[identity] += 1
    return counts


def _default_repeat_label(index: int) -> str:
    if index < 26:
        return chr(ord("A") + index)
    return f"R{index + 1}"
