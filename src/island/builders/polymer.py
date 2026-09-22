"""Construction of finite, linear homopolymers from PSMILES."""

from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

from island.chemistry.psmiles import RepeatUnit, parse_psmiles
from island.chemistry.rdkit_adapter import from_rdkit
from island.core import MolecularSystem
from island.exceptions import EmbeddingError, InvalidRepeatUnitError, PolymerBuildError

_CHAIN_ID_PROPERTY = "_island_chain_id"
_REPEAT_INDEX_PROPERTY = "_island_repeat_unit_index"
_SOURCE_INDEX_PROPERTY = "_island_source_repeat_atom_index"
_GENERATED_HYDROGEN_PROPERTY = "_island_generated_hydrogen"


def build_linear_polymer(
    psmiles: str,
    *,
    dp: int,
    generate_3d: bool = True,
    random_seed: int = 2026,
    add_hydrogens: bool = True,
    chain_id: str = "A",
) -> MolecularSystem:
    """Build a finite, head-to-tail linear homopolymer.

    Hydrogens are added after the complete heavy-atom graph is sanitized. When
    ``generate_3d`` is false, deterministic 2D coordinates are created because
    every :class:`MolecularSystem` requires coordinates for all sites.
    """
    _validate_options(dp, random_seed, chain_id)
    repeat_unit = parse_psmiles(psmiles)
    _validate_attachment_chemistry(repeat_unit)

    polymer, unit_atom_maps = _assemble_heavy_atom_graph(repeat_unit, dp, chain_id)
    try:
        Chem.SanitizeMol(polymer)
    except Exception as error:
        raise PolymerBuildError(
            "The assembled head-to-tail polymer graph is chemically invalid"
        ) from error

    head_atom_index = unit_atom_maps[0][repeat_unit.head.neighbor_atom_index]
    tail_atom_index = unit_atom_maps[-1][repeat_unit.tail.neighbor_atom_index]
    if add_hydrogens:
        polymer = Chem.AddHs(polymer)
        _annotate_generated_hydrogens(polymer)

    if generate_3d:
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = random_seed
        try:
            embedding_status = AllChem.EmbedMolecule(polymer, parameters)
        except Exception as error:
            raise EmbeddingError(
                f"RDKit failed to embed linear polymer with DP={dp}"
            ) from error
        if embedding_status != 0:
            raise EmbeddingError(f"RDKit failed to embed linear polymer with DP={dp}")
        coordinate_kind = "rdkit_etkdg_v3_initial_conformation"
    else:
        rdDepictor.Compute2DCoords(polymer)
        coordinate_kind = "rdkit_2d"

    conversion = from_rdkit(polymer)
    system = conversion.system
    _transfer_atom_provenance(system, polymer, conversion.rdkit_index_to_site_id)
    system.metadata["polymer"] = {
        "architecture": "linear",
        "polymer_type": "homopolymer",
        "source_psmiles": psmiles,
        "degree_of_polymerization": dp,
        "number_of_repeat_units": dp,
        "number_of_inter_repeat_unit_bonds": dp - 1,
        "chain_id": chain_id,
        "head_site_id": conversion.rdkit_index_to_site_id[head_atom_index],
        "tail_site_id": conversion.rdkit_index_to_site_id[tail_atom_index],
        "coordinates": coordinate_kind,
    }
    return system


def _validate_options(dp: int, random_seed: int, chain_id: str) -> None:
    if not isinstance(dp, int) or isinstance(dp, bool) or dp < 1:
        raise PolymerBuildError("Degree of polymerization must be an integer >= 1")
    if not isinstance(random_seed, int) or isinstance(random_seed, bool):
        raise PolymerBuildError("random_seed must be an integer")
    if not isinstance(chain_id, str) or not chain_id:
        raise PolymerBuildError("chain_id must be a non-empty string")


def _validate_attachment_chemistry(repeat_unit: RepeatUnit) -> None:
    for point in repeat_unit.attachment_points:
        if point.bond_type != Chem.BondType.SINGLE or point.bond_order != 1.0:
            raise InvalidRepeatUnitError(
                "Phase 3 supports single-bond polymer attachment points only; "
                f"the {point.role} attachment has order {point.bond_order}"
            )


def _assemble_heavy_atom_graph(
    repeat_unit: RepeatUnit, dp: int, chain_id: str
) -> tuple[Chem.Mol, list[dict[int, int]]]:
    editable = Chem.RWMol()
    source = repeat_unit.mol
    real_atom_indices = [
        atom.GetIdx() for atom in source.GetAtoms() if atom.GetAtomicNum() != 0
    ]
    unit_atom_maps: list[dict[int, int]] = []

    for repeat_index in range(dp):
        atom_map: dict[int, int] = {}
        for source_index in real_atom_indices:
            atom = Chem.Atom(source.GetAtomWithIdx(source_index))
            atom.SetAtomMapNum(0)
            atom.SetProp(_CHAIN_ID_PROPERTY, chain_id)
            atom.SetIntProp(_REPEAT_INDEX_PROPERTY, repeat_index)
            atom.SetIntProp(_SOURCE_INDEX_PROPERTY, source_index)
            atom_map[source_index] = editable.AddAtom(atom)
        for bond in source.GetBonds():
            begin = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if begin not in atom_map or end not in atom_map:
                continue
            editable.AddBond(atom_map[begin], atom_map[end], bond.GetBondType())
            copied = editable.GetBondBetweenAtoms(atom_map[begin], atom_map[end])
            copied.SetIsAromatic(bond.GetIsAromatic())
        unit_atom_maps.append(atom_map)

    for repeat_index in range(dp - 1):
        tail = unit_atom_maps[repeat_index][repeat_unit.tail.neighbor_atom_index]
        head = unit_atom_maps[repeat_index + 1][repeat_unit.head.neighbor_atom_index]
        if editable.GetBondBetweenAtoms(tail, head) is not None:
            raise PolymerBuildError("Polymerization would create a duplicate bond")
        editable.AddBond(tail, head, Chem.BondType.SINGLE)

    return editable.GetMol(), unit_atom_maps


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
        atom.SetBoolProp(_GENERATED_HYDROGEN_PROPERTY, True)


def _transfer_atom_provenance(
    system: MolecularSystem,
    mol: Chem.Mol,
    rdkit_index_to_site_id: dict[int, int],
) -> None:
    for atom in mol.GetAtoms():
        metadata: dict[str, Any] = {
            "chain_id": atom.GetProp(_CHAIN_ID_PROPERTY),
            "repeat_unit_index": atom.GetIntProp(_REPEAT_INDEX_PROPERTY),
        }
        if atom.HasProp(_SOURCE_INDEX_PROPERTY):
            metadata["source_repeat_atom_index"] = atom.GetIntProp(
                _SOURCE_INDEX_PROPERTY
            )
        if atom.HasProp(_GENERATED_HYDROGEN_PROPERTY):
            metadata["generated_hydrogen"] = atom.GetBoolProp(
                _GENERATED_HYDROGEN_PROPERTY
            )
        site_id = rdkit_index_to_site_id[atom.GetIdx()]
        system.topology.get_site(site_id).metadata.update(metadata)
