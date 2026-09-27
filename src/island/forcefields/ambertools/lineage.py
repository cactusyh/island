"""Explicit stable-ID lineage and chemistry checks for MOL2 exchange."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path

from island.core import AtomSite, MolecularSystem
from island.exceptions import AmberToolsInputError

ELEMENTS = {"H", "C", "N", "O", "F", "Cl", "Br", "S"}
ATOMIC_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9,
                  "S": 16, "Cl": 17, "Br": 35}


@dataclass(frozen=True)
class PreparedMolecule:
    names: dict[str, int]
    mol2_text: str
    expected_cip: dict[int, str]
    formal_charge: int


def _base36(number: int) -> str:
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    result = ""
    for _ in range(3):
        number, remainder = divmod(number, 36)
        result = digits[remainder] + result
    if number:
        raise AmberToolsInputError("At most 46,656 uniquely named atoms are supported")
    return result


def _mol2_type(system: MolecularSystem, site_id: int) -> str:
    site = system.topology.sites[site_id]
    element = site.element
    if element in {"H", "F", "Cl", "Br"}:
        return element
    if site.metadata.get("aromatic", False):
        if element not in {"C", "N"}:
            raise AmberToolsInputError(
                f"Aromatic atom {site_id} of element {element} lacks a supported MOL2 type"
            )
        return f"{element}.ar"
    orders = [system.topology.bonds[tuple(sorted((site_id, other)))].order
              for other in system.topology.neighbors(site_id)]
    if element == "O":
        return "O.2" if 2 in orders else "O.3"
    if element == "C":
        return "C.2" if 2 in orders else "C.1" if 3 in orders else "C.3"
    if element == "N":
        return "N.4" if site.formal_charge > 0 and len(orders) == 4 else (
            "N.2" if 2 in orders else "N.3"
        )
    if element == "S":
        return "S.3"
    raise AmberToolsInputError(f"Unsupported element {element!r} at site {site_id}")


def _coordinate_cip(
    system: MolecularSystem, positions: dict[int, tuple[float, float, float]],
    expected: dict[int, str], *, stage: str,
) -> None:
    """Infer all explicitly assigned tetrahedral labels from independent 3D data."""
    if not expected:
        return
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    converted = system_to_rdkit_graph(system)
    mol = Chem.Mol(converted.mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for site_id, index in converted.site_id_to_rdkit_index.items():
        conf.SetAtomPosition(index, positions[site_id])
    mol.AddConformer(conf, assignId=True)
    Chem.RemoveStereochemistry(mol)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    for site_id, label in expected.items():
        atom = mol.GetAtomWithIdx(converted.site_id_to_rdkit_index[site_id])
        actual = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if actual != label:
            raise AmberToolsInputError(
                f"{stage}: site {site_id} stereochemistry expected {label}, "
                f"observed {actual}; no external parameterization was accepted"
            )


def prepare_input(
    system: MolecularSystem, charges: dict[int, float] | None,
    *, max_atoms: int = 100,
) -> PreparedMolecule:
    """Validate one finite, closed-shell, explicit-H molecule and serialize MOL2."""
    if not isinstance(system, MolecularSystem) or system.representation != "atomistic":
        raise AmberToolsInputError("AmberTools requires an atomistic MolecularSystem")
    try:
        system.validate()
    except Exception as error:
        raise AmberToolsInputError(
            f"Authoritative graph/coordinates failed validation: {error}"
        ) from error
    topology = system.topology
    if len(topology.connected_components()) != 1:
        raise AmberToolsInputError("Only one connected finite molecule is supported")
    if not topology.sites or len(topology.sites) > max_atoms:
        raise AmberToolsInputError(f"Supported molecule size is 1..{max_atoms} atoms")
    if system.box is not None and any(system.box.periodic):
        raise AmberToolsInputError("Periodic boxes are outside this backend's scope")
    coordinate_source = str(system.metadata.get("coordinate_source", ""))
    if "2d" in coordinate_source.lower():
        raise AmberToolsInputError(
            "Initial coordinates are labeled 2D; supply genuine finite 3D coordinates"
        )
    for site_id, site in topology.sites.items():
        if not isinstance(site, AtomSite) or site.element not in ELEMENTS:
            raise AmberToolsInputError(
                f"Site {site_id} is not a supported H/C/N/O/F/Cl/Br/S atom"
            )
        if site.atomic_number != ATOMIC_NUMBERS[site.element]:
            raise AmberToolsInputError(
                f"Site {site_id} requires an atomic_number matching {site.element}"
            )
    total_electrons = sum(site.atomic_number for site in topology.sites.values()) - sum(
        site.formal_charge for site in topology.sites.values()
    )
    if total_electrons % 2:
        raise AmberToolsInputError(
            "Only even-electron singlet molecules are supported; an odd-electron "
            "formal state requires an explicit spin model"
        )
    try:
        from rdkit import Chem

        from island.chemistry.rdkit_graph import system_to_rdkit_graph

        converted = system_to_rdkit_graph(system)
    except ImportError as error:
        raise AmberToolsInputError(
            "RDKit chemistry validation is required; install island[chemistry]"
        ) from error
    except Exception as error:
        raise AmberToolsInputError(f"Chemical valence/graph validation failed: {error}") from error
    expected_cip: dict[int, str] = {}
    Chem.AssignStereochemistry(converted.mol, cleanIt=True, force=True)
    for site_id, index in converted.site_id_to_rdkit_index.items():
        atom = converted.mol.GetAtomWithIdx(index)
        if atom.GetNumRadicalElectrons() != 0:
            raise AmberToolsInputError(f"Open-shell site {site_id} is unsupported")
        missing_h = atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
        if missing_h:
            raise AmberToolsInputError(
                f"Site {site_id} requires {missing_h} hydrogen(s) as actual sites"
            )
        stored = topology.sites[site_id].metadata.get("cip_label")
        graph_cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if stored is not None and stored not in {"R", "S"}:
            raise AmberToolsInputError(f"Unsupported CIP label at site {site_id}: {stored}")
        if stored and graph_cip and stored != graph_cip:
            raise AmberToolsInputError(
                f"Site {site_id} stored CIP {stored} conflicts with chemical graph {graph_cip}"
            )
        if stored or graph_cip:
            expected_cip[site_id] = stored or graph_cip
    positions = {
        site_id: tuple(float(value) for value in system.coordinates.get(site_id))
        for site_id in sorted(topology.sites)
    }
    from numpy import array, cross, dot

    for center in topology.sites:
        neighbors = sorted(topology.neighbors(center))
        if len(neighbors) != 4:
            continue
        relative = [array(positions[site_id]) - array(positions[neighbors[3]])
                    for site_id in neighbors[:3]]
        volume6 = abs(float(dot(relative[0], cross(relative[1], relative[2]))))
        if volume6 < 1e-3:
            raise AmberToolsInputError(
                f"Site {center} has planar/degenerate tetrahedral geometry; "
                "provide genuine 3D input"
            )
    for bond in topology.bonds.values():
        from math import dist

        distance = dist(positions[bond.site1], positions[bond.site2])
        if not 0.55 <= distance <= 2.4:
            raise AmberToolsInputError(
                f"Bond {bond.key} has implausible input length {distance:.3f} Å"
            )
    _coordinate_cip(system, positions, expected_cip, stage="input coordinates")
    names = {f"A{_base36(index)}": site_id
             for index, site_id in enumerate(sorted(topology.sites))}
    indexes = {site_id: index + 1 for index, site_id in enumerate(sorted(topology.sites))}
    lines = [
        "@<TRIPOS>MOLECULE", "ISLAND", str(len(indexes)) + " " + str(len(topology.bonds))
        + " 1 0 0", "SMALL", "USER_CHARGES", "",
        "@<TRIPOS>ATOM",
    ]
    for name, site_id in names.items():
        x, y, z = positions[site_id]
        charge = (charges[site_id] if charges is not None else
                  float(topology.sites[site_id].formal_charge))
        lines.append(
            f"{indexes[site_id]:>7} {name:<4} {x:>12.6f} {y:>12.6f} "
            f"{z:>12.6f} {_mol2_type(system, site_id):<6} 1 MOL {charge:.10f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond_index, bond_key in enumerate(sorted(topology.bonds), 1):
        bond = topology.bonds[bond_key]
        token = "ar" if bond.aromatic else str(int(bond.order))
        lines.append(
            f"{bond_index:>6} {indexes[bond.site1]:>6} "
            f"{indexes[bond.site2]:>6} {token}"
        )
    lines.extend(("@<TRIPOS>SUBSTRUCTURE", "1 MOL 1 GROUP 0 **** 0 ROOT", ""))
    formal_charge = sum(site.formal_charge for site in topology.sites.values())
    return PreparedMolecule(names, "\n".join(lines), expected_cip, formal_charge)


def parse_and_validate_mol2(
    path: Path, system: MolecularSystem, names: dict[str, int],
    expected_cip: dict[int, str], *, stage: str,
) -> tuple[dict[int, int], dict[int, float], dict[int, tuple[float, float, float]]]:
    """Verify exact generated names, elements, bond semantics, and assigned CIP."""
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    sections: dict[str, list[str]] = {}
    current = ""
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            current = line.removeprefix("@<TRIPOS>").strip()
            sections[current] = []
        elif current and line.strip() and not line.startswith("#"):
            sections[current].append(line)
    atom_lines = sections.get("ATOM", [])
    bond_lines = sections.get("BOND", [])
    if len(atom_lines) != len(names) or len(bond_lines) != len(system.topology.bonds):
        raise AmberToolsInputError(f"{stage}: MOL2 atom/bond inventory changed")
    output_index: dict[int, int] = {}
    positions: dict[int, tuple[float, float, float]] = {}
    charges: dict[int, float] = {}
    seen_names = set()
    for line in atom_lines:
        fields = line.split()
        if len(fields) < 9:
            raise AmberToolsInputError(f"{stage}: malformed MOL2 atom line: {line}")
        number, name = int(fields[0]), fields[1]
        if name in seen_names or name not in names or number in output_index:
            raise AmberToolsInputError(
                f"{stage}: missing, duplicated, truncated, or unknown atom name {name}"
            )
        seen_names.add(name)
        site_id = names[name]
        element = system.topology.sites[site_id].element
        token = fields[5].lower()
        inferred = ("Cl" if token.startswith("cl") else "Br" if token.startswith("br")
                    else token[0].upper())
        if inferred != element:
            raise AmberToolsInputError(
                f"{stage}: atom {name} element {inferred} differs from site {site_id} {element}"
            )
        xyz = tuple(float(value) for value in fields[2:5])
        charge = float(fields[8])
        if len(xyz) != 3 or not all(isfinite(value) for value in (*xyz, charge)):
            raise AmberToolsInputError(f"{stage}: nonfinite atom data for {name}")
        output_index[number] = site_id
        positions[site_id] = xyz
        charges[site_id] = charge
    if seen_names != set(names):
        raise AmberToolsInputError(f"{stage}: source-name lineage is incomplete")
    converted = system_to_rdkit_graph(system)
    comparison = Chem.RWMol()
    site_to_new: dict[int, int] = {}
    for site_id in sorted(system.topology.sites):
        site = system.topology.sites[site_id]
        atom = Chem.Atom(site.element)
        atom.SetFormalCharge(site.formal_charge)
        atom.SetNoImplicit(True)
        site_to_new[site_id] = comparison.AddAtom(atom)
    seen_bonds: set[tuple[int, int]] = set()
    for line in bond_lines:
        fields = line.split()
        if len(fields) < 4:
            raise AmberToolsInputError(f"{stage}: malformed MOL2 bond line: {line}")
        left, right = output_index[int(fields[1])], output_index[int(fields[2])]
        key = tuple(sorted((left, right)))
        if left == right or key in seen_bonds or key not in system.topology.bonds:
            raise AmberToolsInputError(f"{stage}: changed/duplicate bond {key}")
        seen_bonds.add(key)
        token = fields[3].lower()
        types = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE,
                 "3": Chem.BondType.TRIPLE, "ar": Chem.BondType.AROMATIC}
        if token not in types:
            raise AmberToolsInputError(f"{stage}: unsupported MOL2 bond type {token}")
        comparison.AddBond(site_to_new[left], site_to_new[right], types[token])
        if token == "ar":
            bond = comparison.GetBondBetweenAtoms(site_to_new[left], site_to_new[right])
            bond.SetIsAromatic(True)
            comparison.GetAtomWithIdx(site_to_new[left]).SetIsAromatic(True)
            comparison.GetAtomWithIdx(site_to_new[right]).SetIsAromatic(True)
    try:
        new_mol = comparison.GetMol()
        Chem.SanitizeMol(new_mol)
    except Exception as error:
        raise AmberToolsInputError(
            f"{stage}: MOL2 bond-order/aromaticity representation is invalid: {error}"
        ) from error
    for key in sorted(system.topology.bonds):
        old = converted.mol.GetBondBetweenAtoms(
            converted.site_id_to_rdkit_index[key[0]],
            converted.site_id_to_rdkit_index[key[1]],
        )
        new = new_mol.GetBondBetweenAtoms(site_to_new[key[0]], site_to_new[key[1]])
        if (old.GetIsAromatic() != new.GetIsAromatic()
            or (not old.GetIsAromatic()
                and old.GetBondTypeAsDouble() != new.GetBondTypeAsDouble())):
            raise AmberToolsInputError(
                f"{stage}: bond semantics changed for stable sites {key}"
            )
    _coordinate_cip(system, positions, expected_cip, stage=stage)
    return output_index, charges, positions
