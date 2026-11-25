"""Coordinate-free preflight and evidence-qualified EMC failure localization."""

import re
from collections import Counter

from rdkit import Chem

from island.forcefields._emc_graph_input import chemical_input
from island.forcefields.emc import FIELDS, _field_sources


def atom_environments(mol, sites, selected=None):
    """Stable IDs, chemical state and radius-two fragments; never inferred parameters."""
    selected = set(sites if selected is None else selected)
    rows = []
    for atom, site in zip(mol.GetAtoms(), sites, strict=True):
        if site not in selected:
            continue
        index = atom.GetIdx()
        bonds = list(Chem.FindAtomEnvironmentOfRadiusN(mol, 2, index))
        atoms = {index}
        for i in bonds:
            b = mol.GetBondWithIdx(i)
            atoms.update((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        rows.append(
            {
                "site_id": site,
                "element": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
                "aromatic": atom.GetIsAromatic(),
                "isotope": atom.GetIsotope(),
                "neighbors": [sites[a.GetIdx()] for a in atom.GetNeighbors()],
                "fragment": Chem.MolFragmentToSmiles(
                    mol, atomsToUse=sorted(atoms), bondsToUse=bonds
                ),
            }
        )
    return rows


def preflight_assignment(graph, family, *, emc_root=None, field_dir=None):
    """Check graph/domain constraints without predicting full parameter coverage."""
    from island.forcefields.emc_graph import _field_elements

    if family not in FIELDS:
        raise ValueError(f"Unknown EMC field {family!r}")
    field = FIELDS[family]
    _, _, directory, hashes = _field_sources(field, emc_root, field_dir)
    mol, sites, _ = chemical_input(graph)
    original_count = len(sites)
    mol = Chem.AddHs(mol)
    sites += list(
        range(max(sites) + 1, max(sites) + 1 + mol.GetNumAtoms() - original_count)
    )
    elements = set(_field_elements(directory, field).values())
    missing = [
        sites[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() not in elements
    ]
    isotope = [sites[a.GetIdx()] for a in mol.GetAtoms() if a.GetIsotope()]
    blockers = []
    if missing:
        blockers.append({"code": "element_absent_from_source", "sites": missing})
    if field.endswith("-ua") and isotope:
        blockers.append({"code": "united_atom_isotope_unsupported", "sites": isotope})
    fragments = Chem.GetMolFrags(mol)
    charges = [
        sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in part) for part in fragments
    ]
    return {
        "field": field,
        "field_sha256": hashes,
        "ready_for_assignment": not blockers,
        "parameter_completeness": "unmeasured",
        "md": "unmeasured",
        "blockers": blockers,
        "elements": dict(
            sorted(Counter(a.GetSymbol() for a in mol.GetAtoms()).items())
        ),
        "atoms": mol.GetNumAtoms(),
        "bonds": mol.GetNumBonds(),
        "formal_charge": sum(charges),
        "component_charges": charges,
        "warnings": ["charged_components_require_source_supported_ion_parameters"]
        if any(charges)
        else [],
        "specified_stereocenters": sum(
            s.specified == Chem.StereoSpecified.Specified
            for s in Chem.FindPotentialStereo(mol)
        ),
        "environments": atom_environments(mol, sites),
    }


def assignment_failure(error, graph=None):
    """Preserve every reported category; sites only when explicitly mapped by EMC."""
    message = str(error)
    patterns = {
        "atom_typing": r"(?:typing|type).*?(?:fail|not found|undefined)|cannot type|no atom type|no typing|unassigned|unrecognized atom|no rule found|Missing rules",
        "charge_increment": r"increment pair|bond increment",
        "torsion": r"torsion|dihedral",
        "improper": r"improper",
        "angle": r"angle",
        "bond": r"no (?:\w+ )*bond coefficients|bond parameters",
        "charge_balance": r"charge.*(?:disagree|neutral|total|mismatch)|total.*charge",
        "geometry": r"ETKDG|embed|coordinates|3D|packing|construct|stereochemistry",
        "graph": r"closed-shell|radical|dummy|connectivity|Unresolved input|valence|invalid.*SMILES|chemically invalid",
        "mapping": r"mismatch|unmapped|changed the size|bond endpoints|component membership",
        "runtime": r"No such file|not installed|timeout|timed out|not found.*executable",
    }
    categories = [
        name
        for name, pattern in patterns.items()
        if re.search(pattern, message, re.IGNORECASE)
    ]
    source_sites = sorted(set(getattr(error, "source_sites", ())))
    result = {
        "error_type": type(error).__name__,
        "stage": getattr(error, "assignment_stage", "unreported"),
        "elapsed_seconds": getattr(error, "elapsed_seconds", None),
        "assignment_attempts": getattr(error, "assignment_attempts", []),
        "categories": categories or ["unclassified_backend_failure"],
        "diagnostic": message,
        "source_sites": source_sites,
        "backend_sites": re.findall(r"\{group, site\} = \{(\w+), (\d+)\}", message),
        "localization": "backend_reported_sites"
        if source_sites
        else "backend_did_not_identify_sites",
        "environments": [],
        "parameter_tuples": re.findall(
            r"(?:increment pair|coefficients found for)[^\n]*", message, re.IGNORECASE
        ),
    }
    if graph is not None and source_sites:
        mol, sites, _ = chemical_input(graph)
        n = len(sites)
        mol = Chem.AddHs(mol)
        sites += list(range(max(sites) + 1, max(sites) + 1 + mol.GetNumAtoms() - n))
        result["environments"] = atom_environments(mol, sites, source_sites)
    return result


def main():
    """Inspect a capped repeat before spending time on embedding or MD."""
    import argparse
    import json
    from pathlib import Path

    from island.forcefields.assignment_core import PolymerSpec, normalize

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("psmiles")
    parser.add_argument("--field", choices=sorted(FIELDS), required=True)
    parser.add_argument("--dp", type=int, default=4)
    parser.add_argument("--field-dir", type=Path)
    parser.add_argument("--head-end-group")
    parser.add_argument("--tail-end-group")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    graph = normalize(
        PolymerSpec(
            (("A", args.psmiles),),
            args.dp,
            end_groups=(args.head_end_group, args.tail_end_group),
        )
    )
    result = preflight_assignment(graph, args.field, field_dir=args.field_dir)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    if not result["ready_for_assignment"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
