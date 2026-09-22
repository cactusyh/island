"""Chemical validation helpers for local-template coordinate generation."""

from collections.abc import Mapping

from rdkit import Chem

from island.chemistry import to_rdkit
from island.core import Coordinates, MolecularSystem
from island.exceptions import UnsupportedConformationError


def authoritative_cip_assignments(system: MolecularSystem) -> dict[int, str]:
    """Capture explicitly assigned final-graph R/S labels by stable site ID."""
    converted = to_rdkit(system)
    Chem.AssignStereochemistry(converted.mol, cleanIt=True, force=True)
    return {
        site_id: converted.mol.GetAtomWithIdx(atom_index).GetProp("_CIPCode")
        for site_id, atom_index in converted.site_id_to_rdkit_index.items()
        if converted.mol.GetAtomWithIdx(atom_index).HasProp("_CIPCode")
        and converted.mol.GetAtomWithIdx(atom_index).GetProp("_CIPCode") in {"R", "S"}
    }


def validate_explicit_hydrogens(system: MolecularSystem) -> None:
    """Reject omitted required hydrogens while allowing hydrogen-free chemistry."""
    converted = to_rdkit(system)
    missing: list[tuple[int, int]] = []
    for site_id, atom_index in converted.site_id_to_rdkit_index.items():
        atom = converted.mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() == 1:
            continue
        omitted = atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
        if omitted:
            missing.append((site_id, omitted))
    if missing:
        details = ", ".join(
            f"site {site_id}: {count}" for site_id, count in missing[:8]
        )
        raise UnsupportedConformationError(
            "Local-template generation requires chemically required hydrogens "
            f"as explicit sites; omitted hydrogen counts: {details}"
        )


def transfer_fragment_chirality(
    source: Chem.Mol,
    fragment: Chem.Mol,
    template_to_source: Mapping[int, int],
    cap_to_source: Mapping[int, int],
    template_to_site: Mapping[int, int],
    repeat_index: int,
) -> None:
    """Transfer tetrahedral parity through explicit neighbor-index mappings.

    CIP labels are intentionally not copied: temporary isotopic context caps may
    change CIP priorities. Only the source atom's neighbor-order parity is
    preserved.
    """
    neighbor_mapping = {**template_to_source, **cap_to_source}
    for template_index, source_index in template_to_source.items():
        source_atom = source.GetAtomWithIdx(source_index)
        if source_atom.GetChiralTag() not in {
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        }:
            continue
        template_atom = fragment.GetAtomWithIdx(template_index)
        source_neighbors = [
            neighbor.GetIdx() for neighbor in source_atom.GetNeighbors()
        ]
        try:
            mapped_neighbors = [
                neighbor_mapping[neighbor.GetIdx()]
                for neighbor in template_atom.GetNeighbors()
            ]
        except KeyError as error:
            site_id = template_to_site[template_index]
            raise UnsupportedConformationError(
                "Cannot map all tetrahedral neighbors for local-template site "
                f"{site_id} in repeat {repeat_index}"
            ) from error
        if len(mapped_neighbors) != len(source_neighbors) or set(
            mapped_neighbors
        ) != set(source_neighbors):
            site_id = template_to_site[template_index]
            raise UnsupportedConformationError(
                "Attachment context cannot represent tetrahedral neighborhood for "
                f"site {site_id} in repeat {repeat_index}"
            )
        permutation = [source_neighbors.index(index) for index in mapped_neighbors]
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(len(permutation))
            for right in range(left + 1, len(permutation))
        )
        template_atom.SetChiralTag(source_atom.GetChiralTag())
        if inversions % 2:
            template_atom.InvertChirality()


def validate_coordinate_stereochemistry(
    system: MolecularSystem,
    positions: Mapping[int, object],
    expected_by_site: Mapping[int, str],
) -> dict[str, object]:
    """Infer final CIP from coordinates and compare every constrained site."""
    tacticity = system.metadata.get("polymer", {}).get("stereochemical_sequence")
    if not expected_by_site and tacticity is None:
        return {
            "status": "not_applicable",
            "validated_site_ids": [],
            "reason": "no explicitly assigned tetrahedral stereocenters",
        }
    candidate = system.copy()
    candidate.coordinates = Coordinates(positions)
    converted = to_rdkit(candidate)
    Chem.RemoveStereochemistry(converted.mol)
    Chem.AssignStereochemistryFrom3D(converted.mol, replaceExistingTags=True)
    observed_by_site = {
        site_id: (
            converted.mol.GetAtomWithIdx(atom_index).GetProp("_CIPCode")
            if converted.mol.GetAtomWithIdx(atom_index).HasProp("_CIPCode")
            else None
        )
        for site_id, atom_index in converted.site_id_to_rdkit_index.items()
        if site_id in expected_by_site
    }
    mismatches = [
        (site_id, expected, observed_by_site.get(site_id))
        for site_id, expected in expected_by_site.items()
        if observed_by_site.get(site_id) != expected
    ]
    if mismatches:
        details = "; ".join(
            "site_id={site}, repeat_index={repeat}, expected={expected}, "
            "observed={observed}".format(
                site=site_id,
                repeat=system.topology.get_site(site_id).metadata.get(
                    "repeat_unit_index"
                ),
                expected=expected,
                observed=observed,
            )
            for site_id, expected, observed in mismatches
        )
        raise UnsupportedConformationError(
            "Local-template stereochemistry validation failed: " + details
        )
    if tacticity is not None:
        controlled = sorted(
            (
                site.metadata["repeat_unit_index"],
                observed_by_site.get(site_id),
            )
            for site_id, site in system.topology.sites.items()
            if site.metadata.get("controllable_stereocenter")
        )
        observed_sequence = [state for _, state in controlled]
        if observed_sequence != list(tacticity):
            raise UnsupportedConformationError(
                "Local-template tacticity validation failed: expected "
                f"{list(tacticity)}, observed {observed_sequence}"
            )
    return {
        "status": "validated",
        "validated_site_ids": sorted(expected_by_site),
        "expected_cip_by_site": dict(sorted(expected_by_site.items())),
        "observed_cip_by_site": dict(sorted(observed_by_site.items())),
    }
