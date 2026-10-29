"""Versioned extension graph: chemical validation without a fixed element list."""

from copy import deepcopy
from math import isfinite

from island import AtomSite
from island.forcefields.external_pdms import digest, require


def chemical_metadata(metadata):
    """Keep chemical construction identity while allowing documented coordinate updates."""
    result = deepcopy(metadata)
    for key in (
        "coordinate_source",
        "coordinate_generation",
        "coordinate_history",
        "local_minimization",
        "dynamics",
        "production_validated",
        "simulation_readiness",
    ):
        result.pop(key, None)
    if isinstance(result.get("polymer"), dict):
        result["polymer"].pop("coordinates", None)
        result["polymer"].pop("coordinate_generation", None)
    return result


def extension_graph(system):
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    system.validate()
    require(
        system.metadata.get("production_validated", False) is False
        and system.metadata.get("simulation_readiness", "not_established")
        == "not_established",
        "promoted_input_readiness",
    )
    require(
        system.representation == "atomistic" and system.box is None,
        "finite_nonperiodic_atomistic_graph_required",
    )
    require(bool(system.topology.sites), "empty_graph")
    table = Chem.GetPeriodicTable()
    for i, atom in system.topology.sites.items():
        require(type(i) is int and type(atom) is AtomSite and atom.id == i, "atom_id")
        require(
            type(atom.atomic_number) is int
            and 1 <= atom.atomic_number <= 118
            and table.GetElementSymbol(atom.atomic_number) == atom.element,
            "element_atomic_number_mismatch:" + str(i),
        )
        require(
            type(atom.mass) in (float, int) and isfinite(atom.mass) and atom.mass > 0,
            "invalid_mass:" + str(i),
        )
        require(
            type(atom.formal_charge) is int
            and atom.formal_charge == 0
            and not atom.metadata.get("radical_electrons", 0)
            and not atom.metadata.get("isotope", 0),
            "v1_requires_neutral_closed_shell_nonisotopic_sites:" + str(i),
        )
    converted = system_to_rdkit_graph(system)
    require(
        all(
            a.GetNumImplicitHs() == 0
            and a.GetNumExplicitHs() == 0
            and a.GetNumRadicalElectrons() == 0
            for a in converted.mol.GetAtoms()
        ),
        "explicit_hydrogen_completion_required",
    )
    components = [sorted(c) for c in system.topology.connected_components()]
    components.sort(key=lambda c: c[0])
    payload = {
        "schema": "island_external_opls_final_graph_v1",
        "atom_order": list(system.topology.sites),
        "sites": [
            {
                "id": i,
                "element": a.element,
                "atomic_number": a.atomic_number,
                "mass": a.mass,
                "formal_charge": a.formal_charge,
                "metadata": deepcopy(a.metadata),
            }
            for i, a in system.topology.sites.items()
        ],
        "bonds": [
            {"sites": list(k), "order": b.order, "aromatic": b.aromatic}
            for k, b in sorted(system.topology.bonds.items())
        ],
        "impropers": [list(p) for p in system.topology.impropers],
        "components": components,
        "molecule_membership": {
            str(i): str(system.topology.sites[i].metadata.get("chain_id", n))
            for n, c in enumerate(components)
            for i in c
        },
        "construction_provenance": chemical_metadata(system.metadata),
    }
    return {**payload, "identity": digest(payload)}
