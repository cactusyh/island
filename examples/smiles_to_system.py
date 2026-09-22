"""Create an ISLAND MolecularSystem from a SMILES string."""

from island.chemistry import from_smiles

system = from_smiles("CCO", add_hydrogens=True, generate_3d=True)
print(f"sites: {system.number_of_sites}")
print(f"bonds: {system.number_of_bonds}")
for site_id, site in system.topology.sites.items():
    print(site_id, site.element, system.coordinates.get(site_id))
