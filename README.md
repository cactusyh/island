# ISLAND

**ISLAND** (Integrated Simulation Library for Atomistic and Network Dynamics) is
intended to become a multiscale polymer molecular simulation framework.

The package provides its foundational molecular data model: stable-ID sites,
chemical topology, separate coordinates, simulation boxes, graph-level reaction
transformations, and a force-field abstraction. An optional RDKit chemistry adapter
can convert molecules and create small atomistic systems from SMILES. GAFF/GAFF2/PCFF/CVFF/OPLS
parameterization, crosslinking chemistry, coarse-graining workflows, and LAMMPS
integration are intentionally not implemented yet.

```python
from island import AtomSite, Coordinates, MolecularSystem, Topology

topology = Topology()
topology.add_site(AtomSite(id=10, name="C1", mass=12.011, element="C", atomic_number=6))
topology.add_site(AtomSite(id=20, name="C2", mass=12.011, element="C", atomic_number=6))
topology.add_bond(10, 20, order=1)

coordinates = Coordinates({10: [0.0, 0.0, 0.0], 20: [1.54, 0.0, 0.0]})
system = MolecularSystem(topology=topology, coordinates=coordinates)
system.validate()
```

## Optional chemistry adapter

Install the `chemistry` extra to convert RDKit molecules or build a small atomistic
system from SMILES:

```python
from island.chemistry import from_smiles

system = from_smiles("CCO", add_hydrogens=True, generate_3d=True)
print(system.number_of_sites, system.number_of_bonds)
```

RDKit is an adapter, not ISLAND's authoritative representation. Polymer building,
force-field parameterization, LAMMPS workflows, and real crosslink chemistry remain
deferred.

## Development

Install development dependencies, then run `pytest` and optionally `ruff check .`.
