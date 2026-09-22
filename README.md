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

RDKit is an adapter, not ISLAND's authoritative representation. Force-field parameterization, LAMMPS workflows, and real crosslink chemistry remain
deferred.

## Sequence-driven linear polymer building

ISLAND supports finite linear homopolymers and sequence-defined copolymers from
bifunctional PSMILES repeat units. Explicit labels `[*:1]` (head) and `[*:2]` (tail) are recommended for
asymmetric repeat units; two unlabeled dummy atoms use their original parse order.

```python
from island.builders import build_linear_polymer

system = build_linear_polymer("[*:1]CC[*:2]", dp=5, random_seed=2026)
print(system.number_of_sites)
print(system.metadata["polymer"])
```

All polymer chemistry uses an explicit ordered sequence. For example:

```python
from island.builders import build_polymer_from_sequence

system = build_polymer_from_sequence(
    repeat_units={"A": "[*:1]CC[*:2]", "B": "[*:1]CO[*:2]"},
    sequence=["A", "A", "B", "A", "B"],
)
```

Dummy atoms are removed during graph assembly and finite ends receive hydrogens
through normal RDKit valence handling. ETKDGv3 coordinates are deterministic
initial molecular conformations only; they are not production-equilibrated chains.
Convenience APIs generate alternating, multiblock, and reproducible random
copolymer sequences. In every API, DP is the total number of repeat units, not
the number of pairs or blocks. Random integer compositions use largest-remainder
allocation and a local seeded RNG. Tacticity control, custom end groups, branching, cyclic polymers, force fields, packing, LAMMPS, crosslinking,
and production equilibration remain unimplemented.

## Chemical and parameter identity

`AtomSite.formal_charge` stores integer chemical formal charge. Future force-field
atom types, partial charges, and nonbonded parameters belong to
`ParameterizedSystem.site_assignments`, not to `Site` or `AtomSite`.
`BeadSite.bead_type` is retained because it identifies a coarse-grained
representation site rather than an atomistic force-field assignment.

## Development

Install development dependencies, then run `pytest` and optionally `ruff check .`.
