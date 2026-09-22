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
allocation and a local seeded RNG. Custom end groups, branching, cyclic polymers, force fields, packing, LAMMPS, crosslinking,
and production equilibration remain unimplemented.

## Stereochemistry and tacticity

Phase 3.6A supports homopolymer repeat units containing exactly one explicitly
assigned, backbone-relevant tetrahedral stereocenter. The input PSMILES CIP state
is the reference state. Isotactic chains repeat that state; syndiotactic chains
alternate it with the opposite state; atactic chains use a locally seeded shuffle
of a requested opposite-state fraction. The default `atactic_fraction=0.5` uses
nearest-integer allocation via `floor(dp * fraction + 0.5)`.

```python
from island.builders import build_linear_polymer

system = build_linear_polymer(
    "[*:1]N[C@H](F)C[*:2]",
    dp=6,
    tacticity="syndiotactic",
    stereo_seed=2026,
)
```

ISLAND assigns and verifies R/S states on the complete polymer graph, after dummy
replacement, so final substituent priorities are respected. Graph chiral tags are
inverted chemically; coordinates are not mirrored. Achiral, unassigned,
multicenter, and mixed-repeat tacticity requests are rejected in this phase.
Stereochemical sequence generation is independent of monomer sequence generation
and ETKDG coordinate generation.

The per-atom provenance field `repeat_unit_stereochemical_state` records the
assignment of the containing repeat unit. Only an atom marked
`controllable_stereocenter=True` carries that center's intrinsic `cip_label` and
`chiral_tag`.

## Self-avoiding random-walk coordinates

Phase 3.6B adds non-destructive, coordinate-only generation for one finite linear
atomistic chain:

```python
from island.conformations import generate_polymer_conformation

result = generate_polymer_conformation(system, method="random_walk", seed=2026)
conformed_system = result.apply_to(system)  # returns a copy by default
```

The generator treats each repeat unit as a rigid local geometry and places units
in repeat-index order. It samples torsions with a local seeded RNG, rejects
geometric clashes, and can roll back previously placed units when retries are
exhausted. The default fixed-distance criterion excludes bonded 1-2 and 1-3 pairs;
1-4 pairs are checked. An optional elemental van-der-Waals-radius policy is also
available and does not use force-field parameters.

Random-walk output is an initial conformation, not an equilibrated structure.
Steric criteria are geometric rather than energetic. Force-field minimization, MD
relaxation, multiple-chain packing, and explicit bond-through-ring intersection
checking remain future work. Ring-containing repeat units are moved rigidly so
their internal geometry is preserved. Chemical topology, stable IDs, provenance,
and tacticity are not changed by conformation generation.

## Chemical and parameter identity

`AtomSite.formal_charge` stores integer chemical formal charge. Future force-field
atom types, partial charges, and nonbonded parameters belong to
`ParameterizedSystem.site_assignments`, not to `Site` or `AtomSite`.
`BeadSite.bead_type` is retained because it identifies a coarse-grained
representation site rather than an atomistic force-field assignment.

## Development

Install development dependencies, then run `pytest` and optionally `ruff check .`.
