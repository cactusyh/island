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

The generator requires valid **three-dimensional source geometry**. Build with
`generate_3d=True` (the builder default). Phase 3.6B1 rejects RDKit 2D depictions,
degenerate tetrahedral neighborhoods, and periodic systems. Rotating a planar
repeat unit cannot create valid tetrahedral geometry.

Each repeat unit is moved rigidly in repeat-index order. The full source frame
sets the outgoing connection direction, and the next unit rotates about that
existing inter-repeat bond. This preserves bond lengths, all local bond angles,
and the source's geometric handedness, including attachment-centered chirality.
The sampler returns a rotation **increment** relative to the source frame, not a
measured absolute four-atom dihedral. Bond-angle resampling options are rejected.
Retries and rollback retain the same local seeded RNG and bounded failure behavior.

The default fixed-distance criterion excludes bonded 1-2 and 1-3 pairs, derived
directly from chemical bonds; 1-4 pairs are checked. Incremental steric checks
compare newly placed units against accepted units; internal rigid-unit distances
are retained, not repaired. An optional elemental van-der-Waals-radius policy
is independent of force-field parameters. The reported radius of gyration uses
uniform site weights, not mass weights; distances are in angstroms.

Random-walk output is an initial conformation, not an equilibrated structure.
Input geometry and its absolute chemical stereochemistry must already be valid;
the geometric guards do not constitute complete chemical or energetic validation.
The tests independently remove stored chiral tags and reassign CIP from coordinates.
Chemical topology, stable IDs, provenance, and tacticity metadata are unchanged.

Whole-chain ETKDG may fail for longer chains. The DP=50 walker smoke test supplies
an explicit 3D ETKDG conformer using random-coordinate initialization; it is **not**
a guarantee of end-to-end long-chain building. Scalable local 3D templates remain
future work, along with force-field minimization, MD relaxation, multiple-chain
packing, and explicit bond-through-ring intersection checking.

See [Phase 3.6B1 review and next-stage contract](docs/phase_3_6b1.md).

## Local-template 3D polymer construction

For long chains, explicitly select short-fragment ETKDG templates followed by
incremental self-avoiding assembly:

```python
system = build_linear_polymer(
    "[*:1]CC[*:2]",
    dp=50,
    coordinate_method="local_templates",
    template_seed=2026,
    assembly_seed=2026,
)
```

This method never embeds the full polymer. Stable IDs and the completed chemical
graph remain authoritative; explicit methyl-like context caps provide verifiable
head/tail frames for short local templates and are discarded after coordinate
transfer. Inter-repeat bond lengths use elemental covalent radii, and tacticity is
checked from the assembled 3D coordinates after stored chiral tags are removed.
The result is an initial conformation only. Force-field minimization and MD remain
necessary future steps for energetically meaningful structures. See
`docs/phase_3_6c.md` for the algorithm, diagnostics, supported scope, and the
boundary with the planned Phase 4 `AtomTypingEngine`. Phase 3.6C1 additionally
validates every explicitly assigned final-graph stereocenter by stable site ID,
keeps coordinate provenance synchronized through `ConformationResult.apply_to()`,
and enforces chemically required explicit hydrogens in both builder and direct
generator entry points; see `docs/phase_3_6c1.md`.

## Graph-based atom typing

Phase 4A provides deterministic, coordinate-free SMARTS typing with external
stable-site assignments:

```python
from island.forcefields import (
    RDKitSmartsAtomTypingEngine,
    island_demo_v1_ruleset,
)

result = RDKitSmartsAtomTypingEngine().type_system(
    system, island_demo_v1_ruleset()
)
print(result.complete)
```

Rules use exactly one SMARTS target marked `:1` and explicit transitive override
relationships; declaration order is never precedence. Results include complete
per-site match and resolution diagnostics plus deterministic chemistry/ruleset
signatures for reuse checks. Coordinates do not participate in typing. The
`island_demo_v1` `demo_*` labels are illustrative only: they are not a production
force field and provide no charges or numerical parameters. See
[Phase 4A](docs/phase_4a.md) and `examples/type_polymer.py`.

## Typed numerical parameter assignment

Phase 4B adds deterministic exact-type parameter assignment after compatible,
complete atom typing:

```python
from island.forcefields import (
    ParameterAssignmentEngine,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)

parameters = ParameterAssignmentEngine().assign(
    system,
    typing_result,
    island_demo_v1_ruleset(),
    island_demo_parameters_v1(),
)
print(parameters.coverage)
print(parameters.charges_status)          # "unassigned"
print(parameters.simulation_readiness)    # "not_established"
```

Supported forms are per-type LJ 12-6, harmonic bonds, harmonic angles, and
multi-term proper periodic torsions. Matching uses exact type tuples with complete
reversal symmetry; missing and conflicting records produce stable-site diagnostics.
Inventories are regenerated from bonds rather than trusted topology caches, and
coordinates do not affect assignment compatibility.

`island_demo_parameters_v1` contains **SYNTHETIC SOFTWARE-TEST PARAMETERS — NOT
FOR SCIENTIFIC SIMULATION** and covers only the demonstrated explicit-hydrogen PE
subset. Complete coverage does not assign charges, establish production validity,
or make a system MD-ready. See [Phase 4B](docs/phase_4b.md) and
`examples/assign_demo_parameters.py`.

Phase 4B1 additionally separates input compatibility from parameter-result
integrity, validates complete selected record content before snapshot creation,
supports safe reconstruction/deep copying, and permits explicit zero LJ epsilon and
torsion amplitudes without treating missing records as zero. Coordinates remain in
angstroms while parameter lengths use nm; future evaluators must convert units
explicitly. See [Phase 4B1](docs/phase_4b1.md).

## Charges and nonbonded policy

Phase 4C adds independent graph-based partial-charge results and explicit LJ/scaling
policies. Callers may provide exact stable-site charges or use a versioned exact
atom-type table. Each connected component is checked against authoritative formal
charge, so disconnected errors cannot cancel.

`NonbondedPolicy` explicitly selects Lorentz–Berthelot or geometric LJ 12-6 mixing
and independent LJ/Coulomb weights for shortest-path 1–2, 1–3, and 1–4 pairs.
Composition through `ParameterizedSystem.from_components()` validates the unchanged
Phase 4B result, charge result, and policy before making an owned snapshot.

All included charges and settings remain **SYNTHETIC SOFTWARE-TEST DATA — NOT FOR
SCIENTIFIC SIMULATION**. Production validation is absent and simulation readiness
is not established. Coordinates remain angstrom while parameter lengths remain nm.
See [Phase 4C](docs/phase_4c.md) and `examples/compose_synthetic_forcefield.py`.

Phase 4C1 validates charge-result structure before arithmetic or snapshot creation
and applies one shared structural check to precomputed typing results used by both
parameter and charge assignment. See [Phase 4C1](docs/phase_4c1.md).

## Resolved Amber topology import

Phase 4D1 adds `import_amber_prmtop(system, path, source_atom_to_site,
source=...)` for a restricted, already-parameterized fixed-charge Amber prmtop.
Install the optional parser with `pip install '.[amber]'`. Source atom indices are
zero-based and must map bijectively to existing stable site IDs. The importer
checks source elements, bonds, resolved terms, exclusions, 1–4 scaling, and LJ
coefficients before returning a signed `ImportedAmberResult`.
`ParameterizedSystem.from_amber_import(system, result)` creates an owned snapshot.
No SMARTS rules or automatic GAFF/GAFF2 parameterization are implied.

See [Phase 4D1](docs/phase_4d1.md) and `examples/import_amber_topology.py`.
The example topology is **synthetic**, and a successful import is neither
scientific validation nor simulation readiness.
Phase 4D1.1 rejects active Amber 12-6-4 terms, represents audited zero-LJ sites,
and normalizes snapshot family keys. A pinned external phenol topology tests
conversion, but its exact GAFF/GAFF2 provenance is not established. See
[Phase 4D1.1](docs/phase_4d1_1.md).

## Chemical and parameter identity

`AtomSite.formal_charge` stores integer chemical formal charge. Force-field atom
types and numerical assignments belong to `ParameterizedSystem`, not to `Site` or
`AtomSite`; future partial charges must remain in that parameter layer as well.
`BeadSite.bead_type` is retained because it identifies a coarse-grained
representation site rather than an atomistic force-field assignment.

## Development

Install development dependencies, then run `pytest` and optionally `ruff check .`.
