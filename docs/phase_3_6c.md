# Phase 3.6C: local-template initial polymer conformations

Phase 3.6C separates scalable initial-coordinate generation from polymer chemistry.
The final `MolecularSystem` graph, stable site IDs, sequence, atom provenance, and
stereochemical sequence remain authoritative. RDKit molecules are temporary local
geometry adapters only.

## Public builder API

```python
from island.builders import build_linear_polymer

system = build_linear_polymer(
    "[*:1]CC[*:2]",
    dp=50,
    coordinate_method="local_templates",
    template_seed=2026,
    assembly_seed=2026,
)
```

The default `coordinate_method="etkdg"` remains backward compatible. Selecting
`local_templates` together with `generate_3d=False` is an error: local-template
construction produces genuine 3D coordinates, not a relabeled 2D depiction.

## Algorithm and mappings

1. ISLAND first constructs the complete finite chemical graph exactly as in Phase
   3.6A/3.6B1. This determines terminal hydrogens, internal hydrogens, formal
   charges, stable IDs, and repeat provenance.
2. Each repeat is extracted from that final graph by stable site ID. Atoms carry
   an explicit `_island_final_site_id` property; coordinates are mapped back only
   after this property is verified.
3. Missing neighboring repeat environments are represented by temporary
   methyl-like carbon context caps. These caps define explicit head and tail frame
   directions and are never copied into the final `MolecularSystem`.
4. ETKDGv3 embeds only these short capped fragments. It is never called on the
   DP-length polymer in this method.
5. Repeat templates are assembled in sequence order. The existing self-avoidance,
   torsion sampler, retry, and rollback policies control placement. Each accepted
   repeat moves as a proper rigid transform, preserving its internal bond lengths,
   bond angles, ring geometry, and handedness.
6. Inter-repeat single-bond lengths use the sum of RDKit elemental covalent radii.
   The policy is public and replaceable; no universal 1.5 angstrom constant or
   force-field bond parameter is used.
7. For tactic polymers, stored graph chiral tags are removed after assembly and
   RDKit assigns stereochemistry from the generated 3D coordinates. Controlled
   centers must match the authoritative `stereochemical_sequence` or generation
   fails.

Terminal and internal repeats are separate template environments because their
explicit hydrogen counts differ. Template and assembly seeds are independent and
recorded under `system.metadata["polymer"]["coordinate_generation"]`.

## Supported scope and diagnostics

The implementation supports one finite linear atomistic chain, single
non-aromatic inter-repeat bonds, sequence-defined repeat types, explicit
hydrogens, ring-containing rigid repeats, and Phase 3.6A homopolymer tacticity.
Diagnostics record embedding size, template and assembly seeds, placement
attempts, rejected trials, rollbacks, minimum non-excluded distance, attachment
frame policy, bond-length policy, and that the full polymer was not embedded.

These coordinates are an initial geometric conformation, not an energy minimum or
an equilibrated structure. Geometric clash checks do not replace atom typing,
force-field parameter assignment, minimization, or molecular dynamics.

## Boundaries and deferred work

There is no force field, GAFF/GAFF2, PCFF, CVFF, OPLS, MMFF optimization,
LAMMPS, MD, packing, periodic cell, crosslinking, or multiple-chain construction
in this phase. Explicit bond-through-ring intersection testing and a physical
rotational-isomeric-state torsion model remain deferred.

Phase 4's planned `AtomTypingEngine` begins only after this coordinate layer:

```text
MolecularSystem
-> optional initial Coordinates
-> AtomTypingEngine
-> ParameterAssignmentEngine
-> ParameterizedSystem
```

Atom typing and force-field parameters must not become prerequisites for local
template generation.
