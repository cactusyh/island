# Phase 3.6B1: geometric correctness before force-field integration

## Reviewed baseline

Repository commit: `a7e9ec9ebf3b49dab2c1d8e9b5f4223a06c61518`
(`Phase 3.6B — Polymer Conformation Engine`). All 106 baseline tests passed.

Independent coordinate checks found two scientific correctness gaps:

1. The conformation tests and example constructed source coordinates with
   `generate_3d=False`, which calls RDKit's 2D depiction generator. For
   `[*:1]N[C@H](F)C[*:2]`, DP=6, syndiotactic, the four-neighbor tetrahedral
   determinants were zero before and after rigid placement. After removing
   graph stereochemistry and assigning it from the conformer, none of the six
   intended stereocenters had a coordinate-derived CIP label. The old test
   restored stored chiral tags through `to_rdkit`, so it could still pass.
2. Source units were positioned using their incoming-to-outgoing chord and a
   sampled cone direction, not the full attachment geometry. For genuinely 3D
   `[*:1]C[C@H](F)[*:2]` and `[*:1][C@H](F)C[*:2]` chains, DP=6, seeds 31 and
   2026, coordinate-derived CIP labels changed at attachment centers. Preserving
   distances inside a unit did not preserve the external substituent direction.

## Correction

The source is now required to contain valid local 3D geometry. Known 2D depiction
inputs, zero-length bonds, non-finite coordinates, and degenerate tetrahedral
neighborhoods are rejected with actionable errors. Aromatic planar units remain
supported. Periodic input and non-single inter-repeat bonds are rejected because
the current engine does not support their geometry or torsional constraints.

Let `R_prev` be the accepted rigid rotation of the previous unit. The new anchor
is placed using `R_prev` applied to the **original inter-repeat bond vector**.
The new unit's rotation is `R_bond(delta) @ R_prev`, where `R_bond` is a proper
rotation about that transformed bond. Both endpoint neighborhoods therefore
retain their source bond lengths, bond angles, and handedness. Rejected proposals
do not mutate accepted geometry; rollback discards the corresponding rotations.

This remains a coordinate-only operation. It uses no force-field energy, charges,
atom typing, reflection, topology edits, or chemical stereo reassignment.

Additional corrections:

- Steric exclusions are computed from graph distances in the bond graph, so stale
  angle/dihedral caches do not change 1-3 or 1-4 behavior.
- Built-in distance policies reject non-finite thresholds. Non-finite sampled
  rotations are rejected explicitly.
- Diagnostics identify angstrom units and uniform-site-weighted radius of gyration.
- `accepted_rotation_increments_degrees` replaces the misleading
  `accepted_torsions_degrees` key. These are not measured absolute dihedrals.
- `bond_angle_degrees` and `direction_jitter_degrees` default to `None`;
  explicitly requesting angle resampling raises `UnsupportedConformationError`.
  Changing source geometry requires a separate, chemically informed operation.

## Validation contract

Tests check coordinate-derived stereochemistry after clearing stored tags, not
only metadata round trips. They cover internal and attachment-centered chirality,
all graph-derived local bond angles, every chemical bond length, rigid source
transformations, 2D rejection, false 3D metadata, aromatic planarity, and policy
validation. Existing deterministic seeds, non-destructive behavior, ring rigidity,
retry exhaustion, successful rollback, and DP=50 walker tests remain covered.
Core and conformation imports are tested with RDKit blocked.

Validation on Python 3.12 / RDKit 2026.03.6: **128 tests passed** (including
33 conformation tests); Ruff and `python -m pip check` passed. The updated PE
DP=20 example completed with 122 sites, 19 placement attempts, no rejected trials
or rollbacks, and a minimum non-excluded distance of approximately 1.08457 A.

The geometric degeneracy guard is not a full absolute-CIP or energetic validator.
The engine preserves source stereochemistry; callers must supply a chemically
correct source. Internal rigid-unit clashes are not repaired. Rotations about
chemically restricted single bonds are not sampled from a physical torsion model.
Explicit ring-threading checks remain deferred.

## Long-chain boundary

With RDKit 2026.03.6, the default builder ETKDG setup failed for PE DP=50, seed
2026, in the review environment. ETKDG with `useRandomCoords=True` succeeded and
is used explicitly to prepare the genuine 3D source in the DP=50 walker test.
This validates the walker with supplied geometry, not a scalable local-template
builder. The production builder's embedding behavior has not been changed.

A future local-template stage should create chemically validated 3D repeat
fragments with explicit attachment frames, preserve terminal environments and
final-graph stereochemistry, and avoid whole-chain embedding. It must keep the
coordinate layer separate from the authoritative chemical graph and force fields.

## Next planned stage: Phase 4

The project direction remains:

`MolecularSystem -> AtomTypingEngine -> ParameterAssignmentEngine -> ParameterizedSystem`

Phase 4 should extend the existing `ForceFieldBackend` and `ParameterizedSystem`
instead of introducing competing molecular representations. Its minimum contract:

- Atom typing returns assignments keyed by stable site ID, with rule/backend
  identity and traceable matches, unresolved sites, and ambiguity diagnostics.
- Chemical formal charges stay on `AtomSite`; force-field types, partial charges,
  and parameters stay in assignment objects.
- Parameter lookup is separate from atom typing. Missing bonded, nonbonded, or
  required cross terms must be reported explicitly; no silent zero/default fill.
- Schemas record functional forms, units, force-field version/provenance, mixing
  rules, and 1-4 scaling so later GAFF/GAFF2 and class-II backends are not conflated.
- Assignment validity is tied to the chemical graph. Reaction edits invalidate
  affected assignments; local retyping can be added later without coupling the
  reaction engine to a particular force field.
- Demonstration rules and parameters are explicitly test-only. Phase 4 does not
  claim production GAFF/GAFF2/PCFF/CVFF/OPLS support or simulation readiness.

Reference implementations to inspect during Phase 4 are Foyer for typing rules
and precedence, RadonPy/pysimm for backend organization, and HTPolyNet/Polymatic
for future reaction/retyping boundaries. These are planned design references,
not dependencies introduced or implementations copied in this change.
