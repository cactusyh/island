# Phase 4E2.1: result integrity and coordinate provenance

## Repository state and reproduced failures

Phase 4E2 remains unmerged. This correction continues on
`codex/phase-4e2-local-minimization`, starting at reviewed commit
`bab9ffda9cf0c0edd670bb6d456643ffa357b4a0`. Its main base remains
`a2dece53fbd7e4bb73682deae979ba8cb472592c`. Main is not modified or merged.

Two focused regressions were added and executed **before** changing production
code; both failed:

1. Removing one returned coordinate and corresponding force using
   `dataclasses.replace()`, then recomputing the reduced coordinate fingerprint,
   was accepted by `to_system()`. The resulting system failed its own validation.
2. A real DP=3 PE local-template builder followed by analytical minimization left
   polymer current-coordinate labels pointing at `local_templates_etkdg_incremental`
   while the root label said `local_minimization`. The previous source was lost.

The first boundary checked a few fingerprints/status flags without validating
site inventories and actual record contents. The second updated only root
metadata while builders used the polymer namespace; conformation application
also did not retire root minimization provenance.

## Integrity contract

`MinimizationResult.validate_integrity()` returns normally only for structurally
and numerically consistent records. Malformed results raise
`InvalidMinimizationResultError`, a subclass of `MinimizationInputError`.
Construction errors for malformed mapping/vector containers are translated to
this exception as well. Derived result properties and `deepcopy()` validate
before using their contents, avoiding raw attribute/type/reduction errors.
Valid `dataclasses.replace()` reconstruction and `deepcopy()` remain supported.

Validation checks:

- Nonempty initial and returned mappings, integer stable IDs (not booleans),
  finite three-component numeric vectors, and identical site coverage.
- Available records are `EvaluationResult` instances of the supported concrete
  type with finite energies/components/forces, explicit supported units, exact
  force coverage, and matching coordinate fingerprints. Named components sum to
  total energy within `1e-8 kJ/mol` absolute / `1e-10` relative rounding tolerance.
- Initial/final model and parameter fingerprints agree. Backend identification,
  immutable scalar settings, units, and readiness declarations are well formed.
- Options, integer budgets/counts, optimizer strings/flags, termination reason,
  returned-state label, stereo outcome, and verification flag are valid.
- History contains the initial entry and consecutively numbered accepted
  iterations, with increasing evaluator counts, finite force metrics, and a
  consistent initial record. An accepted returned frame matches the last
  accepted coordinate fingerprint. A verified final evaluation requires an
  additional evaluator call after the recorded accepted history.
- Convergence agrees with `force_converged`, requires verified final evaluation,
  and satisfies the unchanged maximum atomic force norm and energy-nonincrease
  criteria. Budget-limit reasons must match their declared limits.

The result owns coordinate mappings/tuples and history tuples. Supported nested
records contain immutable scalars and mappings; validation rejects mutable or
malformed nested content before `deepcopy()` may return the original object.
No result signature or claim of cryptographic authenticity is introduced.
Consistency validation cannot prove that caller-supplied energies or forces were
actually computed, nor reconstruct intermediate coordinates from history hashes.

### Diagnostics, verification, convergence, and application permission

These remain distinct:

- An initial evaluation failure is a valid diagnostic with no energy records,
  one attempted evaluator call, no accepted history, and unchanged coordinates.
  It cannot be applied as an evaluated geometry.
- An exhausted one-call budget retains the initial evaluated state as diagnostic
  output, with `final_evaluation_verified=False`. It can be applied only with
  explicit `allow_unconverged=True` if its contents and geometry are admissible.
- A verified final evaluation does not imply force convergence. Energy stopping,
  line-search failure, and budget exhaustion retain their actual statuses.
- `allow_unconverged=True` bypasses only the requirement for convergence. It does
  not bypass integrity, topology compatibility, or stereochemistry checks.

`to_system()` first validates the record and the target system, checks exact
coordinate coverage against the target topology and the existing system identity,
and independently checks assigned tetrahedral stereo on initial and returned
coordinates against the target graph. A stored `stereochemistry="passed"` string
is not evidence for reconstructed geometry. Inverted or degenerate assigned
centers are rejected, including diagnostic application. The shared Cartesian
validator removes stored tags before assigning CIP from geometry. Missing RDKit
still raises `MinimizationUnavailableError` when required stereo checks cannot run.
The completed copied system is validated before it is returned; neither success
nor failure mutates the caller's system.

## Current versus historical coordinate provenance

A small helper in `island.core.coordinate_provenance` is shared by minimization
and conformation application. Coordinate fingerprint serialization is unchanged.

After either operation:

- Root `coordinate_source` and `coordinate_generation` describe the current
  coordinates. If a polymer namespace exists, its `coordinates` label and
  `coordinate_generation` record are synchronized with those root fields.
- Prior coordinate labels, generation diagnostics, and any current
  `local_minimization` record are copied into a flat chronological
  `coordinate_history` list. Each entry records the actual superseded frame's
  coordinate fingerprint, root records, and polymer records where present.
- Old generation diagnostics are replaced wholesale in current fields: an old
  minimum-distance measurement or assembly rejection count is not presented as
  a measurement of minimized coordinates.
- Previous source resolution uses an existing root source first, then the
  polymer generation source / polymer coordinates label, or the non-polymer
  generation source. A system with no prior source retains `None`; no source is
  invented. Historical records themselves are preserved verbatim.
- Current minimization records include actual initial/returned coordinate
  fingerprints, model/parameter identity, previous source, termination status,
  convergence and final-verification flags. Applied diagnostics use
  `local_minimization_diagnostic`; converged output uses `local_minimization`.
- Repeated minimization appends the preceding record to history. Subsequent
  conformation application removes current `local_minimization` and preserves it
  in history, so it no longer describes the current frame.

Chemical topology, IDs, repeat metadata, assigned stereo, charge values,
parameter provenance, and historical AmberTools preparation records/signatures
are unchanged. The helper edits only coordinate-provenance fields. Original
preparation must still be validated against its original preparation system.

New metadata fields change newly generated system identities where appropriate;
historical signed preparation records are neither migrated nor re-signed. Existing
chemical-snapshot tests exclude the explicitly coordinate-related root fields,
while still comparing all chemical/repeat provenance. Dedicated tests now check
current and historical coordinate records directly.

## Regression and numerical verification

Coverage includes the exact missing-coordinate reconstruction, unknown and extra
IDs, missing forces, malformed/nonfinite vectors, wrong fingerprints/identities,
malformed nested records, false convergence/verification claims, valid failure
and one-call diagnostics, caller-data isolation, reconstructed inverted/planar
centers, missing RDKit, non-polymer provenance, repeated minimization, and later
conformation application. A real OpenMM minimization/application test verifies
that synthetic validated preparation records and their signatures remain intact,
and every successfully applied system passes `MolecularSystem.validate()`.

Verification environment: Python 3.11.16, SciPy 1.17.1, OpenMM 8.6.1 Reference,
NumPy 2.4.6, ParmEd 4.3.1, RDKit 2026.03.6, in the existing separate Phase 4E2
validation environment. No base environment changes were needed.

Verification commands and results:

- Complete ordinary `pytest -q -ra`: **540 passed, 3 skipped**. The three
  default skips are live AmberTools regeneration, archived single-point checks,
  and archived minimization checks; all require explicit opt-in.
- Archived minimization acceptance, explicitly enabled separately: **1 passed**,
  covering all five retained cases at their unchanged declared tolerances.
- Ruff: passed. `python -m pip check`: no broken requirements.
- `examples/minimize_geometry.py`, `build_local_template_polymer.py`,
  `generate_random_walk.py`, and `build_linear_polymer.py`: all passed.
- Ordinary tests include blocked-optional-dependency import checks and explicit
  missing-RDKit application checks. No required optional validation was
  unavailable; live AmberTools and the archived single-point test were left
  opt-in, as this patch requires the archived minimization test.

The retained Phase 4D2.1 archive is available and its existing archived
minimization acceptance test was rerun with the original **0.1
kJ/(mol*angstrom)** force tolerance, **500 iterations**, **2000 evaluations**,
and original source-reference comparison tolerances. All five cases passed:
GAFF/GAFF2 AM1-BCC phenol, GAFF2 provided-charge phenol, capped PE DP=3, and
assigned Cl/Br. Historical hashes and charge provenance were verified by the
existing loader. No archive or charge values were changed, tolerances loosened,
or preparation records re-signed. Live AmberTools regeneration is not required
for this patch and was not rerun.

The minimization example and relevant builder/conformation examples were run.
The minimization example retained energy **67.827126 → 0 kJ/mol**, fmax
**336.920582 → 0 kJ/(mol*angstrom)**, one accepted iteration and four evaluator
calls.

## Limits

This patch adds no optimization algorithm, physical model, or preparation scope.
Structural consistency and coordinate stereo validation do not establish
computational authenticity, absence of bond-through-ring intersections, global
minimality, thermal equilibration, or scientific suitability.
`production_validated=False` and `simulation_readiness="not_established"` remain
unchanged. MD, periodic systems, packing, new force fields, MLIP, and expansion of
the AmberTools preparation scope remain outside this patch.
