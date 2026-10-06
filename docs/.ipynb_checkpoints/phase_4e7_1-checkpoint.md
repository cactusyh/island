# Phase 4E7.1: workflow consistency and valid boundary publication

## Base and reproduced defects

Phase 4E7 remains unmerged: fetched `origin/main` was `606d439`, while the
reviewed feature tip was `d75dc1cfe3ecd932f0ce9fac3b476470039e9f07`.
This correction stays on `codex/phase-4e7-single-chain-workflow`. Main and
unrelated checkpoint directories were not changed.

Six parameterized regressions failed before the correction:

1. Both +5e-9 and -5e-9 kJ/mol injected into the first fresh startup potential
   and one component passed checkpoint verification and publication, but later
   workflow reads rejected exact frame inequality. This is controlled boundary
   injection, not observed OpenMM nondeterminism.
2. Independently valid segments/checkpoints using twice the recorded velocities
   or a different velocity seed were accepted despite disagreeing with the saved
   `VelocityInitialization`. All replacement files had correct checksums.
3. A rechecksummed manifest changing temperature from 300 to 777 K passed status
   and frame reading in both paused and completed runs. Completed resume returned
   early without validating the immutable bundle/config relationship.

The causes were separate per-file validation without setup-to-trajectory binding,
whole-frame equality at a numerically verified restart boundary, and bundle
validation limited to propagation paths.

## Shared consistency contract

`workflow_status()`, `read_workflow_frames()` and `resume_workflow()` now use the
same manifest consistency validation. Completed workflows may return without
propagation only after these checks. Every production manifest publication also
validates the candidate against this contract **before** atomic replacement.
An invalid candidate can leave unreferenced immutable files, but the prior
manifest bytes and checkpoint remain authoritative. Recovery continues to ignore
orphans; no historical records are silently repaired or re-signed.

Boundary checks require exact frame type, stable-ID positions and synchronized
velocities, step/time, units, coordinate/velocity fingerprints and DOF.
Kinetic energy and derived temperature use the existing frame-arithmetic
tolerances (`atol=1e-10`, `rtol=1e-12`). Evaluation coordinate/calculation fingerprints also match
exactly. Existing `verify_final()` supplies exact model, parameter, backend and
settings identity checks and the existing numerical comparisons for potential
energy, all forces and components: **atol 1e-8, rtol 1e-10**. Each frame must first
pass its own integrity checks, including finite values and energy arithmetic.
Total energy and energy deviation may reflect the independently recomputed
potential; they are not required to be bitwise equal. No tolerance increased.

The combined reader deterministically retains the **earlier accepted frame** at
a compatible shared boundary. The later segment's startup evaluation remains
unchanged in its original signed record. It is not substituted into the earlier
record. Discrepancies outside tolerance, altered positions/velocities, model
changes or conflicting timing remain errors.

## Binding durable setup and trajectory

Whenever a bundle is present, validation reconstructs and validates its existing
preparation/import signatures, minimum, starting system and initialization. It
checks the original trajectory against those durable records:

- Original positions and coordinate fingerprint equal the minimized starting
  system; velocities and velocity fingerprint equal the saved initialization.
  Velocities are never regenerated for validation.
- Initialization masses agree with authoritative starting-system masses. Each
  segment has that same stable-ID mass inventory and physical configuration.
- The minimum's parameter identity equals the reconstructed imported content
  signature. The trajectory's original evaluation agrees with the minimum under
  existing verification rules, including physical-model/backend identities.
- The initialization's system fingerprint uses its established unnormalized
  `system_identity()` definition. Checkpoint compatibility uses the existing
  coordinate-provenance-normalizing `compatibility()` definition. These distinct
  hashes are never compared as if interchangeable.
- The first accepted frame agrees with the original trajectory state. Every
  later segment preserves the complete original trajectory identity, physical
  settings and environment and extends the prior accepted lineage. Existing
  parent-checksum and selected-checkpoint equality checks still apply.
- Ineligible diagnostic segments cannot enter the accepted segment sequence.

This is internal structural/numerical consistency, not proof that arbitrary
caller-supplied energies were computed or that rechecksummed data are authentic.
Checkpoint schemas, original signatures, parameter content, integration, budgets,
charge policies, physical settings and readiness flags are unchanged.

## Inspection and dependencies

Missing bundles remain legitimate during incomplete or failed setup. Those runs
remain inspectable through manifest/file checks without pretending that they
have a resumable prepared trajectory.

For a published bundle, status and frame reads now perform deeper reconstruction:
ParmEd is required to validate the original Amber import, and RDKit is needed
when assigned stereochemistry requires coordinate verification. No AmberTools
executable, minimization run, velocity initialization, OpenMM Context or new force
evaluation is invoked by consistency validation. SciPy optimization and OpenMM
execution remain requirements of their respective execution operations, not
these read-only checks. Imports remain lazy. The existing conservative checkpoint
environment policy still governs actual continuation.

## Verification

The focused regressions cover both tolerated offset signs, larger energy and
force disagreements, changed boundary coordinates/velocities/model, two valid
but mismatched initializations, paused/completed config disagreement, transactional
candidate rejection preserving prior bytes, deterministic boundary deduplication,
unchanged historical segment bytes and completed inspection with scientific
execution functions prohibited. Existing tests retain failure/setup inspection,
relocation, separate-process continuation, ownership, budgets, optional imports
and exact evaluator accounting.

Actual acceptance was rerun using the unchanged Phase 4E7 script/configuration:
PE DP=3, 20 explicit sites, GAFF2, 300 K, friction 5/ps, timestep 0.1 fs, 60 steps,
three 20-step segments, recording every ten steps, velocity seed 78123 and
thermostat seed 99181. Numerical comparison remained `atol=1e-10, rtol=1e-12`.

| Actual execution | Charge tolerance / residual (e) | Minimized fmax | Minimum calls | Split / uninterrupted calls |
| --- | --- | ---: | ---: | ---: |
| Live PSMILES → GAFF2/AM1-BCC | 1e-4 / -1.999967073326725e-6 | 0.0752321217956569 | 141 | 66 / 62 |
| Checksummed archived PE, provided charges | historical 1e-4 / 0 | 0.08006897945685111 | 139 | 66 / 62 |

Both completed 0.006 ps after directory relocation and separate-process resumes.
Maximum coordinate, full-step velocity, force and energy differences against
uninterrupted execution were **zero**; complete PCG64 states matched exactly.
Minimization used the original fmax 0.1 kJ/(mol*angstrom) criterion and independent
verification. Neither charge method nor any tolerance was changed after results.
Archived parameters were not regenerated. Generated live artifacts and complete
reports are retained outside the repository in the user validation directory
`island-validation/phase4e7-1/` (`live-report.json`, `archive-report.json` and
corresponding relocated bundles).

Environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1 Reference, ParmEd 4.3.1,
RDKit 2026.03.6, SciPy 1.17.1, AmberTools 24.8. Controlled offsets and injected
publication/engine failures are software regressions, separate from actual live
and archived numerical acceptance.

Delivery checks:

- Complete ordinary suite: **774 passed, 8 skipped** (87.62 s), including the
  14 new consistency regressions. The skips are the opt-in five-case AmberTools
  generator and seven archived acceptance modules.
- Separately enabled archived workflow and checkpoint tests: **2 passed**
  (22.59 s). Other archived modules and the older five-case live regeneration
  test were not enabled in this focused phase.
- Dedicated live and archived PE workflow scripts: both passed, as detailed above.
- Ruff and `python -m pip check`: passed.
- Cross-process checkpoint example and workflow inspection, status and completed
  resume examples: passed. Checkpoint example split differences were zero for
  NVE and Langevin, with matching RNG states (21 split / 17 uninterrupted calls).
- A valid retained Phase 4E7 live bundle also passed the corrected status reader;
  there is no schema migration or rewrite of older valid records.

## Remaining limitations

Deep inspection now incurs parameter reconstruction and requires its optional
dependencies. It performs no dynamics and does not establish computational
authenticity. Existing POSIX publication/ownership constraints, manual stale-lock
recovery, planned-boundary persistence, same-compatible-host continuation,
100-site/nonperiodic scope and BAOAB finite-timestep limitations remain.
Completion is not equilibration or scientific suitability.
`production_validated=False` and `simulation_readiness="not_established"` remain.
