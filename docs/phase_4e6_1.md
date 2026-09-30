# Phase 4E6.1: transactional continuation startup

This correction is based on reviewed Phase 4E6 commit
`8c2630af19c678bab06241d4bc20ccaf832b69dd`. At implementation time,
`origin/main` remained at `12fc65f1049bc8a4dea24416b55335dc5024e6a8`;
Phase 4E6 was unmerged. The fix stays on
`codex/phase-4e6-dynamics-checkpoints` and depends on that feature.

## Reproduction and root cause

A two-site analytical system with zero velocities and `Harmonic(0)` produces
a valid one-step NVE checkpoint. On resuming, an otherwise identical evaluator
adds either +5e-9 or -5e-9 kJ/mol to its potential energy and component.
These differences pass the existing independent verification tolerance but
violate an original NVE energy guard of either zero or 1e-9 kJ/mol.

Before this fix, all four boundary regressions raised
`InvalidDynamicsCheckpointError("NVE original energy guard violated")`.
Startup had already installed the candidate evaluation and cumulative maximum
before checking the guard. Restoring the saved frame left the rejected maximum
in the diagnostic payload. Integrity validation correctly rejected that payload.
Larger offsets of +/-5e-7 kJ/mol were already rejected by independent energy
verification; all four corresponding regression controls passed before the fix.

## Corrected contract

Startup constructs the candidate frame and candidate cumulative maximum locally.
It commits the frame, model, initial energy reference, cumulative maximum and
verified flag only after all startup checks succeed. The failure message can
include the rejected signed deviation; accepted statistics never include it.

A rejected startup returns a structurally valid `DynamicsSegment` with
`termination_reason="startup_verification_failed"` and `startup_verified=False`.
The saved synchronized coordinates, velocities, evaluation, trajectory origin
and cumulative accepted maximum remain unchanged. Exactly one local evaluator
call is counted and added to the prior cumulative count. No step or thermostat
draw is consumed, and no propagation or final verification occurs.

Such a diagnostic cannot create a checkpoint or apply coordinates, including
with `allow_incomplete=True`. The original checkpoint can still resume normally.
Changed model, parameter, backend or settings identities continue to raise
`DynamicsCheckpointCompatibilityError` through the existing typed boundary.
No checkpoint schema, public API, guard, verification tolerance, original energy
reference, or integrity-validation rule changed.

## Verification

The eight new parameterized analytical cases check both offset signs, both
guards, small and larger disagreement paths, diagnostic integrity, accepted
state/origin/maximum preservation, exact counters, checkpoint/application
rejection, input and global RNG non-mutation, and subsequent successful resume.
The existing checkpoint tests also cover larger force disagreements, Langevin
startup failure retaining its complete RNG state, and exact uninterrupted versus
three-segment continuation. The targeted file passes all **40 tests**.

Executed in the existing validation environment (Python 3.11.16, NumPy 2.4.6,
OpenMM 8.6.1 Reference, ParmEd 4.3.1, RDKit 2026.03.6, SciPy 1.17.1):

- Complete ordinary suite: **734 passed, 7 skipped** (61.86 s). Skips are the
  opt-in live AmberTools test and six archived acceptance modules; checkpoint
  archive acceptance was separately enabled below.
- Ruff and `python -m pip check`: passed.
- `examples/checkpoint_continuation.py`: passed in separate processes for NVE
  and Langevin, three five-step segments versus fifteen uninterrupted steps.
  Maximum state differences were zero and RNG states matched. Split runs used
  21 calls versus 17 uninterrupted calls.
- `examples/run_nve.py`: completed 20 steps with 22 calls; maximum absolute
  energy deviation 0.00013617064024629144 kJ/mol.
- `examples/run_langevin.py`: completed 200 steps with 202 calls, independent
  final verification passed, and 1,200 thermostat normal draws were recorded.
- Separately enabled `tests/test_checkpoint_archive.py`: **1 passed** (13.17 s).
  This reran the checksummed capped PE DP=3 and assigned Cl/Br references for
  NVE and BAOAB, 60 steps versus three saved/reloaded 20-step segments at 0.1 fs,
  with unchanged `atol=1e-10`, `rtol=1e-12` and original physical settings.
  No parameterization was rerun or historical record modified.

Analytical offset injection is distinct from actual OpenMM execution in the
ordinary tests, NVE example and archived acceptance. Other opt-in acceptance
suites and live AmberTools execution were not enabled for this corrective phase.

## Limitations

The existing finite nonperiodic scope, conservative environment compatibility,
planned boundary checkpointing, and BAOAB finite-timestep limitations remain.
A valid diagnostic records a rejected startup; it is not a verified restart.
No equilibration or production-readiness claim follows from these tests.
`production_validated=False` and `simulation_readiness="not_established"`
remain unchanged.
