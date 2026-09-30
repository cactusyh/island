# Phase 4E6: versioned dynamics boundaries and cross-process continuation

## Base and scope

Branch `codex/phase-4e6-dynamics-checkpoints` starts from updated `origin/main`
`12fc65f` (merged PR #13). Its tree matches reviewed Phase 4E5 `89a7aaa` and
includes prerequisite `c33671c` through PR #12. No unmerged prerequisite remains.
Main and unrelated checkpoint directories are not modified.

This adds planned segment-boundary persistence for the existing finite,
nonperiodic, unconstrained velocity-Verlet and BAOAB models. It does not save
asynchronously inside a step. The original `run_nve()` and `run_langevin()`
result contracts remain unchanged. Their identical kick/drift expressions were
extracted into shared private kernels used by both original and segment drivers.
There is no new parameter assignment, force construction, or parameterization.

## Public API

All APIs are exported from `island.dynamics`:

- `run_dynamics_segment(system, evaluator, velocities, options,
  velocity_unit="angstrom/ps")` starts a trajectory using existing
  `DynamicsOptions` or `LangevinOptions`.
- `DynamicsSegmentOptions(steps, max_evaluations, max_frames,
  recording_interval=1)` contains execution budgets only.
- `resume_dynamics(checkpoint, system, evaluator, segment_options)` resumes with
  saved physical settings and new segment budgets.
- `create_dynamics_checkpoint(segment)` produces a `DynamicsCheckpoint` only
  from an eligible `DynamicsSegment`.
- `save_dynamics_checkpoint(checkpoint, path, replace=False)` and
  `load_dynamics_checkpoint(path)` persist strict data-only JSON.

```python
from island.dynamics import (
    LangevinOptions, DynamicsSegmentOptions, initialize_velocities,
    run_dynamics_segment, create_dynamics_checkpoint, save_dynamics_checkpoint,
    load_dynamics_checkpoint, resume_dynamics,
)

initial = initialize_velocities(system, temperature_kelvin=300, seed=78123)
with evaluator.open_session() as session:
    first = run_dynamics_segment(
        system, session, initial.velocities,
        LangevinOptions(0.1, 300, 5, 99181, 20, 22, 3, 10),
    )
save_dynamics_checkpoint(create_dynamics_checkpoint(first), "boundary.json")

# In a NEW process: rebuild the external system and evaluator first.
checkpoint = load_dynamics_checkpoint("boundary.json")
with reconstructed_evaluator.open_session() as session:
    next_segment = resume_dynamics(
        checkpoint, reconstructed_system, session,
        DynamicsSegmentOptions(20, 22, 3, 10),
    )
```

The external MolecularSystem must have the same authoritative chemical topology,
stable IDs, masses, assigned stereo, parameter/chemical metadata and provenance.
The caller must reconstruct a compatible bound potential from the same validated
parameter inputs. The checkpoint is not a parameter archive and never contains a
live Context, Integrator, evaluator, executable, pickle or arbitrary Python object.
The reconstructed system's current coordinates may differ: saved coordinates are
explicitly supplied to evaluation without changing caller/default coordinates.

## Schema and owned records

The checkpoint envelope has exactly `payload` and `sha256` fields. The payload
schema is **`island_dynamics_checkpoint_v1`**. SHA-256 covers canonical UTF-8 JSON
of the entire payload (sorted object keys, compact separators, finite JSON
numbers). This is content integrity, not computational authenticity or protection
against a caller who fabricates and rechecksums data.

The payload includes:

| Field | Content |
| --- | --- |
| `units` | angstrom, angstrom/ps, dalton, kJ/mol, ps |
| `physical` | versioned integrator, timestep; original NVE guard or temperature/friction/thermostat seed |
| `environment` | Python/NumPy identities, NumPy build hash, host/OS/CPU identifiers, byte order |
| `system_fingerprint` | chemical system with only established coordinate bookkeeping normalized |
| `masses` | ascending stable-ID/mass rows |
| `origin` | original synchronized state, original total energy, original trajectory fingerprint |
| `lineage` | segment index, start/end absolute steps, local calls/budgets, parent checksum, termination |
| `state` | synchronized accepted coordinates/velocities, absolute step/time, U/K/E, forces/components, coordinate/velocity hashes |
| `final_evaluation` | independently verified evaluation at that boundary |
| `rng` | complete PCG64 state plus algorithm and NumPy version; null for NVE |
| `counters` | cumulative evaluator calls, random steps and scalar Gaussian draws |
| `max_abs_energy_deviation` | cumulative maximum relative to the original trajectory energy |

Backend name/version/platform/settings and physical-model/parameter fingerprints
are carried by the accepted evaluations and compared with the original state.
The trajectory fingerprint binds original state, physical settings, system,
masses and execution environment. Parent checksums identify earlier boundaries;
loading the current file does not require every earlier file or authenticate a
claimed history.

`DynamicsCheckpoint` and `DynamicsSegment` own immutable canonical JSON strings.
Their `payload` accessors return fresh dictionaries. `checkpoint.state` and
`segment.final_state` return owned immutable full-step frame models;
`checkpoint.rng_state` returns a copy. Valid records support `deepcopy` and
`dataclasses.replace`. Constructors and `validate_integrity()` reject malformed
records with `InvalidDynamicsCheckpointError`.

Segment records use `island_dynamics_segment_v1`, additionally holding their
bounded absolute frames and runtime diagnostics. Existing DynamicsResult and
LangevinResult invariants have not been weakened to accommodate nonzero starts.
Old results do not contain captured PCG64 state and cannot be promoted into
resumable checkpoints.

## Actual RNG state and eligibility

The driver reads `rng.bit_generator.state` from its private generator at the
boundary. All PCG64 fields are retained: `bit_generator`, 128-bit `state` and
`inc`, `has_uint32`, and `uinteger`. Load validates their exact types/ranges,
requires an odd increment and verifies assignment round trips. Resume assigns
that state to a new private PCG64 instance. It never reconstructs it from a seed,
replays old steps, advances by `normal_draws`, or assumes Gaussian draws consume
one generator output each.

The established `(N,3)` draw shape and ascending stable-ID/xyz order remain
unchanged, including T=0 and gamma=0. Random-step and scalar-normal counters are
cumulative; eligible boundaries have one random step per accepted BAOAB step.
A failed attempted trial may consume another set of normals in a diagnostic
segment, but that segment cannot create a resumable checkpoint. Global NumPy RNG
state is not touched.

V1 permits checkpoint creation only after:

1. requested segment completion, or an evaluation/frame budget stop before the
   next trial;
2. successful independent startup and final verification.

Initial/startup failures, failed trials, stereo/geometry/energy-guard failures,
and failed or unavailable final verification remain diagnostics and are rejected
by checkpoint creation. A one-call budget cannot produce a checkpoint. A
verified zero-step budget boundary is valid and consumes no thermostat noise.

## Continuation and accounting

Physical settings are immutable across continuation. Resume accepts only
`DynamicsSegmentOptions`, so changing timestep, NVE guard, temperature, friction
or thermostat seed is not an execution-budget option. Mass, force-field/charge,
nonbonded, model and backend changes are rejected through compatibility checks.

Frames carry **absolute accepted steps** and `time_ps = absolute_step*dt_ps`.
Each segment stores its starting boundary, relative interval checkpoints and
final accepted boundary without duplicates. Two adjacent segments therefore
share a boundary; their frame collections should not be concatenated blindly.

Local `segment.completed_steps` is the number advanced by that segment;
`segment.final_state.step` is absolute. `segment.evaluations` is local. The
payload carries cumulative counters and all earlier segment counters. Every
initial/fresh startup, trial, failed evaluator invocation and final verification
counts. Successful M-step segments take **M+2** calls. Three five-step segments
use 21 calls versus 17 for a single fifteen-step segment. Startup verification
provides the initial force, so it is not followed by an uncounted duplicate call.

With a newly opened OpenMM session, each segment constructs three Contexts:
the reusable session and separate fresh startup/final Contexts. Context counts
are separate from logical calls. Evaluations never step an OpenMM Integrator.
A suspect session is invalidated by the existing implementation, while the
parent fresh path remains available to verify the last accepted state.

NVE always checks `abs(E_trial - E_original)` against the **original** guard,
including after several restarts. It does not reset the reference at a segment
boundary. Cumulative maximum deviation includes unsaved accepted steps and
previous segments. Langevin reports the same cumulative energy diagnostic but
has no NVE conservation guard. BAOAB endpoint temperatures retain the 3N DOF
convention and finite-timestep interpretation from Phase 4E5.

On failure the returned segment keeps the last synchronized accepted position
and full-step velocity, plus attempted step, stage, message, verification flags
and actual call counts. Failed startup compares against the saved boundary but
never propagates its saved forces. The saved boundary remains diagnostic data.

## Compatibility and independent verification

Before propagation, load/validation checks schema, checksum, exact fields/types,
finite coordinates/velocities/numerical content, units, inventories, masses,
frame fingerprints, kinetic/total/temperature arithmetic, original identity,
lineage, times, budgets/counters and PCG64 state.

Resume checks the external system with both the checkpoint compatibility identity
and the original bound evaluator's `validate_system()`. Only these established
coordinate metadata fields are normalized: root current source/generation,
coordinate history, minimization and dynamics records, polymer coordinate
source/generation, and the established false/not-established readiness labels.
Chemical, repeat, atom, parameter and historical preparation metadata are retained.
This permits legitimate application of prior coordinate results without ignoring
chemical or potential changes.

V1 conservatively requires exact equality of recorded Python implementation/version,
NumPy version/build configuration hash, OS/platform, host label, machine/CPU
signature and byte order. Backend version, platform and settings must also match.
The Linux CPU signature uses static model/vendor/stepping/microcode/flags, not
varying frequency. There is no override to silently accept a different environment.
These descriptors are compatibility checks, not hardware attestation or a promise
of cross-platform bitwise identity. Reference-platform OpenMM is the tested path.

Fresh startup evaluation compares saved energy, components and forces with
`rtol=1e-10, atol=1e-8` in public units, using the existing independent verification
routine. Fresh forces, never blindly saved forces, drive the next step. Model,
parameter or backend-setting mismatch raises
`DynamicsCheckpointCompatibilityError`; when discovered by startup evaluation,
the exception's `evaluations` attribute is one. A same-identity numerical
verification failure returns a one-call `startup_verification_failed` diagnostic.
Missing optional dependencies retain the existing domain behavior.

## Persistence and application safety

Strict JSON parsing rejects duplicate object keys, duplicate stable-ID rows,
noninteger IDs, unsupported schemas, NaN/Infinity, malformed mappings, extra fields
and inconsistent records. Inventories use explicit integer-ID rows rather than
JSON object keys. Decimal float representations round-trip to the same Python
binary floats; PCG64 integers remain arbitrary-precision integers. No executable
object deserialization is used.

Saving validates first and writes, flushes and fsyncs a temporary file in the same
directory. By default an exclusive hard-link publication refuses an existing
name. Explicit `replace=True` uses `os.replace` atomically. Failures before
publication leave the prior file usable; temporary files are cleaned up on ordinary failures. Abrupt process termination
can leave an unreferenced temporary file; it is not a published checkpoint. IO errors
raise `DynamicsCheckpointIOError`. This guarantees planned process-boundary
publication, not transactional recovery from hardware/filesystem failure.
Parent directories must already exist.

`segment.to_system(system, allow_incomplete=False)` validates compatibility and
assigned coordinate stereochemistry and returns a validated copy. Incomplete
output requires explicit opt-in; failed startup cannot be applied. The shared
coordinate-provenance helper records absolute step/time, original trajectory and
segment identities while preserving earlier records as history. Charges,
parameter assignments and historical preparation signatures are unchanged.

**MolecularSystem coordinates alone are not a restart. Reusing a seed is not
restoring the random stream.** External system/parameter inputs and the saved
checkpoint are required. A completed segment establishes neither equilibration
nor production readiness. Only the last successfully saved boundary survives a
later process interruption; work after that boundary must be repeated.

## Executed acceptance

Analytical unequal-mass, noncontiguous-ID tests compare both old drivers and
uninterrupted/two-/three-segment execution. Coordinates, synchronized velocities,
energies, forces, absolute time, original reference, cumulative deviation and
PCG64 state agree **exactly** in the tested environment. Ordering reversal,
T=0, gamma=0, budget boundaries and retained NVE guarding pass.

`examples/checkpoint_continuation.py` launches three separate Python worker
processes for each integrator. Every worker reconstructs its external analytical
system/evaluator, loads only the previous saved file and exits after saving its
own boundary. The controller compares with an uninterrupted reference: zero
state differences, identical PCG64 state, 15 accepted steps and 0.0075 ps,
21 split calls versus 17 uninterrupted calls. This example needs no optional
scientific backend.

Actual synthetic Amber OpenMM tests reconstruct evaluators and sessions between
segments. Predetermined split/uninterrupted tolerances are absolute `1e-10`,
relative `1e-12` in the respective public coordinate, velocity, force and energy
units. Session-failure injection verifies safe independent final evaluation after
invalidation. Persistence tests cover corrupt/truncated JSON, duplicate keys/IDs,
nonfinite values, missing/malformed RNG, exact 128-bit state, default no-overwrite,
explicit replacement and interrupted replacement preserving the old file.

### Retained PE and Cl/Br references

The existing checksummed archive was available. No parameters were regenerated.
After valid minimization with the existing 0.1 force criterion, the declared
acceptance ran both NVE and BAOAB on capped PE DP=3 and assigned Cl/Br:
60 steps at 0.1 fs (0.006 ps), compared with three saved/reloaded 20-step segments.
Frames are recorded every ten steps. NVE retains the existing explicit
`0.2*sin(site+1+0.7*axis)` velocities and 0.01 kJ/mol guard. BAOAB uses explicit
300 K initialization seed 78123, thermostat seed 99181 and friction 5/ps.

All four comparisons completed and independently verified. Maximum coordinate,
velocity, force and energy differences were **zero**, and PCG64 states and
cumulative energy deviations matched exactly. Each comparison used 62 calls for
the uninterrupted run and 66 across the three segments, with 12 total Contexts
for those four runs. Assigned Cl/Br stereo passed. Source checksums, original
preparation signatures, full checkpoint payloads and observed runtime are in
[the acceptance report](references/phase_4e6/acceptance.json).

```sh
python examples/checkpoint_continuation.py
python scripts/validate_checkpoint_references.py \
  --manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" \
  --output /new/path/continuation.json
ISLAND_AMBERTOOLS_REFERENCE_MANIFEST=/path/to/references.json \
  pytest -q tests/test_checkpoint_archive.py
```

Delivery verification:

- Complete ordinary suite: **726 passed, 7 skipped**, 40.95 s. Skips are the
  opt-in live AmberTools test and six archived acceptance modules.
- Archived checkpoint continuation separately enabled: **1 passed**, 8.24 s,
  covering the four PE/Cl-Br × NVE/BAOAB comparisons above. Other archived suites
  were not re-enabled in this phase; no live AmberTools execution was performed.
- Ruff and `python -m pip check`: passed.
- Checkpoint cross-process, NVE, Langevin, session and velocity-initialization
  examples: all passed.
- Actual real OpenMM evaluations and process-exit tests are distinct from
  targeted injected backend/atomic-write failures.

Recorded four-run acceptance runtimes were 2.755/1.698 s (PE NVE/BAOAB) and
0.630/0.799 s (Cl-Br NVE/BAOAB), including binding, contexts, IO and comparisons.
These are observed wall times, not performance guarantees.

Environment: Python 3.11.16,
NumPy 2.4.6, OpenMM 8.6.1 Reference/double, ParmEd 4.3.1, RDKit 2026.03.6,
SciPy 1.17.1 for optional minimization setup. Analytical checkpoint imports and
execution are tested with OpenMM, RDKit, ParmEd and SciPy blocked together.

No live AmberTools execution, asynchronous mid-step checkpoint, periodic cell,
packing, constraint, COM removal, force-field family, MLIP, or expansion of the
AmberTools preparation scope is introduced. BAOAB finite-timestep bias and lack
of equilibration/production validation remain. Every checkpoint retains
`production_validated=False`, `simulation_readiness="not_established"`.
