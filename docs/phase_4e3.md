# Phase 4E3: bounded nonperiodic NVE dynamics

## Base and public API

Phase 4E2.1 was merged through PR #10. Branch
`codex/phase-4e3-nve-dynamics` starts from `origin/main` at `69a11a7`, whose
content matches reviewed `a35248f`. Main is not modified or merged here.

`island.dynamics` exports `DynamicsOptions`, `DynamicsFrame`, `DynamicsResult`,
and `run_nve(system, evaluator, velocities, options, *, velocity_unit="angstrom/ps")`.
The core integrator needs only NumPy. OpenMM remains optional through the existing
bound `PotentialEvaluator`; no parameter parsing, assignment, force construction,
or AmberTools execution occurs inside dynamics. SciPy is not required for NVE.

```python
from island.dynamics import DynamicsOptions, run_nve

options = DynamicsOptions(
    timestep_fs=0.1, steps=20, max_evaluations=22, max_frames=5,
    recording_interval=5, max_energy_deviation=0.01,
)
# Caller supplies every stable ID, including hydrogens; no velocity initialization.
velocities = {site: (0.1, 0.0, -0.1) for site in system.topology.sites}
result = run_nve(system, bound_evaluator, velocities, options)
print(result.completed_steps, result.requested_steps, result.termination_reason)
print(result.final_state.total_energy, result.max_abs_energy_deviation)
if result.completed:
    moved_system = result.to_system(system)
    final_velocities = result.final_state.velocities
```

The system must be nonempty, atomistic, finite, and have no simulation cell.
Masses must be positive finite daltons. Coordinates, velocities, masses, and
forces cover exactly the authoritative stable IDs. Cartesian order is ascending
stable ID, independently of insertion order. The OpenMM entry point checks the
original system against its bound graph/masses/imported parameters before any
propagation. The implementation owns a system copy; it fixes graph, masses,
assigned stereo and model/parameter identity for the run.

## Units and synchronized velocity Verlet

| Quantity | Public unit |
| --- | --- |
| Position | angstrom |
| Velocity | angstrom/ps |
| Mass | dalton |
| Force | kJ/(mol*angstrom) |
| Potential, kinetic, total energy | kJ/mol |
| Frame time | ps |
| Timestep option | fs (`timestep_fs`) |

The timestep is converted once: `dt_ps = timestep_fs * 0.001`.
Using these conventions, `a = 100*F/m` in angstrom/ps² and
`K = 0.005*sum_i(m_i*dot(v_i,v_i))` in kJ/mol. The factors follow from
angstrom/nm length conversion squared (100 in acceleration, 0.01 in kinetic
energy), with the kinetic one-half factor. No extra force conversion is applied
to the evaluator output.

Each step computes:

```text
v_half = v_n + 0.5*dt*a_n
x_next = x_n + dt*v_half
F_next = evaluator(x_next)
v_next = v_half + 0.5*dt*(100*F_next/m)
```

Only after validation is the full `(x_next, v_next)` state committed. U, K and
U+K are reported at the same physical time. The explicit velocity-Verlet
construction follows the official
[OpenMM CustomIntegrator documentation](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html).
Ordinary OpenMM Verlet uses leapfrog velocity timing and is not used as an
identical-velocity reference.

## Bounds, recording and failure semantics

Timestep, requested steps, evaluator budget and stored-frame budget are explicit;
positive finite values/integers are required. `recording_interval` defaults to 1.
`max_energy_deviation` defaults to 1.0 kJ/mol and must be finite and nonnegative.
It is an **absolute** guard on `abs(E_step - E_initial)`, including when initial
energy is zero or negative. NVE energy is not required to decrease.

Every actual evaluator invocation is counted. A successful N-step run uses
**N+2** calls: initial, N new-coordinate forces, and one independent final check.
The previous accepted force is reused. Failed evaluator calls count too. One
call is reserved for final verification before starting a trial. A one-call
budget can return only the evaluated initial diagnostic state, unverified at the
final boundary. An initial evaluator failure returns an unevaluated initial
frame with finite velocities/K and absent potential/total energy; it cannot be
applied as an evaluated structure.

Stored frames are precisely step zero, accepted multiples of the recording
interval, and the last accepted step (once, even when it is an interval sample).
If another accepted step would exceed that schedule's frame budget, execution
stops before its evaluator call. The interval is not changed and old frames are
not discarded. Every trial is checked even if it will not be stored.

Each trial checks finite positions/full-step velocities/energies/forces, exact
coverage, supported units, matching coordinate fingerprint, component sums,
constant model/parameter/backend/settings, and assigned tetrahedral stereo using
the existing shared Cartesian validator. Stored tags are not treated as proof
of coordinate handedness. Missing required stereo dependencies raise
`DynamicsUnavailableError`; malformed inputs raise `DynamicsInputError`.

Recoverable failures return the **last fully accepted state**, never a
minimum-energy state or a half-updated velocity. Results retain `attempted_step`,
`failure_stage`, message and scalar failure details. Reasons include
`maximum_evaluations`, `maximum_frames`, `evaluation_failed`, `invalid_geometry`,
`stereochemistry_changed`, `energy_guard_exceeded` and `final_evaluation_failed`.
An energy-guard failure records the rejected total energy and signed deviation.
There is no adaptive timestep, velocity rescaling, retry, or automatic minimization.

The independent final check re-evaluates the returned coordinates and compares
energy/components and forces with `rtol=1e-10, atol=1e-8` in their public units.
It never replaces the accepted state's forces or synchronized velocities with
values from a different evaluation. Its record is stored separately. Failure
keeps the accepted trajectory and reports incomplete status, including when all
requested propagation steps had already completed. If an earlier failure also
exists, its primary diagnostics remain and `final_check_error` records the extra
final-check failure.

## Ownership, integrity and provenance

Frames own immutable coordinate and velocity mappings, explicit fingerprints,
kinetic/total energy, signed deviation and evaluation record; `fmax` exposes the
maximum atomic force norm. The result owns masses, options, bounded frames,
completed-step count, actual evaluator accounting, verification record/flags,
maximum absolute accepted energy deviation (including unsaved steps), and stereo
outcome. Backend version, platform and settings are carried by evaluation records.
`dynamics_fingerprint` includes initial coordinates **and velocities**, masses,
system/potential/parameter identities, timestep, budgets, recording/guard options,
and the versioned full-step integrator identifier.

`validate_integrity()` checks actual data before trusting these declarations:
coverage, types, finite content, units, hashes, mass-dependent kinetic arithmetic,
U+K, deviations, frame ordering/times/schedule, model identity, options/counts,
termination/completion consistency, and independent final-record consistency.
Malformed reconstructed results raise `InvalidDynamicsResultError`, a subclass
of `DynamicsInputError`. Valid `dataclasses.replace()` and `deepcopy()` are
supported; immutable nested content is validated before sharing it on deepcopy.
Standalone frame validation checks its own content; the result additionally binds
it to masses, initial energy, options and other frames.

This is consistency validation, not cryptographic authenticity or proof that
caller-supplied forces were genuinely calculated. A sparse stored trajectory
cannot prove all intermediate states from its fingerprints; the engine validates
every trial and tracks the maximum accepted deviation during execution.

`to_system()` validates the result and target identity/coverage/masses, rechecks
stereo for stored frames against the target graph, and returns a validated copy.
An incomplete trajectory requires `allow_incomplete=True`; that choice cannot
bypass malformed records or invalid stereo. A `MolecularSystem` copy contains
coordinates **only**, not a complete restart. Velocities remain explicitly
available on the dynamics result.

The shared coordinate-provenance helper synchronizes root/polymer current fields
and archives superseded generation/minimization/dynamics records. Successful
application uses `nve_dynamics`; diagnostic application uses
`nve_dynamics_diagnostic`. Subsequent minimization or conformation application
retires the current `dynamics` record into history. Chemical/repeat metadata,
charges, parameter records, historical AmberTools preparation and signatures are
not changed or re-signed.

## Independent numerical acceptance

Analytical tests include constant-velocity motion, the closed harmonic solution,
unequal masses, one-step position/velocity algebra, kinetic units, time reversal,
zero-total-energy guarding, and second-order timestep convergence. At the same
0.1 ps end time, observed maximum position errors were:

| dt (fs) | Steps | Error (angstrom) | Maximum energy deviation (kJ/mol) |
| ---: | ---: | ---: | ---: |
| 1.0 | 100 | 3.12307653e-6 | 1.99115004e-5 |
| 0.5 | 200 | 7.80762226e-7 | 4.97785594e-6 |
| 0.25 | 400 | 1.95190125e-7 | 1.24446279e-6 |

Position error ratios were 4.000035 and 4.000009. Time reversal is checked with
absolute tolerances 2e-14 angstrom and 2e-13 angstrom/ps. Force-free positions use
2e-15 angstrom; explicit one-step formulas use 1e-15.

Real OpenMM tests construct path A from ISLAND's resolved records and path B from
the original prmtop using OpenMM's Amber reader. Path B uses a three-operation
`CustomIntegrator` and kinetic expression `m*v*v/2`, Reference platform,
NoCutoff, no constraints, no COM removal, no switching/dispersion correction,
no thermostat/barostat, and no automatic state-update operations. Tests cover
proper and ordered improper models and a real singular trial rejection.

Both paths use the **same authoritative graph masses**, explicitly assigned to
the independent reference particles. Original-prmtop masses are retained in the
report: their maximum differences from graph masses are 0.001 Da in phenol/PE
and 0.004 Da in the Cl/Br case. No mass is redistributed and no charge or potential
parameter is changed. This declared alignment avoids comparing different mass
models and is independent of force-record conversion.

Predetermined checkpoint tolerances:

- Positions: absolute 1e-8 angstrom, relative 0.
- Full-step velocities: absolute 1e-6 angstrom/ps, relative 0.
- Potential: absolute 2e-7 kJ/mol, relative 2e-10.
- Kinetic: absolute 1e-8 kJ/mol, relative 2e-10.
- Forces: absolute 2e-6 kJ/(mol*angstrom), relative 2e-9.

### Retained real cases

The checksummed Phase 4D2.1 archive was loaded via an explicit configurable path
using the existing signature-preserving loader. Each original molecule received
the established deterministic 0.025-angstrom perturbation, then valid local
minimization (fmax 0.1, 500 iterations, 2000 calls). NVE then used the following
settings, declared before outcomes were examined:

- dt 0.1 fs; 20 steps; duration 0.002 ps.
- `v(site, axis) = 0.2*sin(site + 1 + 0.7*axis)` angstrom/ps.
- No velocity rescaling or COM removal; all sites have explicit nonzero motion.
- Absolute energy guard 0.01 kJ/mol; 22 calls; five frames at 0, 5, 10, 15, 20.

All five trajectories completed with verified final states and passed independent
checkpoint comparisons:

| Case | Steps / calls | Max absolute ΔE (kJ/mol) | Runtime (s) |
| --- | ---: | ---: | ---: |
| Phenol GAFF AM1-BCC | 20 / 22 | 1.47325e-6 | 0.569 |
| Phenol GAFF2 AM1-BCC | 20 / 22 | 1.50409e-6 | 0.537 |
| Phenol GAFF2 provided | 20 / 22 | 1.49200e-6 | 0.531 |
| Capped PE DP=3 provided | 20 / 22 | 2.06965e-6 | 0.822 |
| Assigned Cl/Br provided | 20 / 22 | 9.98065e-8 | 0.452 |

The Cl/Br trajectory retained assigned tetrahedral stereo; other cases have no
assigned centers. Largest observed independent errors were 9.61e-13 angstrom
(position), 9.51e-10 angstrom/ps (velocity), 2.11e-8 kJ/mol (potential),
1.08e-11 kJ/mol (kinetic), and 6.80e-9 kJ/(mol*angstrom) (force).

Runtime covers NVE including initial and final evaluator calls; it excludes
minimization and the independent reference. The first phenol run took 0.569 s
for 22 calls (about 0.026 s/call including Python bookkeeping). This is adequate
for the declared short acceptance run; the evaluator's per-call private Context
is retained. No shared Context or alternate force construction was introduced.

[Acceptance manifest](references/phase_4e3/acceptance.json) includes checkpoints,
velocities, both mass inventories, historical source checksums/signatures,
identities, numerical errors, timings and charge residuals. Phenol AM1-BCC retains
its 0.002 e tolerance and approximately -0.001 e residual, with no neutralization
or new attribution of the residual. The older external phenol fixture's unknown
exact force-field/charge provenance remains unchanged.

```sh
python scripts/validate_dynamics_references.py \
  --preparation-manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" \
  --output /new/path/nve-acceptance.json
pytest -q tests/test_dynamics_archive.py
```

The script refuses to overwrite an existing report. Incomplete minimization or
NVE outcomes are reported as failures, without changing tolerances.

## Verification and limitations

Environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1 Reference, ParmEd 4.3.1,
RDKit 2026.03.6; SciPy 1.17.1 is used only for optional minimization setup.
The existing separate Phase 4E2 validation environment was reused unchanged.

Final verification:

- Complete ordinary suite: **614 passed, 4 skipped**. Skips are opt-in live
  AmberTools generation and the single-point, minimization, and NVE archive tests.
- Archived NVE test, explicitly enabled separately: **1 passed**, covering all
  five actual retained cases. The acceptance script also completed all five and
  retained the reproducible numerical manifest.
- Ruff and `python -m pip check`: passed.
- `run_nve.py`, `minimize_geometry.py`, `evaluate_singlepoint.py`,
  `build_local_template_polymer.py`, `generate_random_walk.py`, and
  `build_linear_polymer.py`: all passed.

Ordinary tests include actual OpenMM execution; archived tests use retained real
AmberTools outputs.
No live AmberTools generation or parameterization was needed or executed in
this phase. Optional-import tests block SciPy, OpenMM, RDKit and ParmEd together
and still execute/apply analytical NVE successfully.

`examples/run_nve.py` requires optional OpenMM/ParmEd, but no AmberTools. Its
synthetic stretched bond has initial total energy 67.83132995 kJ/mol and final
67.83119378 after 20 steps / 22 calls, maximum deviation 1.36171e-4 kJ/mol. Other
executed examples cover single points, minimization, local-template building,
random-walk conformation and linear-polymer construction.

Numerical agreement and bounded energy error establish short-trajectory
integration consistency only. They do not establish thermal equilibration,
scientific suitability, absence of bond-through-ring intersections, production
readiness, or long-chain scalability. Finite-precision velocity Verlet does not
conserve total energy exactly. Timestep and energy guard are caller choices.
`production_validated=False` and `simulation_readiness="not_established"` remain
fixed. NVT/NPT, pressure, periodicity, packing, constraints/restraints, trajectory
formats, export, new force fields, MLIP, long-chain charge transfer and expansion
of the existing AmberTools preparation limit remain outside this phase.
