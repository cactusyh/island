# Phase 4E5: explicit thermal initialization and isolated Langevin dynamics

## Repository and prerequisite

Phase 4E4 is unmerged at implementation time. The prerequisite correction is
commit `c33671c` on `codex/phase-4e4-evaluation-sessions`, following reviewed
`9588b42962b26a10cec3eea07463452370f762a3`. This feature branch,
`codex/phase-4e5-langevin-dynamics`, starts at the corrected tip. Merge Phase 4E4
including its correction before Phase 4E5. Main and checkpoint files are untouched.

The regression first failed before the correction. Opposite finite coordinates
of magnitude 1e308 overflow during geometry subtraction. The shared validator now
raises `EvaluationInputError` chained from the numerical exception in both fresh
and session paths, before position installation. The regression checks that no
Context positions were set and subsequent valid session evaluation succeeds.
Post-installation backend failures still invalidate sessions. All 24 session
regressions passed after this separate correction.

## Public APIs

```python
from island.dynamics import (
    initialize_velocities, LangevinOptions, run_langevin,
)

initialization = initialize_velocities(system, temperature_kelvin=300, seed=78123)
options = LangevinOptions(
    timestep_fs=0.1, temperature_kelvin=300, friction_per_ps=5,
    thermostat_seed=99181, steps=200, max_evaluations=202,
    max_frames=5, recording_interval=50,
)
with evaluator.open_session() as session:
    result = run_langevin(system, session, initialization.velocities, options)
if result.completed:
    moved = result.to_system(system)
```

`VelocityInitialization`, `LangevinOptions`, `LangevinFrame`, and
`LangevinResult` are distinct thermal contracts. Existing NVE options and
velocity-Verlet propagation are unchanged. The result envelope and content
validators are shared, with explicit integrator/frame type checks. There is no
new parameter parser or OpenMM force construction.

The initializer needs only core ISLAND/NumPy. The analytical Langevin interface
also works without OpenMM, RDKit, ParmEd or SciPy. Bound OpenMM evaluators and
sessions retain authoritative graph, mass, stable-ID, stereochemical and parameter
validation before energy calls. Assigned tetrahedral validation still requires
RDKit and raises `DynamicsUnavailableError` when it is missing.

## Units, degrees of freedom and thermal sampling

| Quantity | Public unit/convention |
| --- | --- |
| Coordinates | angstrom |
| Full-step velocity | angstrom/ps |
| Mass | dalton |
| Force | kJ/(mol*angstrom) |
| Energy | kJ/mol |
| Time | ps; `timestep_fs` converted once by 0.001 |
| Temperature | kelvin |
| Friction | 1/ps |
| Degrees of freedom | all 3N unconstrained Cartesian components |

The molar gas constant is `R=0.00831446261815324 kJ/(mol*K)`, from exact SI
Boltzmann and Avogadro constants. Cartesian components are independent Gaussians
with mean zero and variance `100*R*T/m`. The initializer reports masses and
mass/system/velocity identities, target and instantaneous temperature, kinetic
energy, DOF, velocity seed, RNG algorithm and NumPy version.

No translation or rotation is removed; no sample is rescaled. An individual draw
need not have its requested temperature or exact mean kinetic energy. T=0 is
supported and produces zero velocities. Temperatures must be finite and
nonnegative. Seeds must be Python integers in `[0, 2**128)`; booleans are rejected.
All IDs must be integers, masses positive and finite, and inventories complete.
Malformed inputs raise `DynamicsInputError` with chained causes where applicable.

For every synchronized frame:

```
a = 100*F/m                         # angstrom/ps^2
K = 0.005*sum_i(m_i * dot(v_i,v_i)) # kJ/mol
T_inst = 2*K/(3*N*R)               # kelvin, including translation and rotation
E_total = U + K
```

The engine does not convert evaluator forces again.

## BAOAB full-step scheme

The versioned scheme is `baoab_full_step_v1`. With dt in ps:

```
v_B = v + (dt/2)*100*F(x)/m
x_A = x + (dt/2)*v_B
c = exp(-gamma*dt)
s = sqrt(-expm1(-2*gamma*dt) * 100*R*T/m)
v_O = c*v_B + s*normal
x_next = x_A + (dt/2)*v_O
F_next = evaluate(x_next)
v_next = v_O + (dt/2)*100*F_next/m
```

The product gamma*dt is computed before forming its doubled exponent, avoiding
premature overflow when gamma is large and dt tiny. `expm1` avoids cancellation
for tiny products. Gamma=0 is the deterministic velocity-Verlet limit, within
floating-point drift-order rounding. It does not thermalize. The stored velocity
is the endpoint after the final B operation, not a leapfrog or stochastic-substep
velocity. U, K, and temperature refer to that same endpoint.

The splitting follows Leimkuhler and Matthews:
[Rational Construction of Stochastic Numerical Methods for Molecular Sampling](https://arxiv.org/abs/1203.5428)
and [Robust and efficient configurational molecular sampling via Langevin Dynamics](https://arxiv.org/abs/1304.3269).
Finite-timestep endpoint observables have a discretized stationary distribution.
For a scalar harmonic well, this implementation's stationary variances are
`<x²>=R*T/k` and `<v²>=(100*R*T/m)*(1-omega²*dt²/4)`, with
`omega²=100*k/m`, in the stable regime. Thus the endpoint kinetic temperature
can be biased even when harmonic configurational variance is exact. General
molecular observables need their own timestep and sampling assessment.

## Random streams, execution bounds and failure behavior

Both operations explicitly create a local
`numpy.random.Generator(numpy.random.PCG64(seed))`. No global NumPy RNG is changed.
Each stream draws standard normals in ascending stable-ID, then x/y/z order.
Initialization's `velocity_seed` and propagation's `thermostat_seed` are separate;
calling the initializer does not advance the thermostat stream.

Every begun BAOAB trial consumes exactly 3N normal variates, including T=0 or
gamma=0. `random_steps` includes an unsuccessful trial, and `normal_draws` counts
scalar Gaussian variates, not raw PCG64 words. Evaluator and frame capacities are
checked before beginning a trial or drawing noise. There are no retries, hidden
rescalings or timestep changes. Failed noise is not rewound. Seed, NumPy version,
RNG identity and integrator/options belong to the dynamics fingerprint along with
initial coordinates/velocities, masses and physical-model identity.

Identical inputs/seeds reproduce on the tested software/platform configuration.
There is no cross-platform bitwise guarantee, and these records are not a restart
format. A seed and endpoint coordinates do not constitute RNG or velocity state.

Initial evaluation, accepted trials, failed evaluation attempts and final checks
all count as actual calls. Previous accepted forces are reused. Successful N-step
runs take N+2 calls: a new session uses one reusable Context plus one fresh Context
for final verification. The latter remains available after session invalidation.
The exact same budget and last-accepted-state rules used by NVE are retained.
Storage contains step zero, interval checkpoints and the last accepted state,
without duplicates. When another state would exceed the frame budget, no trial
or noise is attempted.

Every trial validates finite coordinates, velocities, energy/components/forces,
coverage, coordinate fingerprint, physical and backend identities and assigned
stereochemistry, even if it will not be saved. A failed trial is never committed:
the result retains the last synchronized endpoint, with attempted step, failure
stage, failed-call count and normal-draw accounting. Final verification is separate
from completion. A one-call budget can return a valid unverified diagnostic.

There is **no NVE conservation guard**. Energy deviations remain reported as
observations; thermal exchange may increase or decrease U+K. Completion means
requested propagation plus independent final verification, not convergence,
equilibration or canonical-ensemble validation of that trajectory.

## Ownership, integrity and coordinate provenance

Records own immutable mappings; valid records support `dataclasses.replace` and
`deepcopy`. `validate_integrity()` checks types, inventories, finite scalar/vector
content, units, coordinate/velocity identities, mass-dependent kinetic energy,
3N temperature arithmetic, RNG identity/draw accounting, times and frame schedule,
call counts, termination and final-verification consistency. Malformed records
raise `InvalidDynamicsResultError`. This is structural and numerical consistency,
not authentication of caller-supplied calculations.

`to_system()` checks result integrity and target compatibility, rechecks assigned
coordinate stereo and returns a validated copy. Incomplete diagnostics require
`allow_incomplete=True`; this never permits malformed content or invalid stereo.
The current source is `langevin_dynamics`, or `langevin_dynamics_diagnostic` for
explicitly applied incomplete output. The shared provenance helper updates
polymer and ordinary namespaces and retires previous minimization, conformation
and dynamics records to history. Later minimization/conformation application
retires the Langevin record. Preparation records, historical signatures, charges
and parameter provenance are unchanged. MolecularSystem alone is not a restart;
velocities and RNG metadata remain in the explicit result/initialization records.

## Numerical evidence and reproducibility

### Deterministic and contract tests

Ordinary tests cover gamma=0/NVE agreement (absolute 3e-15, relative 1e-14), a
prescribed-noise one-step calculation with unequal masses, zero temperature,
tiny gamma*dt, stable-ID/insertion-order independence, isolated RNGs, exact call
and frame budgets, partial-trial and final-check failures, mutable input isolation,
reconstructed-result rejection and stereo/provenance application. Real OpenMM
checks include session failure after position installation and fresh verification
of the unchanged accepted state. Failure injection supplements actual backend
execution; it does not replace it.

### Independent ensemble statistics

[Statistical report](references/phase_4e5/statistics.json) was generated by
`scripts/validate_langevin_statistics.py` with settings declared before execution:
2048 independent analytical atoms, half mass 2 and half mass 8 dalton, providing
3072 independent Cartesian samples per mass group. T=300 K; velocity seed 52117,
thermostat seed 67231. All acceptance thresholds are six standard errors.
No seeds or thresholds were changed after examining results.

- Gaussian initialization mean/variance and per-atom chi-square(3) kinetic-energy
  mean/variance passed. Largest Gaussian standardized discrepancy was 1.219.
- Force-free OU: 20 steps at 20 fs, gamma=3/ps, starting velocities all one.
  Finite-time mean `exp(-gamma*t)` and variance
  `(100*R*T/m)*(1-exp(-2*gamma*t))` passed; maximum discrepancy 1.558 standard errors.
- Harmonic confinement k=m: gamma=5/ps; independently derived 2x2 transition and
  discrete Lyapunov covariance, checked against the closed formula above.
  After 200 burn-in steps, a single endpoint from each independent replica was
  used. Timesteps 100 and 50 fs give kinetic-variance factors 0.75 and 0.9375;
  observed coordinate/velocity moments and cross-covariance passed (largest
  standardized discrepancy 2.140). These intentionally large analytical-test
  timesteps are not suggested for molecules.

No correlated time samples are treated as independent. The ensemble report takes
about 19.2 s in the recorded environment. A preliminary execution was stopped to
avoid repeatedly validating the same immutable endpoint inside a comprehension;
this bookkeeping correction did not change sampling, seeds or criteria.

### Independent OpenMM and five archived molecules

The reference path constructs its forces from the original prmtop using OpenMM's
Amber reader and a separate explicit BAOAB CustomIntegrator, with no COM removal,
constraints, implicit solvent or cutoffs. Authoritative masses are installed in
both paths, with original prmtop masses retained in the report. The reference
receives prescribed Gaussian increments via a per-DOF variable at each step;
it never assumes NumPy and OpenMM seeds generate the same values. The production
engine has no public noise-file/replay feature; the validation script replaces
only its private normal-draw function, and ordinary tests also compare the
unpatched production RNG path. Force evaluation is real and unmocked.

Predetermined tolerances remain those of the independent Phase 4E3 reference:
positions 1e-8 angstrom, velocities 1e-6 angstrom/ps, potential 2e-7 kJ/mol
(rtol 2e-10), kinetic 1e-8 kJ/mol (rtol 2e-10), forces
2e-6 kJ/(mol*angstrom) (rtol 2e-9). The four-site synthetic test passes as well.

The [archived acceptance report](references/phase_4e5/acceptance.json) covers all
five checksummed references after valid minimization at the existing force
threshold 0.1 kJ/(mol*angstrom), 500 iterations and 2000 calls. Declared settings:
300 K initialization/thermostat, gamma=5/ps, dt=0.1 fs, 200 steps (0.02 ps),
velocity seed 78123, thermostat seed 99181, checkpoints every 50 steps. Every case
completed with 202 calls, two Contexts, a fresh final check and valid stereo.

| Reference | Initial T (K) | Final T (K) | Final fmax (kJ/(mol*angstrom)) |
| --- | ---: | ---: | ---: |
| Phenol GAFF AM1-BCC | 433.86 | 232.36 | 266.58 |
| Phenol GAFF2 AM1-BCC | 433.86 | 259.92 | 240.45 |
| Phenol GAFF2 provided | 433.86 | 261.41 | 239.71 |
| Capped PE DP=3 provided | 394.86 | 236.25 | 239.14 |
| Assigned Cl/Br provided | 310.84 | 172.43 | 153.27 |

Maximum independent discrepancies across cases: position 3.51e-11 angstrom,
velocity 5.16e-9 angstrom/ps, potential 2.19e-8 kJ/mol, kinetic 1.05e-9 kJ/mol,
force 3.54e-8 kJ/(mol*angstrom). All pass unchanged tolerances. Energy changes
of several kJ/mol are bath-exchange diagnostics, not conservation failures.
Short molecular runs are execution checks, not equilibration evidence.
AM1-BCC historical 0.002 e tolerance and approximately -0.001 e charge residual
remain unchanged; exact residuals and source checksums are in the report.

```sh
python examples/initialize_velocities.py
python examples/run_langevin.py
python scripts/validate_langevin_statistics.py --output /new/path/statistics.json
python scripts/validate_langevin_references.py \
  --preparation-manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" \
  --output /new/path/acceptance.json
ISLAND_AMBERTOOLS_REFERENCE_MANIFEST=/path/to/references.json \
  pytest -q tests/test_langevin_archive.py
```

Scripts refuse to overwrite reports. Archives are configured explicitly; no
server path is embedded in a public API. No AmberTools parameterization is run.

## Verification and limitations

Environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1 Reference/double,
ParmEd 4.3.1, RDKit 2026.03.6, SciPy 1.17.1. No base environment was modified.
Executed delivery checks:

- Ordinary suite: **689 passed, 6 skipped** in 36.05 s. The six opt-in skips are
  live AmberTools and archived single-point, minimization, NVE, session and Langevin
  acceptance. No ordinary test needs AmberTools execution or external archives.
- Separately enabled archive suite: **5 passed** in 107.97 s, covering all five
  archived acceptance modules (including the Phase 4E4 1000-step reference case).
- Statistical ensemble script: passed with the recorded seeds and thresholds.
- Ruff: passed. `python -m pip check`: no broken requirements.
- Eight examples passed: thermal initialization, Langevin, session, single-point,
  minimization, NVE, local-template polymer and random-walk conformation.
- No live AmberTools execution or parameter regeneration was performed.

Recorded archived Langevin runtimes were 1.313, 1.522, 0.720, 2.000 and 0.271 s
in table order, including session initialization/final verification/cleanup.
These measurements overlapped other validation jobs and are not speed guarantees.
The ordinary synthetic two-site example completed 200 steps (0.02 ps), 202 calls,
with initial/final instantaneous temperatures 454.83/413.66 K and a verified final
state. Actual backend execution and numerical references are separate from the
explicit failure-injection tests.

This supports thermostatted isolated finite molecules. It establishes neither
bulk NVT/density/pressure nor equilibration, absence of topological crossings,
production suitability or long-chain scaling. There is no COM/rotation removal,
periodic cell, packing, constraint, restraint, restart format, new force field,
MLIP integration or expansion of the AmberTools 100-site preparation scope.
The larger analytical statistical ensemble is a test of independent oscillators,
not an expansion of that preparation scope. Scientific flags remain
`production_validated=False`, `simulation_readiness="not_established"`.
