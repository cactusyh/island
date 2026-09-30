# Phase 4E4: explicitly owned OpenMM evaluation sessions

## Repository state and scope

Phase 4E3 was merged through PR #11. Branch
`codex/phase-4e4-evaluation-sessions` starts at `6db3984` from `origin/main`,
whose content matches reviewed `728cda7`. Main is not modified or merged here.
Existing partial session work on this branch was inspected and completed.

This changes resource lifetime, not the physical model, force construction,
minimization algorithm, NVE integrator, units, or acceptance thresholds.
No parameterization is performed within an evaluation or run.

## Explicit API and lifecycle

```python
from island.dynamics import DynamicsOptions, run_nve
from island.minimization import minimize_geometry

with evaluator.open_session() as session:
    point = session.evaluate(coordinates, coordinate_unit="angstrom")
    minimum = minimize_geometry(system, session)
    if minimum.converged:
        starting = minimum.to_system(system)
        trajectory = run_nve(
            starting, session, velocities,
            DynamicsOptions(0.1, 20, 22, 5, recording_interval=5,
                            max_energy_deviation=0.01),
        )
```

`OpenMMEvaluationSession` implements the existing `PotentialEvaluator` contract.
`open_session()` immediately creates one private System, Integrator and Context.
The shared private resource bundle deserializes the evaluator's existing owned
XML; it does not parse prmtop or duplicate force construction. Ordinary
`evaluator.evaluate()` remains a fresh resource bundle on every call. Existing
callers do not silently opt into reuse.

`close()` deterministically releases Context, Integrator, then System. It is
idempotent. Context-manager exit closes on both normal and exceptional exit.
The public `closed` property exposes lifecycle status without exposing resources.
`session.evaluate()` and `validate_system()` after closure raise
`EvaluationSessionClosedError`; neither allocates a replacement Context.
Sessions cannot be copied or deep-copied; open a separate session instead.

### Concurrency policy

A session belongs to its creating Python thread. Evaluation, verification,
binding validation, entry, and close from another thread raise
`EvaluationSessionBusyError`. A nonblocking lock rejects overlapping or reentrant
operations on the owner thread. Nested entry of the same context manager is also
rejected. Different sessions own independent resources and can be interleaved;
there is no global cache or pool. Always use explicit close or a `with` block;
resource lifetime is not defined in terms of garbage-collection timing.

### Validation and failure behavior

Every evaluation validates complete integer-ID coverage, finite N x 3 angstrom
coordinates and singular geometry before touching Context positions. Each valid
call explicitly installs the entire coordinate frame, even when its fingerprint
matches a previous frame. Results are not cached. The unchanged extraction path
checks finite potential, every component, and all forces and returns an owned
immutable `EvaluationResult`. Positions, velocities and parent defaults are never
modified on the caller's system. The Context Integrator is never stepped.

An input error leaves the reusable Context available. Any backend exception after
state installation begins invalidates and closes the session; future ordinary
session evaluations fail closed. Earlier result objects remain unchanged.
Failed Context construction releases partially constructed resources. Actual
OpenMM tests cover these paths; targeted injected failures supplement real
numerical execution.

## Compatibility and independent final verification

`OpenMMBoundPotential` is the shared compatibility marker implemented by both
the original evaluator and its session. Dynamics and minimization check it and
call `validate_system()` before any energy call. The session delegates to the
original bound evaluator's existing content validation: actual imported parameter
integrity, chemical graph, IDs, masses, and assigned stereochemical metadata.
Regressions change graph/masses with unchanged IDs and reject both NVE and
minimization before their first evaluation. Analytical evaluators continue using
`PotentialEvaluator` without requiring OpenMM, ParmEd, RDKit, or SciPy imports.

Both engines dispatch their final independent call through `evaluate_fresh()`
for bound OpenMM potentials. On a session this explicitly delegates to its parent
fresh-context path, including after session invalidation/closure. This is the
sole intentional exception to the closed-session evaluation rule: it performs an
explicit separate verification calculation, never reopens or repairs the session.
It remains subject to owner-thread/non-reentrancy checks. The parent evaluator
itself remains independently usable after session failure.

The call is charged once to the existing run budget. For successful N-step NVE:

| Execution | Actual evaluator calls | Context constructions |
| --- | ---: | ---: |
| Fresh calls | N+2 | N+2 |
| One new session | N+2 | 2 |

The session counts N+1 evaluations (initial plus steps); final verification uses
one fresh Context. A one-call run budget cannot perform final verification and
remains diagnostic. Failed trial calls still count. The original last-accepted
state behavior is preserved even when backend failure invalidates the reusable
Context. Minimization uses the same separation: a new session plus one fresh
final verification, with unchanged logical evaluation budget/convergence rules.

Lifecycle and timing diagnostics are separate from physical settings. Model,
parameter and evaluation fingerprints do not change solely because a Context was
reused. Backend version, platform and precision reporting remain unchanged;
there is no result-schema or physical-model version change.

## Baseline and benchmark procedure

Initial fresh-path profiling (50 evaluations, after warmup) measured:

| Fixture | Sites | Binding (s) | 50 calls (s) | Context constructor cumulative (s) |
| --- | ---: | ---: | ---: | ---: |
| Synthetic Amber | 4 | 0.00274 | 0.98191 | 0.623 |
| Capped PE DP=3 | 20 | 0.03770 | 1.67076 | 0.783 |

These initial figures include profiler overhead and are motivation, not targets.
The reproducible benchmark uses three unprofiled measured repetitions after a
complete unmeasured warmup per workload/mode. Both modes use identical coordinates,
velocities, Reference platform/double precision, validation, component reporting
and full run bookkeeping. It reports medians and every sample for:

- Model binding separately from evaluation.
- Resource initialization (including deserialization/System/Integrator) and
  Context constructor time separately; the latter is a subset of initialization.
- Repeated evaluation or NVE work duration.
- Resource cleanup separately; fresh cleanup occurs inside each call.
- End-to-end duration including session opening, final checks and closing.
- Actual Context construction counts, distinct from logical evaluator calls.

Instrumentation wraps real operations without substituting numerical results.
A separate cProfile run identifies remaining costs; its timings are not included
in the benchmark medians. No wall-clock assertion is part of ordinary tests.

Measured medians (seconds, three repetitions, including initialization and cleanup):

| Fixture | Workload | Fresh | Session | Ratio | Contexts fresh/session |
| --- | --- | ---: | ---: | ---: | ---: |
| Synthetic, 4 sites | 50 evaluations | 0.99614 | 0.02696 | 36.95x | 50/1 |
| Synthetic, 4 sites | 20-step NVE | 0.48895 | 0.05701 | 8.58x | 22/2 |
| PE DP=3, 20 sites | 50 evaluations | 1.55042 | 0.21534 | 7.20x | 50/1 |
| PE DP=3, 20 sites | 20-step NVE | 0.67523 | 0.20813 | 3.24x | 22/2 |

Binding medians were 0.00212 s (synthetic) and 0.01783 s (PE). The manifest
separates resource initialization, Context construction, querying, cleanup and
all individual samples. Timings vary with host load; these are measured outcomes,
not portable speedup guarantees. The PE session profile spent 0.140 s in singular
geometry checks and 0.104 s in bound-system validation out of a profiled 0.415 s
NVE run. These checks now account for much of the remaining work; none was removed
or weakened. Profiling overhead is excluded from the table.

## Numerical acceptance

Fresh/session single points compare total energy, every energy component, and
all forces with absolute tolerance `1e-10` in their respective public units and
relative tolerance `1e-12`. A -> B -> A uses a non-rigid perturbed B frame.
Interleaved sessions and reversed mapping insertion order pass the same checks.

Minimization acceptance uses the unchanged force threshold 0.1 kJ/(mol*angstrom),
500-iteration/2000-call bounds. Both fresh and session runs converge; their final
energies/components/forces are re-evaluated through independent fresh calls using
the above tolerances. Final fmax uses absolute `1e-8`, relative `1e-10`.
No claim of identical optimizer history across environments is required.

NVE fresh/session comparisons cover all stored coordinates, synchronized
velocities, potential/components/forces, kinetic and total energy, termination,
and evaluator counts. Pair tolerances remain `atol=1e-10, rtol=1e-12`.

The independent source-prmtop path remains OpenMM's Amber reader plus explicit
synchronized velocity-Verlet `CustomIntegrator`, using the same authoritative
masses and physical settings as Phase 4E3. Original source masses are retained in
the report; there is no mass repartitioning. The unchanged independent tolerances
are 1e-8 angstrom for positions, 1e-6 angstrom/ps for velocities,
2e-7 kJ/mol (rtol 2e-10) for potential, 1e-8 kJ/mol (rtol 2e-10) for kinetic,
and 2e-6 kJ/(mol*angstrom) (rtol 2e-9) for forces.

### Five archived short trajectories and a longer case

All source files and historical preparation signatures are checked using the
existing archive loader. Initial geometries receive the established seeded
0.025-angstrom perturbation and valid local minimization. No historical record
is re-signed. The five cases are GAFF/GAFF2 AM1-BCC phenol, GAFF2 provided-charge
phenol, capped PE DP=3, and the assigned Cl/Br stereocenter.

Declared before execution:

- Short runs: 20 steps at 0.1 fs (0.002 ps), 22 evaluator calls, five frames.
- Longer GAFF2 AM1-BCC phenol: 1000 steps at 0.1 fs (0.1 ps), 1002 calls,
  11 frames at steps 0, 100, ..., 1000.
- Velocities: `0.2*sin(site_id + 1 + 0.7*axis)` angstrom/ps.
- Absolute energy guard: 0.01 kJ/mol for both lengths.
- No rescaling, adaptive timestep, COM removal or thermostat.

All five short runs completed with two Contexts each. Maximum absolute energy
deviations were 1.47325e-6, 1.50409e-6, 1.49200e-6, 2.06965e-6, and
9.98065e-8 kJ/mol, respectively. Assigned Cl/Br stereo passed.

The 1000-step trajectory completed in 1.478 s with 1002 logical calls and two
Contexts. Fresh/session stored position and velocity differences were zero.
Maximum absolute energy deviation was **2.15702e-6 kJ/mol**. Maximum independent
errors: positions 1.60e-11 angstrom, full-step velocities 1.48e-9 angstrom/ps,
potential 4.52e-10 kJ/mol, kinetic 1.22e-11 kJ/mol, forces
1.67e-8 kJ/(mol*angstrom), all within the predeclared tolerances.
No tolerance was relaxed to obtain these outcomes.

[Acceptance and benchmark manifest](references/phase_4e4/acceptance.json) contains
source hashes, charge provenance, initial states, mass inventories, checkpoint
fingerprints, quantities/errors, Context counts, minimization outcomes, repeated
timing samples and environment details. AM1-BCC phenol retains its 0.002 e
historical tolerance and approximately -0.001 e residual; charges are unchanged.

```sh
python scripts/validate_evaluation_sessions.py \
  --preparation-manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" \
  --output /new/path/sessions.json
pytest -q tests/test_evaluation_sessions_archive.py
python examples/evaluation_session.py
```

The script requires the repository's synthetic fixture and dev/scientific test
dependencies; it refuses to overwrite an existing report. Omit the manifest to
benchmark the synthetic fixture only, with the archive omission explicitly
recorded. No server archive path is hard-coded into a public API.

## Verification and limits

Environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1 Reference/double,
ParmEd 4.3.1, RDKit 2026.03.6, SciPy 1.17.1 for minimization. The existing isolated
validation environment was reused without altering base scientific environments.

Executed verification:

- Ordinary suite: **637 passed, 5 skipped** in 35.65 s. Skips are explicit opt-in
  live AmberTools, archived single-point, minimization, NVE and session acceptance.
- Session archived acceptance separately enabled: **1 passed** in 60.39 s,
  exercising all five cases and the 1000-step trajectory. The benchmark/acceptance
  script also completed independently and produced the retained manifest.
- Ruff: passed. `python -m pip check`: no broken requirements.
- Runnable examples passed: `evaluate_singlepoint.py`, `minimize_geometry.py`,
  `run_nve.py`, `evaluation_session.py`, `build_local_template_polymer.py`,
  and `generate_random_walk.py`.

Numerical and resource-lifecycle tests execute real OpenMM. Mocks are limited to timing/count
instrumentation, forbidden Integrator stepping, and injected failures. Existing
optional-import tests still block OpenMM, ParmEd, RDKit and SciPy together.
No live AmberTools parameterization or regeneration was run for this phase.

Sessions are explicit, thread-owned and non-reentrant, not a worker pool or
restart mechanism. Their physical scope is unchanged: finite nonperiodic,
unconstrained supported Amber models, with the same 100-site preparation limit.
The speedup does not establish thermal equilibration, scientific suitability,
production readiness or long-chain scalability. No NVT/NPT, automatic velocities,
packing, periodic cells, restraints, force-field family or MLIP functionality was
added. `production_validated=False` and
`simulation_readiness="not_established"` remain fixed.

### Numerical input boundary correction

Extreme finite coordinates whose geometry arithmetic overflows now raise
`EvaluationInputError` in both fresh and session paths, chained from the original
numerical exception. Rejection precedes Context position installation and leaves
an existing session usable. Backend failures still invalidate the session.
The regression uses opposite coordinates of magnitude 1e308, checks no positions
were installed, and then successfully evaluates an ordinary frame.
