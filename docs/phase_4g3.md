# Phase 4G3 — reusable OPLS sessions and bounded local minimization

## Base and preserved contracts

Branch `codex/phase-4g3-oplsaa-minimization` starts from reviewed, still-unmerged
4G2.1 `c040f175e51167ee0234cbbad083155a1f001e6a`. Fetch showed main at
`b61c96628d47e93417e809c175e2f77e5a7486e8`, without 4G2/4G2.1. This branch
therefore depends on that feature tip. No automatic merge, main modification,
force-push, historical-artifact rewrite or change to Phase 4F5 was made.

The [4G2 source and forms](phase_4g2.md) are unchanged: Foyer revision
`dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`, XML SHA256
`c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`,
native library charges, six-coefficient RB torsions, geometric sigma/epsilon,
shortest-path exclusions and separate 0.5 LJ/Coulomb 1–4 scaling. Historical
parameter, model and evaluation identity algorithms and settings remain intact.
No impropers, guessed parameters, charge normalization or new source were added.

## Public path

```python
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa
from island.minimization import MinimizationOptions, minimize_geometry

source = load_oplsaa_source('/installed/oplsaa.xml')
parameters = parameterize_oplsaa(system, source)
evaluator = OPLSSinglePointEvaluator(system, parameters, source)
options = MinimizationOptions(
    force_tolerance=0.1, max_iterations=5000, max_evaluations=10000,
)
with evaluator.open_session() as session:
    minimum = minimize_geometry(system, session, options)
if minimum.converged and minimum.final_evaluation_verified:
    optimized = minimum.to_system(system)  # validated copy, never an in-place edit
    parameters.validate_integrity(optimized, source)
```

The public session class is
`island.evaluation.oplsaa_session.OPLSEvaluationSession`; opening it through the
validated evaluator is preferred. Retained `OPLSParameterizationResult` records
can be validated and bound directly without invoking Foyer again. OpenMM is
loaded only when needed, SciPy by the existing minimizer, and Foyer only for
new typing/assignment or the independent acceptance reference. The optional
backend remains isolated from core imports.

Run the complete PSMILES example:

```sh
python examples/oplsaa_minimize.py --xml /installed/oplsaa.xml --dp 3
```

It builds an explicit-H local-template PE chain with template/assembly seed
2026, parameterizes, minimizes with a session, reports the force criterion and
independent-verification flag, and applies only a converged result. An unmet
criterion exits nonzero. A converged local minimum is not equilibration or proof
of scientific suitability.

## Session ownership and lifecycle

Each session exclusively owns a deserialized System, unstepped Integrator and
Reference Context. It reuses the existing private `_OpenMMResources` holder.
The OPLS fresh and reusable paths share the exact same geometry checks and
energy/force/result construction. Every evaluation installs the complete frame;
there is no result cache. Identical inputs retain parameter/model/coordinate/
calculation identities, components and forces across both execution paths.

- Explicit `close()` is idempotent; `with` cleans up on success and exceptions.
- Sessions cannot be copied or deep-copied. Separate sessions have distinct
  mutable resources and can be used in interleaved sequences.
- A session belongs to its creating thread. Cross-thread, overlapping calls
  and nested context-manager entry raise `EvaluationSessionBusyError`.
- Closed-session evaluation and binding validation raise
  `EvaluationSessionClosedError`; they do not allocate a replacement Context.
- `evaluate_fresh()` is an explicit exception to that closed-state rule: like
  the Amber session, it remains available after closure/invalidation for
  independent recovery verification through the parent evaluator. It creates
  a new Context and never accesses the session Context or a cached result.
- Geometry (including finite overflow checks and singularity checks) is validated
  before backend mutation. Such input errors leave a healthy session usable.
- Backend failure after installation invalidates and disposes of the session.
  An active failure is preserved if cleanup also raises; the secondary error
  is attached as an exception note. The shared resource constructor gained the
  same protection for failed initialization, without changing normal Amber use.

Caller systems, source records, default coordinates and prior EvaluationResults
remain unchanged. Resources are private; there is no global pool or cache.

## Existing minimizer contract

No optimizer, convergence rule or minimization result schema was introduced.
`OpenMMBoundPotential` validation rejects changed chemical graph, stable IDs,
elements, masses and assigned stereochemistry before the first energy call.
Compatible coordinate changes reuse the bound potential. Sessions delegate this
check to their validated parent.

The engine counts every actual evaluation, including failures and independent
final verification. It reserves final-check capacity when the budget permits;
a one-call budget can only yield an unverified diagnostic. Cached optimizer
requests do not count as additional evaluator calls. A successful session-backed
minimum uses one reusable Context plus a separate final Context, independent
of its iteration count. Context counts are distinct from evaluator calls.

Force convergence uses the maximum atomic force-vector norm in
kJ/(mol*angstrom), not SciPy's success flag. RMS force is
`sqrt(sum_i |F_i|^2 / N)`. Coordinates are angstrom and analytical gradients are
minus force. Existing energy-increase/change tolerances, line-search bounds,
stereo checks, finite-value checks and result integrity/application rules are
unchanged. Unconverged results cannot be applied without explicit opt-in;
malformed/incompatible results never become admissible through that opt-in.

## Declared real-source acceptance

The command has separate exclusive declaration and execution steps:

```sh
python scripts/validate_oplsaa_minimization.py \
  --retained /path/to/phase4g2-energy-final --xml /installed/oplsaa.xml \
  --output /new/phase4g3-run
# Inspect declaration.json before execution:
python scripts/validate_oplsaa_minimization.py \
  --retained /path/to/phase4g2-energy-final --xml /installed/oplsaa.xml \
  --output /new/phase4g3-run --execute-declared
```

The four required systems and parameter envelopes were reused from retained
4G2 output. Parameter identities were checked against that run's outcomes, and
complete parameter integrity was verified against each authoritative system
and pinned XML. The new declaration seals the consumed system/parameter file
SHA256 values, prior outcome-file hash, original coordinate fingerprints and
parameter identities **before** optimization. No new construction, typing or
parameter assignment was run for those four cases. Source hashes were checked
again after execution.

Settings: Reference platform; force tolerance 0.1; 5000 iterations; 10000
energy/force calls; existing `max_line_search_steps=20`,
`energy_change_tolerance=1e-12`, `energy_increase_tolerance=1e-8` kJ/mol.
No retries, seed changes or budget extensions. Each case retains its complete
MinimizationResult, optimization history and accepted copied system. A failed
required case would retain its diagnostics and make CLI acceptance nonzero.

Independent final reference: the 4G2 `Foyer.apply -> ParmEd.createSystem`
geometric-mixing path, with no ISLAND conversion or pair-list helpers. At the
exact returned coordinates it compares total energy, bond/angle/RB/combined
nonbonded components and all stable-site forces. Predetermined absolute tolerance
is 1e-5 kJ/mol for energy and 1e-5 kJ/(mol*A) for force, both with relative
2e-10. The reference constructs a third, separately counted Context.

### Actual results

All four cases returned `force_converged` with independently verified final
states and validated copies. Units: energy kJ/mol; force kJ/(mol*A).

| Case | Initial U | Final U | Final fmax | Final RMS | Iterations | Evaluations |
|---|---:|---:|---:|---:|---:|---:|
| Butane | 49.95322514 | 12.02375143 | 0.09506875 | 0.06336087 | 82 | 89 |
| Ethanol | 17.85531760 | 3.59023158 | 0.09033426 | 0.06555005 | 39 | 42 |
| PE DP3 | 65.93767568 | 15.90352677 | 0.08044012 | 0.05070752 | 126 | 132 |
| PEO DP3 | 143.41163733 | 28.34549540 | 0.08788877 | 0.05236517 | 249 | 263 |

Each used **two minimization Contexts**, plus one independent Foyer/ParmEd
reference Context. Maximum final-reference energy difference was
5.33e-15 kJ/mol; maximum force difference was 1.47e-14 kJ/(mol*A). All component
checks passed. Timings from case start through minimization were respectively
2.276, 0.562, 1.190 and 1.347 s. These include retained-input reconstruction,
validation, evaluator binding and first-use imports; they are not optimizer-only
benchmarks or a promised speedup.

PS DP3 was separately built with the predeclared 4G2 local-template recipe and
rejected by native charge validation at -0.23 e, before minimization. No PS
charge or source repair was attempted. PE50 convergence was not requested or run.

Machine-readable declaration and outcomes are in `phase_4g3_acceptance.json`;
full new artifacts were retained in `island-validation/phase4g3-minimization`.
The latter is an execution location, not a hard-coded API path. Historical
4G2 artifacts are unchanged.

## Verification and limitations

Synthetic session tests use the labelled H2 harmonic software fixture with
real OpenMM. They exercise A→B→A frames, interleaved sessions, immutable results,
closed/thread/reentrant/copy rules, no Integrator steps, pre-mutation invalid
geometry rejection, mismatched binding, actual Context counts, final-call
reservation, backend failure recovery, cleanup failure chaining and application
of converged versus diagnostic results. They do not count as real-source
scientific evidence.

The real acceptance used the pinned Foyer stack, OpenMM 8.6.1 Reference,
ParmEd 4.3.1 and NumPy 2.4.6; its exact SciPy and dependency versions are recorded
in the declaration. Complete ordinary tests, Ruff, pip check and affected
examples are recorded in final verification below.

A finite local minimum is neither an equilibrium ensemble nor an accuracy
assessment. The pinned source's absent explicit impropers and PS neutrality
failure remain. No dynamics, periodic systems, packing, PCFF, crosslinking or
new workflow/checkpoint schemas were implemented. Both readiness flags remain
`production_validated=False`, `simulation_readiness="not_established"`.

### Final verification

- **1,081 ordinary tests passed; 10 expected opt-in skips.** No Amber/archive
  generation or QM acceptance was needed for this phase.
- **18 new OPLS session tests passed** (real OpenMM, labelled synthetic inputs).
  The combined OPLS/Amber session, binding and parameter test run passed 90 tests
  with one opt-in skip before the two initialization-cleanup regressions were
  added; those two are included in the final full suite and new-session count.
- **23 parameter tests passed in the pinned Foyer environment**, including the
  actual-source stable-ID/insertion-order permutation case.
- Ruff and pip check passed; pip check ran in both the ordinary and isolated
  pinned-source environments.
- New OPLS minimization example, existing OPLS single-point example, and Amber
  minimization/session examples passed. The new example reproduced the accepted PE DP3 result (126 iterations,
  132 calls, fmax 0.08044012) and returned a validated copied system.
- Real-source minimization acceptance: four required passes, PS expected rejection;
  no failed case was retried or replaced, and all original consumed hashes
  remained unchanged.
