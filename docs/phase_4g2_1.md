# Phase 4G2.1 — authoritative OPLS binding at consumer entry points

## Repository and reproduced defect

Correction on the existing, unmerged branch
`codex/phase-4g2-oplsaa-parameters-energy`, based on reviewed
`11ec193658c3f2ebdddc0e05a9e346053c81cbcc`. Fetch confirmed main remains
`b61c96628d47e93417e809c175e2f77e5a7486e8`. No merge, main modification,
force-push, historical-record rewrite or unrelated checkpoint change occurred.

The explicitly synthetic H2 source fixture in `test_oplsaa_parameters.py` was
used before changing production code. Both sites were copied and changed to
F (atomic number 9, mass 18.998) while retaining IDs and coordinates. Direct
`validate_system()` raised `EvaluationInputError`, but `minimize_geometry()`
returned `force_converged`, `converged=True`, and
`final_evaluation_verified=True`. This is a software binding failure; neither
the toy H2 spring nor the resulting F2 calculation is scientific OPLS data.

All 27 initial regression cases failed before correction: 20 early-binding
checks (five public routes × graph/element/mass/stereo changes), one valid
checkpoint resume binding check, and six compatible-coordinate/fresh-path
checks. Some malformed stereo inputs already encountered later stereo
validation; they still bypassed the required bound-system check. The tests
specifically require that check before any energy invocation.

## Root cause and minimal correction

Consumers recognize `OpenMMBoundPotential`, not merely the presence of a
`validate_system()` method. OPLS implemented that method but omitted the marker
base class, so it followed the unrestricted analytical-evaluator route.

`OPLSSinglePointEvaluator` now inherits the established marker and implements
`evaluate_fresh()` by calling its existing `evaluate()`. Every such evaluation
constructs and destroys its own OpenMM System/Integrator/Context from the owned
serialized model. No session, cache, integrator step or alternative force path
was added. A real OpenMM construction-count test verifies two contexts for an
ordinary call followed by a fresh call, with identical returned results.

No consumer needed a backend-specific conditional:

| Entry point | Existing contract now applied to OPLS |
|---|---|
| `minimize_geometry()` | Authoritative binding before optimization; fresh final verification |
| `run_nve()` | Binding before initial evaluation/propagation; fresh final verification |
| `run_langevin()` | Binding before initial evaluation/noise; fresh final verification |
| `run_dynamics_segment()` | Binding before execution; fresh startup and final checks, for NVE and BAOAB |
| `resume_dynamics()` | Checkpoint/system compatibility, then bound validation before evaluation; fresh startup and final checks |

Incompatible calls raise the existing `MinimizationInputError` or
`DynamicsInputError`, with the binding failure chained by existing consumer
code. Generic analytical evaluators and Amber marker/session behavior remain
unchanged. Coordinate-only changes remain permitted without reparameterization.

The resume regression first constructs a valid F2 checkpoint with an explicitly
analytical test potential numerically equal to the toy spring. This is legitimate
for an analytical evaluator and preserves exactly the energy/model identities
needed for the test. No checkpoint fields, checksums, masses or compatibility
identities are patched. Resuming with the *bound H2* OPLS evaluator therefore
passes external checkpoint/system compatibility and reaches the intended bound
check, which rejects before its first energy call.

## Verification

- Focused binding suite: **28 passed**, including the additional actual Context
  construction test. The original 27 failing regression cases now pass.
- Complete ordinary suite: **1,063 passed, 10 skipped**. Skips are the existing
  opt-in checks; the OPLS real-source test is executed separately below.
- Pinned Foyer environment: **23 parameter/evaluation tests passed**, including
  the real-source ethanol permutation test with noncontiguous IDs, reversed
  insertion order and stale interaction caches.
- Existing `examples/oplsaa_singlepoint.py`: passed with the original pinned XML.
- Ruff and pip check: passed; pip check ran in both the ordinary and isolated
  pinned Foyer environments (OpenMM 8.6.1 Reference).

The tests exercise actual public entry points and real OpenMM evaluations with
labelled synthetic parameters. Rejection tests count binding calls and forbid
energy calls, preserve caller systems and prior results, and check that a
compatible changed coordinate frame keeps model/parameter identities. Fresh
verification is counted explicitly, including both checks on new segments and
resume. The ordinary suite retains Amber and analytical behavior checks.

## Unchanged scientific scope

No numerical model, source pin, native charge, mixing rule, exclusion, schema,
signature, calculation setting or tolerance changed. Phase 4G2's independent
numerical evidence is reused; no QM or full acceptance matrix was rerun.

Contract compliance allows existing generic consumers to enforce their existing
safety and verification semantics. It does **not** establish scientific OPLS
optimization, dynamics or checkpoint-continuation validation. The unchanged PS
native-charge failure and missing explicit source impropers remain limitations.
No reusable OPLS session or production dynamics capability is claimed.

`production_validated=False` and `simulation_readiness="not_established"`.
