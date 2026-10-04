# Phase 4H6 — reusable PCFF sessions and bounded minimization

## Base and scope

Branch `codex/phase-4h6-pcff-minimization` starts at updated `origin/main`
`b799fcff07c821117941c034cd8ef4a31aeb5cfd`. Main's tree exactly matches reviewed
H5 `e602ebfb0f15927c3e97a01aa7ce46f3f8ac1b0a`:
`6407dc7d53ad64dd8007a280a0464f206df1cf79`. No prerequisite merge was needed.

The H4/H5 physical model, source pins, native charges, explicit pair policy,
all coefficients/dependencies, units, electrostatic constant, and identity
formulas are unchanged. H1–H5 records are not rewritten. This phase adds
resource reuse and acceptance through the existing generic minimizer. It does
not add an optimizer, dynamics, new chemistry, periodicity, or readiness claims.

## API and lifecycle

```python
from island.evaluation import PCFFSinglePointEvaluator
from island.minimization import MinimizationOptions, minimize_geometry

evaluator = PCFFSinglePointEvaluator(system, validated_model_specification)
options = MinimizationOptions(
    force_tolerance=0.1, max_iterations=5000, max_evaluations=10000,
)
with evaluator.open_session() as session:
    minimum = minimize_geometry(system, session, options)
if minimum.converged and minimum.final_evaluation_verified:
    optimized = minimum.to_system(system)  # validated copy
```

`PCFFEvaluationSession` follows the established OPLS lifecycle. It exclusively
owns one deserialized System, Integrator, and Reference Context. It never steps
the integrator, caches an energy/result, pools Contexts, or shares mutable
backend state with another session. Every evaluation replaces all positions.
Fresh and session paths now share `_coordinates()` and `_evaluate_context()`;
result components and H5 fingerprints are computed by the same code.

- Sessions belong to their creating thread. Cross-thread, concurrent/reentrant
  operations and repeated context-manager entry raise
  `EvaluationSessionBusyError`.
- Copy/deepcopy raises `EvaluationInputError`.
- Context-manager exit disposes resources; `close()` is idempotent.
- `evaluate()`/`validate_system()` after closure raise
  `EvaluationSessionClosedError` without allocating a Context.
- Bad coverage, booleans, nonfinite/overflow coordinates, wrong units, and
  singular geometry are rejected before backend mutation. A healthy session
  remains usable after these input failures.
- A backend failure during installation/extraction closes and invalidates the
  session. Secondary cleanup failures are attached to the original exception.
- `evaluate_fresh()` always delegates to the parent's independent Context path.
  It intentionally remains available after closure/invalidation for recovery
  verification; thread/overlap restrictions still apply.

The evaluator continues to own its model and default coordinates. Session use
cannot mutate the caller system or earlier `EvaluationResult` objects. No FRC
parsing or H3/H4 record reconstruction occurs within optimizer evaluations.
OpenMM is loaded only for backend use; SciPy only by the generic minimizer.

## Minimization contracts

Both evaluator and session implement `OpenMMBoundPotential`; the existing
engine calls binding validation before its first energy evaluation. Graph,
element, mass, stable-ID, and assigned-stereo changes are rejected. Coordinate
replacement alone remains compatible.

No generic minimization code changed. Force convergence is the **maximum atomic
force-vector norm**, not SciPy's success flag or a small energy change. Every
invocation, including failed calls and independent final verification, counts
against the budget. The engine reserves the final call and retains its existing
best-valid-state/history behavior. Result integrity and coordinate application
remain non-destructive.

An unverified/unconverged diagnostic cannot be applied through default
`to_system()`. The established explicit `allow_unconverged=True` diagnostic path
is preserved: it labels coordinates `local_minimization_diagnostic`, retaining
`converged=False` and the verification flag. It does not turn a diagnostic into
a verified minimum. Tests cover this distinction, backend recovery, insufficient
budgets, and failed independent verification.

## Frozen real-source experiment

All H5 evidence-manifest files were hash-checked before and after execution.
The four starting frames are the retained `p0-f0/coordinates.json` records;
external chemical systems are reconstructed from the graph in each validated
H5 parameter record. No construction, embedding, typing, charge assignment,
parameter assignment, or starting-coordinate replacement was run for acceptance.
The reconstructed starting systems are explicitly persisted. They contain the
retained authoritative chemical graph; absent builder metadata is not invented.
Original H5 signatures and files remain unchanged.

A declaration was published before optimization in
`../island-validation/phase4h6-minimization/declaration.json`, binding input
hashes, source/executable identity, versions, options, comparisons, and no retries.

Settings for every case:

- Reference, finite/nonperiodic, unchanged native charges;
- caller-selected LJ=(0,0,1), Coulomb=(0,0,1), not universal PCFF defaults;
- force tolerance 0.1 kJ/(mol·Å), maximum 5,000 iterations/10,000 evaluations;
- existing line-search limit 20, energy-change tolerance 1e-12,
  allowed energy increase 1e-8 kJ/mol;
- LAMMPS comparisons: energy/force atol=1e-5 in kJ/mol and kJ/(mol·Å),
  rtol=2e-10; no tolerance or budget adjustment after results.

The independent reference reuses H5's raw converter data and reconstructs the
AA equilibrium-role overrides from its own angle inventory. Raw data remain
unchanged. See [H5](phase_4h5.md) for the ABC/CBD/ABD converter discrepancy and
ABC/ABD/CBD executable convention. H6 records the explicit overrides, LAMMPS
inputs, logs, full final coordinates/forces, and checksums. It does not generate
the reference from ISLAND's compiled term list.

Pinned LAMMPS revision: `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
Executable SHA256:
`ee945ead72f4997905ea42636bc987f3cc4a2a9cf89ab2af30d8e02c4c1cfebc`.
FRC SHA256:
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

## Actual outcomes

All four terminated `force_converged`, with fresh verification true, independent
LAMMPS maximum force below 0.1, and decreasing total potential energy.

| Case | Initial → final energy (kJ/mol) | Final fmax | Final RMS force | Iterations / calls |
|---|---:|---:|---:|---:|
| Butane | 15.800644807 → −23.097559198 | 0.07603609 | 0.04711809 | 104 / 109 |
| Ethanol | −12.089209366 → −31.964244819 | 0.04243823 | 0.02480792 | 51 / 56 |
| PE DP3 | 6.847369714 → −38.687915562 | 0.07347116 | 0.04509578 | 143 / 152 |
| PEO DP3 | 71.380601230 → −12.235700679 | 0.08904442 | 0.04665491 | 301 / 311 |

Forces are kJ/(mol·Å). RMS is the root mean square of atomic vector norms.
LAMMPS recomputed both metrics independently. All total/grouped energies and
all force components passed, with maxima:

| Case | Energy/component error (kJ/mol) | Force error (kJ/(mol·Å)) |
|---|---:|---:|
| Butane | 1.24e-14 | 6.29e-13 |
| Ethanol | 2.49e-14 | 6.13e-13 |
| PE DP3 | 2.16e-14 | 1.99e-12 |
| PEO DP3 | 7.64e-14 | 2.28e-12 |

Each minimization constructed exactly **two Contexts**: one session plus one
fresh final verification. There were zero extra diagnostic evaluator calls.
LAMMPS is a separate operation, not part of that call count.

| Case | Bind/validate/compile (s) | Session construction (s) | Optimization including fresh final check (s) |
|---|---:|---:|---:|
| Butane | 2.819 | 0.122 | 1.392 |
| Ethanol | 2.670 | 0.120 | 0.283 |
| PE DP3 | 3.594 | 0.127 | 1.051 |
| PEO DP3 | 3.626 | 0.121 | 2.186 |

These are single-run measurements, not a claimed speedup. LAMMPS checks took
0.026–0.029 s per case. Full results, optimization histories, identity links,
commands, hashes, and versions are in the retained directory and
[evidence/phase_4h6.json](evidence/phase_4h6.json).

## Reproduction and verification

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY scripts/validate_pcff_minimization.py --source "$FRC" \
  --inputs ../island-validation/phase4h5-whole-system-final \
  --lammps ../island-validation/phase4h4-upstream/build/lmp \
  --output /new/exclusive/acceptance-directory
$PY examples/pcff_minimize.py --source "$FRC" --lj 0 0 1 --coulomb 0 0 1
```

The acceptance command refuses existing output directories, missing/changed
retained inputs, wrong model policy, or a different reference executable. It
never regenerates missing acceptance inputs. Failed minima retain diagnostic
records; the aggregate exit status is nonzero if any required gate fails.
The example independently demonstrates construction through validated copied
application; its run is separate from retained-input acceptance.

Regression fixtures are explicitly synthetic, executed with real OpenMM:
alternating nonrigid frames, interleaved sessions, unchanged fingerprints and
components, copying/thread/overlap/closed rejection, input errors before backend
mutation, Context counts, backend invalidation/fresh recovery, verification
failure, budgets, graph/mass/stereo binding, caller ownership, diagnostic and
successful application, cleanup exception chaining, lazy imports, and no
integrator stepping. Historical H4/H5 tests remain unchanged. No generic-engine
defect was found; a test assumption about opt-in diagnostic application was
corrected to match the existing contract.

Final verification: **1,242 passed, 10 opt-in skipped** (189.27 s); **21 new
focused session tests passed**. Ruff and `python -m pip check` passed. Both
`pcff_minimize.py` and the unchanged H5 `pcff_singlepoint.py` example ran
successfully. All 172 H5 artifact hashes were rechecked unchanged. Outcomes
and log hashes are recorded in the evidence manifest. Existing unrelated archive/QM opt-ins are distinct from the
executed pinned-source PCFF acceptance. PE DP50's H5 single-point evidence is
preserved; no PE DP50 minimization was attempted.

## Next gate

The bounded minimization/resource/numerical gates passed. Bounded PCFF dynamics
validation can begin next, with its own integration, lifecycle, continuation,
and independent-reference acceptance. No dynamics claim is made here.
A converged local minimum is not equilibration or scientific validation.
`production_validated=False`, `simulation_readiness="not_established"`.
