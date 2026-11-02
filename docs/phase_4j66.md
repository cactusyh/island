# Phase 4J66 — HMDS constrained feasibility and block identification

The unchanged 58-parameter HMDS model **can satisfy all four frozen training
limits simultaneously**. J65's least-squares result was an objective tradeoff,
not a demonstration of model/data/bounds incompatibility. A constrained trial
was executed, but **zero candidates qualify**: its optimum approaches a bound,
solver accuracy remains explicitly qualified, and block-resampling robustness
fails. The three reserved validation configurations remain untouched.

This study executed **zero new QM calculations and zero research MD steps**.
It uses the exact 27-atom, neutral-singlet, explicit-H HMDS graph
`d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`.
Si–H remains stopped; no capped-DP2 calculation was initiated.

## Base, reconstruction and frozen contract

Fetched `origin/main` is `712e588a8ca9b3c0d3ba03114a51064d3a2a83bb`, and
contains required J65 `1b03f835718599b7f6e53e3331400518112882ec` (PR88).
The J66 branch starts at that fetched main. Main was neither merged into the
branch nor pushed. Unrelated untracked files were preserved.

J65's native archive, producer identities, exact protocol reuse, decimal
energies, ASCII gradients, unit conversions and complete fit were replayed
in a relocated directory before new optimization. Reconstruction agrees with
all four reported J65 metrics to better than 1e−9 in their physical units:
energy RMSE0.0104172274, force RMS0.3751669382, maximum energy error0.0343103913,
and maximum local vector-force error4.2905539501. The 22 configurations, all
58 physical columns, fixed MM contributions, 1804×59 weighted matrix, parameter
names/units/bounds and all64 original block draws are bound in
[reconstruction.json](evidence/phase_4j66/reconstruction.json). No producer
was recreated for a reused observation. `parameter-provenance.json` exposes
every source-query/site binding, units, equilibrium dependency, unchanged charge
identity and equation authority for the58 physical variables. All historical files remain unchanged.

The [scientific plan](evidence/phase_4j66/plan/scientific-plan.json) was frozen
before solving. The target remains PBE0/def2-TZVP, DF/def2-universal-jkfit,
unpruned250×1454 grid, SCF1e−12/1e−10, integral threshold1e−12, without dispersion.
J65's original least-squares fit and failed gates are unchanged.

## Executed conic problem

There are58 unchanged physical coefficients and one global energy-only
intercept. No source equilibria, K3/K4, charges, nonbonded settings, bounds,
observations or acceptance limits changed. With residuals MM−QM, the problem is

```
minimize rho >= 0
subject to original physical parameter bounds
  norm(rE)       <= 0.1 * sqrt(22) * rho          [kcal/mol]
  abs(rE[i])    <= 0.3 * rho                    [kcal/mol]
  norm(rF)      <= 0.5 * sqrt(3*22*27) * rho    [kcal/mol/angstrom]
  norm(rF[i,j]) <= 2.0 * rho                    [kcal/mol/angstrom]
```

The local norm always includes all three components of one atom's error.
The global force RMS includes all1782 components. No atom or configuration was
removed. The intercept never affects forces or shifts configurations separately.

The direct Clarabel0.11.1 implementation uses column norms of the actual weighted
matrix, including the intercept. Residual cones are divided by their original
physical limits, and bound inequalities by their widths. NumPy2.4.6 and
SciPy1.17.1 are pinned in a separate solver environment. Clarabel and independent
SCS3.2.9 backends are lazy optional research imports; the ordinary package/API
has no new mandatory dependency. One worker/BLAS thread,120seconds per solve,
4GiB and6hours campaign limits were declared. Canonical sparse matrices, full
primal/dual/slack vectors, settings, scaling, statuses, residuals, gaps and
runtimes are retained. SCS is a distinct solver, not another Clarabel frontend.

The Clarabel rho run ends `AlmostSolved`, with rho0.99540328094 and absolute
gap8.41e−11. Independent physical-array and complete-MM evaluation establish a
**feasible witness**, without claiming exact optimality from that status:

| Error | Witness | Original limit |
| --- | ---: | ---: |
| Energy RMSE, kcal/mol | 0.0995403281 | 0.1 |
| Force RMS, kcal/mol/Å | 0.4977016405 | 0.5 |
| Maximum absolute energy, kcal/mol | 0.2986209843 | 0.3 |
| Maximum local vector force, kcal/mol/Å | 1.9908065619 | 2.0 |

SCS returns rho0.9954030934. Its raw coefficients violate one bound by3.79e−7;
that raw result is not called a physical witness. Explicit projection onto the
original bounds, followed by independent reevaluation, supplies a second
feasible witness. No performance threshold was relaxed. A further rho run with
the inherited1e−6 bound-activity margin removed from both sides of each interval
also yields a feasible witness with **no original active bound**. Bounds are
therefore not an obstacle to training feasibility. A diagnostic without bounds
gives rho0.9924977144; those outside-bound values are not eligible coefficients.
The force-only diagnostic is explicitly partial and is not a complete-fit claim.

The independently checked dual/support bound includes every stationarity
residual. A feasible primal upper bound bounds rho, and an energy row then
bounds the otherwise unbounded intercept. The dual is projected into its
self-dual cones; minimizing its affine Lagrangian over this derived finite box
and subtracting a conservative arithmetic allowance produces the retained
numerical lower bound. No unbounded-variable residual is ignored. This gives
the interval0.9954030738–0.9954042810 for the rho optimum and, separately,
certificates for failed resamples.

## Constrained least squares and numerical qualifications

The second problem minimizes exactly J65's original objective,
`mean(rE²)/0.1² + mean(rF²)/0.5²`, subject to all four limits and original bounds.
No regularization or duplicate gradient-derived curvature measurement is added.
Full nominal rank59 makes this objective strictly convex. The deterministic
scaling and solution selection are recorded.

The initial direct quadratic formulation produced a feasible trial but a less
useful absolute objective gap because it included a large additive constant.
Its complete original state remains retained. An equivalent conic epigraph,
`norm([2*weighted_residual, t-1]) <= t+1`, minimizes `t` without subtracting that
constant. This amendment changes neither the objective nor physical problem.
Its raw `AlmostSolved` result exceeds the local limit by6.85e−11. The strict gate
correctly fails. Before any bootstrap solve, an additive numerical-selection
contract fixed the first feasible convex blend toward the verified rho witness
from a finite alpha sequence. The selected alpha1e−8 increases the objective by
2.68e−9 and passes the original limits. Raw failures/statuses remain visible;
this is a verified near-optimal feasible trial, not an unqualified exact optimum.

| Metric | J65 nominal | J66 constrained trial |
| --- | ---: | ---: |
| Energy RMSE, kcal/mol | 0.01041723 | 0.07839563 |
| Force RMS, kcal/mol/Å | 0.37516694 | 0.49999999996 |
| Maximum energy error, kcal/mol | 0.03431039 | 0.25442484 |
| Maximum local vector force, kcal/mol/Å | 4.29055395 | 1.99999999997 |
| Original weighted objective | 0.57385279 | 1.61458750 |

The constrained trial redistributes residuals: primary force RMS rises from
0.17201748 to0.35652584 and joint force RMS from0.72020482 to0.81720481, while
the joint local maximum falls below2. All per-configuration signed energies,
Cartesian errors, all atom-vector norms, bound distances and coefficient changes
are retained in `constrained-epigraph.json`. The central tradeoff concerns local
maxima versus the average objective; it does not imply improved physical accuracy.

`middle_bond-torsion_3:hc,c3,sio,osi:block0:n1` approaches its−100kcal/mol/Å
bound. This alone prevents inherited candidate eligibility, even though a
separately verified interior feasible vector exists. The constrained trial's
Si–C/C–H quadratic corrections are18.38588277 and30.39856170kcal/mol/Å².
Complete-model unknown-null path curvatures are411.48745677 and781.43960720,
against the selected QM protocol's425.81581623 and762.51145542. These correlated
training diagnostics are not extra independent fitting observations or bare
source bond constants. No equilibrium or higher-order row was released.

## Block robustness and actionable information requirement

Every exact J65 block draw is retained, including repetitions and omitted blocks.
J65 has48 bound-active draws and3 condition failures; all3 overlap with bound
failures (replicates5,51,59), leaving16 eligible. Nominal full rank does not
establish robustness under these changes in information.

J66 executes64 rho problems and56 conditional constrained fits on those same
selections. Eight sampled-weight problems have independently checked rho lower
bounds1.00687556–1.16033599 and are incompatible with their own reweighted global
limits. This does **not** contradict feasibility for the original22-configuration
weighting. Their certificates, support-error allowances and limiting arrays are
retained. No failed draw was replaced.

There are50 bound-active fits,3 conditioning failures, and30 conditional fits
with `AlmostSolved` status. Categories overlap. Only2/64 qualify under the frozen
accurate-status, rank/condition, bound and physical checks:3.125%, far below80%.
Even ignoring the status failures would leave only6/64 without the other failures;
the uncertainty conclusion is not created solely by conservative status handling.
Survivor-conditional energy/force prediction RMS95 values are0.09579360 and
0.50265629, with force also above its0.5 limit. Two-survivor intervals and
correlations are explicitly unqualified; they are not credible physical
confidence intervals. Absolute intervals and prediction errors are retained,
including coefficients near zero. No basis/model/physical uncertainty is included.

The most recurrent J66 bound direction is the `hc,c3,sio,osi` middle-bond first
harmonic (38draws); the corresponding `c3,sio,c3,hc` direction occurs12times.
The weakest nominal scaled combination mixes the energy intercept, angle-angle
coefficients and third-harmonic end/middle terms. Removing CH, joint-040 or
joint-044 raises J65 bound-failure fractions from approximately63–67% when present
to88–92% when omitted. These correlated inclusion associations are not isolated
causal experiments. Replicates5,51,59 lose several joint/angular constraints
and reach scaled conditions15386,10339,17372 respectively. The complete
multiplicities, weak vectors, conditional bound-combination sensitivities and
energy/local-force prediction effects identify which observations constrain them.

The bounded geometry-only proposal selects `joint-032`, `joint-002`, `joint-024`
and `joint-034` from the existing graph-controlled pool using information in
measured weak directions, not QM outcomes. It includes explicit actions,
coordinates and parameter sensitivities. Proper graph-aware alignment excludes
all existing training configurations and the reserved validation geometries,
covering methyl/H/whole-end symmetries. Its diagnostic augmented condition drops
from1556.02 to864.84. This is targeted redundancy to investigate fragile
combinations, not proof that four new observations will pass uncertainty.
**None of these proposed QM calculations was executed.** A subsequent phase
must freeze a bounded sequential training-information experiment and assess
its block robustness before activating the still-reserved validation set.
Adding easy points must not dilute the existing local errors; future trials
must retain the original22 configurations and their limits as a separate
protected subset. These proposals are deformations of the existing center,
not independently sampled conformational basins. A separate geometry-only
block stress check retains each old multiset and appends all four proposed
feature blocks; its conditions are descriptive information checks, not a new
bootstrap distribution or uncertainty qualification. In this guaranteed-addition
stress check, the worst retained-draw condition decreases from17372.05 to1742.43;
future resampling must also test omission of those new observations.

## Verification, replay and scope

Synthetic conic tests are explicitly separate from native QM evidence. They
cover feasible recovery, a least-squares local failure with a feasible constrained
solution, checked incompatible constraints, inaccurate statuses, vector norms,
RMS normalization, scaling, intercept handling, bounds/interior restrictions,
rank/conditioning, failed-draw retention, graph/protocol/matrix/units/bounds/order
mutations, validation nonactivation and fresh-process offline relocation.

The supported full suite and J54–J65 regression tests, historical receipt checks,
Ruff, compileall, applicable environment pip checks and diff/manifest checks are
recorded under `validation/`. The terminal suite has5265unique cases:
5204passed and61skipped, including25J66 cases and770J54–J65 regressions.
It passed after two explicitly retained interrupted preliminary layouts;
focused reruns add no distinct coverage. Historical receipts are checked against their
explicit original implementation snapshots; J63's old broad module glob is not
allowed to redefine historical implementation by finding later modules. J66
uses a fixed explicit implementation-file list. Terminal JUnit coverage counts
each test once; focused reruns are retained without adding overlapping coverage.
Original analysis/evaluator and numerical-selection failures are retained
separately from corrections. The initial evaluator-interface exception occurred
after a rho solve whose vector had not yet been captured; its exact problem
and trace remain, and that lost vector is not invented. Subsequent raw solver
states are captured before independent evaluation. A staged diff check also caught whitespace in captured third-party CLI help;
the exact raw diagnostics were archived and inventoried, and the retry passed
without editing those bytes. See `validation/attempt-history.json`.

Offline replay needs the repository evidence and NumPy/SciPy, not a QM executable
or network. A separate environment with `clarabel==0.11.1`, `scs==3.2.9`,
`numpy==2.4.6`, `scipy==1.17.1` executes the synthetic solver tests and optional
new optimization. Existing native observations retain their original producers.

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/replay_phase_4j66.py
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/validate_phase_4j66_receipt.py docs/evidence/phase_4j66
```

Software verification, training feasibility and candidate qualification are
separate outcomes. Training feasibility is established; candidate qualification,
untouched validation, physical validation and chemical transfer are not.
Historical failed/unestablished basis qualifications remain unchanged.
No native-PCFF equivalence, complete-basis accuracy, condensed-phase accuracy,
Si–H transfer, capped-DP2 coverage or complete PDMS parameter assignment follows.
The smallest methyl-domain extension beyond HMDS remains a separately audited
capped-DP2 interior graph/charge/source inventory and data for
`sio(osi,osi,c3,c3)`, O–Si–O and Si–O–Si–O; none is supplied by this pilot.

The28 confirmations remain user-attested; original author artifacts remain
`not_yet_supplied`. Ordinary `forcefield="pcff"` behavior is unchanged. There is
no automatic registration, complete PCFF-PDMS profile or production MD.
`production_validated=False`; `simulation_readiness="not_established"`.
