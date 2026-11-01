# Phase 4J65 — HMDS cross-site controls and protocol-specific coupled fitting

J65 executed **27 native gradients** and the complete **58-parameter trial fit**.
All numerical eligibility checks passed. The selected objective meets the global
energy/force RMSE targets, but fails the frozen maximum local force criterion
and block-uncertainty qualification. **Zero candidates are accepted**, and the
three reserved validation configurations remain untouched. Historical failed
basis qualifications remain unchanged.

This experiment concerns the exact 27-atom, explicit-H, neutral-singlet HMDS
graph `d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`.
The Si–H investigation remains stopped. The ordinary API remains
`forcefield="pcff"`; no research candidate is automatically assigned.

## Base and reconstruction

Fetched main is `76db72c0ab4c83f6af711d897b33761efbe12e6a`, containing J64
`71777122f12179f6c5ce2ca6828d55d92a2fe635` (PR87). This exact fetched main is
the base of `codex/phase-4j65-hmds-protocol-specific-joint-fit`. Main is neither
merged into this branch nor pushed by this task. Unrelated local files remain
untouched.

Before new QM, the J64 receipt validator reconstructed its retained native
gradients, path responses, conditional corrections and 58-column design. Its
receipt SHA-256 remains
`92688de1ef4302b6321d39083ec6a7b989d8a99847674a6633d6c22674cc3120`.
Historical source rows, charges, native files, identities and failed scientific
gates remain unchanged. The 28 convention confirmations remain user-attested;
original author artifacts are `not_yet_supplied`.

## Frozen experiment

The [scientific plan](evidence/phase_4j65/plan/scientific-plan.json) SHA-256 is
`93043b3bf43ca1c5352d94d2a4822104e6b0ebd7202ba96e697bf0ebb8c9af88`.
It fixes the protocols, exact requests, controls, selection, bounds, objective,
uncertainty and activation rules before new results. Analysis source bindings
were additionally frozen before any new training observation.

The fitting target is **PBE0/def2-TZVP**, DF/def2-universal-jkfit, unpruned
250×1454 quadrature, energy/density convergence 1e−12/1e−10, integral threshold
1e−12, without dispersion. This is a cost-informed method-specific research
target. Successful reproduction would not establish complete-basis accuracy,
physical validation or chemical transfer. J64's complete-force-fit authorization
remains false; J65 uses a distinct experiment contract.

The initial cap is **32 counted native attempts**, including startup failures,
interruptions and retries. Thirty exact requests are prepared: five required
QZVP controls; ten missing primary training points; five conditional joint
points; seven numerical controls; and three conditional validation points.
Two attempts remain unallocated. There is no implicit Hessian, large scan or
additional parameter family.

There is one worker and four threads. The cgroup has 32 GiB and eight CPU
equivalents; its anonymous and reclaimable memory are inventoried separately.
QZVP uses 8 GiB configured memory, a separately justified 12 GiB process allowance
and 7,200-second wall limit. The allocation uses J64's measured approximately
9.5 GiB peak RSS and 682-second center timing, while allowing for the larger
unpruned grid. QZVP runs separately. TZVP uses 4 GiB and 1,800 seconds per job.
Scratch is capped at 64 GiB, and campaign wall time at 36,000 seconds.

The unchanged J63 native capture and J64 bounded process/resource wrapper retain
exact requests, native output, full gradients, derived forces, software/settings,
hashes, process commands, return codes, stdout/stderr and failed attempts. Resume
requires an exact verified completed request; failures require additive retry
identities and remaining budget. Scratch is external to Git.

## Prediction-only QZVP controls

The four stretching requests use Si–C(6–7) and C–H(7–19), each ±0.003 Å, at the
exact retained center. Basis, auxiliary, grid 150×974 and SCF1e−10/1e−8 match the
J64 QZVP calibration. Before execution, the frozen corrections predict total
gradient-derived curvatures of **421.14141338** and **754.92913144 kcal/mol/Å²**.
Their verified fixed total-MM curvatures are 372.13327596 and 720.51116272,
respectively. All 56 unknown coupling second-path derivatives vanish on these
graph-fragment stretches. No unknown coefficients are defaulted to zero, and
the two corrections are not refitted on these observations.

HMDS graph-aware proper alignment includes methyl/H/whole-end permutations.
The new and original displaced paths pass equivalence at 0.001 Å but fail at
0.000001 Å. They therefore supply near-symmetry cross-site controls in the same
conformer, with small local geometric differences; they do not provide an
independent conformational basin or independent methyl-environment transfer test.

The fifth control changes only the QZVP center integration grid to 250×1454.
Its geometry, SCF thresholds, auxiliary basis and other settings match the
retained QZVP150×974 center exactly. Atom-wise vector-force differences are
reported separately from the historical common-grid orbital-basis contrast.
Neither this point nor a successful SCF establishes a QZVP stationary point or
all-path basis convergence.

## Executed control observations

The matched QZVP150×974/250×1454 center contrast changes the Cartesian-force
RMS by **0.000337882 kcal/mol/Å**, with maximum local vector difference
**0.002807818** at oxygen stable ID5. The electronic-energy change is
−0.000107342 kcal/mol. Independent ASCII-gradient reconstruction agrees within
6×10−14 kcal/mol/Å. Configured memory is 8 GiB and measured peak RSS8,228,444 KiB.
This small center-only grid response is distinct from the retained0.535341
TZVP/QZVP common-grid force RMS shift; it supplies neither all-path basis
qualification nor a QZVP stationary point.

For C–H(7–19), the observed QZVP gradient-derived total path curvature is
**754.94580536 kcal/mol/Å²**, versus the preregistered **754.92913144** prediction.
The difference is 0.01667392 kcal/mol/Å². The J64 correction remains unchanged.
This is a same-conformer, near-symmetry cross-site check, not an independently
sampled chemical environment.

All five required controls completed without a native failure or retry. At
Si–C(6–7), the observed total curvature is **420.92938284 kcal/mol/Å²**, versus
**421.14141338** predicted, a difference of−0.21203054 (−0.0503466%). The C–H
relative difference is+0.00220867%. These numerical prediction errors preserve
the conditional nature of the original corrections; they neither fit new
parameters nor establish complete-force, basis or transfer validation. The
full native values, independent chain-rule reconstruction and atom-wise grid
contrast are in `targeted-controls.json`.

## Exact-protocol reuse and numerical eligibility

Seven of the 17 primary geometries have exactly matching J64 observations:
the bitwise exact center, Si–C/C–H ±0.003 Å, and siloxane ±0.02 rad. Their original
native job identities, producer metadata, result checksums and parent archive
references are preserved. They acquire no new execution record. The ten remaining
primary observations require new dense-TZVP calculations; QZVP/TZVPP or weaker
SCF/grid values are never pooled into this fitting dataset.

All observations must pass native convergence, graph/order/charge/multiplicity,
units, hash and mapping checks. Native decimal energies are subtracted before
conversion; independent line and ASCII parsers reconstruct energies and forces.
Inherited energy/gradient path-curvature checks gate joint-data activation.

The numerical controls compare center and joint-040 at 250×1454 versus 300×1454
with identical SCF, plus a center-only SCF1e−13/1e−11 comparison. Four additional
target-protocol gradients at joint-040 ±0.001 and ±0.0005 Å along a frozen unit
Cartesian direction check energy/gradient consistency. All full Cartesian
responses remain available. Floors, relative criteria and the direction are
frozen in the plan. These are bounded numerical checks for the declared target;
historical orbital-basis failures remain separate.

## Primary numerical results

All 17 primary observations were assembled from 10 new dense-TZVP calculations
and 7 exact-protocol J64 executions. All 8 preregistered energy/gradient
path-curvature checks pass. The largest ratio of error to its frozen allowance
is 0.0472013 (siloxane0.04 rad); no tolerance was changed.

| Path | Step | Energy curvature | Gradient curvature |
| --- | ---: | ---: | ---: |
| Si–C | 0.003 Å | 425.785288 | 425.815816 |
| Si–C | 0.006 Å | 425.754610 | 425.859659 |
| C–H | 0.003 Å | 762.527587 | 762.511455 |
| C–H | 0.006 Å | 762.576599 | 762.596060 |
| Methyl rotation | 0.02 rad | 6.640408 | 6.636792 |
| Methyl rotation | 0.04 rad | 6.630316 | 6.625211 |
| Siloxane rotation | 0.02 rad | 1.434301 | 1.432301 |
| Siloxane rotation | 0.04 rad | 1.427494 | 1.425204 |

Curvature units are kcal/mol/Å² for stretches and kcal/mol/rad² for rotations.
The rotational quantities include the curved-path gradient/acceleration term.
They are total molecular responses, not isolated coefficients. Passing these
checks activates the 5 frozen joint geometries, without promoting historical
basis or physical-validation claims.

## Complete coupled trial

The Hamiltonian is the existing
`island_research_class2_hmds_j64_quadratic_v1`: all verified fixed bonded and
nonbonded contributions, 56 explicit missing coupling coefficients, and additive
Si–C/C–H quadratic corrections. Source equilibria, K3/K4, charges, torsions,
pair conventions and other fixed rows remain unchanged. Its relationship is

```
QM = verified fixed MM + explicit 58-parameter contribution + discrepancy
```

The complete physical evaluator requires every missing coefficient explicitly.
The trial is linear in its 58 physical variables, with one additional unbounded
training energy intercept. The exact objective is

```
mean(energy residual²) / (0.1 kcal/mol)²
  + mean(all Cartesian force-component residual²) / (0.5 kcal/mol/Å)²
```

The intercept is fitted once on training; validation energies retain it without
per-conformer shifts. Source baselines and offsets are included in the actual
energy/force Jacobian. Column scaling includes the offset. On the 22 selected
geometries its geometry-only augmented condition number is approximately 1556,
distinct from J64's987 screening metric. The frozen experiment screen requires
rank 59 and scaled condition≤10000, with all singular values and correlations
reported. That screen is not a physical identifiability certificate.

Bounds and units are explicit for every parameter. Quadratic corrections keep
total K2 between0.5 and 2 times its source value. Coupling envelopes are broad
local research bounds comparable to source mixed-stiffness scales; harmonic
bond/angle-torsion amplitudes are capped at 100 in their stated units. These
limits are not a global stability certificate. Any active bound prevents
candidate eligibility. The optimizer is bounded TRF linear least squares,
tolerance1e−12, maximum 1000 iterations, without regularization.

The unchanged target gates are energy RMSE≤0.1 kcal/mol and Cartesian force
RMS≤0.5 kcal/mol/Å. Additional frozen maxima are0.3 kcal/mol per energy and
2 kcal/mol/Å per local vector-force error. Every configuration and every atom
is reported. Selected-protocol curvature estimates remain diagnostics, avoiding
duplicate statistical use of the same gradients. QZVP correction values are
not imposed on the TZVP fit.

Uncertainty uses64 deterministic whole-configuration/path-block bootstrap
replicates. Each primary path's signs and steps travel together. Every
rank/conditioning/bound failure is retained without replacement resampling.
At least80% of replicates must qualify, and frozen energy/force prediction-spread
limits must pass. Absolute parameter intervals cover near-zero coefficients;
correlations and prediction spreads are reported. This small-block analysis
does not include model-discrepancy, basis or physical uncertainty.

The three reserved configurations remain untouched unless numerical, solver,
training, identifiability, bound and uncertainty gates all pass. Leakage is
screened against all planned configurations before execution and rechecked
against actual observations before validation activation. They are fresh
deformed configurations, not independently sampled conformational basins.

## Executed numerical precision and compute

All 27 counted native attempts succeeded, without startup failure, interruption
or retry. They comprise5 QZVP controls,10 missing primary observations,5 joint
observations and 7 precision controls. Together with 7 reused J64 executions,
they provide22 consistent training geometries and separate sensitivity controls.
The three validation requests were never launched. No new Si–H calculation ran.

The 250×1454/300×1454 force RMS differences are1.56880e−5 at the center and
1.79451e−5 kcal/mol/Å at joint-040. Maximum local vector differences are
9.19841e−5 and 1.14627e−4, respectively. Their relative deformation-energy
contrast differs by 7.13478e−11 kcal/mol. The 1e−13/1e−11 SCF center control
ultimately met its native criteria within the frozen wall limit after energy
changes fluctuated near double-precision resolution. Its force RMS difference
is 1.24683e−5 and maximum local difference8.23754e−5 kcal/mol/Å; the native-energy
change is 8.56048e−10 kcal/mol. No threshold was relaxed.

At joint-040, the directional Cartesian-gradient contraction is−1.61007927
kcal/mol/Å. Central energy differences give−1.61108939 at 0.001 Å and−1.60931713
at 0.0005 Å. Errors−0.00101012 and+0.00076213 pass the frozen
`0.01 + 0.001*abs(reference)` criterion. Full signed Cartesian responses,
atom-wise comparisons, decimal energies and checks remain in
`fitting-numerical-checks.json`. These tests qualify internal consistency for
this specific target and bounded geometry selection; they do not qualify
complete-basis or physical accuracy.

## Trial outcome and local residuals

The actual 1804×59 weighted matrix includes all 58 physical variables and the
training-only intercept. Its rank is 59 and column-scaled condition 1556.0212,
passing the separately frozen10000 experiment screen. The strongest column
cosine is 0.9998867, demonstrating substantial correlation despite full rank.
The inherited1000 geometry-screen flag is false; it is not substituted for the
actual fitting criterion. All singular values and the weakest combination are
retained with the actual matrix. Nominal algebraic rank is distinct from
resampling robustness and scientific identifiability.

The unconstrained optimum lies inside every declared bound, so the bounded
solver returns it directly (status 3, no TRF iterations required). Its scaled
optimality is 2.87e−12. The single training intercept is 175.94134173 kcal/mol;
it is a nuisance energy reference, not a physical coefficient.

| Training group | Configurations | Energy RMSE (kcal/mol) | Cartesian force RMS (kcal/mol/Å) | Maximum local vector-force error |
| --- | ---: | ---: | ---: | ---: |
| Primary paths/center | 17 | 0.00177440 | 0.17201748 | 0.61499038 |
| Joint deformations | 5 | 0.02160503 | 0.72020482 | 4.29055395 |
| All training | 22 | 0.01041723 | 0.37516694 | 4.29055395 |

The global energyRMSE≤0.1 and forceRMS≤0.5 gates pass. Maximum absolute energy
error0.0343104 passes0.3. The local vector-force maximum fails2.0. The largest
Cartesian component error is 3.72634. Every joint configuration exceeds the
local-vector gate, so the global RMS cannot be used to conceal those failures:

| Configuration | Energy error | Force RMS | Maximum local error | Stable atom ID |
| --- | ---: | ---: | ---: | ---: |
| joint-040 | +0.0158561 | 0.847586 | 3.735275 | 3 (C) |
| joint-028 | −0.0343104 | 0.857111 | 4.290554 | 10 (H) |
| joint-038 | −0.0227996 | 0.546552 | 2.168066 | 2 (Si) |
| joint-044 | +0.0193485 | 0.656474 | 3.403272 | 5 (O) |
| joint-020 | −0.0033289 | 0.640903 | 2.299572 | 11 (H) |

Energy errors are MM−QM in kcal/mol; force errors are in kcal/mol/Å.
`trial-fit.json` retains all 22 configuration errors, all 27 atom-wise vector
norms per configuration, full signed force residuals and every trial value.
The joint geometries account for83.7549% of force residual sum of squares.

The joint corrections areΔK2(Si–C)=13.17871585 andΔK2(C–H)=−1.19404104
kcal/mol/Å². They differ from the selected-protocol conditional curvature
estimates25.55006250 and 20.93448581. No QZVP values were imposed. Independent
complete-model force contractions, supplying every one of the 58 parameters,
give unknown-null stretch curvatures401.073123 and 718.254402, versus target
425.815816 and 762.511455. The deficits are24.742693 and 44.257054 kcal/mol/Å².
They agree with verified fixed curvature plus twice the fitted correction to
better than1.5e−11. This diagnostic uses correlated training gradients and was
not added as a duplicate statistical measurement or an after-the-fact fitting
constraint. These trial corrections do not replace J64's conditional candidates.

The weighted objective and independently reconstructed global unbounded
least-squares lower bound are both0.57385279. **This does not demonstrate model
infeasibility**: a lower bound below2 cannot exclude simultaneous RMSE passes,
and those RMSEs already pass. The chosen soft objective fails the local maximum.
No alternate-weight fit, minimax fit or constrained local-force feasibility
calculation was performed. Source-baseline incompatibility and objective
tradeoff remain hypotheses, not established sole causes. No parameter family
is automatically added and no campaign is enlarged.

## Uncertainty and validation decision

All 64 deterministic whole-block replicates are retained, without replacement
resampling. They all have nominal rank 59;3 fail the conditioning screen and 48
hit parameter bounds, with overlap. Only 16/64 (25%) qualify, below the frozen80%
requirement. The nominal fit itself has no active bound. This distinguishes
conditional algebraic rank from robustness under block removal/reweighting.

Conditional prediction RMS95 among those16 survivors is 0.0857265 kcal/mol for
energy and 0.4920204 kcal/mol/Å for force. Those aggregate spreads pass their
limits, but do not cure the 75% failure fraction. Joint-specific force spreads
range0.870796–1.574859 kcal/mol/Å. Absolute parameter intervals, covariance and
correlations are retained. For example, the near-zero AAT trial value−0.08542
has absolute standard deviation5.45091 and conditional interval
[−12.87599,3.48647] in its recorded units. ΔK2 Si–C has conditional interval
[−10.10343,25.45841], and C–H[−49.10538,16.97850] kcal/mol/Å². These are
unqualified survivor-conditional, small-block development intervals, not
physical confidence intervals or uncertainty that includes basis/model error.

The local-force and bootstrap gates fail. Consequently the three reserved
configurations remain **untouched**, and no fresh-validation errors exist.
There are58 trial coefficients, **zero accepted method-specific candidates**,
zero force-qualified/transferable parameters and no automatic registration.
The fit executed successfully; reproduction of the declared protocol is **not
accepted**, while physical/transfer validation and basis qualification remain
unestablished for separate reasons.

The [next-data decision](evidence/phase_4j65/next-data-requirements.json) requires
**zero new QM initially**: preregister a development-only convex feasibility/
tradeoff audit on the existing 22 training geometries with the same58 variables,
bounds and training intercept, enforcing the unchanged RMSE and local-force
limits. It should determine whether the local failure is avoidable within this
model before a new objective or parameter family is selected. Audit the
bound-active block directions at the same time. Only a demonstrated remaining
capacity or identifiability requirement justifies a separately frozen targeted
design. The three validation configurations must stay reserved through that
selection. More basis calculations or more samples of the same design are not
the default remedy for this result. Capped-DP2 accounting/data remain separate.

## Reproduction

Execution stages use the existing pinned interpreters:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/phase_4j65.py controls \
  --work ../island-validation/phase4j65-qm \
  --python ../island-validation/phase4j36-native-env/bin/python
```

The retained coordinator and command records show dependent primary, joint,
precision and conditional-validation stages. For a new numerical rerun, copy
the frozen evidence and requests/plan/jobs into a new empty external work
directory; keep original ledgers and outputs in the archive. Do not overwrite
historical evidence or generate a replacement plan under its old identity.

The dependent execution commands, using the same `--work` and `--python`
arguments above, are `phase_4j65.py primary`,
`analyze_phase_4j65.py primary --work ...`, `phase_4j65.py joints`,
`phase_4j65.py precision`, and `analyze_phase_4j65.py trial --work ...`.
Analyze the controls before primary execution. The joint command enforces the
recorded primary authorization. Do not execute validation unless the trial
explicitly authorizes it and `analyze_phase_4j65.py validation-preflight` has
checked actual observed-geometry leakage. That stage was not activated here.
Exact commands, return codes and coordinator logs are retained under
`validation/`; native subprocess commands are additionally hash-bound inside
the native archive. Reissuing the controls command against the completed work
was verified to leave the attempt ledger byte-for-byte unchanged.

Offline replay reconstructs native intake, protocol contrasts and the trial in
a fresh relocated directory without invoking QM:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j65.py replay
```

## Software verification and preservation

The supported suite covers **5240 unique cases:5179 passed,61 skipped**.
Two disjoint base shards cover5215 cases; the final 25 J65 tests all pass,
including native reconstruction in a fresh relocated process. There are745
applicable J54–J64 regression cases. Earlier focused reruns are preserved and
are not counted again as terminal coverage. Tests cover exact-protocol reuse
and producer preservation, cross-site/grid contrasts, manufactured recovery of
all 58 nonzero parameters plus offset, bounds/rank loss, offset handling,
leakage, uncertainty failures and absolute near-zero uncertainty, native hash
mutations, startup diagnostics/caps and relocation.

Ruff over `src tests scripts`, compileall, all 5 applicable environment pip
checks, offline native/fit replay, verified resume, receipt/manifest and diff
checks pass. Analysis/source bindings remain unchanged. All 7526 files from the
exact fetched base are verified byte-for-byte; no historical identity, source
row, charge or native output is rewritten. Native scratch remains outside Git.

Historical J54–J62 and J64 receipt checks pass in the exact-base checkout.
The first J63 receipt check there failed because its broad implementation glob
includes later J64 modules. The original failure is retained; its retry passes
in J63's preserved delivery snapshot
`3fd458b33fe92e9bd206d82e2b7d9cdbb7cb75e3`. This is a snapshot-selection
correction, not a historical byte change. Early lint failures and their
corrections remain separately visible. The final checks are successful terminal
runs, not a claim that every earlier software command passed. In contrast,
all 27 native attempts in this campaign completed without a retry. Immutable
subprocess logs retain original whitespace under narrowly scoped attributes.

The native archive, normalized observations, actual weighted matrix, all trial
coefficients and bootstrap outcomes, independent checks, exact requests,
settings, process exit statuses and resource measurements are retained. Total
native elapsed time is 16080.60 seconds; peak process RSS is 9,989,252 KiB. No large
scratch arrays are committed. `training-capacity-diagnostics.json` is reproduced
by the retained `validation/postfit_diagnostics.py` helper; it performs complete
58-value force contractions and raw-matrix checks without new fitting or QM.

## Scope

Parameter completeness, backend correctness, target reproduction, orbital-basis
sensitivity, optional cross-engine equivalence, and scientific/chemical transfer
validation remain independent. A passing target fit could support only a bounded
method-specific HMDS research candidate, not native-PCFF equivalence, a physically
validated parameter set, condensed-phase accuracy or complete PDMS support.

Capped-DP2 interior `sio(osi,osi,c3,c3)`, O–Si–O and Si–O–Si–O environments
still need separate exact graph/source/charge accounting and reference data.
No Si–H values transfer to Si–C. No complete profile or production MD is created.

`production_validated=False`, `simulation_readiness="not_established"`.
