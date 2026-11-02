# Phase 4J68 — executed HMDS model-span and baseline diagnosis

The original 58-variable affine HMDS model cannot satisfy the four frozen training limits on the existing data even after removing every original physical coefficient bound. Two independent backends give checked lower bounds above one. A shared Si–O stretching polynomial reduces the minimax error substantially, but remains incompatible. **Accepted coefficients: zero.** The three reserved validation configurations remain untouched.

Fetched `origin/main` was `fbf093d7a75262a652192a57981cee0872c7c372`, which does not contain required J67 `873e1cc7d5d20f0d6f7c7e5345c3841e0892011c`. The J68 branch, `codex/phase-4j68-hmds-model-span-baseline-diagnosis`, starts directly from that exact prerequisite. No merge or push of main is involved. Unrelated working files are preserved.

## Reconstruction and frozen scope

Strict J67 receipt replay completed before the diagnostic solves (`validation/J67-parent-reconstruction.log`). This reconstructs the 27-atom neutral-singlet explicit-H HMDS graph, original 22 observations, four J67 native gradients and their original producer records, decimal native energies and complete Cartesian gradients, conversions, fixed MM terms, 58 physical columns, energy-only intercept, canonical constraints, nominal and independent SCS results, J66 witnesses, and all retained bootstrap draws. No historical observation receives a replacement execution identity. `dataset-and-producers.json` records the mapping and unchanged parameter specifications.

The reference protocol remains PBE0/def2-TZVP, DF/def2-universal-jkfit, unpruned 250×1454 integration, SCF energy/density thresholds 1e-12/1e-10, integral threshold 1e-12, no dispersion, retained Psi4 1.9.1. These are reconstructed observations, not new QM results. J66's original-22 feasibility, J67's expanded-data incompatibility and J67's disclosed proposal-sensitivity timing deviation are preserved. Parent receipts, native archives, source rows, charges and implementation snapshots are unchanged.

The signed plan was frozen before new solves: SHA256 `147b5785f64fd859cdacff207065c484b442da9f5b77a4efc0c3b26f0f7d0775`. It allocates at most 12 scientific solver calls, including failures and interruptions; eight were used. One worker and one BLAS thread used the isolated retained solver overlay: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, Clarabel 0.11.1 and SCS 3.2.9. Limits were 4 GiB address space per worker, 180 seconds per solver, 210 seconds per process and 3,600 seconds for the diagnostic campaign. The resource check found a 32 GiB cgroup, about 9.4 GiB then available, eight CPU quota cores and sufficient disk. Actual solver runtimes were approximately 0.53–4.92 seconds, with recorded worker peak RSS around 122–128 MiB. The recovered reporting process has its own explicitly identified RSS measurement.

Scientific campaign counts are **zero new QM attempts, zero MD steps, zero new bootstrap draws and zero validation executions**. Existing native integration fixtures exercised by the supported software regression suite are software checks, not HMDS reference observations or an MD research campaign. No Psi4/native execution interface is provided by the J68 runner.

## Problem and independent certificate

For both protected original22 and full26, the uniform-rho problem retains:

- energy RMSE ≤0.1 kcal/mol ×rho;
- Cartesian-component force RMS ≤0.5 kcal/mol/Å ×rho;
- each absolute configuration energy error ≤0.3 kcal/mol ×rho;
- each atom's three-component vector-force error ≤2.0 kcal/mol/Å ×rho.

The intercept is a single training energy offset. Its energy column is one and its force column is zero. It is unbounded and included in every rank, scaling and certificate calculation. No observation or atom is removed. Original equilibrium dependencies, K3/K4, charges, LJ96/Coulomb conventions and fixed Class-II terms remain unchanged. The bound-free diagnostic removes only the physical coefficient box. Its canonical matrices contain no coefficient-box rows. Shape-only placeholders required by the historical canonical helper for newly added variables enter neither constraints nor original-bound reporting. Its verifier reports `performance_limits_satisfied`, `original_bounds_satisfied` and both subset performances separately.

The historical finite-box certificate is used only for the bounded controls. For the bound-free problem let `z` include all column-scaled physical coefficients and the scaled intercept, `M=W/scale`, and `y` be the full-data normalized target. The energy and force norm constraints imply `||Mz-y|| ≤ sqrt(2) rho`. A direct original-unit check on **both subsets** supplies a feasible upper bound U. A checked approximate left inverse B supplies

```
epsilon >= ||B M - I||_2, epsilon < 1
sigma_lower = (1-epsilon)/||B||_F
R = (||y|| + sqrt(2) U)/sigma_lower
```

The implementation computes the defect in extended precision, bounds spectral norms by Frobenius norms and includes conservative division, accumulation and norm-error allowances. Thus every optimum with rho≤U lies inside the derived computational ball `||z||≤R`; this is not an imposed physical bound. Rank-deficient cases without this argument return an inconclusive certificate and never ignore unbounded null directions.

The raw dual is projected into the self-dual nonnegative/Lorentz cones. With `g=q+A.T lambda`, the independently checked lower bound is

```
-b.lambda - R ||g_z|| + min(0, U*g_rho) - arithmetic_allowance.
```

Every coefficient and intercept stationarity residual is included. For the original model the checked smallest-singular-value lower bound is 0.0025719444 and R≈654,093.45. The left-inverse defect is 1.69e-13 or less and its allowance is 2.99e-9. The final support arithmetic allowance is about 1.15e-5 for Clarabel. Raw stationarity vectors, cone margins, left inverse, scaling, primal/dual/slack vectors and canonical matrices are retained. No solver status alone establishes incompatibility; the margins below remain well above the frozen 1e-5 classification margin. SCS and Clarabel are independent backends.

## Executed comparisons

All rows enforce protected22 and full26. The two variants are separate bound-free diagnostic models; neither is combined with the other.

| Model | Clarabel rho | Checked lower, Clarabel | Checked lower, SCS | Checked upper, Clarabel |
|---|---:|---:|---:|---:|
| Original bounded | 1.196374166 | 1.196373926 | 1.196372588 | 1.196375166 |
| Original bound-free | 1.196230308 | 1.196218605 | 1.195918879 | 1.196231309 |
| Add shared bridge-angle polynomial | 1.192869063 | 1.192856945 | 1.189856970 | 1.192870063 |
| Add shared Si–O polynomial | 1.023051982 | 1.023040834 | 1.021157220 | 1.023052982 |

The bounded controls reproduce J67's Clarabel and SCS rho values and lower bounds at serialized precision. Clarabel reports `AlmostSolved` for the bounded and angular rows, and `Solved` for the original bound-free and Si–O rows. Inaccurate statuses are retained. Classification uses the checked bounds rather than promoting these status strings to candidate qualification.

| Clarabel model/subset | Energy RMSE | Force RMS | Maximum absolute energy | Maximum local force |
|---|---:|---:|---:|---:|
| Original bound-free/full26 | 0.11962303 | 0.59811515 | 0.34051247 | 2.39246062 |
| Original bound-free/protected22 | 0.08735399 | 0.51890663 | 0.27857086 | 2.39246062 |
| Si–O amendment/full26 | 0.10230520 | 0.51152599 | 0.22121092 | 2.04610396 |
| Si–O amendment/protected22 | 0.07431065 | 0.43366606 | 0.22121092 | 2.04610396 |
| Si–O amendment/new4 (reported separately) | 0.19406159 | 0.81634433 | 0.21893648 | 2.04610396 |

Energy units are kcal/mol and force units kcal/mol/Å. New4 is a diagnostic subset, with no separate acceptance threshold substituted for the required full/protected limits. Full residual vectors, configuration ordering, local errors and original-unit slacks are retained in the solver records and comparison summary.

Removing bounds improves rho by only 0.000143858. The original bound-free Clarabel vector places two middle-bond torsion coefficients below their original −100 kcal/mol/Å bounds by 35.5145 and 39.1377 kcal/mol/Å. Every outside-bound value, original interval, signed violation, normalized width fraction, name and unit is recorded in `original-bound-violations.json`. No such vector is an eligible force-field set. Minimum bound relaxation is **not activated**: even complete removal of the box cannot make this original span feasible.

## Projection, localization and hypotheses

The unbounded weighted least-squares projection retains all singular values, coefficients, intercept and target component outside the 59-column span. The scaled rank is 59 and condition is 864.8371. Weighted residual norm is 0.98523617; global force RMS is 0.47989851 but maximum local force error is 4.71205273. This projection floor alone does not prove incompatibility with the four distinct limits: its simple combined-norm lower bound is only 0.69667. The minimax certificates supply that conclusion. Full rank does not establish robust physical identification.

`fixed-family-localization.json` independently differentiates every actual fixed family, including source Si–O quartic bonds, Si–O–Si quartic angle, other bonds/angles, fixed bond-bond/bond-angle/angle-angle terms, source torsions, LJ96 and Coulomb. Family energies and Cartesian forces close against the original complete fixed evaluator. The explicit 58-coefficient adjustable contributions are grouped separately. No missing coefficient is silently zeroed. The audit retains forces at every atom and residual-direction projections at O5, Si2 and Si6.

For example at joint-002/Si2, component norms sum to 66.49 kcal/mol/Å while the norm of their total is only 5.39. At O5 the Si–O bond family contributes a 19.58 vector norm while its residual-direction projection is −19.47. The residual is 2.39246. These cancellations prevent identification of a cause from the largest component alone.

Before any variant solve, the two hypotheses were frozen in `plan/baseline-hypotheses.json`, SHA256 `40dc0193e635883c7e86ccfc04f61cce3333c04da6b1cb34b37917509b68085c`. Both use the scalar potential

```
Delta E = a sum_s(q_s-q0_s) + b sum_s(q_s-q0_s)^2
Delta F = -sum_s [a + 2b(q_s-q0_s)] grad(q_s).
```

There is no extra 1/2. The bridge-angle variant has site (Si2,O5,Si6), source type `sio,osi,sio`, reference 2.740620711236616 radians; a and b have kcal/mol/radian and kcal/mol/radian² units. The Si–O variant shares a and b across (Si2,O5) and (O5,Si6), source types `sio,osi` under reversal, references 1.6562 Å; units are kcal/mol/Å and kcal/mol/Å². Source rows remain `quartic_angle:cff91:2365` and `quartic_bond:cff91:1775`. Their original quartic conventions, equilibria and higher powers are unchanged.

These are additive research polynomials about fixed source references. A linear term can exert force at the source equilibrium; it is not a fully propagated equilibrium revision across Class-II dependencies. No per-atom force offset, per-configuration variable, arbitrary residual column or separate nonconservative force correction is introduced.

Both pairs add independent energy/force directions. Their column-normalized norms outside the old span are (0.04903,0.08292) for angle and (0.06241,0.47271) for Si–O stretch. Weighted projection capture norms are 0.13548 and 0.29946. Independent scalar finite differences and Jet derivatives check both columns at all 26 geometries. The amended scaled ranks are 61, with conditions 926.97 and 881.73. These information measures supported the predeclared tests, not parameter acceptance.

The Si–O diagnostic selects a≈10.41315 kcal/mol/Å and b≈18.77721 kcal/mol/Å², while refitting the other 58 free coefficients. These values describe a failed training diagnostic, not source PCFF rows, a stationary-point calibration or accepted coefficients. Its approximately 88.25% reduction in excess rho does not establish a unique physical cause. Significant coefficient cancellation remains: the maximum sum of individual adjustable-column force norms is about 151.51 versus maximum norm of their sum 7.09 kcal/mol/Å.

Limiting amended-model local residuals include joint-032/H10, joint-002/Si2/O5/Si6, joint-020/C4, joint-040/C3 and joint-028/H10, each near 2.04610. The original model also limits bridge O5 at several configurations. Full ranked atom lists and all residual vectors are preserved; `remaining-residual-localization.json` decomposes the amended model at these limiting Si/O/C/H atoms and includes their original weighted-projection residuals. The broad remaining pattern does not justify a unique additional parameter family or automatic combination of these variants.

## Stopping decision and retained failures

Stop this allocation. The existing adjustable span with its fixed baseline is insufficient under the frozen constraints. The original bounds are not the sole obstruction. The tested angular amendment barely helps; the shared Si–O amendment is quantitatively relevant but insufficient even with unconstrained coefficients.

The next dependency is a source Si–O baseline/equilibrium and Class-II dependency audit tied to the retained residual directions, including the limiting carbon/hydrogen sites. A separately frozen conservative-potential revision must demonstrate relevant independent columns and pass the **existing** protected22/full26 problem before further data acquisition is considered. The current evidence does not support widening bounds, combining variants automatically, adding random configurations, another basis ladder, or claiming a unique physical explanation. More observations cannot remove the demonstrated incompatibility of these existing constraints.

Eight research solver calls were launched, including one process that failed during postprocessing after successfully retaining its raw conic solution. The J65/J67 prediction helper correctly requires exactly 58 coefficients; the amended model has 60. The failure log and raw solution are preserved. A J68 generic report and complete amended-MM check recovered that exact state without a ninth solver call. The failed process's missing peak-RSS measurement is disclosed separately from recovery-process resource use. Raw solver vectors are never moved inward, clipped or repaired in these bound-free solves.

Preliminary software failures are retained separately: fresh-process argument handling, NumPy scalar comparison in replay, and the archive extractor's requirement for a nonexistent destination were corrected. Successful terminal checks do not erase those earlier logs. Synthetic solver tests are software controls and never substitute for the real HMDS observations. No new uncertainty qualification or bootstrap campaign is claimed.

## Verification and replay

Implementation is additive and research-only, reusing existing canonicalization, native evidence, complete-MM evaluation, graph mapping, manifests and receipt/JUnit infrastructure. Solvers remain optional and lazy. The ordinary `forcefield="pcff"` behavior is unchanged.

The terminal supported full suite and J54–J67 regressions, Ruff, compileall, environment-specific pip checks, historical byte checks and strict receipts are recorded under `validation/`. `collection.json` and `validation/junit-files.json` define unique terminal node coverage. Focused reruns and preliminary failures contribute no additional unique coverage. The first full-suite worker collected the old J68 relocation test body before its correction, retaining one failure. Terminal selection keeps every non-J68 case byte-for-byte from that worker and uses the final corrected 69-case J68 file once. The original six-shard logs/XML, failure and return code remain in an inventory-checked archive; the strict validator checks the exclusion/replacement node sets and rejects unrelated hidden failures or rewritten historical cases. Terminal unique coverage is 5,414 cases: 5,353 passed and 61 skipped, including 69 J68 cases. All 8,416 historical files are byte-preserved. The receipt records these exact counts. Explicit implementation inventories prevent later research modules from changing this receipt's scope.

Using the retained analysis environment (no solver package or network needed for replay):

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j68.py

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j68.py --deep-parent

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/validate_phase_4j68_receipt.py \
  docs/evidence/phase_4j68
```

Replay verifies canonical matrices, complete physical vectors, both subset constraints, certificate residual allowances and computational radii, fixed-family closure, derivatives, independent added span and retained failures. Compressed solver states have an exact file inventory. The fresh-process relocation test forbids importing Psi4, Clarabel or SCS. Re-signed scientific declarations and graph/protocol/matrix/bounds/units/order mutations are tested independently of matching outer digests.

Software verification, training feasibility, diagnostic hypothesis testing, parameter completeness, target reproduction, basis sensitivity, optional cross-engine equivalence, chemical transfer and scientific validation remain separate. No research candidate is accepted or registered. This exact-HMDS study establishes neither native PCFF equivalence nor complete-basis accuracy, PDMS/capped-DP2 transfer, condensed-phase accuracy or production readiness. Si–H remains stopped; no capped-DP2 or production campaign is initiated. The 28 confirmations remain user-attested; original author artifacts remain `not_yet_supplied`.

```
accepted_coefficients = 0
automatic_assignment = False
production_validated = False
simulation_readiness = "not_established"
```
