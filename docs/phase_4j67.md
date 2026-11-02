# Phase 4J67 — targeted HMDS information augmentation

This phase executes the four exact J66 proposals and assesses the unchanged
58-parameter Hamiltonian against the expanded method-specific HMDS dataset.
Four native attempts produced four verified gradients. The expanded, unchanged
model/data/bounds problem is numerically incompatible with the simultaneous
training limits. No nominal constrained least-squares trial was activated; zero
coefficients were accepted and all three reserved validation configurations
remain untouched. This preserves J66’s separate original-22 feasibility result.

## Frozen scope and reconstruction

The branch starts from fetched `origin/main`
`fbf093d7a75262a652192a57981cee0872c7c372`, which contains required J66
`6e1e3ae7dd20e4d153c6e23dc9107379a475350f` and J65
`1b03f835718599b7f6e53e3331400518112882ec`. Main is neither merged nor pushed.
Unrelated untracked files are preserved.

The target is the exact neutral-singlet, explicit-H, 27-atom HMDS graph
`d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`,
with stable atom order 1–27. J66 native observations, the 22-configuration
matrix, nominal/interior witnesses, constrained trial and all 64 bootstrap
outcomes were replayed before new execution. J66 remains a training-feasibility
result with zero accepted candidates. Its failed gates and all historical
archives, source rows, charges and implementation snapshots remain unchanged.

An additive parent audit checks receipt/manifest schemas and scientific
semantics, including feasibility versus qualification, parameter completeness,
validation/transfer/assignment declarations, counts, readiness and author
artifact status. Newly executed J67 tests mutate these declarations, rebuild
valid identities and manifest links, and require rejection. Explicit receipt
and manifest field inventories also reject unknown scientific override claims.
Each current field-mutation test first accepts an actual valid baseline receipt,
then checks the correctly re-signed mutation. These are new
mutation controls; they are not claims about tests performed during J66.

The exact stored coordinates, actions and geometry identities of `joint-032`,
`joint-002`, `joint-024` and `joint-034` are frozen in
[the scientific plan](evidence/phase_4j67/plan/scientific-plan.json).
Graph-aware proper alignment excludes duplicates and equivalence with all
existing observations and reserved validation, covering methyl/H/whole-end
permutations rather than the nine-atom shortcut. They are center deformations,
not independently sampled conformational basins. No replacement geometry was
selected.

The additional independent replay of all four proposed weak-direction
sensitivities agrees exactly. Its timing is explicit in
[proposal-sensitivity-reconstruction.json](evidence/phase_4j67/proposal-sensitivity-reconstruction.json):
this extra check ran after the first gradient completed. The exact coordinates,
parent bindings, leakage checks, requests, analysis rules and new bootstrap
already had been frozen before any new gradient; none changed. Thus the
additional numerical sensitivity check was later than the requested initial
reconstruction sequence. No new QM values entered that geometry-only check.

## Native experiment and data intake

The unchanged fitting target is PBE0/def2-TZVP, DF/def2-universal-jkfit,
unpruned 250×1454 grid, energy/density SCF thresholds 1e−12/1e−10,
integral threshold 1e−12, without dispersion. Effective native settings are
written explicitly; stale nested defaults in the historical plan are not
used as settings. Native requests exactly match the existing training protocol.
The pinned engine is Psi4 1.9.1; the orbital basis has 399 functions in 175
shells. Basis qualifications from J63–J65 remain unchanged.

The cap is eight counted native attempts, including startups, failures,
interruptions and retries: four training gradients, at most three conditional
validation gradients and at most one additive failed-attempt retry. There is
one native worker, four threads, 4 GiB configured memory, 1800 seconds per job,
64 GiB scratch and a 18000 second campaign limit. Retained dense-TZVP training
measurements of 470–568 seconds justify the timeout. Current cgroup memory,
anonymous/kernel/cache distinctions, CPU quota and free disk were checked
before execution and before each job. No optimization, Hessian, additional
scan, basis ladder, Si–H or capped-DP2 job is authorized here.

Native requests, command lines, return codes, stdout/stderr, full Cartesian
gradients, energies, resource use and hashes are retained. Independent decimal
energy and ASCII-gradient conversion checks establish signs and units without
subtracting large converted electronic-energy offsets. Only successfully
verified matching observations enter the dataset. The original 22 retain their
original producer metadata; reused observations acquire no new execution
identity. The four new records are separate producers.

## Protected constrained problem

The full dataset has 26 configurations, 58 unchanged physical columns and one
energy-only intercept. Fixed source equilibria, K3/K4, charges, exclusions,
pair conventions and parameter bounds are unchanged. Every complete MM check
supplies every physical coefficient; no missing coupling is silently zeroed.

The conic feasibility problem minimizes nonnegative rho with the original
bounds and the following limits on **both** original 22 and full 26:

- Energy RMSE≤0.1kcal/mol; maximum absolute energy error≤0.3kcal/mol.
- Cartesian force RMS≤0.5kcal/mol/Å; maximum local atom-vector error≤2kcal/mol/Å.

Rho multiplies these limits only for feasibility diagnosis. Acceptance always
uses the original limits. Force RMS includes every Cartesian component; local
errors are three-component Euclidean norms. The protected constraints prevent
adding easier configurations from diluting an existing RMS failure.

If a witness is independently verified, the second problem minimizes the
original objective, `mean(rE²)/0.1² + mean(rF²)/0.5²`, on the full 26, subject
to both subsets' constraints and original bounds. The epigraph and variable
scaling reuse J66. Protected cones are transformed to the same 59-column
scaling; they add no second objective or intercept. No regularization,
parameter release or duplicate gradient-derived curvature measurement is added.

Clarabel 0.11.1 uses the inherited settings: 300 iterations, 120 seconds,
absolute/relative gap tolerance1e−9 and feasibility tolerance1e−10.
`AlmostSolved` remains unqualified. Qualification additionally checks primal/
dual residuals≤1e−8 and absolute gap≤1e−7. Numerical tolerances never widen
physical limits. Raw primal/dual/slack states, canonical matrices, column
scales, gaps, bound distances and physical slacks are retained. Independent
original-unit arrays and the complete MM evaluator check the results.

The frozen selection rule permits only the first strictly feasible convex
blend toward a verified rho witness in the retained finite sequence. The raw
vector, selected vector, fraction and objective change remain visible.
Original-bound activity of the raw vector cannot be erased by this selection;
a tiny inward move cannot manufacture candidate eligibility.

## Genuine expanded bootstrap and validation gate

The new 64 draws use seed 4067 and 14 whole blocks: the ten inherited blocks
plus one for each new configuration. Each draw samples 14 blocks with
replacement. Each primary path's signs/steps remain together; each new block
can be omitted. The frozen draws are in
[expanded-bootstrap-draws.json](evidence/phase_4j67/plan/expanded-bootstrap-draws.json).

Every replicate optimizes and constrains only its sampled observations with
multiplicity. Protected22/full 26 predictions are evaluated afterward, without
refitting. Failed, infeasible, inaccurate, ill-conditioned and bound-active
draws remain; there is no replacement resampling. The original J65/J66 draws
and the J66 guaranteed-addition geometry stress check remain separate
historical diagnostics. That stress check is not bootstrap qualification.

Eligibility retains rank 59, scaled condition≤10000, no original-bound activity,
qualified solver accuracy and verified sampled physical constraints. Overall
qualification also requires≥80% eligible and prediction RMS95≤0.1kcal/mol
energy and≤0.5kcal/mol/Å force. Survivor intervals are absolute, explicitly
conditional and unqualified when the overall fraction fails. Descriptive
intervals including failed fits are separately labelled and cannot qualify
a candidate. Atom-wise prediction effects and all overlapping failures are
reported.

All nominal physical, numerical, bound, identifiability and uncertainty gates
must pass before the three reserved configurations may execute. Activation
binds the fixed candidate/intercept, trial and dataset hashes, repeats leakage
checks and checks the remaining attempt budget. Validation uses the frozen
training intercept, with no refitting. Consumption is recorded before startup,
including failed attempts. A failed prerequisite leaves all three untouched.

## Executed outcome

The four requests completed without a failed attempt or retry. There were no
validation attempts, geometry optimizations, Hessians, Si–H/DP2 calculations or
research MD steps. The count is **4 of the maximum 8 native attempts**, with
**4 successful new gradients**.

| Exact configuration | Native energy (hartree) | Relative energy (kcal/mol) | Process time (s) |
| --- | ---: | ---: | ---: |
| joint-032 | −893.3788898790668327 | 1.8539767285 | 813.098 |
| joint-002 | −893.3792263549414656 | 1.6428349294 | 642.239 |
| joint-024 | −893.3795351887379184 | 1.4490387962 | 869.922 |
| joint-034 | −893.3788768531811684 | 1.8621505951 | 527.419 |

The relative-energy origin is the retained dense-grid center. Exact result
hashes, original/new producers, normalized samples and full converted forces
are in [fitting-dataset.json](evidence/phase_4j67/fitting-dataset.json).
[native-outputs.tar.xz](evidence/phase_4j67/native-outputs.tar.xz) contains
all requests, native outputs, full gradients, ledgers and resource records;
its unpacked inventory is independently checked during offline replay.

The actual weighted matrix has rank **59** and scaled condition **864.8371**,
compared with J66’s **1556.0212**. Better conditioning does not establish that
the model can reproduce the new observations. The frozen J66 constrained
vector predicts new-four force RMS **1.09850 kcal/mol/Å** and maximum local
force error **4.85219 kcal/mol/Å**.

The new protected minimax solve gives rho **1.1963741659**. Clarabel terminates
`AlmostSolved`, which remains unqualified under the inherited accuracy rule.
An independent cone projection and support calculation, accounting for every
stationarity residual over the derived finite variable box, gives the checked
interval **[1.1963739257, 1.1963751660]**. Its arithmetic allowance is
2.32×10⁻⁷. Both the checked lower bound and its margin above 1 are retained;
the conclusion does not depend on accepting the solver status.

A genuinely different solver, SCS 3.2.9, gives rho **1.1963741924** and an
independently reconstructed lower bound **1.1963725878**. Its numerical status
and residuals remain separate from the Clarabel qualification policy. Both
lower bounds exclude rho≤1 by approximately 0.19637. This demonstrates
incompatibility of this exact combined model/data/bounds problem, without
asserting a unique physical cause.

The diagnostic minimax vector has the following errors. It is neither a
passing witness nor a constrained least-squares trial.

| Subset | Energy RMSE (kcal/mol) | Cartesian force RMS (kcal/mol/Å) | Maximum absolute energy error (kcal/mol) | Maximum local atom-vector error (kcal/mol/Å) |
| --- | ---: | ---: | ---: | ---: |
| Protected original 22 | 0.0873311 | 0.519061 | 0.278582 | 2.392748 |
| New four | 0.2260271 | 0.918722 | 0.340158 | 2.392748 |
| Full 26 | 0.1196374 | 0.598187 | 0.340158 | 2.392748 |
| Frozen limits | 0.1 | 0.5 | 0.3 | 2.0 |

The new observations are also reported individually, using the same single
intercept as all other configurations.

| Configuration | Signed energy error (kcal/mol) | Force RMS (kcal/mol/Å) | Maximum atom-vector error (kcal/mol/Å) | Worst stable atom ID |
| --- | ---: | ---: | ---: | ---: |
| joint-032 | 0.210736 | 0.847548 | 2.392748 | 5 |
| joint-002 | −0.340158 | 1.038626 | 2.392748 | 5 |
| joint-024 | −0.183871 | 0.847588 | 2.392748 | 5 |
| joint-034 | 0.102116 | 0.927745 | 2.392748 | 5 |

Every original-unit constraint, each configuration and each atom-vector
residual, bound distance, dual, slack, solver setting and canonical problem is
retained in `nominal-rho*`. Complete MM evaluation with all 58 physical
coefficients agrees with the physical arrays to **2.85×10⁻¹⁴ kcal/mol** in
energy and **3.56×10⁻¹⁵ kcal/mol/Å** in force. Both original 22 and full 26
checks fail. The nominal constrained least-squares stage was therefore
**not activated**, as frozen; this is not a failed least-squares optimizer.
No curvature targets or parameter families were added.

## Expanded-data robustness and limiting directions

All **64 genuine expanded block draws** ran their sampled feasibility problem.
Of these, **31** have checked incompatible sampled constraints; their
least-squares stages remained inactive. The other **33** ran constrained
least squares. All outcomes remain in
[bootstrap-solves.tar.xz](evidence/phase_4j67/bootstrap-solves.tar.xz), with
355 inventory-checked files and the unchanged draw-to-outcome mapping.

Only **1/64 (1.5625%)** is eligible, versus J66’s **2/64 (3.125%)**. The
required fraction remains 80%. All sampled matrices have rank 59; their
conditions range from 659.27 to 2426.32, below 10000. This is not robust
parameter identification.

The overlapping failure categories are:

| Failure | Draws |
| --- | ---: |
| No verified sampled feasibility, each with checked incompatibility | 31 |
| Unqualified constrained solver accuracy (`AlmostSolved`) | 32 |
| Original-bound activity | 25 |
| Rank or conditioning failure | 0 |

All 25 bound-active fits also fail solver accuracy; 7 more fail solver accuracy
alone. The remaining 1 fit has `Solved` status and passes the residual/gap,
physical, rank and bound checks. Ignoring solver accuracy would leave at most
8/64 eligible, still far below 80%. Numerical-status failures therefore do not
explain away the statistical failure. No outcome was replaced or a new seed
selected.

The new blocks were omitted in 18 (`joint-032`), 19 (`joint-002`), 22
(`joint-024`) and 27 (`joint-034`) draws. Replicate 46 is the sole eligible draw;
it omits `joint-002`, `joint-024` and `joint-034`. These are correlated block
inclusion associations, not isolated causal experiments. Its predictions on
the original 22 and all 26 were evaluated afterward without refitting or
constraints from omitted observations. All per-draw protected/full prediction
errors and atom-wise effects are retained. For the 31 draws without an
activated fit, separately labelled posthoc minimax predictions cover all 26
observations; they do not enter fitted uncertainty or qualify a candidate.

Recurring original-bound activity affects middle-bond/torsion terms:
`hc,c3,sio,osi:block0:n1` (11 fits), `c3,sio,c3,hc:block0:n1` (10), and
`hc,c3,sio,osi:block0:n2` (2), plus end-bond/torsion
`c3,sio,osi,sio:block1:n1` (5). Counts overlap. In the nominal minimax vector,
`middle_bond-torsion_3:hc,c3,sio,osi:block0:n1` is near its original −100 bound.
No inward selection erased raw bound activity.

There is only one eligible survivor, so no statistically meaningful qualifying
survivor interval or RMS95 is reported. Separately labelled **descriptive**
intervals/correlations use all 33 executed fits, including inaccurate and
bound-active outcomes. They exclude the 31 unexecuted fits and do not qualify
a candidate. For example, their absolute 2.5–97.5% ranges for the two K2
corrections are −10.0919 to 37.4546 (Si–C) and −13.5925 to 29.5737 (C–H), in
kcal/mol/Å² under the retained correction convention. The problematic
middle-bond/torsion n1 term spans essentially its entire [−100,100] range.
Corresponding end-bond/torsion n1 terms for `c3,sio,c3,hc` and
`hc,c3,sio,osi` correlate at 0.99958. The weakest scaled nominal combination
mixes the energy intercept, angle-angle terms and several end-bond/torsion
terms. Its unit effects are reported for every observation; the largest local
force effect is 0.02631 kcal/mol/Å at `joint-024`.

Descriptive full 26 prediction spread about the minimax vector is
energy RMS95 **0.14554 kcal/mol**, force RMS95 **0.62245 kcal/mol/Å**.
The largest descriptive atom-vector95 effect is **7.48587 kcal/mol/Å**,
O5 in `joint-044`. These include failed fits and are sensitivity diagnostics,
not accepted uncertainty estimates or validation errors. The full absolute
intervals, covariance/correlation matrices, omission associations and
atom-wise effects are in
[block-identifiability-diagnosis.json](evidence/phase_4j67/block-identifiability-diagnosis.json).

## Bounded decision

The allocation stops after the exact four observations. **Zero candidates and
zero coefficients are accepted.** The three validation configurations remain
unexecuted and unconsumed, confirmed by an actual gated runner call that
launched no native worker. The hash-bound nonactivation record preserves the
failed training, numerical, bound and uncertainty gates. Identifiability at
nominal matrix rank/conditioning passes; it cannot override the other gates.

The limiting minimax local residuals involve O5/Si2 in `joint-034`, O5/Si6/Si2
in `joint-002`, O5 in `joint-032` and `joint-024`, and retained `joint-040`
C3 and `joint-028` O5/H10. They reach approximately 2.392748 kcal/mol/Å.
Original-bound activity and weak correlated directions identify where to
investigate, but do not distinguish a bounds-only obstruction from the
functional form/fixed baseline.

[The bounded next decision](evidence/phase_4j67/next-step-decision.json)
proposes a separately frozen **zero-QM model-span feasibility diagnostic** on
the same 26 observations, retaining the same four thresholds, protected22,
units and single intercept, with physical bounds removed only for diagnosis.
It has not been executed in J67. A feasible outside-bound vector would
localize the required violations without authorizing those coefficients;
a checked incompatible bound-free problem would motivate one specifically
localized fixed-baseline or functional-form hypothesis. An inconclusive solve
would remain numerical. Additional observations cannot remove the present
incompatible constraints, so another generic sampling campaign or basis
ladder is not prescribed.

## Verification and replay

The supported terminal suite has **5,345 unique cases: 5,284 passed,
61 skipped, zero failures/errors**. Its six file-disjoint historical/general
shards contribute 5,265 cases; the separate J67 terminal file contributes 80
passing cases. The union contains 795 J54–J66 regression cases. Test collection
and terminal JUnit nodes are checked for an exact unique multiset match.
The intermediate 78-case passing J67 run is separately retained; it and
other overlapping preliminary/focused runs add no distinct coverage. Ruff,
compileall, all applicable environment pip checks and every J54–J66 historical
receipt validation pass. Final diff checks and the strict J67 receipt validator
also verify historical byte preservation and all artifact inventories.

Terminal counts and failed/retried checks are recorded under
`evidence/phase_4j67/validation/`. The original failed terminal test, preliminary
suite interruption, capture streams and metadata have separate archives/logs;
successful corrections do not overwrite them. A preliminary full-suite launcher had an incorrect optional solver
bridge path; its interruption and artifacts are retained separately from the
corrected terminal run. Pip-check failures from applying that bridge to
unrelated interpreters are retained, with environment-specific retries. The
analysis bridge's missing declared CFFI dependency was installed only in that
isolated bridge; the native environment and numerical solver versions were
unchanged. An initial preflight exception from an unavailable aggregate cgroup
field was corrected before plan freezing or native execution. An initial
terminal J67 run passed 77 tests and failed 1 synthetic protected-force test:
its exactly zero parameter column used the historical 1e−30 scaling floor,
making the conservative certificate arithmetic allowance uninformative. The
corrected fixture keeps an identifiable energy column and the same
parameter-independent force residual, isolating protected force-RMS dilution.
The original failed log and JUnit bytes, including their padded blank lines,
are retained in `validation/terminal-attempt1.tar.xz` with a checked inventory;
the corrected terminal run remains directly readable. No actual matrix,
physical limit, coefficient or solver result changed. An additive diagnostic
collection initially passed an existing directory to the safe unpacker; its
rejected attempt and corrected child-directory collection are retained. This
collection added all-draw prediction reporting, with no optimizer or native
retry.

The executed analysis environment is Linux x86_64, Python 3.11.16,
NumPy 2.4.6/SciPy 1.17.1, with an isolated lazy solver bridge containing
Clarabel 0.11.1, SCS 3.2.9, CFFI 2.1.1 and pycparser 3.12. The native environment
is separate: Python 3.11.14, Psi4 1.9.1, NumPy 2.3.4, SciPy 1.17.1,
QCEngine 0.33.0 and QCElemental 0.29.0. Both execution environments pass
`pip check`; all applicable retained project environments also pass their
checks without unrelated solver overlays. The cgroup limit is 32 GiB with
8 CPUs; native jobs used the smaller frozen four-thread/4 GiB allocation.
Neither an arbitrary NumPy/BLAS combination nor an unpinned QM engine is
claimed to reproduce historical serialized floating-point identities.

Offline replay requires the pinned NumPy/SciPy analysis environment and the
repository evidence, without Psi4, solver execution or network. Optional new
optimization additionally uses the isolated Clarabel/SCS bridge. Explicit
implementation inventories prevent later modules from changing historical
receipts through a broad glob. The 61 supported-suite skips retain their
reasons: 45 require explicitly supplied author sources/engine docs; the rest
require opt-in external engines, physical qualification or archived reference
fixtures. Skipped controls are not claimed as executed scientific validation.

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/replay_phase_4j67.py
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/validate_phase_4j67_receipt.py docs/evidence/phase_4j67
```

Software verification, training feasibility, candidate qualification, untouched
validation, basis sensitivity, chemical transfer and scientific validation are
separate outcomes. This is exact-HMDS, protocol-specific research. It does not
establish native PCFF equivalence, complete-basis accuracy, condensed-phase
accuracy, capped-DP2 interior coverage, methyl-PDMS transfer or complete PDMS
assignment. The smallest domain extension remains a separately audited
capped-DP2 interior graph/charge/source inventory and data for
`sio(osi,osi,c3,c3)`, O–Si–O and Si–O–Si–O; none is initiated here.

The 28 confirmations remain user-attested; original author artifacts remain
`not_yet_supplied`. Ordinary `forcefield="pcff"` behavior is unchanged. No trial
is automatically registered, no source row is overwritten and no production
MD is run. `production_validated=False`;
`simulation_readiness="not_established"`.
