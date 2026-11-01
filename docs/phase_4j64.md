# Phase 4J64 — HMDS reference precision and conditional stretching calibration

The matched dense-grid control closes the retained siloxane numerical-consistency
question at its declared tolerance. Two unknown-null stretching observables
support conditional quadratic corrections, while full-force basis adequacy and
scientific candidate validation remain unestablished. There are 17 successful
new native gradients and 18 counted process attempts, including the preserved
startup failure. No force-qualified parameters are accepted.

This phase studies the exact neutral-singlet, explicit-H, 27-atom HMDS control,
Me3Si–O–SiMe3, graph
`d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`.
The Si–H investigation remains stopped. The ordinary API is `forcefield="pcff"`.
No complete PDMS profile, automatic parameter assignment or production MD is
created by this experiment.

## Base, reconstruction and frozen plan

Fetched main and exact branch base are
`3fd458b33fe92e9bd206d82e2b7d9cdbb7cb75e3`, which contains J63 prerequisite
`c10e76f05316c05e66d1d4c0df6581af775572cd` (merged by PR86).
The branch is `codex/phase-4j64-hmds-precision-baseline-calibration`.
Main is neither merged nor pushed. Unrelated local work remains untouched.

Before new QM, the existing J63 replay reconstructed all 67 native gradients,
the stationary candidate, path responses, complete 56-column mapping and
unknown-null curvature checks. Historical inputs, source rows, charges,
identities and failed gates remain unchanged.

The [frozen plan](evidence/phase_4j64/plan/scientific-plan.json) SHA-256 is
`6c324b0e8fdc9dae7955324f2d20b9cd3758bda4b5fe807073eaad407d160b1e`.
All process attempts, including failures and interruptions, count against 24.
There are 21 possible native requests: required dense-grid triple and QZVP
center; eight matched stretch grid/SCF controls; four conditional QZVP stretch
points; and five conditional soft-torsion step/rigid-rotation controls. No joint
design or reserved-validation QM is authorized. Activation rules and precision
thresholds were frozen before execution.

The pinned Psi4 1.9.1 environment runs PBE0 without dispersion, density fitting
with `def2-universal-jkfit`, integral threshold 1e−12, and four threads. Orbital
bases have 399/175 (TZVP),553/213 (TZVPP) and 1,079/343 (QZVP) functions/shells;
the auxiliary has 1,107/359. Historical shell inventories and basis-file hashes
are retained. QZVP is a larger orbital basis; the earlier TZVP/TZVPP contrast
is specifically additional-polarization sensitivity.

The actual cgroup limit is 32 GiB and 8 CPU equivalents, despite the host's much
larger apparent free memory. Preflight records anonymous and reclaimable memory,
disk space and retained native timings. QZVP runs alone with 8 GiB Psi4 memory
and 7,200 seconds per job; TZVP has 4 GiB and 1,800 seconds. A full QZVP three-index DF
tensor would occupy about9.6GiB, so native integral blocking/disk behavior and
actual RSS are retained. Scratch is capped at 64 GiB, one worker at a time,
campaign wall time 36,000 seconds. No full or numerical Hessian is implicit.

The initial coordinator launch failed before starting a gradient because
`/usr/bin/time` was unavailable. Its actual coordinator exit1, traceback,
command, stdout/stderr and unstarted native status are retained. An initial
collection note used conventional unavailable-command status127; the terminal
record distinguishes that note from the actual coordinator exit. An additive
[maintenance amendment](evidence/phase_4j64/plan/startup-maintenance.json)
archives the original sources, replaces that dependency with Linux Python
`getrusage`, retains process-group timeout handling and gives the retry a
separate process identity. Request geometries, settings, acceptance criteria and
24-attempt cap are unchanged. This is not an unqualified clean first attempt.

An additional [exact-center correction](evidence/phase_4j64/plan/exact-center-correction.json)
was frozen after discovering that the zero-angle path construction rounded one
center coordinate by2.220446049250313e−16Å. The successful original calculation
is retained as a separate roundoff control. A new request uses the exact stored
J63 center, so all main grid-only and basis-only contrasts retain **bitwise equal
input coordinates**. The maximum is now 22 native jobs plus the original failed
startup, 23 counted attempts within 24. The correction changes neither an optimized
center nor a numerical acceptance criterion. The original analyzer is archived;
its additive update selects the explicitly identified exact-center job. The
waiting follow-on orchestrator was replaced before launching any follow-on QM;
no native calculation was interrupted.

## Objective-specific reference checks

The dense-grid triple uses the exact retained J63 center and siloxane
±0.02rad geometries, PBE0/def2-TZVP, unpruned 250×1454 quadrature and SCF
1e−12/1e−10. Its grid-only contrast uses the retained 150×974 triple with
identical thresholds. Center, method and auxiliary basis never change in that
contrast. The separate QZVP comparison uses150×974 and SCF1e−10/1e−8,
matching the retained primary TZVP/TZVPP geometries and settings.

Native decimal total energies are read and subtracted before unit conversion
with 50-digit decimal arithmetic. An independent line parser repeats the native
energy reconstruction using 60-digit arithmetic and a separately derived unit
conversion; its additive source revision was frozen before its first execution.
The report retains both
printed final-energy strings and total energies, full ASCII Cartesian gradients,
force conversion/sign checks, path tangents and accelerations, torques and
curvatures. Text-rounding and binary serialization/subtraction bounds are
quantified separately from SCF/integration error. The old convert-before-
subtraction result is preserved for comparison, rather than attributing a
physical discrepancy to cancellation without a numerical bound.

Path rotations still cut the selected acyclic graph bond and rotate the whole
connected fragment. Independent endpoint cross products reconstruct tangents
and accelerations. Rotational second derivatives include `t^T H t + g·a`.
These molecular curvatures are not isolated torsion coefficients. Any conditional
whole-molecule rotation is a proper rigid rotation of the entire triple,
including tangents and accelerations; it tests grid-orientation sensitivity.

The [dense-grid control](evidence/phase_4j64/dense-grid-decision.json) passes the
unchanged siloxane energy/gradient consistency gate. At 0.02 rad, decimal-energy
curvature is 1.43430145 and gradient-derived curvature is 1.43230072
kcal/mol/rad². Their difference is 0.00200073, versus the frozen allowance
`0.02 + 0.02*abs(gradient curvature)`. The matched 150×974/250×1454 contrast
passes all declared energy, Cartesian/local-force, torque and curvature checks.
The gradient-curvature grid shift is 0.01198151 kcal/mol/rad²; the two torque
shifts are 0.00092394 and 0.00044468 kcal/mol/rad.

Printed-decimal rounding contributes at most 3.14e−10 kcal/mol/rad² and the
conservative binary serialization/arithmetic bound is 7.13e−7. Neither explains
the old 0.105309 kcal/mol/rad² discrepancy. The original zero-path center and
the exact center differ in energy by −3.44142e−7 kcal/mol, which would change
the second difference by 0.00172071 kcal/mol/rad². Their native screened-grid
inventories also differ by one point. This is an observed numerical repeatability
control, not a pure arithmetic experiment or proof of a sole mechanism.
Only the exact center enters the main matched contrasts.

The denser grid resolves the retained consistency question at the declared
pilot tolerance. J63's failed gates remain unchanged. Since the preregistered
activation condition is false, the additional five torsion step/rigid-rotation
QM controls are not launched. A single QZVP center cannot establish a QZVP
stationary point, all-path convergence or complete-basis accuracy.

The [QZVP center comparison](evidence/phase_4j64/QZVP-center-decision.json)
completed in 682 seconds. Peak process RSS was 9,924,916 KiB, within the
separately justified 12 GiB process allowance; configured Psi4 working memory
was 8 GiB. The native calculation retained 1,079 orbital functions in 343
shells and the unchanged 1,107-function auxiliary basis. Contrasts use identical
coordinates, grid, SCF thresholds and auxiliary settings.

| Center contrast | Cartesian-force RMS | Maximum local vector difference | Atom |
| --- | ---: | ---: | --- |
| TZVP / QZVP | 0.535341 | 1.441261 | Si, stable ID 6 |
| TZVPP / QZVP | 0.308181 | 1.385143 | Si, stable ID 6 |

Force units are kcal/mol/Å. Silicon ID 2 has nearly the same local differences;
the full signed vectors and all 27 atom-wise norms are retained in the report.
These sizeable shifts do not establish full-force basis adequacy. In particular,
the TZVPP/QZVP shift is not evidence that the earlier TZVP/TZVPP force failure
has disappeared. Stretch-curvature precision is tested independently below.

## Minimal baseline hypothesis

The separate internal research identity is
`island_research_class2_hmds_j64_quadratic_v1`. The original J63 Hamiltonian and
all source rows remain unchanged. The pinned FRC equation at line1648 is

```
E_bond = K2 (r-r0)^2 + K3 (r-r0)^3 + K4 (r-r0)^4
```

There is no additional 1/2. The total bare-bond second derivative is
`2K2 + 6K3(r-r0) + 12K4(r-r0)^2`, but the calibration uses **all** verified
bonded and nonbonded contributions to the molecular path curvature.

The only new variables are shared additive Si–C and C–H quadratic corrections:

```
E_correction = delta_K2_SiC sum_6_SiC (r-r0_source)^2
             + delta_K2_CH  sum_18_CH (r-r0_source)^2
```

Units are kcal/mol/Å². Source equilibrium distances, cubic/quartic coefficients,
angles, charges, torsions, exclusions, LJ mixing and isolated electrostatics
remain fixed. Because no equilibrium changes, dependent cross-term equilibrium
values and the original 56 kernels are unchanged.

A pure graph-fragment stretch changes only one graph bond linearly; all graph
angles and proper torsions stay fixed. Each missing coupling is at most linear
in that displacement. Consequently **all 56 unknown second path derivatives
vanish**. Independent scalar energies and Cartesian-force parameter columns
verify this relation at both retained steps. The two correction sensitivities
are diagonal2,0 and0,2. Thus a conditional estimate is
`delta_K2=(QM_total_curvature−verified_fixed_total_curvature)/2`.
This never equates molecular curvature with a bare source bond coefficient.

A complete energy/force call to the extension still requires all 58 explicit
finite coefficients. Neither the original 56 nor the two corrections have zero
defaults. Manufactured nonzero values are software controls only. Curvature-only
estimates are conditional, protocol-specific research candidates: they are not
force-qualified, transferable, native PCFF rows or validated complete HMDS
parameters. Step/grid/SCF/basis shifts are sensitivity diagnostics, not statistical
confidence intervals or complete-basis uncertainty bounds.

Soft-torsion precision, stretch-curvature precision, calibration completion,
physical validation and parameter completeness are reported independently.
A failed soft-torsion check cannot invalidate a mathematically independent
stretch observable; a passing stretch check cannot authorize a full force fit.

The [calibration record](evidence/phase_4j64/baseline-calibration.json) includes
eight new TZVP stretch-control gradients and four QZVP stretch gradients.
At the retained 0.003 Å step, the verified fixed total curvatures are
374.71569123 for Si–C and 720.64248379 for C–H, in kcal/mol/Å².
The following estimates use gradient-derived molecular path curvatures:

| Protocol | Si–C total QM curvature | C–H total QM curvature | ΔK2 Si–C | ΔK2 C–H |
| --- | ---: | ---: | ---: | ---: |
| TZVP, 150×974, standard SCF | 425.890158 | 762.506606 | 25.587234 | 20.932061 |
| TZVP, 150×974, tight SCF | 425.882090 | 762.515064 | 25.583199 | 20.936290 |
| TZVP, 250×1454, tight SCF | 425.815816 | 762.511455 | 25.550063 | 20.934486 |
| TZVPP, 150×974, standard SCF | 425.513372 | 759.152781 | 25.398840 | 19.255149 |
| QZVP, 150×974, standard SCF | 423.723829 | 755.060453 | 24.504069 | 17.208984 |

All table entries have units kcal/mol/Å². Standard SCF thresholds are
1e−10/1e−8; tight thresholds are 1e−12/1e−10. The TZVPP observations are
retained J63 data, reconstructed without new calculation. Each row is a
separate conditional estimate; protocols are not pooled into a fitting dataset.

For both paths, the matched SCF and grid full precision gates pass. The QZVP
curvature-only gates pass, but the corresponding full energy/Cartesian-force/
path precision gates fail. Every declared energy/gradient curvature-consistency
check passes. The retained 0.003/0.006 Å standard-TZVP correction estimates
change by 0.008789 for Si–C and 0.000876 for C–H. In contrast, the TZVP/QZVP
correction changes are −1.083165 and −3.723077, respectively. These correlated
step and protocol sensitivities are not statistical confidence intervals.

Independent scalar and force-derivative checks include every fixed interaction.
The largest scalar unknown-column second difference is below 1.8e−10, and
the force-chain-rule counterpart is below 6e−13. This verifies the specific
null observables without assigning any of the 56 missing physical coefficients.
Energy-derived finite-step estimates are also retained, with their separate
truncation and numerical sensitivity; they are not silently substituted for
gradient-derived estimates.

The [scope decision](evidence/phase_4j64/candidate-and-scope-decision.json)
accepts calibration completion only for these two conditional total-curvature
observables in the exact HMDS control. There are **zero** force-qualified,
transferable, scientifically validated or automatically assigned candidates.
The soft-torsion precision flag covers the matched TZVP grid/SCF and
energy/gradient-consistency experiment only. It supplies no torsion basis or
global force-reference qualification. Parameter completeness stays incomplete.

## Joint design and chemical scope

The geometry-only design appends two baseline columns to the complete 56-column
energy/Cartesian-force Jacobian. Diagnostic scales remain0.05kcal/mol and
0.5kcal/mol/Å, without adopting physical fitting weights. It rechecks the original
five proposals, then, only if needed, screens at most three extra selections
from the original48-member pool. Minimum-cardinality statements apply only to
the bounded selected subset. None of those proposed QM jobs is executed here.

The freshly replayed 58-column Jacobian has rank 58. Its scaled condition number
is 73,573.747288 on the existing 17 geometries. The same five J63 proposals,
`joint-040`, `joint-028`, `joint-038`, `joint-044` and `joint-020`, reduce it to
987.014857, below the unchanged 1,000 design criterion. Thus two baseline
variables do not require another selected training geometry within this bounded
design. This is a geometry-only result, not physical fitting readiness.

The three reserved validation deformations remain untouched. Leakage checks use
HMDS's graph-aware methyl/H/whole-end permutation handling, not the nine-atom
shortcut. Deformation labels do not establish independently sampled basins.

HMDS addresses `sio(osi,c3,c3,c3)` and `c3(sio,hc,hc,hc)` terminal environments.
Capped-DP2 interior `sio(osi,osi,c3,c3)`, O–Si–O and Si–O–Si–O environments
need a separate exact graph, source/type/charge accounting and reference design.
There is no Si–H transfer, foreign-force-field borrowing or charge normalization.
The 28 convention confirmations remain user-attested; original author artifacts
are `not_yet_supplied`.

The [next-data decision](evidence/phase_4j64/next-data-requirements.json) separates
two immediate questions. A frozen future four-gradient check at Si–C(6–7) and
C–H(7–19), each ±0.003 Å with the retained center, would test the shared
curvature corrections at other graph-screened sites. These are other sites in
the same conformer, not independent conformational basins. Before full-force
fitting, one matched QZVP center gradient at 250×1454, with unchanged geometry,
auxiliary basis and standard SCF thresholds, would isolate QZVP grid-force
sensitivity from the observed orbital-basis contrast. It needs its own resource
preflight and frozen plan. Neither experiment is executed here. Subsequent
joint fitting still needs a qualified single-protocol objective, the five
training proposals and the three untouched validation configurations; it must
not infer force qualification from curvature-only calibration.

## Reproduction and preservation

Offline reconstruction needs the retained native archives and installed numerical
Python dependencies, not a QM executable or external scratch. It unpacks into a
fresh directory and independently reconstructs native values, protocol contrasts,
curvature calibration and the geometry-only design:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  scripts/replay_phase_4j64.py replay
```

For numerical reruns, copy the frozen evidence into an isolated output directory.
Unpack the retained native archive with the existing `bundle.unpack` helper into
a separate temporary source directory; copy only its `scientific-plan.json`,
`jobs.json` and `requests/` into a **new**, empty external work directory. The
archive's job map includes the additive exact-center request. Leave its original
attempt ledgers, results and process records in the archive; a new campaign gets
its own ledgers. Do not overwrite the retained archive. The staged commands are:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/phase_4j64.py dense \
  --evidence /absolute/copied-evidence --work /absolute/new-work \
  --python ../island-validation/phase4j36-native-env/bin/python
# Insert the additive exact-center job after dense, before analyzing its result:
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python -c \
  'from pathlib import Path; import sys; from scripts.phase_4j64 import run_job; sys.exit(run_job(Path.cwd(), *map(Path, sys.argv[1:]), "exact/center"))' \
  /absolute/copied-evidence /absolute/new-work \
  /absolute/phase4j36-native-env/bin/python
# Then qz-center, stretch, qz-stretch as permitted by the frozen native/resource gates.
# Analyze the dense triple before any conditional torsion-control activation.
```

Only exact hash-verified completed jobs resume. Failures and interruptions retain
nonzero statuses and require additive attempts; original outputs are not erased.
The existing J63 native capture/importer, bundle, history and JUnit contracts
are reused. Large scratch files stay outside Git.

Historical receipt validators with broad implementation globs run in an isolated
checkout of the exact base, so additive J64 modules cannot retroactively expand
a historical manifest. The current tree is separately checked byte-for-byte
against every base file. Original failures and retries remain separate from
terminal software coverage and scientific precision decisions.

## Software verification

Terminal coverage is 5,215 unique supported-suite cases: **5,154 passed and
61 skipped**, combining four disjoint base-suite shards and the final 42 J64
tests without counting earlier focused reruns twice. Applicable J54–J63
regressions are included. Focused controls exercise protocol separation,
independent native reconstruction and mutations, path-force signs and curved
paths, rounding bounds, unknown-null scalar/force calibration, correction
energy/force/parameter derivatives, incomplete coefficient rejection, objective
gates, budgets/failure retention, and fresh-process relocation.

Ruff, compileall, all five applicable environment pip checks, exact verified
resume without new QM, receipt/manifest validation and diff checks pass.
Historical J54–J63 receipt validators pass in the isolated exact-base checkout;
all 6,950 base files also match the current tree byte-for-byte. Raw native files,
requests, normalized records, checksums, commands, process statuses and resource
records are retained. Large scratch stays outside Git.

Early lint/helper failures, the missing timing-command startup, and an
evidence-copy helper's Path/string error remain visible alongside their
corrections. Initial staged diff checks identified native trailing whitespace
and then the diagnostic quoting it. A local J64 attribute exempts only those
immutable raw logs; code and other text keep normal whitespace checking.
The terminal suite is a successful final run, not a claim that
every earlier command passed. Native calculations and software checks have
separate attempt records. The final receipt validator reconstructs the numerical
results and checks exact JUnit coverage against terminal test collection.

`production_validated=False`, `simulation_readiness="not_established"`.
Parameter completeness, backend correctness, optional cross-engine comparisons
and scientific validation remain separate. Exact LAMMPS/OpenMM trajectory
replication is not a prerequisite for this study. No universal PCFF/COMPASS,
complete FRC/PDMS coverage, production MD, pressure, Ewald/PPPM or long-NPT claim
is made.
