# Phase 4J59 — physical residual accounting and QM convergence

J59 audits the **9-atom disiloxane Si–H control**, executes a bounded Psi4
convergence study, and makes a numerical **no-fit** decision. It does not create
a PCFF-PDMS profile or run PCFF MD. The ordinary `forcefield="pcff"` interface is
unchanged. Parameter completeness remains `incomplete`, independently of the
optional, `not_claimed` LAMMPS/OpenMM equivalence diagnostic. Historical backend
regression tests remain separately scoped from the J59 control/PDMS candidate
numerical campaign; the no-MD claim applies to that candidate domain.

The exact fetched `origin/main` and branch base is
`d361f7e0a7c16a3a3a3513987971db6127690386`. It contains J58
`361aaa3eabe77df9bdbc7bdbcba5a444807c4a40`, J57
`b56659d13d5e0502ebf57a7c0fc56bdc1fdb5df2` and J56
`24ad7698f64d62aeb08ed893eb8372bdbbf33819`. No prerequisite is assumed merged;
`plan/ancestry.json` records the checks. This branch does not merge main. The
complete base tree is the byte-preservation baseline, including the original
J58 failures, normalization revisions and input bundles.

## Complete accounting, deliberately partial subtraction

`interaction-table.json` enumerates every source assignment, dependency and
potentially applicable remainder, retaining coefficients or explicit missing
values, units, roles, source rows/equations and provenance. There are 74 fixed
records (including 27 typing/mass/charge records), 115 unresolved records and 8
structurally inapplicable Wilson impropers. A separate unresolved physical
applicability record prevents that converter rule from becoming source omission
authority. These are interaction records, not
115 independent transferable parameters.

| Family | Site records | Treatment |
| --- | ---: | --- |
| Quartic bond | 8 | Fixed source coefficients and equilibrium distances |
| Quartic angle | 13 | Fixed source coefficients; angular displacements in radians |
| Adjacent bond–bond | 13 | Fixed, with two source bond equilibria |
| Bond–angle | 13 | Fixed, with bond and angle equilibria |
| Main proper torsion | 6 | Numerical rows exist; source Phi/phase-to-geometry binding remains explicit |
| Angle–angle | 24 | Missing numerical coefficients, source role symmetry retained |
| Bond–bond 1–3 | 6 | Missing numerical coefficients, path reversal retained |
| End/middle bond–torsion | 6 each | Missing numerical coefficients, both end blocks explicit |
| Angle–torsion / angle–angle–torsion | 6 each | Missing coefficients; literal source equations retained |
| Torsion–torsion | 9 | H–Si–O–Si–H paths exist; empty FRC section provides no absence authority |
| Atom nonbonded rows | 9 | r-min/epsilon supplied; pair execution policy unbound |
| LJ9–6 and Coulomb pairs | 36 | Every pair, graph distance, charge product, missing pair weights and policies explicit |
| Wilson improper | 8 | Existing converter's degree-four structural inapplicability, clearly labeled |

The fixed evaluator implements the four fully bound bonded families, without
one-half prefactors or harmonic reduction. It supplies analytical Cartesian
forces and separately returns every fixed row energy. Tests independently use
atan2 angle geometry, Horner scalar polynomials and finite differences at two
steps, plus manufactured small-system controls for every fixed family. Net force
and torque vanish. The source LJ9–6 mixing/kernel and literal scalar plus-cos
source torsion equation also have analytical controls; neither selects an
unbound pair weight or atom-to-Phi mapping. Asking for a complete residual
raises an error listing the actual unbound records.

The 28 confirmations remain unchanged **user-attested author confirmations**:
`author_confirmed_audited`, original author artifacts `not_yet_supplied`.
They remain effective interpretation records. No absent coefficient, author
email, author date or primary document is invented. Their numerical values
remain absent. The 19 directed numerical queries remain `qm_fit_required`.

## Joint model and fitting decision

The model is explicitly

```text
QM = verified fixed source contribution
   + explicit unknown contributions
   + model discrepancy
```

`residual-identifiability.json` records every affected coordinate/dependency,
fixed energy and force subtraction over all 50 historical scan points, and
conditional source-equation nuisance Jacobians. The conditional 31 columns
include the four target columns, 16 additional torsional-coupling columns plus a
separate linear-Phi phase-mapping stress column, a main-torsion mapping sensitivity, a
source-empty torsion–torsion column, and eight lower-bound pair-policy columns.
The latter group pairs by topological separation **only as a diagnostic
hypothesis**. This does not authorize exclusions or tie source pair policies.
The main torsion rows already supply a first phase of zero; it is not an absent
numerical value. The report also evaluates the 30-column alternative linking
AAT `Phi1(0)` to that phase, without the extra phase-mapping stress column.
Neither alternative binds the atom-to-Phi geometry. Linear-Phi evaluations use a continuous local branch anchored to the historical
center, with no claim of source-authorized Phi geometry.

The energy design has rank **22/31**. Adding Cartesian forces gives conditional
rank **31/31**, column-scaled condition **868.13**. This diagnostic concatenates energy in kcal/mol and force
multiplied by 1 Å, then normalizes columns; it is not a selected fitting objective
or an uncertainty estimate. End-bond and middle-bond n=3
columns have cosine correlation **0.99994654**. Thus the historical four-column
full rank is insufficient, and full rank in this conditional expansion still
cannot establish physical identifiability. Pair weights, Coulomb
prefactor/dielectric, cutoff policy, torsion mapping and empty-section
applicability must first be bound. Unknown terms are never subtracted as zero.
The report quantifies per-unit nuisance energy/force sensitivities and unexplained
fixed-subtraction energies/forces; these are **not fitted errors**.

No fit is performed, and no coefficient uncertainty or transferable fitted value
is invented. Before inference, the smallest source information needed is the
numerical hsi–sio–osi–sio end/middle bond–torsion, angle–torsion and
angle–angle–torsion rows with phase/equilibrium dependencies; an explicit
H–Si–O–Si–H torsion–torsion inclusion/absence rule; and exact intramolecular
pair policies. Once bound, separate Si–O–Si torsion and center-angle deformations
can target the nearly correlated columns. The two original grids share one
conformer; a distinct held-out conformer is required. Neither a grid column
split nor corner holdouts supply that validation. Fit weights/bounds have not
been selected; no fit is authorized before complete-model and reference gates.

## Predeclared convergence study and native evidence

`plan/convergence-plan.json` predates execution. It selects the common J58
center and all four sign combinations of each original two-coordinate scan:
±0.08 radian angle displacements and ±0.04 Å bond displacements. Graph, atom
order, neutral singlet and all 27 Cartesian coordinates are preserved at every
comparison. Protocols are kept separate rather than mixed into one fit dataset.
The method/basis/grid/SCF contrasts and accuracy criteria are recorded before
results. PBE0 without dispersion remains a proposed pilot protocol.

The original QZ/PK preflight produced approximately 73 GiB of integral scratch
for two incomplete jobs. Those native attempts are retained in
`native-pk-attempts.tar.xz`; only regenerable external integral caches were
discarded. `plan/df-amendment.json` changes later work to explicit DF with
`def2-universal-jkfit`, adds three PK/DF TZVP controls and preserves the original
plan. `plan/cpu-amendment.json` retains the two completed two-thread DF preflight
jobs, then permits two simultaneous eight-thread jobs, 4 GiB each, with an
additional paired thread-count control. Limits remain 1800 seconds/job and
14400 seconds campaign wall time. No shared scientific environment is modified.
Executed preflight adapter sources are compressed and hash-bound in `engine/`;
the current adapter is additive. Other QM engine adapters remain unexecuted.

The study compares PBE0/def2-TZVP and def2-QZVP at nine common geometries,
99×590 versus 150×974 unpruned Treutler grids, SCF energy/density thresholds
1e-10/1e-8 versus 1e-12/1e-10, B3LYP versus PBE0 at the same QZ basis/grid,
and PK versus DF separately. Relative deformation energies, full and
center-subtracted gradients, and energy/gradient mixed scan secants are retained.
Mixed secants are total-QM path derivatives, not missing cross coefficients.
A separate bounded QZ/dense-grid center optimization diagnoses center effects;
its identity is distinct and its data are never blended with the old-center
scan. The original J58 data are read, not rewritten.

Native inputs/logs, indexed gradient text, lossless arrays and hashes remain
separate from normalized ISLAND records. The settings-aware importer reuses
J58's verified completion/matrix infrastructure and J57's explicit SI conversion
contract. Forces are minus gradients. Native geometry uses Psi4's recorded
Bohr constant; observable conversion retains J57's distinct documented constant.
Independent Decimal/native-text reconstruction uses neither production importer
nor production residual evaluator. Convergence subsets are **not** represented
as complete fitting datasets and do not bypass J57 deformation-coverage gates.

## Scope and reproduction

Disiloxane's `sio/osi/hsi` control does not establish methyl-PDMS transfer.
HMDS requires distinct `sio`–`c3`/`hc` environments, methyl internal coordinates,
charge dependencies and couplings. DP2 adds interior Si with two O neighbors and
different chain/end-group couplings. These require separately bound targets,
reference data, conformers and validation. No cyclic, crosslinked, Si–OH,
Si–O–C or phosphorus scope is inferred.

Use the existing validation interpreters and pinned FRC distribution (paths below
are this session's declared environments):

```bash
export PYTHONPATH=.:src
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
QM_PY=../island-validation/phase4j36-native-env/bin/python
CHECK_PY=../island-validation/phase4g1-env/bin/python
FRC=../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files
WORK=../island-validation/phase4j59-qm-convergence-fast

# Existing completed jobs are verified and skipped; failed outputs never overwritten.
$CHECK_PY -m island.pcff_qm_pilot.convergence execute --work "$WORK" --python "$QM_PY"
$CHECK_PY scripts/phase_4j59.py analyze --distribution "$FRC" --work "$WORK" \
  --j58-work ../island-validation/phase4j58-qm-pilot
$CHECK_PY scripts/phase_4j59.py residual --distribution "$FRC"
$CHECK_PY scripts/phase_4j59.py replay --distribution "$FRC"
$CHECK_PY scripts/validate_phase_4j59_receipt.py docs/evidence/phase_4j59 \
  --distribution "$FRC"
```

For reconstruction without live scratch, extract `native-outputs.tar.xz` to
`$WORK` using `island.pcff_qm_pilot.bundle.unpack`; it contains all exact requests
and completed outputs. Extract J58's immutable archive for the comparison path.
Replay and receipt validation extract their own temporary trees, support a
relocated evidence directory, and do not execute QM. Re-execution requires a new
empty work tree containing the retained request files, the pinned environment
and the recorded adapter. The saved manifest limits execution to declared jobs.
Fresh execution dispatches the two original two-thread requests to the retained
v2 source and the eight-thread requests to v3. Both source identities are checked
against executed native metadata. The DF campaign binds the decompressed v2
source checksum; `engine/preflight-adapter-hashes.json` separately binds the
compressed archive checksum. The dispatch repair changes no QM request,
executed adapter or numerical result.

Scientific validation and production readiness remain `not_established`,
`production_validated=False`, `simulation_readiness="not_established"`.
No PCFF MD, final PCFF profile, pressure/Ewald/PPPM validation, production
equilibrium or long-NPT readiness is claimed.

## Executed numerical results

The retained campaign contains 39 completed energy/full-gradient calculations, including the paired thread repeat, and a separately recorded center optimization. Two initial QZ/PK resource-preflight attempts were interrupted and preserved.

| Contrast | Maximum relative deformation-energy difference (kcal/mol) | Maximum gradient RMS difference (kcal/mol/Å) | Predeclared gate |
| --- | ---: | ---: | --- |
| SVP → TZVP, with separate DF/PK control | 0.88739709 | 4.4145081 | Sensitivity only |
| TZVP → QZVP | 0.10373436 | 0.55301832 | Failed |
| Grid 99×590 → 150×974 | 8.0384896e-05 | 0.0042072374 | Passed |
| Tighter SCF thresholds | 1.5931437e-07 | 1.1675701e-05 | Passed |
| DF → PK | 0.00096894102 | 0.0080915953 | Passed |
| PBE0 → B3LYP | 0.097541696 | 1.8536559 | Sensitivity only |

The TZVP/QZVP objective gate **fails**: the maximum deformation-energy difference exceeds the predeclared 0.05 kcal/mol limit, and the maximum gradient RMS difference exceeds 0.5 kcal/mol/Å. The full signed tables retain every relative energy, gradient difference and total-QM energy/gradient mixed secant. Method sensitivity is not a final adequacy test. The technical thread repeat has energy difference -1.59664e-07 kcal/mol and gradient RMS difference 4.23514e-06 kcal/mol/Å.

The center optimization status is `completed`. At unchanged QZ/dense-grid/DF settings its energy changes by -0.5855077 kcal/mol; force RMS changes from 4.332214 to 0.0003812952 kcal/mol/Å. The new versioned targets/scans are generated separately and **unexecuted**. No new minimum-stability Hessian is claimed.

The concrete fitting decision is **no fit**: verified fixed terms do not close the physical remainder, source/mapping/pair policies remain unbound, and no independent conformer validation is available. The reference gate result is retained without changing tolerances. `fitting-decision.json` links both demonstrated obstacles and the next bounded source/data requirements.

The total-QM mixed energy secants change from 31.96459 to 31.78886
kcal/mol/radian² for the angle scan and from −3.777664 to −3.819778
kcal/mol/Å² for the bond scan on going from TZVP to QZVP. Their smaller
changes do not override the failed energy/gradient gate or identify isolated
MM cross terms. Center-subtracted gradient differences are retained separately.
The QZ optimization changes Si–O distances from approximately 1.6554 to
1.6309 Å. This is a same-protocol center effect, separate from the same-geometry
basis contrast; the revised-center deformation curves have not been executed.

The optimization trajectory is reconstructed from the actual native SCF geometry
tables at their printed 12-decimal precision. Analytic calls in this environment
do not print Input QCSchema blocks. The final full-precision exported geometry,
energy and gradient are independently bound to the native tables and arrays.
Step 20 meets all four active force/displacement criteria. The parser repair
changes analysis only; executed adapters and native numerical evidence are
unchanged.

## Software verification

The supported full-suite collection contains **4900 unique nodes**: **4839 passed and 61 skipped** after one unchanged isolated retry. The focused J59 file passes all **53 tests**. The actual J54–J58 regression subset contains 377 nodes (0 skipped); it is extracted from those executions without double counting.

The first full-suite attempt failed the historical
`test_real_final_graph_source_completeness_and_exact_execution[crosslinked-gaff]`
case with an AmberTools charge-import completeness/consistency diagnostic. Its
unchanged isolated rerun passed. The exact original cause is unconfirmed because
the historical harness no longer retained its temporary native artifacts. The
original failed JUnit/log, successful isolated retry, and explicitly labeled
unique-outcome aggregate are all preserved. The validator accepts only that
exact replacement and does not erase or relabel the original event.

Ruff, compileall, five applicable environment pip checks, J54–J58 receipt
validators, independent native unit reconstruction, relocated fresh-process
replay, staged/combined diff checks and byte-for-byte historical preservation
checks pass. The J59 validator seals the recursive artifact/implementation
manifest and historical baseline, then replays native normalization and source
reconstruction. It is also run from a relocated directory in a fresh process
against the exact receipt SHA. Validation does not execute QM or authorize MD.
