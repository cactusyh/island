# Phase 4J62 — bounded angular/torsional control and HMDS entry

This is a molecular force-field research experiment. It stops the Si–H control
at a bounded milestone and begins the exact methyl-siloxane graph independently.
It neither registers a complete PCFF-PDMS profile nor promotes research trials
into automatic assignment. The ordinary API remains `forcefield="pcff"`.

## Base, immutable inputs and plan

Fetched `origin/main`, and the branch base, is
`1604a19f4f6245a5e35b302467c519aa62e7814d`. It contains J61
`b4cae05e3666d6adb232e0bad711db85292994c7`. The branch is
`codex/phase-4j62-angular-torsional-control-hmds-entry`; main is not merged or
pushed. All 6,208 historical base files, including J17–J61 evidence and source
implementations, are preserved byte for byte.

J61's native gradients, projection, fit and diagnostics were reconstructed
before new QM. The [staged plan](evidence/phase_4j62/plan/scientific-plan.json)
SHA-256 is `2fbb54e771e865d233b004153eb37799acca94f6aa003cfc956522b641f5e5d5`.
Its hypotheses separate objective tradeoff, missing parameter directions,
functional-form/model discrepancy and reference precision. HMDS entry does not
depend on a successful Si–H fit.

The cap is **48 new energy/gradient evaluations**: 16 primary angular/torsional
responses, eight representative precision controls, at most 12 HMDS optimizer
callbacks, six matched HMDS responses, three HMDS basis controls and three
conditional fresh validation points. Numerical derivatives and failed or
interrupted calls count individually. No hidden numerical QM Hessian job or
889-geometry campaign is permitted. Limits are two concurrent processes,
four threads and 4 GiB per process, 20 GiB scratch, 1,200 seconds per gradient,
9,000 seconds for the capped optimizer and 28,800 seconds campaign wall time.
The resource choices use retained native timings. Original failures and retries
are retained separately.

One planning error was corrected additively **before control-response inspection
or fitting**: the source n3 torsion amplitude is signed, −0.13 kcal/mol. The
original correction bound `[-source_n3,5]` excluded the zero extension. The
[amendment](evidence/phase_4j62/plan/signed-torsion-bounds-correction.json)
preserves that original text and declares correction bounds [−5,5] kcal/mol.
No acceptance threshold was changed.

## Existing-fit weighting audit

The [training-only audit](evidence/phase_4j62/training-weighting-tradeoff.json)
reconstructs J61's weighted SS of 2,028.103278: 164.318603 from relative energies,
1,742.929303 from Cartesian forces, 118.106266 from stationarity and 2.749107
from projected curvature. Neither J60 nor J61 validation errors select weights.

| Stationarity penalty multiplier | Stationary force RMS (kcal/mol/Å) | Stretch curvature relative error | Training energy RMSE (kcal/mol) | Training force RMS (kcal/mol/Å) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.542042 | 0.020622 | 0.136648 | 0.856479 |
| 10 | 0.030785 | 0.018596 | 0.195405 | 0.979528 |
| 100 | 0.003319 | 0.018396 | 0.204133 | 0.991064 |

The existing 24-variable model therefore supplies a feasible stationarity and
curvature witness. J61's soft-penalty optimum is not proof of infeasibility.
Stronger stationarity worsens training energy and force errors: objective
tradeoff explains part of the stationarity failure, while substantial residual
model discrepancy remains. The multiplier-1,000 diagnostic was rejected for
scaled condition 3,377.56; its failure is retained. J61's values, gates and
rejected scientific decision remain unchanged.

## Angular/torsional experiment and research models

The reference is the exact J61 stationary candidate at PBE0/cc-pV(5+d)Z,
DF/`def2-universal-jkfit`, unpruned 150 × 974 grid, SCF thresholds 1e−10/1e−8,
integral threshold 1e−12, and no dispersion. Psi4 1.9.1 and the pinned existing
environment are reused. Actual contracted orbital shells and basis files are
recorded in [the inventory](evidence/phase_4j62/engine-environment.json).

The four seed coordinates are Si–O–Si (1,2,3), H–Si–H (4,1,5), H–Si–O
(4,1,2), and silyl torsion (4,1,2,3). Their Cartesian coordinate gradients are
projected against six rigid modes and the rank-eight stretch subspace, then
positive-diagonal QR gives rank four. The original gradients, projection and
coordinate-to-Cartesian mapping are retained. Coordinates are unweighted
Cartesian, not mass weighted. Displacements are ±0.004 and ±0.002 Å along
normalized Cartesian directions.

Every response retains the full 27 × 4 gradient derivative, the unsymmetrized
combined 12 × 12 projection, stretch/nonstretch cross-block reciprocity, energy
second differences and two-step comparisons. Representative central-bend and
torsion directions test QZ versus 5Z; torsion tests grid and SCF precision
separately. Frozen absolute-plus-relative error floors avoid dividing a soft
near-zero torsion by itself. These are molecular total-QM versus total-MM
responses, not isolated source coefficients. Only 12 of 21 internal directions
are sampled; positive curvature there cannot establish a full minimum.

The separate research identity is
`island_research_class2_sih_disiloxane_j62_nested_v1`. Nested models are:

- M0: J61's 24 stretch/baseline/coupling variables;
- M1: M0 plus one shared signed n3 amplitude correction;
- M2: M1 plus H–Si–H and H–Si–O equilibrium-angle shifts;
- M3: M2 plus central Si–O–Si quadratic stiffness correction.

Charges, exclusions, nonbonded settings, torsional phases, central equilibrium
angle, local angle stiffnesses, higher-order source polynomials and all other
source rows remain fixed. Angular equilibrium shifts propagate through quartic
angles, fixed bond-angle and fitted angle-angle, angle-torsion and
angle-angle-torsion terms, including Cartesian and parameter derivatives.
No Si–H trial is treated as an author-confirmed PCFF source row.

Training uses the 22 existing J61 training configurations. The energy/force
weights remain 0.05 kcal/mol and 0.5 kcal/mol/Å; the stationarity block uses the
training-only feasible multiplier 10 and the combined curvature block is treated
as one correlated constraint. Finite-difference gradient pairs are not also
counted as independent statistical observations. There is no regularization.
The actual combined Jacobian, column-scaled rank/condition, correlations,
absolute sigma proxies (also near zero) and active bounds are retained.

Historical energy/force gates remain 0.1 kcal/mol RMSE and 0.5 kcal/mol/Å RMS.
Stationarity must be ≤0.05 kcal/mol/Å, combined curvature error ≤5%, each full
nonstretch response column must satisfy its absolute-plus-relative gate, rank
must be complete, scaled condition <1,000 and no bound may be active.
Whole-configuration bootstrap uncertainty and three fresh points are activated
only for a training-qualified model. The reserved points use seed 4062 and are
checked under all 72 graph automorphisms and proper rigid alignment against
all J60/J61 development geometries. They are new deformations, not independently
sampled conformational basins. A single training-derived energy offset applies
to future validation; validation does not select new families or weights.

All **24 Si–H gradient calls** completed: 16 primary and eight precision controls.
The two-step full-response column RMS differences are 0.00129, 0.00237,
0.000737 and 0.00397 kcal/mol/Å². The representative QZ/5Z differences are
0.16940 and 0.11228 kcal/mol/Å²; torsional grid and tighter-SCF differences
are 0.05190 and 0.000303 kcal/mol/Å². Every frozen precision gate passes,
including cross-block reciprocity and energy/gradient curvature checks.
This supports this bounded response experiment, not full-Hessian or final
method adequacy. See [the actual responses](evidence/phase_4j62/angular-torsional-control.json).

| Model | Variables | Training energy RMSE (kcal/mol) | Training force RMS (kcal/mol/Å) | Stationary force RMS (kcal/mol/Å) | Combined curvature relative error | Scaled Jacobian condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| M0 | 24 | 0.195419 | 0.979502 | 0.030798 | 0.026352 | 31.27 |
| M1 | 25 | 0.154935 | 0.762180 | 0.012684 | 0.028956 | 31.10 |
| M2 | 27 | 0.154027 | 0.744784 | 0.011382 | 0.030685 | 76.71 |
| M3 | 28 | 0.087860 | 0.663636 | 0.008762 | 0.029745 | 79.46 |

All four optimizers converge with complete local rank and no active bounds.
Main-torsion release reduces weighted SS by 39.80%; local equilibrium release
adds 4.10%, and central stiffness release adds 27.04%. These are successive,
training-only comparisons, not independent causal fractions. M3 passes energy,
stationarity and aggregate curvature gates, but fails force RMS and all four
full nonstretch-response gates. Its column RMS errors are 4.11993, 2.65776,
5.51080 and 5.25963 kcal/mol/Å², versus allowed 2.79797, 2.01272, 2.37615
and 2.19677. Aggregate stretch-dominated curvature error conceals these
angular/torsional response mismatches.

The rejected M3 extension values and absolute local sigma proxies are:

| Extension | Trial value | Absolute sigma proxy | Units |
| --- | ---: | ---: | --- |
| Shared n3 amplitude correction | 0.170241 | 0.006030 | kcal/mol |
| H–Si–H equilibrium shift | 0.006286 | 0.002478 | radian |
| H–Si–O equilibrium shift | −0.002356 | 0.002640 | radian |
| Central Si–O–Si K2 correction | −5.684828 | 0.354850 | kcal/mol/radian² |

The near-zero H–Si–O shift is poorly resolved even by this local proxy.
Full parameter correlations and per-configuration energy/Cartesian-force errors
remain in [the trials](evidence/phase_4j62/nested/M3.json). These proxies omit
method/source/model systematics and correlated-measurement effects; they are
not confidence intervals or acceptance evidence. No model qualifies for
bootstrap acceptance or the three reserved fresh QM points, so neither is run.
Historical validation is not reused as fresh validation.

**Zero coefficients are accepted.** The experiment supports a torsional
amplitude contribution and additional angular sensitivity, but these selected
directions do not remove the model discrepancy. Neither reference precision
failure nor local rank deficiency explains the remaining tested response
errors. This does not prove global infeasibility or uniquely establish which
fixed baseline or functional-form assumption is responsible. The Si–H study
ends here; another expansion of this control is not the next dependency for
methyl-siloxane work. See [the diagnosis](evidence/phase_4j62/control-diagnosis.json)
and [the fitting decision](evidence/phase_4j62/fitting-decision.json).

## Exact methyl-domain entry

The audited graph is **27-atom hexamethyldisiloxane**, Me3Si–O–SiMe3:
two Si, one O, six C and 18 explicit H; 26 single bonds; neutral singlet.
Stable IDs and original atom ordering are preserved. The graph identity is
`d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`.
Each `sio` has one `osi` and three `c3` neighbors; each methyl `c3` has one
`sio` and three `hc` neighbors. Initial geometry, valence, separations,
source-specific typing and baseline rows are audited in
[the domain record](evidence/phase_4j62/HMDS-domain-audit.json).

The pinned PCFF candidate bond-increment charges are complete without
normalization: Si +0.52, O −0.44, methyl C −0.259 and H +0.053 e; their sum
is numerical zero. Missing classical couplings do not prevent a gas-phase QM
calculation with known formal charge and multiplicity, but do prevent a complete
MM residual fit. There are 378 missing interaction sites: 78 angle-angle and
60 each of first/third bond-bond, end/middle bond-torsion, angle-torsion and
angle-angle-torsion. Eighteen angle-angle sites are source-complete. The 26
directed missing queries are not asserted to be 26 independent coefficients.
No OPLS/GAFF rows or hsi/Si–H fitted numbers are imported.

The exploratory protocol is PBE0/def2-SVP, DF/`def2-universal-jkfit`, 99 × 590
grid, SCF 1e−10/1e−8, integral threshold 1e−12, no dispersion. The 12-gradient
optimization cap is reached without convergence. Its last successful geometry
is **not stationary**. Matched Si–C (1,2), C–H (1,10) and methyl-torsion
(3,2,1,10) responses use ±0.003 Å; TZVP center/Si–C calculations use the exact
same geometries. This single-step exploratory contrast does not establish
Hessian convergence, final QM adequacy or methyl coefficients.

The retained final Si–C distances are 1.87455–1.87981 Å, C–H distances
1.10113–1.10268 Å and Si–O distances 1.66121–1.66215 Å. The closest nonbonded
pair is H13/H14 at 1.77816 Å. These are graph/geometry compatibility checks,
not a quantum bond-order determination or a minimum proof.

The 21 HMDS gradients reconstruct independently. The final SVP gradient RMS/max
are **3.52643e−4 / 1.20567e−3 hartree/bohr**, above the declared stationary gates.
At the same geometry, the TZVP/SVP center-gradient RMS difference is
**3.84427 kcal/mol/Å**. The Si–C full-Cartesian-response difference is
**6.32551 kcal/mol/Å² (4.66702%)**; relative deformation energies differ by
−0.0157989 / +0.0154107 kcal/mol for the positive/negative displacement.
No final accuracy gate is invented for this exploratory dataset.

| HMDS Cartesian seed direction | Gradient-response diagonal (kcal/mol/Å²) | Energy second difference (kcal/mol/Å²) |
| --- | ---: | ---: |
| Si–C | 1207.66853 | 1205.91420 |
| C–H | 1679.27825 | 1678.55748 |
| Methyl torsion seed | 973.09340 | 960.58895 |

These are total molecular Cartesian responses at an unfinished center. In
particular, the dihedral-gradient seed changes other internal coordinates; its
large response is not a methyl torsion coefficient. The single-step energy/
gradient curvature discrepancies and basis contrast require numerical/reference
controls before inference. HMDS SVP uses 224 functions/112 shells, TZVP uses
399/175; the measured cost supports a bounded TZVP continuation rather than
imposing the nine-atom 621-function 5Z protocol on this larger graph.

## Replay, failures and scientific scope

Native inputs, outputs, analytic gradients, engine settings, atom ordering,
hashes and subprocess diagnostics are retained separately from normalized
ISLAND records. Independent native-text/ASCII conversion checks verify energies,
full gradients, force sign and units. Scratch and wavefunctions stay outside Git.
Point-gradient requests retain some inherited, inactive optimization and
finite-difference options. The actual Optking settings and twelve-call cap are
bound by the separately retained instrumentation; see
[the protocol-context qualification](evidence/phase_4j62/native-protocol-context.json).
No numerical QM Hessian was requested.
The dimension-independent coordinate/response path handles HMDS while historical
nine-atom importers and reconstruction identities remain unchanged.

The campaign executes **45 new gradients**, all retained and independently
reconstructed: 24 Si–H and 21 HMDS. No QM gradient fails; the 12-call HMDS
optimization attempt remains explicitly unfinished. Its subprocess return
code is zero because the bounded optimizer exception is caught and preserved,
not because optimization converged. The native traceback, optimization log,
per-callback inputs and full gradients are retained. Gradient wall times range
from 31.71 to 515.39 seconds; summed per-gradient elapsed time is 7,056.23
seconds across the two-process campaign. Three reserved calls remain unused.

Replay from the tracked native archive without QM or refitting:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  scripts/replay_phase_4j62.py replay
```

The original execution command also resumes completed jobs by exact hash
verification; it does not overwrite failed or interrupted attempts:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/phase_4j62.py \
  --work ../island-validation/phase4j62-qm \
  --python ../island-validation/phase4j36-native-env/bin/python
```

The exact executed coordinator source is archived. A subsequent additive
maintenance record explains its intentional catch-all failure collector to
Ruff and corrects an optimizer output-directory reporting pointer. Original
process records and scientific settings are unchanged; see
[the amendment](evidence/phase_4j62/plan/coordinator-maintenance-amendment.json).

The 28 confirmations remain **user-attested author confirmations**;
`original_author_artifact_status="not_yet_supplied"`. Their interpretation
authority does not invent missing numerical coefficients. See
[the scope and next-data map](evidence/phase_4j62/chemical-scope-and-next-data.json).

The measured basis-dependent gradient makes completion of SVP optimization alone
insufficient. The smallest next numerical dependency is a separately budgeted
PBE0/def2-TZVP stationary-candidate optimization starting from the retained
geometry, followed by two-step Si–C/C–H and methyl/Si–O torsional responses with
matched reference-precision controls. TZVP is a cost-informed candidate, not an
established final protocol. Binding the complete methyl coupled MM dependency
model and independently reserving joint deformations remains necessary for an
identifiable fit, independently of whether this QM optimization can execute.
See [the numerical next-experiment decision](evidence/phase_4j62/HMDS-next-experiment-decision.json).
Capped-DP2 is a separate
37-atom domain: interior silicon with two `osi` and two `c3` neighbors introduces
O–Si–O and Si–O–Si–O environments absent from HMDS end groups.

No Si–H coefficients transfer automatically to Si–C. Cyclic, crosslinked,
Si–OH, Si–O–C and phosphorus systems remain unsupported. No universal PCFF,
universal COMPASS, complete `.frc`/PDMS coverage, pressure validation,
Ewald/PPPM closure, long-NPT readiness or production MD is claimed.
Parameter completeness, backend execution, optional cross-engine equivalence
and scientific validation remain independent under J55. Throughout:
`production_validated=False`, `simulation_readiness="not_established"`.

## Software verification

The supported suite covers **5,058 unique source-collected cases**: 4,997 pass
and 61 skip. Eight file shards cover all 5,001 historical cases; the 57 focused
J62 tests all pass, including the real nine-/27-atom importer, native mutations,
manufactured nonzero parameter recovery, combined Jacobian rank loss,
angle-equilibrium propagation, rigid covariance, failure retention, validation
leakage, relocation and fresh-process replay. This includes 531 J54–J61
regression cases. No original full-suite shard fails or requires a retry.

Ruff, compileall, diff checks, five applicable environment `pip check` commands,
J54–J61 receipt validation and the historical hash check pass. All 6,208
historical base files are unchanged. Exact commands, JUnit coverage, skipped
cases and retained subprocess diagnostics are in
[the verification records](evidence/phase_4j62/validation/checks.json).

Earlier failures are preserved in
[the attempt ledger](evidence/phase_4j62/validation/early-attempts.json): one
focused test incorrectly aggregated a per-row helper by family, the first basis
inventory used the wrong installed filename aliases, and the initial Ruff run
required an explanation of intentional failure capture. Their original logs
and successful corrections remain distinct. The basis-inventory failure did
not consume a QM evaluation. The first staged diff check flagged padding in
retained failure/native stdout. A narrowly scoped local `.gitattributes`, as in
existing evidence directories, preserves those exact bytes while keeping
implementation and documentation whitespace checks active. The original diff
failure and passing retry are both retained. The unfinished HMDS optimizer is
also preserved as an unfinished scientific attempt, independently of passing
software tests.

Validate the receipt and manifest, native reconstruction and historical bytes:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  scripts/validate_phase_4j62_receipt.py docs/evidence/phase_4j62 --repo .
```
