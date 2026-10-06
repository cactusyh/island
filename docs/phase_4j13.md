# Phase 4J13 — executable converter-compatibility tranche

## Scope and repository state

Branch: `codex/phase-4j13-pcff-executable-tranche`. Base:
`19c55bb65b5a8d20c69c4f6e15e1cd828fa3b51e` (merged J12, PR #45).
Its complete tree equals reviewed J12 `99bd1a39af20e640f33d37099190730ca49ee882`.
No prerequisite merge, main modification, or historical record rewrite was performed.

Two retained **PSMILES final graphs** now have independently checked executable
models under a new, explicit converter-compatibility policy: PEO DP3 and thioether
DP2. They are still incomplete under strict FRC-source semantics. No new atom
predicate, increment, physical equation, numerical engine, source version, or
scientific accuracy claim is introduced.

`production_validated=False`; `simulation_readiness="not_established"`.
Full-source completion remains **false**. A converged minimum and four-step
trajectory check do not establish equilibration or polymer-property accuracy.

## Policy and primary evidence

The external source remains LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
No library or compiled executable is bundled.

The unchanged strict policy is `island_pcff_msi_guarded_source_v3`:
exact before wildcard, family/position-specific equivalence, legal reversals,
highest version per identical ordered source pattern, conflict rejection, and
an explicit source row for every active coupling. Historical v1/v2 interpretation
and identity payloads are unchanged.

New opt-in policy: **`island_pcff_msi_non_cp_compatibility_v4`**.
It uses the same guarded matcher and charge convention, with one narrowly
specified model-assembly exception:

* The family is BB13, the source search reports **missing**, and all four supplied
  labels are different from `cp`.
* Both terminal equilibrium bond parameters are assigned. Their numerical values,
  source identities and physical roles remain dependencies.
* The coefficient is zero, with origin **`converter_derived_zero`**, empty
  `source_rows`, and explicit converter evidence. Raw source coverage remains false.
* Assigned rows (including explicit source zero), ambiguous rows, cp at any of the
  four positions, missing dependencies, and other missing families are not replaced.

This is a guarded subset of converter behavior, not a general zero-fill rule.
The actual converter skips non-cp BB13 lookup altogether; this policy only permits
its initialization when the guarded source lookup is missing, preserving source
assignments and conflict rejection.

Audited pinned [GetParameters.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c):

* Lines 365–392 initialize BB13 K to zero and attach `rab`/`rcd`.
* Lines 506–525 perform BB13 lookup only if any supplied proper-torsion type is `cp`.
* Lines 1055–1169 implement exact/wildcard forward/reverse matching; native conflict
  rejection and center-preserving improper roles remain the J12 safeguards.
* Lines 1241 onward map ordinary family equivalences.

The policy stores revision, routine/line evidence, `GetParameters.c` SHA256
`add196c115ed9066aa5eba6da5bd8036bd29a3619395385140276e108c769591`,
and verified msi2lmp executable SHA256
`a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`.
These identify validation evidence; they do not create a runtime executable dependency.

LAMMPS executable SHA256:
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`
(30 Sep 2026, pinned revision build). Compiled matcher/ordinary-equivalence controls
were also rerun from the checksummed pinned C routines, without adopting incidental
file-order choices as physical authority.

## Reproduction, declarations, and reference independence

The pre-edit seven-case rerun reproduced **7 typed, 4 native-charge complete,
1 strict-model complete**. PEO lacked 42 BB13 requests; thioether lacked 27.
The pre-implementation [declaration](evidence/phase_4j13_declaration.json) retained
source/executable/input hashes, geometry, policies, budgets and tolerances.
Case-specific declarations are retained beside the external acceptance records.

No retained starting structure was regenerated:

| Case | Construction | Sites | Converter BB13 zeros | Independently selected unique source rows |
|---|---|---:|---:|---:|
| PEO | `[*:1]CCO[*:2]`, DP3, seed 2026 | 23 | 42 | 69 |
| Thioether | `[*:1]CCS[*:2]`, DP2, seed 2026 | 16 | 27 | 68 |

The graphs include actual chain ends and inter-repeat bonds. Existing graph
predicates resolve C/H/O/S environments (`c2`, `c3`, `hc`, `oc`, `oh`, `ho`, `sc`,
`sh`, `hs`) as applicable. This does not promote these labels globally.

`pcff_j13_reference.py` authors independent bounded graph labels, uses the separate
raw-FRC Decimal reader, and rebuilds inventories from authoritative bonds. It checks
selected source rows, oriented coefficients, native component charges and terminal
BB13 equilibrium dependencies. Independently supplied labels and charges enter
msi2lmp through validation-only CAR/MDF; **no ISLAND coefficient list enters the
oracle**. Actual converter BB13 coefficients are checked as zero with equilibrium
lengths reconstructed from the converter's own bond inventory. No `-ignore`.

The reference retains raw converter output and the established separate H5
angle-angle equilibrium-role corrections derived from its own angle inventory.
Degree-four centers retain all angle-angle couplings; Wilson terms are structurally
inapplicable there. The new policy does not authorize nonzero Wilson equilibria or
invent torsion–torsion parameters. All other H4/H5 forms, units, geometric roles,
special pairs, electrostatic constant and force definitions remain unchanged.

Source increments independently selected for PEO are lines **493, 505, 519, 810**;
for thioether **493, 505, 524, 794** (source namespace `cff91_auto`; namespace does
not by itself imply automatic-equivalence lookup). Both components are neutral.
The [coverage receipt](evidence/phase_4j13_coverage.json) lists every selected
bond/angle/torsion/cross-term/nonbonded row by family, full inventories and identities.
External `independent-selection.json` retains per-interaction sites and dependencies.

## Numerical and workflow acceptance

OpenMM Reference, finite nonperiodic systems, explicit LJ=(0,0,1) and
Coulomb=(0,0,1). Initial and asymmetric coordinates were compared; perturbation is
`0.015*sin(arange(3*N).reshape(N,3)+0.3)` Å. LAMMPS comparisons also cover the
minimum and retained dynamics steps 2 and 4.

Declared comparison tolerances: energy `1e-5 kJ/mol`, force
`1e-5 kJ/(mol*angstrom)`, rtol `2e-10`. Finite differences use displacements
`1e-4`, `1e-5`, `1e-6` Å; production forces are analytic, never finite differences.

| Measured result | PEO DP3 | Thioether DP2 |
|---|---:|---:|
| Max total/grouped energy difference, kJ/mol | 2.416e-13 | 6.395e-14 |
| Max force component difference, kJ/(mol·Å) | 2.229e-12 | 1.187e-12 |
| FD force error at 1e-4 Å | 5.205e-5 | 2.869e-5 |
| FD force error at 1e-5 Å | 5.260e-7 | 2.881e-7 |
| FD force error at 1e-6 Å | 4.738e-8 | 3.891e-8 |
| Initial / minimum energy, kJ/mol | 82.363466 / −23.583234 | −24.569147 / −83.757991 |
| Minimum maximum / RMS force, kJ/(mol·Å) | 0.0868384 / 0.0491022 | 0.0637809 / 0.0455887 |
| Iterations / actual evaluations | 201 / 210 | 143 / 150 |

Both minima are `force_converged` with independent fresh verification at the
unchanged 0.1 maximum-force criterion, 5,000 iteration / 10,000 evaluation limits.
No optimizer, line-search or verification tolerance was changed.

Each exact saved initialized state (300 K, velocity seed 78123) ran four BAOAB
steps, dt=0.1 fs, friction=5/ps, thermostat seed=99181. Uninterrupted: 6 evaluations,
5 frames. Split: 2+2 steps, 4 evaluations/3 frames per segment; cumulative 8
evaluations, 5 deduplicated frames. The directory was relocated and continuation
executed in a distinct Python process, without initialization/preparation reruns.

All frame differences were zero on this host (declared atol=1e-10, rtol=1e-12),
with exact full PCG64 state and trajectory-origin equality. This is integration
consistency, not independent integrator validation or a cross-platform equality
promise. Context counts: start/minimize/first segment **5**, uninterrupted **3**,
child continuation **3**, per case. Each minimum uses one reusable Context and
one fresh final Context. LAMMPS operations are separate reference checks.

Both bundles reject rechecksummed policy and model-source tampering. Completed
status, frames, and no-op resume passed again with OpenMM/RDKit/SciPy/Foyer/ParmEd
imports blocked. [Numerical receipt](evidence/phase_4j13_numerical.json).

## Charge and coupling limitations preserved

| Case | Outcome under new policy |
|---|---|
| Pyridinium | `nh+–cp` (two) and `nh+–hn` (one) increments absent under direct, ordinary and automatic searches. Ordinary row 276/automatic row 399 preserve the charged increment label. No `nh` substitution. |
| Neutral amine chain | `na–hn2` missing. `hn2` ordinary row 248 and automatic row 373 do not map to `hn`/`h*`. No silent alias. |
| Closed-shell protonated sidechain | Typed `n4/hn`; `n+–h*` missing after equivalence. The distinct source `h+–n+` row 813 is not authorization to change the historical hydrogen predicate. Exact blocked requests retained. |
| Guanidinium | Raw Decimal total **0.9999 e**, formal +1: residual **−0.0001 e**, fails 1e-12 e component criterion. |
| Chlorinated chain | Compatibility addresses 27 missing BB13 requests; 66 other required couplings remain missing: AA18, BB6, BA6, EBT9, MBT9, AT9, AAT9. Model remains ineligible. |
| Peroxide near-miss | Existing typing/charges complete; other BB/BA/torsional cross terms missing. Public preparation rejects. |
| Sulfoxide near-miss | S/O graph typing unresolved. Public preparation rejects. |

Guanidinium raw charge rows: c+–nr line **528**, endpoint values 0.2653/0.0680,
three bonds; h*–nr line **808**, 0.4068/−0.4068, six bonds.
A separately declared diagnostic converter run received these charges and retained
the exact vector and **0.9999 e** total in its data file. Its console rounds the
printed total; that is not a charge correction. No native model was published.
[Charge-echo receipt](evidence/phase_4j13_charge_echo.json).
Pinned [README limitations](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README)
exclude bond-increment support and auto-equivalence supplementation.
`ReadCarFile.c:192–227` reads supplied `q` and sums it; no corrective convention
was established.

The first protonated test `[*:1]CC[NH2+][*:2]` produced a terminal three-coordinate
N+ radical in the existing builder, correctly rejected. Its graph and diagnostic
are retained. A **separate** charged-sidechain declaration used
`[*:1]C(C[NH3+])C[*:2]`, DP1, seed2026; it is closed-shell but still charge-incomplete.
No original graph or failed receipt was replaced. The first control-report driver
attempt accessed `.assignments` on incomplete typing; its traceback is retained,
and the corrected diagnostic reader uses owned partial payload evidence.

[Gap receipt](evidence/phase_4j13_gaps.json) retains exact stable-site requests,
raw searches, source candidates, components and model failures. Metal/surface,
coordination, isotope and protonation context constraints are not relaxed.

## Coverage and integrity

For the seven retained inspection cases:

| Measure | Strict J12 before | Strict v3 after | Opt-in J13 compatibility |
|---|---:|---:|---:|
| Typed | 7/7 | unchanged | 7/7 |
| Native-charge complete | 4/7 | unchanged | 4/7 |
| Model complete | 1/7 | unchanged | 3/7 |
| Newly independently evaluated/workflow-verified | — | 0 | 2/2 declared positive cases |

Typing remains **86/133 labels with predicates**, **47 unresolved**; this is not
86 globally validated domains. There are zero newly introduced typing labels and
zero newly strict source-complete priority models. Full-source acceptance remains
nonzero. Source inventory is unchanged (133 atom labels, 134 ordinary equivalences,
108 automatic equivalences, 564 increments).

Native policy/model identities differ intentionally for the new policy. Historical
v1/v2/v3 models, source-derived coefficients, wrapper/bundle/checkpoint schemas and
the shared offline fingerprint derivation are untouched. Bundles retain model
interpretation, converter evidence and original preparation coordinates. Runtime
preparation/evaluation/reconstruction does not require CAR/MDF, converter, or LAMMPS.

Historical J8/J9/J10/J12 bundles and completed workflows were reconstructed
read-only alongside the two J13 workflows, with scientific imports blocked.
[Historical receipt](evidence/phase_4j13_historical.json). Snapshot hashes cover
historical J1–J12 evidence and retained files; preservation totals and test results
are in the verification receipt. No unrelated checkpoint directory was modified.

## Reproduction

The environment-local paths below are explicit inputs, not runtime search rules.
Set `PY=../island-validation/phase4e2-env/bin/python` and
`FRC=../island-validation/phase4h1-sources/lammps/pcff.frc`.
Use new output directories; publication protects existing records.

```sh
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" --output NEW/before
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" --declaration docs/evidence/phase_4j13_declaration.json --output NEW/after --require-full-source
# Expected exit 1: compatibility is not full source coverage.
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration docs/evidence/phase_4j13_peo_declaration.json --output NEW/peo
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration docs/evidence/phase_4j13_thioether_declaration.json --output NEW/thioether
$PY scripts/inspect_pcff_tranche_controls.py --source "$FRC" --output NEW/controls
$PY scripts/inspect_pcff_tranche_controls.py --source "$FRC" --declaration docs/evidence/phase_4j13_charge_declaration.json --output NEW/charged-sidechain
$PY scripts/audit_pcff_converter_charge_echo.py --source "$FRC" --output NEW/charge-echo
$PY scripts/audit_pcff_c_matcher.py --source "$FRC" --msi2lmp-source ../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/src --output NEW/c-match --compare-native
$PY scripts/audit_pcff_completion.py --source "$FRC" --declaration docs/evidence/phase_4j11_declaration.json --retained-matrix ../island-validation/phase4j10-declared/parity-matrix --output NEW/full-source.json --require-full-source
# Expected exit 1.
$PY scripts/revalidate_pcff_tranche.py --source "$FRC" --root ../island-validation/phase4j13-declared/peo-vertical --root ../island-validation/phase4j13-declared/thioether-vertical --output NEW/offline.json
$PY examples/pcff_converter_compatibility.py --source "$FRC"
$PY -m pytest -q tests/test_pcff*.py tests/test_prepared_workflow.py tests/test_prepared_bundles.py tests/test_psmiles*.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

Artifacts: `../island-validation/phase4j13-declared/`; `peo-vertical/` and
`thioether-vertical/` contain raw converter files, explicit AA overrides, LAMMPS
commands/outputs, bundle relocation, workflow and child receipts. Scientific
inputs/executables stay external. [Artifact hashes](evidence/phase_4j13_artifacts.json).

The bounded executable compatibility tranche is established. Strict source charge
and interaction gaps remain requirements for future work; neither converter success
nor this tranche authorizes full PCFF coverage, production dynamics, periodic
systems, or automated crosslinking.

## Verification executed

* Ordinary suite: **1,683 passed, 10 skipped**, 337.44 s.
* Focused PCFF/PSMILES/bundle/workflow suite: **540 passed**, 228.93 s;
  final converter/guarded-resolver regressions: **17 passed**, 6.49 s.
* Ruff and both environment `pip check` commands passed.
* Both affected examples executed; compiled C matcher controls, independent raw
  selections, ten LAMMPS configurations, two FD scans/minima/relocated workflows,
  and six offline completed-workflow reconstructions passed.
* **3,300 historical path entries / 2,496 unique files** checked unchanged, including
  the pre-existing aliased paths in historical manifests. [Preservation receipt](evidence/phase_4j13_preservation.json).
* Both full-source commands returned **1**, as required. Negative chemistry is
  retained as rejected preparation, not counted as successful executable coverage.

The peroxide control's initial expectation was overly restrictive: historical v4
already types it and assigns charges. The actual retained outcome is missing cross
terms and public model rejection; no rule was changed to fit the expectation.
The separately fixed reporting/test-fixture errors and original logs remain in the
artifact manifest. [Commands, environment, exit codes and log hashes](evidence/phase_4j13_verification.json).
