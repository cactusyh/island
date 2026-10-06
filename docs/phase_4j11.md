# Phase 4J11 — source obstructions and validation-cache corrections

## Acceptance status

**J11 is not complete. The required full-source gate exits 1.** No source label
or incomplete model is promoted by this change. In particular, the requested
native guanidinium charge gate cannot pass simultaneously with the unchanged
printed increments, the +1 formal charge, and the unchanged 1e-12 e tolerance.
The inspected source and converter establish no correction convention.

Branch: `codex/phase-4j11-pcff-full-source-resolution`.
Base: `980c5b5ad7ea16e8f53ab32bbb6e291bfd8632c8`, as requested. After fetching,
origin/main was `63a3561` (J10 squash); its tree equals that reviewed base.
Main was not modified or merged. No additional graph profile is introduced.

The implementation fixes reproduced cache invalidation defects and adds
reproducible independent source/interaction and executable C-matcher audits.
These are not substitutes for the unmet chemical and numerical requirements.

## Immutable source and declared checks

Source: LAMMPS revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

`phase_4j11_declaration.json` was published before implementation. The separate
performance declaration binds exact historical input hashes and budgets before
measurement. Source files, C extraction/build products, native records, and
workflow/checkpoint artifacts stay outside Git in
`../island-validation/phase4j11-declared`. Prior J1–J10 evidence is not rewritten.

## Native charge obstructions

An independent Decimal reader scans all 564 bond-increment rows, with separate
ordinary `bond` and automatic `bond_increment` equivalence columns. It preserves
all endpoint orientations and candidate rows. Existing retained typing/charge
records are also loaded through their native validators. This invokes no typing
job, parameter assignment, embedding, optimization or force evaluation.

* Pyridinium: `nh+–cp` and `nh+–hn` remain absent through direct, ordinary and
  automatic searches. No increment row contains `nh+`. Its ordinary row 276 and
  automatic row 399 preserve the charged label. Neutral `nh` is not substituted.
* Amines: `na–hn2` remains absent; ordinary rows 272/248 and automatic rows
  395/373 retain `na` and `hn2`. No `hn` substitution is authorized.
* Guanidinium: row 528 (`c+–nr`) contributes 0.2653 + 0.0680 per bond; row 808
  (`h*–nr`) contributes 0.4068 − 0.4068. Three C–N and six N–H bonds give
  **0.9999 e exactly in decimal arithmetic**, residual **−0.0001 e** from +1.
  Both component and whole-system conservation fail the 1e-12 e contract.
  Decimal precision is evidence of printed numbers, not recovered unrounded data.

The [pinned msi2lmp README, lines 140–149](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README#L140)
explicitly excludes bond increments and automatic-equivalence supplementation.
It consumes supplied charges and cannot establish an intended charge residual
or supply a native correction. Consequently no numerical model, minimization or
workflow is published for these failed native inputs.

Required new authority is specific: charged-pyridine increments, an amino-H
increment for the declared labels, and a documented guanidinium charge convention
or more precise source data. None is established by the pinned source. Changing
the required source pin or charge tolerance would be a separate scientific policy.

## Source lookup and executable converter controls

`audit_pcff_c_matcher.py` checksum-verifies and extracts the actual `find_match`
and `match_types` functions (GetParameters.c lines 1055–1169) to an external C
library. It compiles with `cc -shared -fPIC -Wall -Werror`. It uses the pinned
`Forcefield.h` layout and independently parsed raw source patterns.

Exact hashes:

* GetParameters.c: `add196c115ed9066aa5eba6da5bd8036bd29a3619395385140276e108c769591`
* Forcefield.h: `4d6c8fa7772d235a3e0fa5ed7398da11b7ced28794170214d754f79d23c3c139`

The code tries exact matches first, then wildcard rows in file order, checking
forward and reverse per row. For `h_–c=_–s'_`, the original table selects line
1875 (K=60); a separately labelled reversed-table control selects line 1879
(K=40). Reversing query order preserves each selection in this control. The
independent physical candidate check still finds conflicting equally specific
rows 1877 (K=37.5) and 1879 (K=40). The exact `h_–c'_–h_` control selects line
1807 in both orientations/table orders. The unknown-center control finds no row.

This validates leaf-function behavior. It does **not** claim the full converter
supports automatic rows or that first-file-order priority resolves their physical
ambiguity. No runtime precedence, numerical equation, parameter value or
historical resolution policy changes. No `-ignore` is used.

## Interaction-level and type ledger

The fresh full-source row audit still accounts for 5,319 versioned records,
133 atom labels, 134 ordinary and 108 automatic equivalences. It does not equate
presence with assignment or numerical verification.

The independent interaction audit rechecks 162 missing requests in retained
J10 PEO, thioether and halogenated final graphs:

| Case | Missing source requests | Families |
|---|---:|---|
| PEO chain | 42 | BB13 |
| Thioether chain | 27 | BB13 |
| Halogenated chain | 93 | BB13, BB, BA, AA, EBT, MBT, AT, AAT |

Each receipt retains the exact sites, supplied types, native candidate/dependency
record and independent raw-source search. These are **raw source gaps**. Existing
historical non-cp BB13 compatibility zeros remain readable under their original
policies; they are not source rows and do not satisfy J11's strict requirement.
No unrelated coupling is removed or filled with zero.

The only nonzero Wilson equilibrium row is `wilson_out_of_plane:cff91:3231`,
`az ob hb sz`, K=116.0100 kcal/mol/rad², equilibrium 3.8934 degrees. Its signed
orientation/peripheral permutation authority remains unresolved for native
assembly. The empty torsion–torsion section supplies no coefficient and does not
establish universal physical inapplicability.

All 47 unresolved labels are retained with atom row/line, descriptions,
ordinary/automatic columns, direct increment rows, nonbond searches and required
next evidence. Their corrected grouped domains are 16 metal, 2 halogen, 1 ionic,
6 charged-environment, 10 zeolite/surface, and 12 organic alias/specialized
environments:

`Ag Al Au Br Cl Cr Cu Fe K Li Mo Na Ni Pb Pd Pt Sn W az c_a ca+ cg ci cr h* h+
hb hi hoa hos n+ n1 n2 nb nho ni nz oah oas ob oe osh oss p= s- sf sz`.

Metal phase/coordination cannot be inferred from an element alone. The source
aliases and descriptions do not resolve precedence, charged resonance,
coordination/base-charge or surface/isotope conventions. These remain genuine
implementation/interpretation gaps, not declarations of chemical inapplicability.
Typing stays **86/133 predicate labels**, with **47/133 unresolved**. No additional
label or domain gains independent numerical verification in J11.

## Cache correction and measured integrity performance

Five regressions failed before the fix: replacing a rule's `__code__` in place,
mutating positional defaults, keyword defaults or closure data, and changing a
parser function's code all reused earlier successful validation. Function object
identity alone did not detect these changes.

Cache keys now snapshot function code, defaults and closure data (including
wrapped functions and cyclic containers). Class static/class methods are also
included in the trusted context. The parser cache uses the same callable-state
helper. No arbitrary `repr`, user property or serialization method is executed
by this helper. This is process-local invalidation, not authentication against
arbitrary code running in the process. Existing content keys, nested-record
checks, immutable-source guards, cold semantic validation and live graph binding
remain in effect. Persisted parameter/model/evaluation identities do not change.

Declared budgets: short cold load ≤60 s and 20 validations ≤5 s; retained long
record cold validation ≤120 s and 20 validations ≤10 s. Measured:

| Retained record | Sites | Cold load/validation | 20 validations |
|---|---:|---:|---:|
| J10 fused prepared bundle | 34 | 4.234 s | 1.234 s |
| H5 PE DP50 parameter/model records | 302 | 41.859 s | 0.235 s |

Native/source caches were cleared for each cold case; Python modules may already
be imported. Short checks validate the facade and final graph; long checks
validate native model and final graph. These are distinct workloads, not a
speedup comparison. Original hashes and native/facade identities match. The long
case is historical integrity/performance evidence, not new strict-source coverage.

## Reproduction

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
ROOT=../island-validation/phase4j11-declared
$PY scripts/audit_pcff_completion.py --source "$FRC" \
  --declaration docs/evidence/phase_4j11_declaration.json \
  --retained-matrix ../island-validation/phase4j10-declared/parity-matrix \
  --output /new/source-and-interaction-obligations.json --require-full-source
# Exit 1: required native charge/model obligations remain unmet.
$PY scripts/audit_pcff_c_matcher.py --source "$FRC" \
  --msi2lmp-source ../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/src \
  --output /new/c-matcher
$PY scripts/audit_pcff_source_rows.py --source "$FRC" \
  --output /new/source-rows.json --require-full-source
# Exit 1: unchanged full-source acceptance contract.
$PY scripts/benchmark_pcff_validation.py \
  --declaration docs/evidence/phase_4j11_performance_declaration.json \
  --output /new/performance.json
$PY -m pytest -q tests/test_pcff*.py tests/test_prepared_workflow.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

Outputs are exclusive. No new whole-system LAMMPS, FD, minimization or dynamics
matrix is claimed: physical equations are unchanged and no new source-complete
domain has been established. J10's successful independent numerical receipts
remain historical evidence only. The new C checks are assignment-matcher controls.

`production_validated=False`; `simulation_readiness="not_established"`.

## Executed verification and preservation

* Ordinary suite: **1,662 passed, 10 skipped**, 311.94 s, exit 0.
* Focused PCFF/prepared-workflow suite: **476 passed**, 214.23 s, exit 0.
* The five reproduced cache failures pass after correction. Additional controls
  cover nested policy mutations, rehashed changed source bytes, cyclic closure
  ownership, explicit source zero versus absence, and endpoint orientation.
* Ruff and both applicable environment pip checks pass.
* The independent charge audit agrees with retained partial vectors within
  1e-16 e maximum difference; all three native charge cases correctly remain failed.
* Compiled pinned C controls and both declared performance budgets pass.
* A separate Python process reloaded J8/J9/J10 bundles and validated completed
  workflow status, five frames and no-op resume, with external OpenMM, RDKit,
  SciPy, Foyer and ParmEd imports blocked. Native/facade/profile identities match
  the prior receipts. This is offline reconstruction, not new propagation.
* All **2,292** historical files in the pre-execution snapshot are byte-identical
  after execution. The separately declared historical PE DP50 input hashes also
  remain unchanged. Unrelated checkpoint directories were not edited.

The raw-source row command and completion-obligation command both exit **1**
with `--require-full-source`. No new complete domain, independent whole-system
energy/force check, minimized structure or trajectory is claimed.

`phase_4j11_verification.json` records executable outcomes, environment/log hashes,
input preservation and initial harness errors. Initial wrong-directory, C extraction
and import-hook failures are retained separately from corrected software controls.
No chemical fixture, charge, coefficient, seed or numerical tolerance was changed.

Machine-readable deliverables are `phase_4j11_obligations.json` (47 type obligations,
source charge proofs, nonzero Wilson row and 162 exact interaction gaps),
`phase_4j11_c_matcher.json`, `phase_4j11_performance.json`,
`phase_4j11_historical.json`, and `phase_4j11_preservation.json`. The fresh full
5,319-row ledger is external at `source-rows.json`, bound by its hash in the
verification receipt; J10's independently observed row ledger remains unchanged.

**The requested full-source completion is still false.** This commit fixes
validation correctness and provides independently executable evidence explaining
why the currently required source/charge/model gates remain unmet. It does not
claim that the remaining chemical-domain implementation has been completed.

## Corrective source-obligation taxonomy

The diagnostic taxonomy originally grouped uppercase `Br`/`Cl` with `ca+` as
`ionic`. This correction classifies `Br` and `Cl` as `halogen`, preserving the
pinned atom-row descriptions at lines 49/50 (`bromine ion`, `chlorine ion`).
`halogen` describes the chemical family here; it does not assert neutral atoms,
authorize an ionic charge model, or equate uppercase and lowercase labels.

`ca+` remains `ionic` based on its explicit `calcium ion` description (line 85).
Source descriptions explicitly mentioning charged or protonated environments
place `ci`, `h+`, `hi`, `n+`, `n1`, and `ni` in `charged_environment` (lines 87,
100, 104, 118, 119 and 136). This category does not assign a formal charge to
each atom in such an environment. Neither an element name nor a sign in a type
label establishes that classification: `Na` remains metal, and `s-` remains a
specialized environment with the source description `partial double sulfur`.
Names of functional groups alone do not resolve their existing alias ambiguity.

The corrected counts are **16 metal, 2 halogen, 1 ionic, 6 charged-environment,
10 zeolite/surface, and 12 other specialized environments**, totaling the same
47 unresolved labels. The receipt changes only domain counts, domain categories
and their required-evidence descriptions. Source rows, equivalences, numerical
values, charge failures, interaction obligations and unresolved/readiness flags
are unchanged. No runtime typing, charge, parameter, model or cache code changes.

Four new regression cases failed with the previous taxonomy and pass after the
correction. The pinned-source audit was regenerated in the separate
`../island-validation/phase4j11-taxonomy-correction` directory;
`--require-full-source` still exits **1**. The original external receipt and
prior verification receipts are preserved. Detailed checks and hashes are in
`phase_4j11_taxonomy_correction.json`.

Correction verification: **19 focused completion-audit/cache tests passed**
(16.88 s); **1,666 ordinary tests passed, 10 skipped** (323.71 s). Ruff and both
primary/isolated pip checks passed. Separate-process read-only reconstruction of
J8/J9/J10 bundles, completed workflow status, five retained frames per workflow
and completed no-op resume passed with scientific imports blocked. All native,
facade and profile identities match the prior receipt. All **2,338 protected
historical files** are byte-identical; the authorized taxonomy receipt is the
only updated historical evidence file, with its original bytes retained in the
separate correction directory. Readiness and full-source acceptance stay false.
