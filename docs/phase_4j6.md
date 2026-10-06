# Phase 4J6 — authoritative adjudication of missing PCFF interactions

## Scope and repository

Branch: `codex/phase-4j6-pcff-authoritative-resolution`. Base:
`42b9502c467661b20cc5144b8427df31f39513ae` (updated `origin/main`). Its tree
is identical to reviewed J5 `79d04b17b367498ae5568616c991311ad18163ce`.
No prerequisite was merged, main was not changed, and no historical scientific
record was rewritten.

This phase distinguishes a failed source lookup from a genuinely blocking model
term. It adds a recomputed diagnostic contract and executable authority controls.
It **does not establish a new wildcard priority or charge correction**. All six
J5 target models remain incomplete. Full-source completion remains false;
`production_validated=False`, `simulation_readiness="not_established"`.

The external source is unchanged:

- LAMMPS revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- `tools/msi2lmp/frc_files/pcff.frc`.
- SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

## Declaration and reproduction

[The pre-execution declaration](evidence/phase_4j6_declaration.json) fixes the six
checksummed J4 inputs reused by J5, seed 2026, v4 typing, v2 resolution, existing
LJ/Coulomb `(0,0,1)` choices, candidate source rows, asymmetric angle geometry,
reference implementations, numerical tolerances and external timeouts. No QM,
new minimization or propagation is declared. Original failures were reproduced
using `validate_pcff_operational_support.py` before production edits; the full-source
command exits 1. All typing, charge and available model identities match J4/J5.

Detailed new artifacts are external under `../island-validation/phase4j6-declared/`.
J5 receipts remain immutable. The J6 receipts retain exact interaction IDs,
ordered sites/types, candidate rows, equilibrium dependencies and charge searches.
[The gap receipt](evidence/phase_4j6_gaps.json) is not a successful model record.

## 1. Wildcard precedence remains unresolved

The FRC defines conflicting automatic quadratic-angle patterns. The source has
no audited declaration establishing a unique priority between these equally
specific endpoint patterns. We inspected and executed a second pinned converter's
actual matcher instead of inventing a numeric interpretation of `*2`, `*4`, etc.

External LUNAR revision `67dabeda9e6bd3cc8f968aa0c6a88730886138ef`,
`src/all2lmp/ff_functions.py`, SHA256
`aa4487ca98a8a018edf60cfbf42c82a1c97aa1daae8fa1cd75ace2c68859f1f1`:
[`match_angle_wildcards` and `match_3_body`, lines 685–854](https://github.com/CMMRLab/LUNAR/blob/67dabeda9e6bd3cc8f968aa0c6a88730886138ef/src/all2lmp/ff_functions.py#L685).
The module is a user-supplied, checksum-verified validation tool, GPL-3 licensed;
its code is not bundled and it is not a runtime dependency. The validation
adapter independently reads the raw FRC and supplies its tables to the actual
upstream function. It does not claim to run the complete LUNAR converter.

| Angle environment | Forward selected source line | Reversed selected line | K, forward/reversed (kcal/mol/rad²) |
|---|---:|---:|---:|
| hc–c=–s' | 1879 | 1877 | 40 / 37.5 |
| c3–c=2–s' | 1879 | 1890 | 40 / 36.2 |
| c3h–c_1–n3n | 1797 | 1800 | 53.5 / 40 |
| c4h–c_1–n4n | 1798 | 1801 | 53.5 / 40 |

Forward here means the lexicographically canonical supplied-type sequence used
in the receipt, with its center unchanged. Both orientations describe the same
physical angle. Reversing source table iteration did not change these choices,
but reversing the angle did. This fails the required physical reversal invariance
and therefore does not justify adopting that selection behavior. Renaming only
wildcard suffixes to `*99` in a clearly synthetic table control preserved the
selected coefficients. The exact-row control selected line 1807; an unknown
center was rejected. No modified source file was used for molecular assignment.

Pinned msi2lmp
[`find_match` / `match_types`](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L1055)
tries exact matches then wildcard matches, recognizes a leading `*`, and can
select by file order. It supplies no authoritative numeric suffix ordering.
Its [README limitations](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README#L140)
explicitly exclude automatic equivalence supplementation and bond increments.
J5 converter failures remain intact; no `-ignore` or repaired historical output
is used as acceptance.

## 2. Missing cross terms versus the existing BB13 policy

**Correction to J5 explanatory prose, preserved separately:** J5 calls the
bond–bond 1–3 exception a “central-cp” rule. The actual H4/J1 model has always used
`msi2lmp_non_cp_zero_v1`: a policy-derived zero when **none of the four supplied
proper types is `cp`**, with assigned equilibrium dependencies. This phase does
not introduce a new zero or alter that condition.

The pinned converter initializes BB13 K to zero and installs terminal equilibrium
bond lengths at
[`GetParameters.c` lines 363–392](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L363).
It searches BB13 rows only if **any** of the four types is `cp` at
[lines 506–525](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L506).
The identical native condition is extracted into one helper shared by model
construction and diagnostics. Historical signatures and equations do not change.

The independent J5 raw-source oracle was rerun on all six cases. Missing
bond–bond, bond–angle, middle-bond–torsion and other required rows stay missing
under the declared searches. No automatic cross-term family exists to supply
an invented fallback. Some angle–angle and torsional records are present but
have unresolved equilibrium angles. They remain dependency blockers; tetrahedral
angle–angle couplings are not removed. Nonzero Wilson equilibrium and empty
torsion–torsion semantics are not newly resolved here.

| Case | Native charge gate | Structural blockers after policy adjudication | Existing conditional BB13 policy zeros | Complete model |
|---|---|---:|---:|---|
| Thioformaldehyde | pass | 9: 2 ambiguous, 7 missing | 0 | no |
| Thioacetone | pass | 33: 2 ambiguous, 31 missing | 12 | no |
| Pyridinium | fail: 3 absent increments | 0 | 0 | no |
| Guanidinium | fail: −0.0001 e residual | 0 | 12 | no |
| Three-membered lactam | pass | 14: 1 ambiguous, 5 missing, 8 dependencies | 13 | no |
| Four-membered lactam | pass | 20: 1 ambiguous, 7 missing, 12 dependencies | 25 | no |

The 62 conditional zeros were already part of the operational profile. J5's
raw missing-row counts remain correct as **source lookup counts**, but must not
all be called unresolved model terms. “Missing under declared searches” is a
precise claim about the pinned source and supported lookup semantics, not proof
that every commercial PCFF release or unpublished convention lacks that term.
There is no evidence justifying omission of the remaining required couplings.

## 3. Native charge limitations

An independent Decimal audit of all **564** increment records found **no row
containing `nh+`**. The ordinary bond and automatic increment equivalences retain
that label. Pyridinium has two `cp–nh+` bonds and one `nh+–hn` bond without coverage.
Nearby neutral rows 728 (`cp–nh`) and 805 (`h*–nh`) do not authorize substitution.
The source's atom/equivalence declarations alone cannot create the missing
charged increment. Its partially accumulated charges are never published as a
complete assignment.

Guanidinium's only `c+` increment row is line **528**, `c+–nr`, endpoint values
**0.2653, 0.0680**. Its three bonds contribute exactly
`3 × (0.2653 + 0.0680) = 0.9999 e` in decimal arithmetic. Six hydrogen bonds use
line **808**, `h*–nr`, **0.4068, −0.4068**, and add zero molecular charge.
Formal charge is +1; residual is −0.0001 e. This is already present in the printed
source numbers, not a floating-point summation defect. Their four decimal places
are compatible with rounding, but do not prove unrounded values or establish a
scientific correction. The 1e-12 e gate, raw values and failure remain unchanged.

Exact rows, equivalence columns and original endpoint orientations are in
[the charge authority receipt](evidence/phase_4j6_source_audit.json). No new base
charge, alias substitution, residual redistribution or relaxed tolerance is used.
Required next evidence: a traceable charged-pyridine increment convention and
an authoritative guanidinium charge-conservation convention or higher precision
source. Neither is established by the inspected sources.

## Diagnostic API and persistence

```python
from island.forcefields.pcff import (
    assess_pcff_resolution, save_pcff_resolution_assessment,
    load_pcff_resolution_assessment,
)

from island.forcefields.pcff.fallbacks import DOMAIN_POLICY

assessment = assess_pcff_resolution(
    system, validated_typing,
    resolution_policy=DOMAIN_POLICY,
    special_pairs=declared_special_pair_policy,
)
print(assessment.payload["structural_blocker_count"])
save_pcff_resolution_assessment(assessment, "assessment.json")
restored = load_pcff_resolution_assessment(
    "assessment.json", verified_source, system=system,
)
assert restored.identity == assessment.identity
```

Use the exact current resolution policy constant from
`island.forcefields.pcff.fallbacks.DOMAIN_POLICY` in executable code (see commands below).
The owned `island_pcff_resolution_assessment_v1` record binds original typing,
source, resolution policy, special-pair choices, original native identities,
charge diagnostics, source outcomes, conditional policies and dependencies.
Integrity validation recomputes the native inspection and classification. A
rechecksummed contradictory classification, charge gate or policy is rejected
before exclusive publication. Data-only loading requires core/NumPy, not RDKit,
OpenMM, SciPy, ParmEd, Foyer, LAMMPS or LUNAR.

A diagnostic is **not** a native model and cannot be adopted as a prepared force
field. `native_model_complete` is taken from the native validator, never inferred
from a zero structural-blocker count. Rejected source observations remain owned
and inspectable. Saving never repairs a native record.

## Numerical and durability evidence

Seven declared, individually selected source candidate-angle rows were exercised
at a fixed asymmetric non-equilibrium geometry. LAMMPS inputs are authored from
independently parsed raw K and theta; ISLAND's existing compiled quadratic angle
uses radians and K×4.184. The form is `K (theta−theta0)^2`, without a one-half.
The independent handwritten `acos` energy is centrally differenced at
1e-4, 1e-5 and 1e-6 Å. Error decreases initially, then roundoff dominates the
smallest displacement; no monotonic convergence claim is made beyond that.

All seven comparisons pass predeclared `atol=1e-5` kJ/mol (energy) and
kJ/(mol·angstrom) (forces), `rtol=2e-10`. Maximum discrepancies:

- Energy: **1.4210854715202004e-14 kJ/mol**.
- Analytical force: **8.526512829121202e-14 kJ/(mol·angstrom)**.
- Final finite-difference force: **3.633854817053361e-08 kJ/(mol·angstrom)**.

The installed LAMMPS executable reports `30 Sep 2026`,
`patch_30Sep2026-65-ge891a3e10`; SHA256
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`.
Raw inputs/logs/forces, exact commands, times and hashes are retained.
These are **isolated candidate-term checks**, not independently verified complete
molecules or proof of wildcard selection. No incomplete model was evaluated.

All six new diagnostic records were relocated and reconstructed in a separate
Python process with scientific imports blocked, explicit checked source resolution,
byte-identical systems and exact identity checks. This is diagnostic persistence,
not a prepared-bundle acceptance claim. Since all six native model gates fail,
new target bundle, minimization and workflow gates remain unmet; they were not
replaced with synthetic success or rerun with different inputs.

## Reproduction commands

Run from the repository; tools/sources remain external and explicitly supplied:

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
ROOT=../island-validation/phase4j6-declared
DECL=docs/evidence/phase_4j6_declaration.json
$PY scripts/validate_pcff_operational_support.py --source "$FRC" \
  --declaration "$DECL" --output "$ROOT/reproduction" --require-full-source
# Expected exit 1; existing incomplete models remain incomplete.
$PY scripts/validate_pcff_authority.py --source "$FRC" \
  --declaration "$DECL" --reproduction "$ROOT/reproduction" \
  --matcher ../island-validation/phase4h1-sources/LUNAR/ff_functions.py \
  --lammps ../island-validation/phase4j2-tools/build/lmp --output "$ROOT/authority"
$PY scripts/assess_pcff_resolution.py --source "$FRC" \
  --declaration "$DECL" --output "$ROOT/assessment" --require-full-source
# Expected exit 1 even when diagnostic_gate is true.
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

Directories are exclusive; use new output locations to repeat. The diagnostic
CLI's default exit 0 means successful inspection/reconstruction only; its explicit
operational and full-source gates remain false. Term-check success does not
promote an ambiguous molecular model.

## Verification and remaining gates

[The verification receipt](evidence/phase_4j6_verification.json) records tests,
actual external checks, retained artifact hashes, historical offline inspections,
and harness failures separately. Original failures are preserved. Synthetic
software tests cover policy scope (including `cp` in every role), missing
couplings/dependencies, explicit zero versus missing, ownership, strict persistence,
contradictory records and refusal to adopt a diagnostic as a prepared model.

Coverage does not increase: **86/133** labels have graph predicates, **47/133**
remain unimplemented/ambiguous; that is not a global scientific validation count.
For these six targets: typing **6/6**, native charges **4/6**, complete executable
models **0/6**, independent whole-molecule numerical acceptance **0/6** before
and after J6. Isolated-term acceptance is **7/7**. The corrected distinction
between source gaps and existing policy coverage makes the remaining work precise;
it does not complete any of the six domains.

Next actions are specific: establish reversal-invariant wildcard authority;
obtain required thio/lactam cross-term rows and their equilibrium dependencies;
establish a source-backed charged-pyridine increment; establish a separate,
scientifically justified guanidinium conservation convention. Missing couplings
must not become zeros. No merge, production readiness or full-source completion
is claimed.

Final executed checks: **1,604 passed, 10 skipped** in the ordinary suite
(1073.25 s); **73 passed** in the focused PCFF suite (59.71 s), including the
three new diagnostic tests. Ruff and both environment pip checks pass. Six
historical bundles and two completed workflows passed read-only inspection,
frame reading and no-op resume with scientific imports blocked; their identities
and hashes match J5 receipts exactly. **909 historical files** plus the declared
external source/code pins remain unchanged. All eight rerun chemical negative
controls retain their J5 identities. The initial test-only import error and a
misspelled focused-test path are retained in logs; corrected checks pass.

No new complete-domain whole-system calculation, prepared bundle, minimization
or dynamics was run: native charge/model failures prevent those operations.
Diagnostic relocation and existing eligible bundle checks are reported separately.
