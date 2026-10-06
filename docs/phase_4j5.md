# Phase 4J5 — bounded PCFF operational blockers

## Result and scope

**Operational completion is unmet: 0 of the six requested cases is a complete
native model.** This change supplies executable, independently checked diagnosis
of the remaining source gaps; it does not invent the missing physics. No typing,
charge, resolution-policy, equation, coefficient, tolerance, or historical schema
was changed. `production_validated=False`, `simulation_readiness="not_established"`.

Branch: `codex/phase-4j5-pcff-bounded-completion`. Base:
`8765f8dd214ea488cb8bdfd753a5031bcf9df0bb` (updated origin/main); its complete tree
was verified identical to reviewed J4 `2dbdd2f694169ee6eb68df2706807eeaa5030ec9`.
No merge or modification of main was performed.

The six failures were reproduced **before production changes**, using checksummed
J4 input systems. Typing, charge, assignment and model identities remained equal
to their original J4 records. All later artifacts are in the separate
`../island-validation/phase4j5-declared/` directory. The declaration is
[phase_4j5_declaration.json](evidence/phase_4j5_declaration.json).

| Requested case | Typing | Native charge | Structural lookup assigned / ambiguous / missing / inapplicable | Complete model |
|---|---|---|---|---|
| Thioformaldehyde | complete | complete, residual 0 | 11 / 2 / 7 / 0 | no |
| Thioacetone | complete | complete, residual 0 | 119 / 2 / 43 / 8 | no |
| Pyridinium | complete | three missing increments | 246 / 0 / 0 / 0 | no |
| Guanidinium | complete | residual −0.0001 e | 131 / 0 / 12 / 0 | no |
| Three-membered alpha lactam | complete | complete, residual 0 | 129 / 1 / 18 / 4 | no |
| Four-membered beta lactam | complete | complete, residual 0 | 223 / 1 / 32 / 8 | no |

These are **source-query counts**, not model completeness counts. A selected
cross-term coefficient still requires unambiguous equilibrium dependencies.
Policy-derived BB13 treatment is left to the unchanged native model validator.
In particular, the special central-`cp` rule is not applied to other chemistry.
The report does not promote a coefficient lookup to a successful preparation.

The bounded matrix remains 6/6 typed, 4/6 charge complete, 0/6 complete executable
models, 0/6 independently verified whole-system models, and 0/6 new durable
bundles, both before and after this work. The full source still has 86/133 labels
with bounded predicates and 47 unresolved labels. No label is globally certified
by these fixtures. The earlier ten-case J4 matrix remains 10 typed, 8 charge
complete, 4 model complete; its four verified halide cases are unchanged.

## Source and audited paths

Pinned source: LAMMPS `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
The local reference checkout has that revision and identical FRC bytes. A separate
raw-line count confirms 133 distinct atom labels, 134 ordinary equivalences,
108 automatic equivalences and 564 increment rows.

The [pinned FRC](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc)
lines 18–36 declare ordinary and automatic families. Automatic equivalence is
position-specific for the declared lower-order base terms. It does not declare
automatic bond–angle, angle–angle, bond–torsion or other Class II coupling families.
No such families were inferred here.

Selection remains the versioned J2/J3 policy: direct then family ordinary
matching, and only the declared automatic positional routes; highest version per
source pattern, most specific matching pattern, and rejection of conflicting
equally specific oriented candidates. Numbered wildcards are not assigned an
invented numeric priority. The independent oracle agrees with all 977 applicable
structural requests, including failures. Maximum floating-point conversion
roundoff is 4.44e-16; row identities and statuses agree exactly. The maximum
per-site native versus independent Decimal-sum charge difference is 1.11e-16 e.

The independent oracle in `scripts/pcff_j5_reference.py` does not call the
production parser, selector, conversion or inventory functions. It reads raw
rows with Decimal, independently enumerates bonds/angles/propers/center-preserving
angle pairs, and records ordinary/automatic columns, transformations, candidate
rows, legal permutations and endpoint contributions. Charge summation uses the
original decimal values. This is independent implementation evidence for the
**declared** interpretation, not proof that an undocumented commercial wildcard
priority or missing coupling can be recovered.

The pinned [msi2lmp README lines 145–147](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README#L145)
excludes automatic-equivalence supplementation and bond increments. Its externally
supplied labels therefore cannot establish automatic typing or charge fidelity.
[GetParameters.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c)
ordinary bond/angle paths (121–128, 177–185) and matching logic (1055–1175) remain
bounded controls, including its file-order and angle–angle central-role caveats.
The known prior converter failure/correction receipts are unchanged.

## Exact remaining blockers

Full interaction/site IDs, supplied and transformed labels, searched namespaces,
legal permutations, candidate row IDs and equilibrium dependencies are retained
in [phase_4j5_gaps.json](evidence/phase_4j5_gaps.json). The external `inspection.json`
records include source rows; full force-field libraries are not bundled in Git.

* **Thioformaldehyde:** `s'–c=–hc` has conflicting automatic angle candidates
  including rows 1877 (`*2 c=_ h_`, 120°, 37.5 kcal/mol/rad²) and 1879
  (`*4 c=_ s'_`, 120°, 40). Missing bond–bond, bond–angle and center-preserving
  angle–angle rows remain. The row 2218 `s' c= s'` and row 3444 `s' c= s'`
  describe different endpoint chemistry and cannot fill the H/S case.
* **Thioacetone:** `c3–c=2–s'` maps to conflicting automatic angles, including
  rows 1879 and 1890 (`*7 c=_ c_`, 120°, 36.2). Required S-containing bond–bond,
  bond–angle, torsion couplings and angle–angle entries are missing. Its selected
  C=S bond at row 1700 does not imply coverage of those interactions.
* **Pyridinium:** two `cp–nh+` bonds and one `nh+–hn` bond lack increments after
  direct, ordinary bond-equivalence and automatic increment-column searches.
  `nh+` is not replaced with `nh`; `hn`'s equivalent `h*` does not supply the
  missing row. The known-contribution sum is 0, with incomplete coverage and
  formal charge +1. Reporting −1 as a *complete molecular residual* would be
  misleading; the diagnostic explicitly marks the missing coverage.
* **Guanidinium:** source increment row **528**, namespace `cff91_auto`, gives
  `c+ nr: (0.2653, 0.0680)`. Three bonds contribute exactly **0.9999 e**; the N–H
  increments conserve their endpoint sums. The complete component has formal
  charge +1, so the residual is **−0.0001 e**, outside the unchanged 1e-12 e
  criterion. Raw decimal arithmetic confirms this is present in the source,
  not introduced by floating-point summation. No alternative source-backed
  charge path resolves it. Twelve BB13 requests also lack source rows.
* **Lactams:** the carbonyl-centered `c3h–c_1–n3n` / `c4h–c_1–n4n` angles have
  competing automatic rows 1797/1800 and 1798/1801. The source's specific
  `c_ c'_ n_` row does not match a `c3m_` or `c4m_` endpoint. Missing adjacent
  couplings, BB13 and middle-bond–torsion terms remain; assigned angle–torsion
  rows are blocked where their equilibrium angle is ambiguous. Neither ring
  atom type nor dependency is silently replaced by its acyclic counterpart.

The next needed evidence is authoritative wildcard precedence **and** actual
missing coupling data/applicability for these sulfur and ring motifs; pyridinium
needs charge rows, while guanidinium requires an explicitly justified source or
charge-conservation policy correction. A different policy would need a separate
version and evaluation. None is implemented implicitly. Nonzero Wilson equilibrium
and empty torsion–torsion scope retain their earlier unresolved limitations.

## Implementation and API

```python
from island.forcefields.pcff import inspect_pcff_operational_support, special_pair_policy
report = inspect_pcff_operational_support(
    system, validated_typing,
    resolution_policy="island_pcff_positional_fallbacks_v2",
    special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
)
```

The owned `island_pcff_operational_inspection_v1` diagnostic dictionary exposes all
permitted structural searches even after a native charge failure. It is **not**
a parameter/model/prepared record and cannot be adopted as one. Native assignment
and model construction run only when the original charge validator succeeds.
The actual native validators decide model completion. Complete compatible typing,
source identity, graph and explicit policies are validated at entry.

Assignment and inspection now share the original authoritative interaction
request enumeration. Its ordering and all historical contents remain unchanged.
`fallbacks.lookup/resolve(trace=True)` adds attempted paths without changing the
default historical result payload. No equations, source rows, charge treatment or
model/evaluation identities were revised. Inspection requires core/NumPy only;
it does not import RDKit, OpenMM, SciPy, ParmEd, Foyer or launch external programs.
The example accepts an already retained system, avoiding coordinate regeneration.

## Executed numerical and converter controls

Two separately declared **isolated source-row controls**, at asymmetric stretched
coordinates, used source quartic bond row 1700 (`c=2 s'`) and row 1706 (`c_1 n`).
The reference directly supplied independently selected raw source values to
LAMMPS `bond_style class2`. It did not export the native assignment list.

| Row | Energy error kJ/mol | Maximum force error kJ/(mol Å) | FD error at 1e-4 / 1e-5 / 1e-6 Å |
|---|---:|---:|---|
| 1700 | 0 | 2.84e-14 | 3.06e-6 / 2.44e-8 / 6.81e-8 |
| 1706 | 7.11e-15 | 5.68e-14 | 2.23e-6 / 2.46e-8 / 2.06e-8 |

Independent analytical radial forces also agree. The FD decrease reaches the
roundoff floor; monotonic improvement at the smallest step is not assumed.
Declared energy/force atol is 1e-5 in canonical units and rtol=2e-10. Coordinates,
commands, coefficients, all outputs and hashes are under `controls-2/` and in the
separate receipt. These controls do **not** constitute whole-molecule acceptance.

LAMMPS executable SHA256:
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`;
msi2lmp executable SHA256:
`a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`.
LAMMPS reports `30 Sep 2026 - Development - patch_30Sep2026-65-ge891a3e10`;
the retained CMake cache points to the verified checkout and enables CLASS2.
The retained pinned build provenance is reused; no new executable/source was
substituted. Raw LAMMPS logs retain their reported version/build information.

Converter executions without `-ignore` retain exit codes 12 (thioformaldehyde),
13 (thioacetone), 13 (alpha lactam), 13 (beta lactam). The two charge-failed cases
were not given fabricated charges just to run the converter. No complete model
means no admissible new whole-system energy/force comparison, bundle publication,
relocation, minimization or trajectory. Those gates remain **unmet**.

## Reproduction

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
LMP=../island-validation/phase4j2-tools/build/lmp
MSI=../island-validation/phase4h5-tools/msi2lmp
$PY scripts/validate_pcff_operational_support.py --source "$FRC" \
  --declaration docs/evidence/phase_4j5_declaration.json --output new-audit \
  --require-full-source
# Expected exit 1: native models remain incomplete, and full-source is false.
$PY scripts/validate_pcff_j5_controls.py --source "$FRC" --lammps "$LMP" \
  --converter "$MSI" --declaration docs/evidence/phase_4j5_declaration.json \
  --audit new-audit --output new-controls
$PY examples/pcff_operational_inspection.py --source "$FRC" \
  --system ../island-validation/phase4j4-declared/matrix/pyridinium/system.json
$PY -m pytest -q tests/test_pcff_operational.py
$PY -m pytest -q tests/test_pcff*.py tests/test_prepared*.py \
  tests/test_forcefield_preparation.py tests/test_workflow_consistency.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

CLI outputs use established typed JSON encoding to preserve integer mapping keys;
read with `island.workflows.storage.decode`. Output directories are exclusive.
The inspection example's success means inspection finished. The acceptance CLI
separately requires every declared case, verified native charges/model, independent
whole-system verification and a relocated bundle; partial diagnostics cannot pass
that gate. `--require-full-source` remains nonzero.

## Verification and preservation

See [phase_4j5.json](evidence/phase_4j5.json) for independent row/charge checks and
[phase_4j5_negatives.json](evidence/phase_4j5_negatives.json) for eight unchanged
near-miss controls (O2, N2, isotope substitution, sulfoxide, pyridine, neutral
guanidine and non-lactam small rings). Neutral related chemistry may retain its
own valid labels; it must not acquire the charged/lactam/sulfocarbonyl motif labels.

Initial harness failures are retained: `audit-1` compared converted float lists
with exact equality and exposed only rounding at the last bit; `audit-2` uses the
predeclared numerical tolerance with exact row/status checks. `controls-1` finished
its term executions then failed on a diagnostic-entry field name; `controls-2`
corrects that harness key with identical inputs, source, coefficients and budgets.
The initial synthetic policy test incorrectly treated a valid different caller
policy as malformed; its corrected assertion uses an actually contradictory field.
No scientific failure was repaired by changing inputs or tolerances.

Test counts, read-only historical identities/hashes and artifact inventory are
recorded separately in the final verification/preservation receipts. Historical
J3/J3.1/J4 scientific records and the isotope/radical conversion fix are preserved.
Full-source chemical coverage remains the project objective; this phase has not
completed the bounded operational gate or the remaining 47-label domain work.

Read-only child-process verification passed for six historical J3/J4 bundles and
two completed J3 workflows (including status, frame read and completed no-op resume),
with OpenMM, SciPy, RDKit, ParmEd and Foyer imports blocked. The 790-file historical
hash inventory is unchanged. This is compatibility evidence for existing eligible
records, not a newly completed J5-domain bundle experiment.

Final verification: **1,601 passed, 10 skipped** in the ordinary suite (1048.75 s);
**486 passed** in the PCFF/prepared/workflow suite (892.89 s); **6 passed** in the
new diagnostic tests (10.63 s). Ruff and both environment pip checks pass. The
inspection example ran on the retained pyridinium input. The ordinary suite's
skips are not claimed as executed scientific checks. Real pinned-source inspection,
negative controls, isolated LAMMPS/analytical/FD checks and converter failures are
reported separately from synthetic software tests in
[phase_4j5_verification.json](evidence/phase_4j5_verification.json).

The receipt includes exact environment versions and all new artifact hashes,
including initial harness failures. Read-only compatibility results are in
[historical_offline](evidence/phase_4j5_historical_offline.json) and
[j4_offline](evidence/phase_4j5_j4_offline.json); the unchanged 790-file historical
inventory is [preservation](evidence/phase_4j5_preservation.json).
