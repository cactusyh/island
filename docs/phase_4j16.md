# Phase 4J16 — local PCFF distribution and converter boundaries

## Status and repository

Branch: `codex/phase-4j16-pcff-local-distribution`.
Base: `9f659242358ba576cb4752422af8ebc1f276d097`, J15 squash merge PR48.
After fetching origin, the **complete tree** compared equal to reviewed J15
`9533f9940e690c12f5a7eed354f1158a55a27f83`. No unmerged prerequisite remains.
Main and historical scientific records were not changed.

**Partial chemical coverage remains.** This phase implements a hash-bound,
relocatable distribution inspection API and an artifact-backed converter-boundary
acceptance harness. The actual local FRC contains no new numerical coverage.
The auxiliary files do not establish a complete, executable typing language or
supply missing couplings. There is no justified v7 typing change and no newly
executable polymer family. V1–v6 typing and v3/v4 resolver behavior remain intact.

`production_validated=False`; `simulation_readiness="not_established"`.
Audit completion and converter-boundary agreement are separate from executable
polymer acceptance and full-source completion.

## The three actual local files

All files were accessible, read directly, and preserved. Their resolved directory is:

```
/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/tools/msi2lmp/frc_files
```

| File | Bytes | Encoding | Actual content |
|---|---:|---|---|
| `pcff.frc` | 337914 | ASCII | PCFF parameter library, ordinary/automatic equivalences and increments |
| `pcff.rlb` | 13 | ASCII | Exactly `VERSION\nelib\n`; no rules, templates or numerical coefficients |
| `pcff_templates.dat` | 35450 | ASCII | Text pattern/test blocks and a precedence tree; header names `mcff_templates`, March 1995 |

SHA256 values:

```
pcff.frc
 e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c
pcff.rlb
 82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36
pcff_templates.dat
 dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093
```

The local FRC is **byte-for-byte identical** to the historical numerical source
at LAMMPS `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`. Therefore there are no
numerical or equivalence differences to adopt. Its declarations run from v1.0
(July 1991) through v4.0 (October 8, 2013). The semantic inventory still contains
133 atom labels, 134 ordinary equivalences, 108 automatic equivalences and 564
bond increments. Every section/namespace and numerical record ID is inventoried
in the external distribution audit. Empty torsion–torsion sections remain empty.

The checkout is LAMMPS `36f7c81ced8dd456430b13647e48f0c7d7bd8faa`, remote
`https://github.com/lammps/lammps.git`. The latest template-file change is
`0aa77408f8084f2b6827dc4e91b200046679147d` (SVN-import history). These repository
facts do not establish proprietary application semantics. The local
`frc_files/README` describes the distributed files as publicly available; complete
libraries, converter code and binaries remain external regardless. No upstream
library or implementation was copied into ISLAND.

### What the template file establishes—and what it does not

The bounded reader preserves **197 blocks, 133 distinct labels, 197 patterns and
253 atom tests**, including duplicate candidates and source line IDs. It retains
`type oh` without a colon at line 1152 and non-colon template syntax rather than
silently losing those records. Labels `?` and `lp` have no matching FRC atom row;
FRC labels `nn` and `nr` have no template block. Counts alone are not coverage.

The file contains pattern operators, nested branches, aromaticity/hybridization,
allowed/disallowed elements, and planar/nonplanar ring tests. Its precedence
block is lines 1874–1926. The inspector parses that block's parentheses **only as
structure**. Repeated nodes, pattern atom numbering, bracket/operator meaning,
geometry predicates and override/conflict behavior lack a reader specification
in this checkout. Some pattern text also has unmatched delimiters. These are
retained unsupported semantics, not repaired or executed as code.

Relevant development evidence:

| Template type | Line | Bounded observation |
|---|---:|---|
| `na`, `n4` | 571, 588 | Separate amine/cation families, not a global nitrogen alias |
| `nb` | 608, 621 | Aromatic-family candidates; no `nn` template proving an alias |
| `h*`, `hn` | 1296, 1378 | Distinct H candidates with contextual tests |
| `o_2`, `c_2` | 1594, 1644 | Acyl/ester-family pattern evidence |
| `n_2` | 1666 | Explicit carbamate-N comment and nonaromatic N/acyl pattern |
| `hn2` | 1677 | H attached to nonaromatic acyl N with O/N/C/H substituent tests |

The `n_2`/`hn2` records corroborate the bounded J15 carbamate motif. Their apparent
broader acyl scope does not justify global amide/urea overrides, and they do not
provide complete formal-charge, isotope or radical constraints. Parameter
availability is still independent. No operational rule/template interpreter was
invented; a documented dialect reader, override semantics and validated chemical
fixtures are needed before expanding runtime typing from these records.
`elib` is retained as a marker, not described as a recovered semantic version or
an unavailable binary rule library. No named missing companion file was found.

## Reader and converter audit

The [reader audit](evidence/phase_4j16_reader_audit.json) records source/document
hashes and the exact checkout-wide search. No consumer of `pcff.rlb` or
`pcff_templates.dat` was found. The only other `end_precedence` match is embedded
in `compass_published.frc`; it is not a reader.

Pinned local implementation evidence:

- [msi2lmp.c](https://github.com/lammps/lammps/blob/36f7c81ced8dd456430b13647e48f0c7d7bd8faa/tools/msi2lmp/src/msi2lmp.c),
  lines 456–478: CAR, MDF, topology, FRC, parameter matching.
- [README](https://github.com/lammps/lammps/blob/36f7c81ced8dd456430b13647e48f0c7d7bd8faa/tools/msi2lmp/README),
  lines 145–147: no automatic-equivalence supplementation and no bond increments.
- `ReadFrcFile.c`/`SearchAndFill.c`: FRC numerical sections, no auxiliary-file path.
- `GetParameters.c`: ordinary matching and existing non-cp BB13 initialization.
  It is byte-identical to the historical pinned implementation; no matching rule
  or zero policy is changed here.
- `WriteDataFile.c:391`: coordinates written with `%15.9f`, charges `%9.6f`.

`MakeLists.c`, `ReadMdfFile.c` and `SearchAndFill.c` differ from historical
`e891a3e...`; the actual local build is not mislabelled as that executable.
The checked local source was copied to a **new external build directory**, then
compiled unmodified with its Makefile (`gcc 13.3.0`, `-O -Wall -W -g`).
It identifies as msi2lmp 3.9.11 / 6 Sep 2024. Binary SHA256:
`6e30fd753a99bad66479779c732f1e1dc764c540131ec6955f84a76a9c320904`.
The [build receipt](evidence/phase_4j16_converter.json) records every C/header and
Makefile hash, command and build-log hash.

## APIs and ownership

New inspection-only APIs in `island.forcefields.pcff`:

```python
inspect_pcff_templates(raw_bytes)
inspect_pcff_distribution(
    frc_path=..., rlb_path=..., templates_path=...,
    expected_hashes={"frc": ..., "rlb": ..., "templates": ...},
)
save_pcff_distribution_audit(path, result, **explicit_current_paths)
load_pcff_distribution_audit(path, **explicit_current_paths)
```

`PCFFDistributionAudit` owns its data, binds all three byte hashes and parsed
records, and reconstructs semantics during integrity validation. Original paths
are historical provenance. Loading uses only supplied current paths. Rechecksummed
contradictory records, source mismatches, malformed/truncated envelopes and
replacement of an existing publication are rejected. Relocation was exercised
in a child process with OpenMM/RDKit/SciPy/Foyer/ParmEd imports blocked.

Auxiliary inspection is development evidence, **not a new runtime dependency**.
The normal final-graph preparation, bundles, evaluator/session and workflow paths
continue using their existing explicit FRC. No CAR/MDF/converter is required at
runtime. No chemistry, force equation, charge convention, model identity, source
pin, readiness flag or historical schema was changed.

## Frozen matrix and actual boundaries

The [declaration](evidence/phase_4j16_declaration.json) predates production edits.
It binds all 20 frozen J15 systems and two historical PEO/thioether controls by
checksum. Coordinates, seeds, explicit ends, inter-repeat bonds and metadata were
reused without regeneration. Independent graph predicates from the retained
reference readers supplied converter labels; a separate raw Decimal FRC reader
supplied charges. Per-site charge comparison remains **1e-12 e**, no normalization.

| Case(s) | Typing | Charges | Strict v3 | Compatibility v4 | Actual converter |
|---|---|---|---|---|---|
| Primary/secondary carbamate | complete | complete | incomplete | incomplete | exit 17, missing EBT `hc c3 o_2 c_2` |
| Tertiary carbamate | complete | complete | incomplete | incomplete | exit 13, angle `c3 n_2 c3` |
| Aniline | complete | complete | complete | complete | exit 0 |
| N-methylaniline | complete | complete | incomplete | incomplete | exit 12, bond `c3 nn` |
| Diphenylamine | complete | complete | incomplete | incomplete | exit 13, angle `cp nn cp` |
| Mixed amine/urethane/ether/arylamine | complete | complete | incomplete | incomplete | exit 14, BB `c2 oc cp` |
| Urethane DP1/DP2/DP3 | complete | complete | incomplete | incomplete | exit 17, EBT `hc c3 c2 o_2` |
| Arylamine polymer DP2 | complete | complete | incomplete | incomplete | exit 12, bond `c3 nn` |
| Acetamide negative | complete | complete | incomplete | incomplete | exit 13, angle `c3 c_1 n` |
| Pyridinium, guanidinium, charged carbamate/aniline, heteroarylamine | complete | incomplete | incomplete | incomplete | not run: no invented charges |
| Urea, radical/isotope carbamate | incomplete | not reached | not reached | not reached | not run: no invented types |
| Historical PEO/thioether | complete | complete | incomplete | complete | exit 0 / 0 |

All 14 charge-complete cases have raw component sum **0.0000 e**. Guanidinium
retains printed sum **0.9999 e**, formal charge +1, residual **−0.0001 e**;
missing increment cases remain missing. No renaming or charge repair was applied.

Of 14 executed converters, **3 succeed and 11 stop normally on missing parameters**;
8 prerequisite-blocked cases were not exported. No `-ignore` invocation occurred.
The three successful data files match their retained independently constructed,
previously numerically verified reference files byte-for-byte. This checks
coefficient/inventory parity for these controls; it is not new force-model or
trajectory validation. The [reference links](evidence/phase_4j16_reference_links.json)
bind the older declarations, raw outputs, corrections and numerical receipts.

Independent raw-source diagnosis separates two important causes:

- Carbamate EBT `h c o_2 c_2`, urethane EBT `h c c o_2`, and mixed BB `c o cp`
  have no applicable source row after declared ordinary searches. No automatic
  cross-term family exists to fill them. The complete native gap ledger includes
  AT/AAT/BB/BA/AA and equilibrium dependencies, not just the converter's first error.
- `c3 nn` has automatic quadratic-bond row **1274**; `c3 n_2 c3` and `cp nn cp`
  have automatic quadratic-angle rows **1965/1978**. Acetamide's reported angle
  reaches row **1806** after candidates 1787/1799/1802/1806. These paths already
  exist in ISLAND but not msi2lmp. Their existence does not remove other required
  missing couplings. They are not newly fixed native lookup bugs.

The additional files establish **no demonstrable native lookup defect or new
source-backed coupling**, so the numerical engine and chemical predicates were
left unchanged. Required additional data are exact missing cross-term rows and
their equilibrium/permutation conventions, plus documented template language
semantics if expanded automatic typing is desired.

### Before/after scope

For the identical 22-case matrix, both before and after this change:

- Typing: **19/22** complete.
- Native charges: **14/22** complete.
- Strict models: **1/22** complete.
- Compatibility models: **3/22** complete.
- Reused independently verified controls: **3/22**; new numerical cases: **0**.
- Newly executable J15 polymer targets: **0/4**.

These fixture counts do not promote any of the 133 source labels globally.
The frozen full-source obligation audit still lists 47 unresolved labels and
unverified interaction domains (distinct from the current v6 predicate count). The historical full-source CLI exits **1**.

## Gates, failures and preservation

`check_pcff_distribution.py` is the shared current semantic gate for fresh runs
and retained reports. It rechecks the exact mandatory case set, frozen declaration,
source/build hashes, graph identity, native inspection results, independent labels
and charge records, converter invocation/output, graph roles, artifact inventory
and registered independent reference files. No caller success boolean authorizes
acceptance. Missing/duplicate cases, stale inputs, unrelated errors and wrong
chemical rejection stages fail. Full-source and executable-polymer gates remain
separate from audit/converter-boundary gates.

Two new harness failures are retained, not overwritten:

1. `boundaries/` recorded an inspection bug: reading complete-only assignments on
   three incomplete typings. The runner now reads diagnostic payloads without
   pretending they are complete assignments; `boundaries-corrected/` reran the
   bounded converter experiment on the same frozen inputs.
2. That corrected run's `checked.json` retains the failed raw-coordinate assertion:
   maximum observed difference `4.97806463e-10 angstrom`, caused by converter
   `%15.9f` serialization. The separate final checker uses **exact equality to
   that documented decimal representation**, and records raw rounding differences.
   It does not loosen a numerical tolerance, change coordinates or overwrite raw
   output. A focused regression rejects even a one-last-decimal change.

Scientific tolerances remain energy `1e-5 kJ/mol`, force
`1e-5 kJ/(mol*angstrom)`, rtol `2e-10`; finite-difference displacements
`1e-4/1e-5/1e-6 angstrom`; split comparisons atol `1e-10`, rtol `1e-12`.
No newly complete polymer exists, so no new LAMMPS force run, finite differences,
minimization or dynamics were launched. Existing independent numerical and
angle-angle correction evidence is preserved exactly, not regenerated.

## Reproduction

From the repository, the actual compatible environment is:

```sh
PY=../island-validation/phase4e2-env/bin/python
D=../lammps/lammps/tools/msi2lmp/frc_files
OUT=../island-validation/phase4j16-declared
```

Use **new output paths** to reproduce; publication refuses overwrites.

```sh
$PY scripts/audit_pcff_distribution.py \
  --frc "$D/pcff.frc" --rlb "$D/pcff.rlb" --templates "$D/pcff_templates.dat" \
  --declaration docs/evidence/phase_4j16_declaration.json --output NEW/audit.json

# Build receipt identifies the already executed, unmodified local build.
make -C "$OUT/converter-build"
$PY scripts/validate_pcff_distribution.py \
  --frc "$D/pcff.frc" --rlb "$D/pcff.rlb" --templates "$D/pcff_templates.dat" \
  --converter "$OUT/converter-build/msi2lmp.exe" \
  --converter-receipt "$OUT/converter.json" --output NEW/experiment
# Default executable-polymer gate exits 1; --gate converter checks boundaries only.

$PY scripts/check_pcff_distribution.py \
  --input "$OUT/boundaries-corrected" \
  --frc "$D/pcff.frc" --rlb "$D/pcff.rlb" --templates "$D/pcff_templates.dat" \
  --converter "$OUT/converter-build/msi2lmp.exe" \
  --output NEW/checked.json --gate converter

$PY scripts/audit_pcff_completion.py --source "$D/pcff.frc" \
  --declaration docs/evidence/phase_4j11_declaration.json \
  --retained-matrix ../island-validation/phase4j10-declared/parity-matrix \
  --output NEW/full-source.json --require-full-source
# Expected exit 1, unchanged historical full-source requirements.

$PY -m pytest -q tests/test_pcff_distribution.py tests/test_pcff_distribution_gates.py \
  tests/test_pcff_urethane_domains.py tests/test_pcff_urethane_acceptance.py \
  tests/test_pcff_validation_cache.py tests/test_pcff_tranche_controls.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

The example `examples/pcff_local_distribution.py` accepts the same three explicit
paths plus `--frc-sha256`, `--rlb-sha256`, `--templates-sha256`; it reports unsupported
semantics and cannot prepare a model. Inspection/reconstruction requires only
ISLAND's ordinary Python dependencies, not a converter or scientific backend.
The experiment's reference readers additionally use NumPy. The two pip environments
remain separate. Main environment: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1,
RDKit 2026.3.6, OpenMM 8.6.1.

Final executed results, exact artifact hashes, full missing requests and historical
read-only checks are indexed in [the J16 receipt](evidence/phase_4j16.json).

### Executed verification results

- Focused distribution, gate, J15 domain/acceptance, cache and J13 control tests:
  **87 passed** in 12.78 s. These include labelled synthetic grammar and gate fixtures.
- Complete ordinary suite after the final serialization-check correction:
  **1795 passed, 10 skipped** in 347.64 s. The earlier full run is also retained
  (1794 passed, 10 skipped); the additional test is the serialization regression.
- Ruff: **passed**. Both environment `pip check` commands: **no broken requirements**.
- Local source inspection CLI and example: **passed**.
- Final artifact checker `--gate converter`: **exit 0**, audit/boundary gates true;
  executable-polymer/full-source false. Its separate receipt does not overwrite
  the failed assertion in the original run.
- Historical `audit_pcff_completion.py --require-full-source`: **exit 1** as required.
- **14** retained bundle/workflow sets loaded, status/frame reads and completed
  no-op resume passed with scientific imports blocked; **5 frames each**.
- Real distribution-audit relocation and child-process reconstruction: **passed**,
  original paths unnecessary, all three explicit external copies hash-verified.
- **4474 historical file hashes unchanged**, all three user-source hashes unchanged.
- All **20 J15 case outcomes unchanged**; the complete current v6 coverage ledger
  also equals the retained J15 ledger. V6 has **87/133 bounded predicate labels**,
  **0/133 globally certified labels**. The separate frozen J11 obligation audit's
  47 unresolved labels is not a current-profile predicate count.

The 14 read-only roots are recorded exactly in the receipt. Reproduce that check
with `scripts/revalidate_pcff_tranche.py --source "$D/pcff.frc"`, one `--root`
argument per recorded root, and a new `--output`. The roots span J8/J9/J10/J12,
J13 PEO/thioether, seven J14 cases and J15 aniline. Neither that checker nor the
new distribution loader constructs a Context or reruns parameterization.

**Unmet:** authoritative executable auxiliary semantics, new urethane/arylamine
polymer coupling coverage, new-polymer independent energy/force/minimization and
workflow gates, and full-source completion. This phase does not promise executable
polyurethane from the presence of three files.

The existing PSMILES urethane example was also executed against the supplied FRC,
in both strict and compatibility modes:

```sh
$PY examples/pcff_urethane_typing.py --source "$D/pcff.frc"
$PY examples/pcff_urethane_typing.py --source "$D/pcff.frc" --converter-compatibility
```

Both exit **1** with complete charges and incomplete models, as before. Their
[separate receipt](evidence/phase_4j16_examples.json) preserves the commands,
diagnostic counts and output hashes. These example constructions do not replace
any frozen J15 input or count as executable-polymer acceptance.
