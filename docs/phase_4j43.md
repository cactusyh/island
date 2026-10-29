# Phase 4J43 — local PDMS literature asset ingestion

J43 adds an explicit, hash-bound literature ingestion boundary alongside the
unchanged C5 external Class-II bundle schema. No PDMS literature asset, parameter
bundle, independent pin or license was supplied for this phase. The executed
result is **asset_unavailable**, with acceptance **1**, no execution profile and
no numerical backend execution. There are no extracted physical parameter rows.

The fetched and verified base is origin/main at
`fdf8e8058173faa9c013703aa4910b3065820e7c`, the merge of J42-C1. The branch is
`codex/phase-4j43-literature-pdms-parameter-ingestion`. The tested implementation
is `fc77dfcd77c48f7bf6c401d0aebcfe6979da3ffe`. Main is not modified or merged.

`production_validated=False`; `simulation_readiness="not_established"`.

## Input and provenance contract

`island.forcefields.pdms_literature.ingest_literature()` accepts a final graph,
an explicit local source-manifest path and its independent SHA-256. An optional
graph-bound C5 bundle requires its own explicit path and independent SHA-256.
No files are searched for, downloaded or inferred from a directory. Source
dependencies are read only through the manifest's hash-listed relative paths.
Path escapes and escaping symlinks are rejected. Pins are checked before parsing.

The source manifest is ordinary JSON with schema
`island_pdms_literature_source_v1`. It contains exactly these fields:

| Field | Required content |
|---|---|
| `schema` | The versioned manifest schema above |
| `bibliography` | `doi`, `title`, nonempty `authors` list, integer `year` |
| `classification` | Exactly one classification listed below |
| `source_family`, `model_name` | An isolated `external:` family and explicit model name |
| `license` | `identifier`, `authorized_use=True`, `redistribution` (`allowed` or `local_only`), `evidence_file` |
| `files` | Relative filename to exact `sha256`, `family`, `role`, `provenance` |
| `primary_source` | Explicit numerical row-table file from `files` |
| `parameterization` | `method`, `quantum_chemistry` (or explicit null), `condensed_phase_targets` list |
| `functional_form`, `units` | `class2` and the existing C5 `UNITS` record |
| `charge_method` | Explicit source method; no normalization or redistribution |
| `atom_type_definitions` | Source labels mapped to their chemical definitions |
| `reported_validation_properties` | Author-reported properties, preserving their source attribution |

File roles are `paper`, `supplement`, `raw_table`, `row_table`, `code`, `license`
or `benchmark`. Every file belongs to the same external source family. License
evidence must be a nonempty hash-bound file with role `license`. A local-use-only
record does not authorize redistribution of its source files. Source metadata is
a supplied provenance declaration; a checksum does not independently establish
authorship or scientific validity.

Every candidate gets exactly one classification:

- `pdms_specific_external`;
- `compass_class2_external`;
- `pcff_extension_external`;
- `opls_compatible_external`;
- `incomplete_literature_source`;
- `unusable_or_unlicensed`.

Incomplete extraction is recorded as `incomplete_literature_source`;
unauthorized or prohibited use is `unusable_or_unlicensed`. A classification
does not convert an OPLS functional form to Class-II or confer native coverage.
An incompatible functional form or unit system remains blocked.

The primary row table uses the existing C5
`island_external_class2_source_rows_v1` schema: `schema`, `source_family`, `rows`,
and `provenance`. Each row retains C5's exact `family`, ordered `types`, decimal
string `values`, `units`, `provenance`, and `citations`. Each citation contains
`file`, one-based `line`, exact `raw` text, and zero-based numeric `columns`.
A table location must therefore be supplied as a reproducible text row and
line, with its table provenance retained in the row. PDF bytes alone cannot
supply that evidence. No numerical values are read from figures or neighboring
rows, and no PDF text is silently generated or trusted.

`reconstruct_rows()` independently reads the cited raw text, preserves extracted
decimal strings, and compares its results with the unchanged C5 source-row
reader. The immutable `LiteratureSourceRecord` retains the pinned manifest,
source rows, extracted citations and per-row diagnostics. Each extracted
citation records its source file, one-based line, raw text, numeric columns and
decimal values. Explicit source zeros can be reconstructed; invented or added
zeros fail the row check. These records provide data evidence only.

## Bundle and execution boundary

The C5 schema remains `island_external_pdms_class2_bundle_v1`. Its external
profile namespace remains `island_pcff_external_pdms_class2_v1`. No profile is
registered or instantiated in this phase. `profile_candidate` in an ingestion
report names the requested namespace and carries no execution authority.

The explicit bundle must independently pass C5's graph, atom order, component,
mass, charge, units, interaction coverage, equilibrium dependency, mixing,
1–4, source-row and provenance checks. Its family, files, primary row table and
complete source-row set must agree exactly with the supplied literature record;
every assigned atom type needs a source definition. Incomplete or unusable
source classifications cannot qualify a bundle. Rows are never combined across
PCFF, COMPASS, OPLS or another family. There is no automatic matching or inference
of physical parameters from an incomplete paper.

A complete data import would be recorded as `data_status="data_valid"` and
`status="scientifically_unqualified"`. The existing J42 external adapter remains
the preparation boundary through `ForceFieldRequest` and `prepare_forcefield`.
It does not authorize numerical execution. No complete supplied asset was
available to exercise that conditional data-valid branch or establish a new
energy/force adapter. No LAMMPS/OpenMM construction, component comparison,
finite-difference force, minimization, NVT/NPT, checkpoint, RNG or equilibrium
result is claimed. Positive qualification still requires all those gates with
an independently evidenced real bundle.

## Executed graph fixtures and blocked evidence

The five fixtures reuse the validated C5 final graphs. The repeat is
`[*:1]O[Si](C)(C)[*:2]`, with explicit hydrogens and methyl/hydroxyl terminal
groups for the PSMILES builds. The topology and valence are validated after
construction; the string alone does not authorize a parameter assignment.

| Fixture | Atoms | Missing interactions | Unresolved dependency edges | Status |
|---|---:|---:|---:|---|
| `[SiH3]O[SiH3]` | 9 | 116 | 167 | `asset_unavailable` |
| Hexamethyldisiloxane | 27 | 656 | 977 | `asset_unavailable` |
| Explicit-H capped repeat | 16 | 319 | 469 | `asset_unavailable` |
| PSMILES PDMS DP2 | 26 | 576 | 852 | `asset_unavailable` |
| PSMILES PDMS DP4 | 46 | 1,090 | 1,618 | `asset_unavailable` |

Each ingestion record preserves all exact interaction keys and dependencies,
not just these counts. The literature inventory is empty, and the extracted-row
file explicitly records that no user-supplied source exists. All evaluator
attempts raise `ExternalPDMSBlockedError` with the retained preparation result.
Relocation and a fresh process reconstruct all five cases; all 18 study files
remain byte-identical.

The focused tests use a nonphysical `parser_sentinel` row solely to exercise
decimal/citation reconstruction and evidenced-zero rejection. It is never a
PDMS coefficient, literature source, executable fixture or positive qualification
result. Tests cover independent pins, malformed sources, license and provenance
gaps, PDF-only incompleteness, mixed families, units, escaped paths, citation
mutations, missing bundle records, changed graphs, forged authority, receipt
checksums and fresh-process reconstruction of relocated incomplete assets.

## Verification and historical preservation

The receipt records focused tests, the J22–J42-C1 regression subset of the
supported full-suite execution, Ruff, compileall, five pip checks, diff checks,
receipt/manifest verification and relocation. Exact commands, environments,
JUnit results and node IDs are archived. The focused suite passes **42 tests**;
the supported full suite passes **3,467 tests**, with **10 existing optional
skips** and no failures or errors. Its exact J22–J42-C1 regression subset passes
**1,457 tests**. The first orchestration attempt lacked
the pytest temporary parent directory; its setup failures and aborted full-run
logs are retained under `initial-harness-run`. The directory was created before
rerunning the affected checks; no implementation change was required.

A subsequent new relocation test stripped the serialized record's final newline
before comparing bytes. That test was corrected to preserve the complete stdout
record. Its 41-pass/one-failure output and the aborted accompanying full attempt
are retained under `initial-newline-test-run`. The runtime ingestion code did
not change. The published validation runs use the corrected test commit above.

The independent source-availability audit confirms zero supplied physical source
rows and retains every fixture's missing dependencies. Its initial script used
the generic study decoder for an immutable record envelope; the corrected audit
uses `LiteratureIngestion`. That initial script and failure are retained under
`initial-source-audit`. Actual literature-row reconstruction remains blocked by
absent input; the parser controls are not presented as a literature audit pass.

Independent Git-blob reconstruction preserves all **2,538 historical docs and
evidence files** and **651 existing source/test/script files** at the base.
No historical file is intentionally changed. The unchanged J42-C1 receipt and
manifest validate with the current runtime. The unchanged J42 validator checks
its original implementation pins against its original J42 detached checkout,
while reconstructing the retained study with the current runtime. Historical
hashes and validators are not rewritten.

Polymer Ewald/PME remains `ewald_unvalidated`; backend pressure remains
`pressure_unavailable`. Cutoff recovery and long NPT remain blocked. There is no
experimental agreement or thermodynamic equilibrium claim. Native PCFF, native
OPLS-AA, published COMPASS, GAFF/GAFF2, RESP/Psi4 and universal msi2lmp coverage
are not promoted by this external literature boundary.

## Reproduction

```bash
git fetch origin
git rev-parse origin/main
# Verified implementation base: fdf8e8058173faa9c013703aa4910b3065820e7c

export PYTHONPATH=.:src
export PYTHONDONTWRITEBYTECODE=1 OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
PYTHON=../island-validation/phase4e2-env/bin/python

# Exit 1 is the expected blocked scientific/software acceptance result.
"$PYTHON" scripts/phase_4j43.py --output /tmp/island-j43-study
cp -a /tmp/island-j43-study /tmp/island-j43-relocated
"$PYTHON" scripts/phase_4j43.py --reconstruct /tmp/island-j43-relocated

# Only with explicitly supplied, independently pinned assets:
"$PYTHON" scripts/phase_4j43.py --output /tmp/island-j43-supplied \
  --source /explicit/local/source-manifest.json --source-sha256 "$SOURCE_PIN" \
  --bundle /explicit/local/graph-bound-c5-bundle.json --bundle-sha256 "$BUNDLE_PIN"
# Different final graphs require separately bound bundles. No complete asset
# was supplied for the published run; this conditional example is not evidence.

"$PYTHON" -m pytest -q tests/test_phase_4j43.py
"$PYTHON" -m ruff check .
"$PYTHON" -m compileall -q src scripts tests
git diff --check
"$PYTHON" scripts/validate_phase_4j43_receipt.py docs/evidence/phase_4j43
"$PYTHON" scripts/validate_phase_4j42_c1_receipt.py docs/evidence/phase_4j42_c1
```

Full-suite, regression, five-environment pip, historical snapshot and independent
byte-audit commands are retained in `validation.tar.xz`. The full suite uses the
same supported optional local assets as J42-C1. The regression report is an exact
subset of that execution, not a second claimed test run. No source SHA-256 is
fabricated for absent literature; actual supplied file pins would appear in the
source record and its extraction evidence.
