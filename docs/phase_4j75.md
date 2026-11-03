# Phase 4J75 — authoritative PCFF source-path correction

J75 corrects the source-root declaration and performs read-only intake. The
result is **parameter_incomplete**: all **378** native PCFF requests remain
unresolved (**294 analog-only, 60 convention-blocked, 24 source-missing**).
No execution profile, parameter bundle, automatic assignment path, converter
result or static evaluator result is created.

## Base and scope

Branch `codex/phase-4j75-authoritative-pcff-source-path-correction` starts
directly at J74 `7536e143a563b8111d0f4a6c66c296e5456fd4fd`.
Fetched `origin/main` is `c118ddc06a7bd2532bea93c5e8d6ce17727726ef`;
it is already in the base ancestry. Main is untouched. Historical commits
and unrelated untracked files are preserved.

The J72-C1 kernel/process fence was installed and verified before the first
J75 source-directory inspection. J74's additive J73/J74 Python restrictions
are also required. Source inspection and every verification command use this
fence. Read-only Python implements the requested `test -d`, `realpath` and
recursive `find -type f` semantics; no shell utility bypass is needed.
Missing or uninstalled fences block intake before source metadata access.

J75 has **zero new QM attempts, MD steps, converter attempts, fitting
executions, evaluator executions, parameter-generation executions and reserved
validation executions**. Synthetic process tests create marker-only fake files;
no real scientific executable, integrator or converter is invoked.

## Exact root identity

The only accepted source spelling is:

```
/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/island_ffref/pcff
```

The directory exists, is readable, is inside this checkout and resolves to
that exact path. There are no symlink components. The path report retains
device, inode, permissions and `/proc/self/mountinfo` evidence. There is no
separate mountpoint at the source root; no independent bind-mount provenance
is claimed. A symlink root or child is rejected rather than followed as an
alternative source. An absent exact root produces `source_path_unavailable`
and empty source inventories and parameter facts.

The root manifest was written **before format parsing**. It freezes the
complete list of 27 physical files, their hashes, sizes, permissions and
readability. A second inventory after parsing matched the first. Source
files were neither modified nor moved. The incorrect `/island_ffref/pcff`
path was never inspected or used by J75.

J74's recorded checkout root is this same absolute directory. All 27 current
files match its frozen hashes and sizes. The comparison is between current
bytes at one physical root and historical identities; it is not a comparison
of independently obtained author artifacts. J74 did not retain an inode, so
no historical inode equality is invented. J74's blocked receipt and initial
inspection sequencing disclosure remain unchanged.

## Parsing, provenance and request matching

The capture retains **1,012 literal parameter facts**: 358 from `pcff.frc`,
72 from `cff91.frc`, 103 from `compass_published.frc`, 136 from benchmark
implementation text and 343 hash-matched J71 literature records. Twenty-one
PDF copies yielded text on all 385 pages, including duplicate copies.
`.ipynb_checkpoints/ma00107a006-checkpoint.pdf` is truncated and failed
parsing. This is text extraction, with no OCR, graphical inference or new
claim of manual full-paper reading. J71 numeric records and locations are
inherited unchanged; old literature decisions are not re-adjudicated.

The ZIP's 32 members are read as data in memory, with member and container
hashes. Static topology text supplies source type tuples and charge
declarations only. No coordinates are evaluated, no input script is executed,
and no converter is called. FRC extraction retains relevant HMDS/Si source
rows and counts unrelated rows; it does not transfer unrelated parameters.
Raw numeric lexemes, equilibria, signed coefficients, literal equation/unit
headers, source order, source versions and locations remain separate from
unresolved convention fields. Analytic force metadata states the source
derivative relationship; no force calculation is performed. Missing source
equations and conventions are explicit unresolved evidence, never invented.

Seven categories remain separate: native PCFF/CFF91 distribution evidence,
published COMPASS, external PDMS benchmark, Huang PDMS, Smith/Borodin/Smith
PDMS, Sun/Rigby polysiloxane, and analog/secondary/contextual literature.
The original model subtype is also retained, including COMPASS, OPLS and
united-atom benchmark subdirectories. The family reports concern the exact
J73 native PCFF requests; they do not assert completeness of another model.
Authors' names alone do not classify a model. Sun's CFF93 aromatic-polymer
papers, the Rigby PEO paper, mixed benchmark literature and PCFF/COMPASS
comparisons remain contextual. Individual inherited records use their own
recorded model identity, keeping mixed supplements separate at row level.

`pcff.frc` has SHA-256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
Its retained upstream identity is LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`, path
`tools/msi2lmp/frc_files/pcff.frc`, with literal release history through
4.0 (8 October 2013). This establishes distribution identity, not new
author-primary evidence or physical transferability. Other copied assets
retain unverified provenance where appropriate. License metadata does not
infer a redistribution grant. No PDF, FRC, ZIP or other source asset is added
to Git.

Every one of the 378 records binds its original J70/J72/J73/J74 identity,
ordered atom types, stable sites, graph neighborhood and required evidence.
Queries use exact types and explicit ordinary FRC equivalences; no wildcard
or automatic fallback is allowed. Neither native FRC supplies a missing
row. Candidate comparisons report type, graph, function, units, phase,
multiplicity, order, equilibrium, exclusions, charge and LJ conventions
separately. External benchmark entries are family candidates, with numeric
type IDs never interpreted as native PCFF types. An unresolved field cannot
become a successful comparison.

The 28 user-attested confirmations remain provenance metadata;
`author_confirmed_audited=0` and original author artifacts remain
`not_yet_supplied`. Primary artifact, exact-row, convention and executable
qualification flags remain separate and false. The 16 historical equal
duplicate equivalences and two unequal unresolved COMPASS assignments are
preserved, with no averaging or selection by residual improvement.

## Offline integrity and verification

The explicit four-file J75 implementation inventory is in the manifest;
historical receipt inventories are not broadened. A pinned capture-freeze
hash binds source facts, path identities, original request matrices and
historical source inputs. Re-signing a fabricated source, changing units or
equations, removing a request, hiding duplicate candidates or activating a
profile fails reconstruction. The base inventory checks **9,696 historical
files** byte-for-byte, including all J17–J74 evidence, J73 requests, the J74
blocked result, J72's historical violation and the C1 quarantine.

Terminal software results and exact unique node identities are recorded in
`validation/collection.json`, `focused.xml.xz` and `checks.json`. Focused
J73/J74/J75 tests are software evidence. Live science, converter and
evaluator controls are explicit scope skips. Broad historical scientific
suites are not launched. Historical receipt replay and byte preservation
are recorded separately from newly executed test coverage.
The terminal run covers **381 unique tests: 355 passed, 26 skipped**.
The separate relocated-replay rerun overlaps this set and adds no unique
coverage. Earlier test runs are retained without adding their counts.

Ruff, compileall, supported-environment pip check, diff checks, historical
receipt validators, fresh-process offline replay, relocated replay and source
byte rechecks are fenced. The separate PDF parser environment retains its
existing pip-check failures: torch lacks fsspec and requires sympy 1.13.1
instead of installed 1.14.0. No dependencies are changed. All attempt logs,
including any failed verification attempts and corrections, are retained;
focused reruns are not counted as new terminal coverage.

The preliminary test runs and comparison reports are retained under
`validation/attempts/` and the fenced-run/log inventories. One test run failed
because a new relocation test lacked a `json` import; Ruff also rejected that
version. The source-family audit then corrected overly broad author-name and
mixed-citation classification. Its first recapture stopped on a rendered
citation string; the parser now handles both rendered strings and structured
citation metadata. The successful recapture has exactly the same immutable
root manifest and 1,012 raw facts. The superseded capture is losslessly
retained, with no historical phase or parameter value changed.

## Replay

From a checkout, using the supported interpreter (replace each `/tmp` output
with a fresh directory):

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j75-replay -- -m scripts.phase_4j75_source_path replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j75-validate -- -m scripts.validate_phase_4j75_receipt
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j75-tests -- -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j73.py tests/test_phase_4j74.py tests/test_phase_4j75.py
```

Replay and strict receipt validation require no source directory, pypdf,
network or scientific software. `originals` is a separate optional command
that rechecks only the exact authoritative directory; it cannot accept a
substitute. Relocation copies metadata, historical evidence and explicit
implementation dependencies into a fresh temporary directory; it does not
copy source assets or resolve a new source path.

## Scientific limits and dependency

J72's **four unintended external SQM attempts** and **46 verified / 48
possible MD steps** remain historical violations:
`J72_zero_execution_scope_met=false`. J75 zero-execution counts do not erase
that history. Quarantine artifacts remain excluded from parameter decisions.

Source-root correction establishes which local bytes were parsed. It does
not complete native source rows, resolve native Phi/AAT/dependency authority,
qualify a converter, reproduce QM targets or establish chemical transfer.
The next dependency remains exact authoritative native rows and complete
conventions, with original author/publisher evidence. No new scientific
execution is authorized by this receipt. `accepted_coefficients=0`,
`automatic_assignment=false`, `production_validated=false`, and
`simulation_readiness="not_established"`. Ordinary `forcefield="pcff"`
behavior is unchanged. J68's separate incompatibility remains unchanged.
