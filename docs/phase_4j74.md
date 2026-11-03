# Phase 4J74 — local PCFF source intake and exact-row resolution

The source-resolution result is **blocked**. All **378** J73 requests remain
unresolved: **294 analog-only, 60 convention-blocked and 24 source-missing**.
The local assets supply useful source facts, but no complete exact native PCFF
row contract for these missing instances. No profile or parameter bundle is
created, and no source is promoted by its filename, numerical similarity or
presence on this server.

## Base and immutable inputs

Origin was fetched before implementation. Fetched `origin/main` was exactly
`c118ddc06a7bd2532bea93c5e8d6ce17727726ef`, the merge of J72-C1. J74 starts
from J73 `dff35a0dfbca2c0517a9a7d5fe8e5032ac2d8375` and includes main through
normal merge commit `a0fc6003dd637dc4920eab03f0b881dccd1dbc33` on
`codex/phase-4j74-pcff-local-source-resolution`. Main itself is untouched.
No historical commit is rewritten or force-pushed.

The frozen base inventory contains **9,539 files**. It includes all J17–J73
evidence, historical receipt identities, the J73 request records, the J72
execution violation and the C1 quarantine. The explicit four-file J74
implementation inventory does not broaden historical receipt globs. Unrelated
untracked files remain outside the commit.

## Actual directory and source capture

The requested absolute `/island_ffref/pcff` directory does not exist in this
environment. The discovered checkout-relative `island_ffref/pcff` copy is
explicitly recorded as the actual read-only intake directory; it is not
silently aliased to the requested path. The missing declared path remains a
blocked-receipt reason.

The inventory covers **27 physical files**, with exact relative paths, sizes,
SHA-256, MIME types, parser identities, extraction status, source family,
citations, versions, copy/authority status and license/redistribution metadata.
The benchmark ZIP also has a **32-member** inventory in the artifact ledger,
with each member's hash, size, family and provenance. All members were read as
data in memory; no archive command or scientific input was executed.

Twenty-one PDF copies yielded text on all pages: **385 parsed pages**, including
byte-identical checkpoint copies. One truncated checkpoint PDF,
`.ipynb_checkpoints/ma00107a006-checkpoint.pdf`, failed parsing and remains an
explicit blocker. PDF processing used pypdf 6.20.0, without OCR or graphical
inference. Text parsing is not a new claim of manual full-paper reading. Exact
hash-matched J71 extraction records retain their original identities and
locations; J74 does not re-adjudicate them.

The final capture retains **1,012 parameter facts**:

| Source | Retained facts | Authority |
|---|---:|---|
| Local `pcff.frc` | 358 | Exact copy of the pinned historical LAMMPS asset; no new author-primary evidence |
| Local `cff91.frc` | 72 | Native-named distribution candidate; upstream provenance not independently established |
| `compass_published.frc` | 103 | COMPASS evidence, never native PCFF |
| PDMS benchmark archive | 136 | External implementation facts, separated by model subdirectory |
| Exact hash-matched J71 PDF records | 343 | Historical literature facts, unchanged and unqualified for missing native rows |

FRC extraction parses all section rows and retains HMDS component-only and
explicitly Si/osi-related numeric records. Unrelated rows are counted rather
than transferred. Raw numerical lexemes, source order, signed coefficients,
equilibria, versions, references and section/line identities are preserved.
Literal equation and unit headers remain distinct from independently closed
conventions. No numerical unit normalization, derivative calculation, energy,
force or pressure evaluation is performed.
Missing-row eligibility flags apply to this J74 intake; they do not revise
historical coverage for already implemented native source terms.

Benchmark topology parsing reads only type/interaction mappings and observed
charge declarations. Source numeric type IDs do not become PCFF types.
Coordinates are never evaluated or committed. Instantiated coefficient rows
are retained; unused universal-library rows are counted without fallback.
Source zeros are literal source facts, never replacements for absent terms.
The family inventory explicitly includes bond, angle, torsion, improper, all
requested Class-II couplings, torsion-torsion, nonbonded and bond increments.
An absent extracted family is not a zero-valued coefficient.

## Source authority and row resolution

The local `pcff.frc` SHA-256 is
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
Its literal release history reaches 4.0, 8 October 2013. Its bytes match the
retained asset at `tools/msi2lmp/frc_files/pcff.frc` in LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`. This establishes which distributed
source was inspected; it does not supply any of the 378 missing rows or prove
author authenticity, transferability or physical validity.

Every row-resolution record binds its J70/J72/J73 identities, ordered atom
types, stable graph sites, local neighborhoods, family, candidate facts and
source references. Native queries record exact and documented ordinary
equivalence lookup, with no auto-equivalence, wildcard or lower-family
fallback. Neither local FRC supplies a missing exact/equivalence row.

External candidates retain independent type, graph, function, unit,
convention and equilibrium-dependency comparisons. Benchmark candidates are
also indexed by family without treating numeric atom types as native PCFF
types. The separate coverage reports cover native PCFF, native COMPASS,
external Huang, external Xiang benchmark and other external/analog sources.
Their unresolved counts refer to the required native PCFF contract, not a
claim about the completeness of a different model.

Huang benchmark commands are external implementations, not the missing
Huang author-primary parameter artifact. Sun/Rigby and Smith/Borodin/Smith
facts keep their historical family, chemical domain and equations. COMPASS,
OPLS, united-atom, PDMS-specific models and other analogs cannot supply native
PCFF authority. All alternatives across local releases are preserved without
averaging, selecting favorable values or inferring conventions.

The 16 historical equal duplicate equivalences and two unresolved unequal
COMPASS assignments remain unchanged. This phase resolves no old conflict.
The seven row states remain independent: analog-only, convention-blocked,
source-missing, primary-artifact-supplied, exact-row-resolved,
convention-resolved and executable-bundle-qualified. The last four remain
false for every request. The 28 user-attested confirmations remain provenance
metadata, not audited author-primary evidence.

## Execution fence, verification and limitations

The C1 preflight installed successfully before the fenced capture. The initial
read-only filesystem listing, FRC header/search inspection and ZIP member
listing preceded that explicit J74 preflight. This sequencing deviation is
retained as `initial_sequence_compliant=False`; no receipt claims otherwise.
It launched no scientific executable. All subsequent capture, PDF/text
parsing, tests, replay and terminal checks use the unchanged C1 launcher plus
the J73 import guard and J74 named-driver rejection layer.

The gate is the supported cooperative verification policy, not a sandbox for
arbitrary hostile Python code. Synthetic probes contain marker-only fake
executables and modules. Live scientific controls are explicit scope skips.
The phase launches no QM/DFT, AmberTools, OpenMM/LAMMPS dynamics, minimization,
converter, fitting, parameter generation, reserved validation or scientific
energy/force/pressure evaluation.

The supported verification environment passes `pip check`. The existing PDF
parser environment reports two PyTorch dependency conflicts: missing `fsspec`
and SymPy 1.14.0 versus the required 1.13.1. That failed check is retained with
its exact diagnostic, not reported as passed. pypdf extraction did not import
or use PyTorch. No dependencies were installed or changed in this intake.

Terminal software coverage and retained correction history are recorded under
`validation/`; overlapping runs do not add terminal unique coverage. The first
lint run found import ordering, an alias, a broad exception boundary and an
unused variable. These were corrected before terminal tests; the failed log
is retained. The first successful source capture is also retained: the final
capture adds literal FRC reference/license-notice metadata without changing
source bytes or adjudicating source values.

The terminal run contains **288 unique test nodes: 269 passed and 19 explicit
scope skips**, with no failures or errors. This is the union of the complete
J73 and J74 focused software files; replayed historical test counts are not
new coverage. The earlier run passed 263 controls with 19 skips. It overlaps
the terminal run and is retained without double-counting; six added controls
verify disclosure of the parser-environment failure and exclusion of original
source assets from the manifest. Ruff, compileall, supported-environment pip,
source-byte checks and offline replay pass. No live scientific control is
reported passed. Strict validation accepts the disclosed blocked source
outcome and retains the parser-environment pip failure.
The first standalone relocated launch selected the J73-only fixture directory
and correctly rejected missing J74 inputs. Its failed output is retained; the
correct J74 fixture is used for the successful independent relocated replay.
The source-inspection stdout containing external README/parameter-file
excerpts remains outside Git for licensing; its path, size and hash are
retained instead. The separately extracted facts and all execution ledgers
are included.

## Offline commands

Each C1 launcher output directory must be fresh. Replay and strict validation
require neither local PDFs/FRC/ZIPs nor pypdf, network or scientific backends:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j74-replay-fresh -- \
  -m scripts.validate_phase_4j74_receipt docs/evidence/phase_4j74 --replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j74-strict-fresh -- \
  -m scripts.validate_phase_4j74_receipt docs/evidence/phase_4j74

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j74-tests-fresh -- \
  -m pytest -p scripts.phase_4j73_scope \
  tests/test_phase_4j73.py tests/test_phase_4j74.py -q
```

Optional source-byte verification reads the explicit local copy without
parsing or executing it:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j74-originals-fresh -- \
  -m scripts.phase_4j74_source_intake docs/evidence/phase_4j74 \
  --verify-originals --source island_ffref/pcff
```

Validation reconstructs every derived report from hash-pinned facts and
historical requests. Re-signing outer identities cannot authenticate a new
source, conceal a conflict, alter units/equations or promote row/profile
status. Relocated replay uses the explicit JSON inputs and scripts only.
Historical scientific regression workflows are not launched in this strictly
read-only allocation; unchanged C1/J73 receipts and full historical-byte
checks verify preservation instead.

## Remaining dependency and capability

Exact authoritative native rows and complete Phi/AAT, equilibrium dependency,
typing, exclusion, charge and LJ conventions are still required. The complete
author/publisher files and provenance chain requested in J73 remain the next
dependency. Source absence alone does not establish a QM-fit requirement, and
J74 authorizes no fitting or scientific execution. J68's model-span result
remains separate and unchanged.

All J74 execution counts are **zero**. Historical J72 still records four
unintended external SQM attempts and 46–48 MD steps, with its C1 quarantine and
exclusion checks intact. J72 is not retroactively marked compliant.

**accepted_coefficients=0; automatic_assignment=False;
production_validated=False; simulation_readiness="not_established".**
No profile, bundle, converter result or MD workflow is created. Ordinary
`forcefield="pcff"` behavior is unchanged. Chemical transfer, scientific
validation and production readiness remain unestablished.
