# Phase 4J72-C1 — execution scope correction and quarantine integrity

The decision remains **`source_conflict_unresolved`**. This is a software and
evidence-integrity correction. It adds no papers, parameter values, scientific
fits, QM, MD, converter runs or execution profiles. Ordinary
`prepare_polymer(..., forcefield="pcff")` behavior is unchanged.

The actual base is required J72 **30737605b8252e2bf0f2265352f0937ee931fd6d**.
Fetched origin/main was **4c9066e4501591211049c44401e82b3110118f58** and did
not contain it. The branch is
`codex/phase-4j72-c1-execution-scope-correction`. No main merge or push is
performed. Unrelated untracked files and mounted literature are preserved.

## Historical violation and corrective scope

J72's initial guard did not propagate to J31 child processes. Four unintended
AmberTools SQM AM1 minimizations completed for CHO AM1-BCC regression fixtures.
Retained journals establish **46 MD steps**, with **up to two additional
uncommitted steps** possible in the interrupted final resume. These executions
were unauthorized. Their outputs do not become scientific evidence because
they occurred in tests. J72's original zero-execution contract, actual-count
receipt, interruption records, source decisions and all logs remain unchanged.

| Namespace | Value |
|---|---:|
| `authorized_new_QM` | 0 |
| `authorized_MD_steps` | 0 |
| `HMDS_new_QM` | 0 |
| `reserved_validation_QM` | 0 |
| `unintended_external_QM_attempts` (J72) | 4 |
| `unintended_MD_steps_verified_minimum` (J72) | 46 |
| `unintended_MD_steps_possible_maximum` (J72) | 48 |
| `zero_execution_scope_met` (historical allocation) | false |
| `C1_new_native_QM_attempts` | 0 |
| `C1_MD_steps` | 0 |
| `C1_scientific_fits` | 0 |
| `C1_converter_attempts` | 0 |

There is no unqualified zero-QM claim. The historical violation remains false
for scope compliance even though this corrective execution adds no science.

The [quarantine manifest](evidence/phase_4j72_c1/quarantine-manifest.json)
independently reads the original archive without extracting or executing its
contents. It verifies all **192 member paths, byte sizes and SHA-256 hashes**,
reconstructs the four completed SQM outputs and the six MD journal comparisons,
and preserves command lines verbatim from retained stdout. Original native
PIDs, execution-time binary hashes and separately captured SQM stderr were not
retained; their absence is explicit. Parent stderr and original native output
files are retained. `quarantine_complete` means that every originally retained
member is verified, not that missing PIDs or uncommitted steps were recovered.
The archive remains under J72, byte-for-byte; C1 does not duplicate or edit it.

## Process fence

The opt-in launcher installs the fence **before Python test startup**. It
requires Linux x86-64, seccomp and parent-child ptrace support. Installation
failure stops the run. The ordinary package API never imports this machinery.

The kernel filter is inherited through fork, clone and exec. It stops every
`execve`/`execveat` before execution. A supervisor checks the actual filename,
requested argv, resolved executable hash and environment against a frozen
allowlist. Only declared Python runtimes, read-only Git inventory commands and
Ruff are allowed. Shells, native scientific executables, alternate executables,
Python startup bypass flags and unapproved binaries fail closed. Nonexistent
PATH candidates retain their normal ENOENT behavior; they are logged separately
and do not count as executed programs. `execveat` is rejected rather than
guessing descriptor/relative-path semantics. The supervisor never replaces a
requested command. A denied raw exec child is killed before its syscall;
Python callers receive `ScopeDenied` before launch. `EXITKILL` kills tracees if
the supervisor exits; a retained seccomp TRACE filter fails exec closed without
a tracer.

Python `sitecustomize` verifies policy and implementation hashes in every
descendant. Subprocess environments retain the policy even when a caller
supplies a reduced environment. Startup checks require an active kernel
supervisor. Import loaders permit read-only `find_spec()` discovery while
rejecting execution of QM backend modules and Python LAMMPS/Sander backends. OpenMM integrator construction,
stepping and integrator deserialization/cloning are blocked, including native
SWIG entry points. Static type registration, source parsing and System XML
remain available. Pure Python dynamics steps are blocked. Real HMDS conic fits
are blocked; small manufactured numerical controls remain software tests.

Every rejection records the unmodified command, executable path/hash and
synthetic-fixture identity, process/parent IDs, environment policy, reason and
UTC timestamp. Allowed execs and startup installations are also retained.
The event audit independently checks all allowed native launches against the
recorded policy. This is a fence for the supported, cooperative test suite,
not a security sandbox against hostile Python modifying its own memory or
loading arbitrary machine code.

Synthetic tests reproduce the original inheritance failure with a marker-only
Python executable and a parent-only monkeypatch. Corrective tests then reject
fake executables named `sqm`, `antechamber`, `sander`, `pmemd`, `lammps` and other
forbidden entry points. Raw `execv`, `posix_spawn` and libc exec probes exercise
the independent kernel gate. Fake module/import and fake integrator tests
write markers if execution occurs; no real scientific executable or integrator
is used by these probes.

## Contamination exclusion

The [exclusion report](evidence/phase_4j72_c1/contamination-exclusion.json)
checks an explicit inventory of literature records, source comparisons,
authority records, row decisions, bundles, converter/static gates, capability
records and historical HMDS observation/producer records. It rejects quarantine
paths and nonempty artifact hashes in scientific inputs, including renamed
hash references. Audit/scope ledgers are separate and intentionally describe
the violation. Empty-log hashes are not mistaken for scientific contamination.

Checks independently join all 467 normalized records to original J71 producer
identities, every candidate authority reference to those identities, and all
378 adjudicated sites to the exact J70 instance identity/site mapping. They
preserve the 294 analog, 60 convention-unresolved and 24 source-missing statuses,
16 equal duplicate equivalences and two unequal unresolved COMPASS assignments.
No source is re-adjudicated. Historical receipts are replayed only for
regression validation; no new document reading or parameter extraction occurs.

The frozen base-tree byte inventory and input hash anchor prevent an attacker
from substituting scientific dependencies and simply re-signing the outer
manifest. The validator reconstructs quarantine contents and scope declarations
and uses exact types: booleans cannot stand in for counts. Tests inject
quarantine evidence into each scientific role, bundle hashes, HMDS counts and
trajectory references, and re-sign receipts/manifests/nested declarations.
Those mutations are rejected. An additional fresh launch installs its own
fence from the relocated repository and replays pristine relocated evidence;
neither the replay nor that launcher's policy requires the original repository
path. These combined checks support
`scientific_evidence_contaminated=False` and
`source_adjudication_excluded_quarantine=True`.

## Verification and retained corrections

The implementation inventory is explicit in
`scripts/validate_phase_4j72_c1_receipt.py`; it includes the fence, supervisor,
launcher, pytest policy, integrity reconstruction, strict validator,
verification driver, focused tests and this document. No historical glob or
implementation inventory is changed.

Terminal coverage is **5,978 unique tests: 5,178 passed, 800 skipped, zero
failures/errors**. The skips comprise 740 C1 execution-scope skips and 60
explicit optional-source/runtime/slow-check skips. All 145 focused C1 controls
pass. Whole-file corrective reruns cover all 29 historical files with initial
failures; overlapping runs add no unique coverage. All **9,071** base files
are preserved byte-for-byte. The 19 J54–J72 receipt validators, seven supported
pip checks, Ruff, compileall, diff and offline replay pass.

Verification results and terminal unique coverage are recorded in the
corrective receipt and `validation/terminal-selection.json`. Scope-dependent
checks are **skipped**, never reported as successful scientific verification.
The historical J72 integrator probe is explicitly skipped and replaced with
fake C1 probes. Focused reruns overlap with the full suite and add no unique
coverage. All original failures/interrupted runs remain in the validation
archive, including the first smoke decoding failure, two initial focused
bookkeeping failures, the overly broad registration/PATH checks, and the
deliberate interruption to fence integrator restore/clone entry points.

Adapter regressions also exposed code that catches a scope rejection and
returns a failed backend record. The pytest policy now converts a rejection
to an explicit scope skip in historical tests before an application handler
can swallow it. Descendant rejection records similarly trigger a scope skip;
the child's original command, return code and full captured stdout/stderr are
retained. C1's synthetic fence probes keep the actual `ScopeDenied` exception
so their assertions must pass. One J55 CLI parameter explicitly requires live
OpenMM runtime qualification and is declared a scope skip. Failed pre-correction
runs remain retained. The Python 3.10 tool environment's first pip check failed
at fence startup on a Python 3.11 datetime alias; no scientific fixture started.
The launcher bootstrap now provides the equivalent UTC alias for Python 3.10,
and the repeated pip check passes. A final fake-module probe also covers the
Python LAMMPS import path; no real LAMMPS module was imported by the suite.

The retained launch/check records document the supported full suite,
J54–J72 receipt regressions, Ruff over supported `src scripts tests`, compileall,
supported environment pip checks, diff checks, fresh-process replay and
historical byte checks. Archived evidence Python files are evidence, not Ruff
targets; the initial overly broad `ruff check .` failure is preserved.

Run focused software verification through the required launcher:

```bash
PYTHONPATH=../island-validation/phase4j66-solvers-only:.:src \
  ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 \
  --output /tmp/island-c1-focused-fresh -- -m pytest \
  -p scripts.phase_4j72_c1_test_policy tests/test_phase_4j72_c1.py -q
```

The output directory must be new, so an earlier attempt is never overwritten.
For the full supported software suite, replace the final command with
`-m scripts.verify_phase_4j72_c1 suite --output /tmp/island-c1-suite-fresh --workers 4`.
That driver reuses seven declared historical fixtures and retains every node.

Offline replay needs only repository code and retained evidence, with no
PDFs, network, converter, solver or Git subprocess:

```bash
PYTHONPATH=.:src python -m scripts.validate_phase_4j72_c1_receipt \
  docs/evidence/phase_4j72_c1 --replay
```

Strict receipt/manifest validation also verifies the entire frozen historical
tree, recorded terminal coverage and supervised execution ledger:

```bash
PYTHONPATH=.:src python -m scripts.validate_phase_4j72_c1_receipt \
  docs/evidence/phase_4j72_c1
```

These read-only commands can also be run through the launcher. Replay accepts
`--repo` and an independently relocated evidence directory; it does not use
execution-time absolute paths to find quarantined scientific inputs.

## Unchanged scientific state

J71's 12 papers and three supplements, 467 records and 1,470 cells remain
historical readings/extractions. C1 adds **zero**. Native PCFF authority and
external/analog provenance remain separate. All **378 rows remain unresolved**;
there is no complete native or external executable bundle. Converter execution
and complete static closure remain gated off; canonical assignment remains
disabled. Missing sources alone do not establish a QM-fit requirement.

The next scientific dependency remains exact authoritative Si/O source rows,
native Phi/AAT/equilibrium-dependency conventions and author-primary evidence.
J68's model-span incompatibility remains independent. Chemical transfer and
scientific validation are unestablished. The 28 confirmations remain
user-attested and original author artifacts remain `not_yet_supplied`.

**accepted_coefficients=0; automatic_assignment=False;
production_validated=False; simulation_readiness="not_established".**
Passing corrective software checks establishes neither compliance of the
historical J72 allocation nor physical/production validity.
