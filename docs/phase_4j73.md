# Phase 4J73 — authoritative PCFF source-intake request

The decision remains **`source_conflict_unresolved`**. J73 creates an evidence
request and an offline validator. It reads retained J69–J72/C1 JSON evidence;
it does not read new papers, re-adjudicate literature, extract or assign new
values, construct a parameter bundle, evaluate a Hamiltonian, or execute a
scientific backend.

The actual base is **73eb64e02ab5d970e02f444e1b41c8183c49034d** on
`codex/phase-4j72-c1-execution-scope-correction`. Origin was fetched; main was
`c118ddc` at that fetch. J73 branches from the exact requested C1 commit as
`codex/phase-4j73-authoritative-pcff-source-intake`. Main is not merged or pushed.
Unrelated untracked files and the external literature directory are preserved.

## Exact requests and remaining dependency

[source-intake.json](evidence/phase_4j73/source-intake.json) contains one
immutable request for each of the **378** unresolved interaction instances:

| Inherited J72 status | Instances | Required evidence |
|---|---:|---|
| Analog only | 294 | Exact same-contract native chemical-environment rows replacing analogy |
| Convention blocked | 60 | Native Phi/AAT equation, phase, units, prefactor, reference order and equilibrium dependencies |
| Source missing | 24 | Original missing Si-centered angle-angle rows and their complete conventions |

The six families remain 78 angle-angle instances and 60 instances in each of
angle-angle-torsion, angle-torsion, bond-bond 1–3, end-bond-torsion and
middle-bond-torsion. The historical 56 symmetry-reduced unknowns do not replace
this complete instance inventory.

Each request binds the original J70 and J72 identities, exact ordered types,
stable sites, full local neighborhoods, graph/order identity, symmetry mapping,
bond/angle/torsion/improper/Class-II context, required coefficient fields and
historical equilibrium dependencies. Historical equations and equilibrium
values are context, explicitly distinguished from missing primary evidence.
No value is inferred or filled in.

The requests also bind shared requirements for separate Coulomb/LJ 1–2, 1–3
and 1–4 exclusion/scaling rules, original directed bond-increment/charge origin,
the native LJ equation and mixing convention. J69's retained isolated research
pair policy does not become native source authority. Improper applicability
must be documented; absence in a declared inventory is not numeric zero
authority.

[evidence-request-report.md](evidence/phase_4j73/evidence-request-report.md)
lists every exact row under its ordered-type/family/status request group. It
specifies which original complete FRC release, author/publisher table,
supplement, typing dictionary, equation definition or source-code function is
needed. Existing alternative citations and page/table/row locations are
retained as unqualified context. Unknown original primary filenames and table
locations are explicitly `not_supplied`; J73 invents none.

The dependency is original authoritative source rows and the matching native
Phi/AAT, equilibrium-dependency, exclusion, charge and LJ conventions. Missing
sources do not by themselves establish a QM-fit requirement. J68's target
incompatibility remains an independent historical result.

## Artifact contract and separate evidence states

[artifact-contract.json](evidence/phase_4j73/artifact-contract.json) requires
the original author/publisher file or exact native FRC release, immutable
filename and SHA-256, citation, version, delivery/provenance chain,
license/redistribution status, and exact page/section/table/equation/source-row
or source-code lines. Every coefficient needs complete units, scalar energy
and analytic force equations, orientation, multiplicity, phase and dependency
metadata, plus exact graph applicability.

A hash verifies bytes; it cannot authenticate an author or establish native
PCFF authority. J73's independently trusted new-primary-artifact registry is
empty. A fabricated file with a valid hash, a re-signed source record or a
user-attested statement cannot register itself. A future original artifact
requires provenance review and a separately reviewed, hash-pinned contract
revision. This phase does not accept new numerical source records.

[status-separation.json](evidence/phase_4j73/status-separation.json) keeps:

- `author_confirmed_reported`: 28 historical user-attested queries, with no
  invented per-instance mapping;
- `author_confirmed_audited`: zero;
- `primary_artifact_supplied`: zero newly trusted exact missing-row artifacts;
- `exact_row_resolved`: zero of 378;
- `convention_resolved`: zero complete exact native instance contracts;
- `executable_bundle_qualified`: false.

The historical literature files are not erased or relabeled by the
missing-row primary-artifact count. No status automatically implies another.
Native PCFF, external literature, analogs and research corrections remain
separate; no averaging, zero-fill, wildcard fallback or cross-family overlay
is accepted.

## No-execution gate and historical integrity

All verification runs through the **unchanged J72-C1 launcher and inherited
kernel/Python execution fence**. J73 adds a verification-only import guard
against OpenMM and scientific modules and rejects named Python scientific
drivers before spawning them. The guard is an explicit pytest plugin before
collection; each J73 CLI also requires it. Native descendants remain governed
by C1. A missing C1 installation stops the CLI with a blocked result before
input access. This is the supported cooperative verification policy, not a
sandbox for hostile Python executing arbitrary machine code or `-c` programs.

Synthetic marker files exercise forbidden executable names and Python fitting,
parameter-generation and reserved-validation drivers. Fake imports allow
read-only `find_spec()` but reject code execution. No real scientific
executable, OpenMM import/integrator, converter or fitting driver is used by
these probes. All live scientific controls are explicit scope skips, never
reported as successful scientific checks.

| Execution namespace | Count/status |
|---|---:|
| J73 authorized new QM / MD | 0 / 0 |
| J73 actual new QM / MD / OpenMM / LAMMPS / AmberTools | all 0 |
| J73 converters / fitting / parameter generation / reserved validation | all 0 |
| HMDS new QM / reserved validation QM | 0 / 0 |
| Historical J72 unintended external SQM attempts | 4 |
| Historical J72 unintended MD steps | verified minimum 46, possible maximum 48 |
| Historical J72 zero-execution scope met | false |

The C1 quarantine manifest, its 192 original members, contamination exclusion
and the J72 violation are hash-bound unchanged inputs. Source records and all
conflict classifications remain untouched, including 16 equal duplicate
equivalences and two unresolved unequal COMPASS assignments. The frozen
**9,380-file base inventory** preserves all J17–J72/C1 evidence and code
byte-for-byte. The explicit J73 implementation inventory does not broaden any
historical receipt glob.

## Offline replay and verification

The validator reconstructs requests from the fixed parent bytes and compares
both schemas and scientific semantics. Re-signing outer identities cannot
hide missing coverage, duplicate identities/conflicts, changed conventions,
source hashes, analog promotion, a profile activation or erased historical
execution counts. No historical adjudication or scientific evaluator is called.

Run the focused software controls:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/island-j73-focused-fresh -- \
  -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j73.py -q
```

Offline replay needs repository scripts and the retained explicit JSON inputs,
with no PDFs, external FRC directory, backend, converter, solver or network:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/island-j73-replay-fresh -- \
  -m scripts.validate_phase_4j73_receipt docs/evidence/phase_4j73 --replay
```

Strict receipt/manifest validation additionally checks the full historical
byte inventory, explicit implementation files, unique terminal test coverage,
verification records and fenced launch ledgers:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/island-j73-strict-fresh -- \
  -m scripts.validate_phase_4j73_receipt docs/evidence/phase_4j73
```

Every launcher output directory must be new so earlier attempts cannot be
overwritten. Replay accepts a relocated repository via `--repo` and a separate
evidence directory. The retained `validation/` records report terminal unique
coverage, all skips, Ruff, compileall, the supported environment pip check,
diff checks, unchanged C1 receipt validation and fresh-process relocation.
No broad historical suite is launched: its live engines, converters and
scientific reconstruction workflows are outside this allocation.

The terminal focused run contains **144 unique test nodes: 135 passed and 9
explicit scientific scope skips**, with no failures or errors. Earlier focused
runs overlap this coverage and add no unique terminal tests. The first run
passed 134 controls with 9 skips; a later run retained one failure in a newly
added historical-source mutation fixture (134 passed, 9 skipped). That fixture
initially omitted the other required parent files, so intake correctly stopped
at a missing input before reaching the intended mutation. Copying all 13 fixed
parent inputs before mutating the graph corrected the fixture. The terminal
run verifies both the re-signed freeze and changed historical source bytes.
An earlier Ruff failure identified one unused test import, which was removed.
The failed runs, original logs and correction history remain in `validation/`.
These corrections did not change historical evidence or scientific outcomes.
The first evidence-packaging attempt also stopped when its completed-run scan
encountered its own active fence. The corrected packager excludes only its
active policy directory and retains that failed attempt and original driver.
A staged diff check exposed trailing whitespace in the original failed pytest
traceback. That log is retained with lossless XZ compression, preserving every
byte; the failed diff check is retained too.

## Unchanged capability

No PCFF profile, parameter bundle, converter result, automatic assignment path
or MD workflow is created. Ordinary `forcefield="pcff"` behavior is unchanged.
All 378 instances remain unresolved. Exact rows alone would still require a
separately authorized conversion, closure and bounded-domain qualification
phase before execution or registration.

**accepted_coefficients=0; automatic_assignment=False;
production_validated=False; simulation_readiness="not_established".**
Chemical transfer and scientific validation remain unestablished; original
author artifacts remain `not_yet_supplied`.
