# J89: first real authorized OPLS-AA qualification gate

Base: `4ed7377a66c02362d69da0840b7701ba2b6a797f`.
Branch: `codex/phase-4j89-first-real-authorized-opls-qualification`.
Recorded decision: `source_bundle_required`.

No real authorized source bundle, explicit request or independently trusted audit
was supplied. The terminal result has zero source files/rows inspected, zero
qualified graphs and zero accepted coefficients. Source, manifest and audit
hashes are absent, rather than invented. The real-source integration test is
explicitly skipped without these inputs. Manufactured controls exercise software
contracts only; they do not establish scientific coverage or authorize a profile.

J89 reuses J81 provenance lineage, J84 pinned source capture and independent
claim auditing, J85 Class-I coverage, J86 exact-scope guards, J87 qualification
and J88 structured-form validation, privacy projection and external-byte replay.
It leaves every historical file unchanged and introduces no other force-field
family, ranking, native registry change or parameter table.

## Required input gate and scope

Supply all six explicit values: paths and external SHA-256 pins for the request,
manifest and independent audit. If a path or pin is absent/invalid, J89 emits a
source-required receipt before opening supplied external files. It also checks
that all three JSON files exist before opening them. Pinned request/fixture and
audit JSON may then be inspected; parameter files remain unopened until exact
scope and independent reviewer/provenance metadata pass. No unpinned source is
opened, no source directory is discovered, and no source is fetched.

Requests reuse `island.j87.opls-qualification-request.v1`. J88 validates the
entire structured `expected_forms` schema before consumers access fields.
J89 requires all eight immutable J87 fixture bytes and final graph identities:
hydrocarbon, aromatic carbon, ether, ester, carbonyl, amide, nitrile and the
halogenated motif. Subsets, altered end groups, changed graphs and new fixtures
cannot qualify. Relative paths may change on relocation; hashes and graph IDs
must remain identical. Graphs are finalized before typing. No graph qualification
implies chemical-family transfer or universal organic-polymer coverage.

The manifest uses the unchanged J84 audited canonical JSON format and must pin
every record source and opaque primary attachment, source role, family/version,
license/authorization and provenance. Vendor originals may be opaque attachments;
only explicitly audited canonical records supply terms. No converter is run.
The independent audit pins the manifest identity, every source file and every
mapped claim, including exact graph-scope claims. Its reviewer, authority basis
and independent-review policy must be explicit; unresolved placeholder identities
block source inspection. Caller-supplied audit pins must come from an
independently trusted review channel. A hash, JSON policy string or self-selected
reviewer name cannot independently authenticate authorship. Production code does
not generate or approve audits.

Qualification reopens and hashes all external files, including opaque originals,
before parsing and again during replay. It compares the independently pinned
expected source version, all six Class-I forms and all conventions: energy,
length and angular units; atom typing and scope; charge origin/neutrality;
sigma/epsilon definitions; mixing; exclusions; long-range/cutoff/tail policy;
separate LJ/Coulomb 1–4 scaling; bond/angle equations; torsion phase,
normalization, ordering and equilibrium references; improper definitions.
Missing forms/terms, ambiguity, unsupported scope, convention mismatch and
source/audit defects block qualification. No unit reinterpretation, conversion,
9–6/12–6 coercion, zero substitution, analog transfer, hidden patch, PolyJarvis
numeric import or fallback exists.

The only decisions are `source_bundle_required`, `parameter_bundle_incomplete`
and `family_scoped_qualification_ready_runtime_blocked`. The final state requires
all eight graphs to pass every source, audit, identity, typing, convention and
coverage gate. Parser or build success is insufficient. `accepted_coefficients`
remains zero even on that path.

## Metadata outputs and replay

`docs/evidence/phase_4j89/` contains the input gate, source manifest metadata,
source/claim hashes and locations, audit references, per-graph coverage,
convention comparisons, missing terms/blockers, sealed receipts and verification.
Coverage distinguishes typing, bonds, angles, torsions, impropers, charges,
nonbonded, 1–4 and topology. Class-II cross terms remain explicitly inapplicable
to this Class-I contract; they are never counted as supplied. Audit pointers bind
record locations to claims in the pinned independent audit.

Raw licensed tables, coefficient arrays, captured table JSON, audit attachments
and proprietary originals remain external. J89 persists only metadata, hashes,
source pointers, claims, blockers and receipts. It copies no source assets into
Git. Replay of an external-source intake requires all external paths/pins and
the input-snapshot pin again. It never substitutes a saved metadata report for
missing originals. A no-source receipt can replay only its empty obligation
report; this produces no scientific coverage. Every report is reconstructed and
compared, so re-signing a changed result does not bypass the gates.
The strict validator also accepts the same six external input arguments and
`--snapshot-sha256` for a sealed external-source result; it cannot validate one
using metadata alone. Default receipt validation uses the no-source result.

Use fresh output and fence directories. The inherited J72-C1 fence must install
before any source IO; installation failure stops the command.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-emit-fence-unique -- \
  -m scripts.phase_4j89_opls emit --directory /tmp/j89-evidence-unique

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-ingest-fence-unique -- \
  -m scripts.phase_4j89_opls ingest --directory /external/j89-intake-unique \
  --request /authorized/request.json --request-sha256 REQUEST_SHA256 \
  --manifest /authorized/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/audit.json --audit-sha256 AUDIT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-external-replay-unique -- \
  -m scripts.phase_4j89_opls replay --directory /external/j89-intake-unique \
  --snapshot-sha256 INPUT_PINS_SHA256 \
  --request /relocated/request.json --request-sha256 REQUEST_SHA256 \
  --manifest /relocated/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/audit.json --audit-sha256 AUDIT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-replay-unique -- \
  -m scripts.phase_4j89_opls replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-validation-unique -- \
  -m scripts.validate_phase_4j89_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j89-tests-unique -- \
  -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j89.py
```

The optional real-source integration test accepts only explicit environment
variables `ISLAND_J89_REAL_INPUTS` and `ISLAND_J89_REAL_INPUTS_SHA256`, identifying
a separately pinned JSON descriptor with the six CLI input keys (`request`,
`request_sha256`, `manifest`, `manifest_sha256`, `audit`, `audit_sha256`). It
reuses the same ingestion path and performs no scientific execution. Without
both values it records a clear skip, not a scientific pass. Synthetic controls
are labeled separately. Relocated controls move request, sources, attachments,
audit and metadata, remove the old directory and rehash originals. Another check
relocates code, metadata and the immutable graph fixtures under the inherited
security runtime. Strict terminal validation binds the tested commit, explicit
full regression inventory, actual check results, fence events and every
historical base-tracked byte. Historical J81–J88 replay/receipts remain unchanged.

QM, MD, converter, physical fitting, parameter generation, scientific energy/force
evaluation, reserved validation and production executions are zero.
`runtime_registration=false`, `automatic_assignment_enabled=false`,
`accepted_coefficients=0`, `production_validated=false`,
`simulation_readiness="not_established"`,
`chemical_family_transfer="not_established"`. Main, historical J72 violations,
C1 quarantine records and unrelated untracked assets remain unchanged.
