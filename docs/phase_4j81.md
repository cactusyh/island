# Phase 4J81 — user-supplied PDMS source ingestion

The delivered decision is **source_ingestion_blocked**. No authoritative bundle
was explicitly supplied. J81 does not search local literature directories,
retry the publisher request, fetch a paper or infer a missing source. There are
zero supplied artifacts, zero new parsed/verified literature rows and twenty
historical candidates still blocked. Convention validation is **not evaluated**.
There is no new QM, fitting, MD, converter, evaluator or reserved-validation run.

Branch: `codex/phase-4j81-user-supplied-pdms-source-ingestion`.
The verified base is `637957fc8167e2f14f983cc364512da3fede9c52`.
Fetched `origin/main` was `93128ac53a09390d6801208466abc55afa67e7f0`.
Main is untouched. Implementation is additive in an explicit subpackage; no
historical file, registry, source hash, value, readiness flag or receipt changes.

## Explicit source bundle

`source-bundle.schema.json` is the machine-readable local manifest schema.
`python -m scripts.phase_4j81_ingest schema` emits the same definition. Its strict
offline validator rejects duplicate JSON keys, unknown fields, missing metadata,
duplicate source hashes/IDs, candidate IDs and conflicting duplicate locations.
Manifest bytes require an independently supplied SHA-256 pin. Paths resolve
relative to the manifest directory; explicit absolute local paths also work.
URLs, traversal, symlinks and implicit discovery are rejected. No HTTP or
publisher-access code exists in this phase.

Each artifact records its path, accepted source type, SHA-256, byte size,
bibliography (title/authors/year and identifier or explicit null), version,
acquisition method, local-use authorization status, separate redistribution
status, declared source authority/family, extraction tool/version/timestamp
and parser version. Timestamp is a timezone-qualified supplied provenance fact;
replay does not replace it with the current time. Hashes establish byte identity,
not authenticity, licensing or transferability. User attestation is metadata.
No original PDF, FRC, export or archive is committed.

Accepted inputs are primary PDFs, original PCFF FRC files, authorized Materials
Studio text/FRC exports, precisely cited parameter tables and authorized ZIP/TAR
archives. The parser preserves raw UTF-8 records, including CRLF, spacing and
comments. It recognizes named FRC interaction sections; proprietary/binary or
unknown export formats remain unparsed. It never runs Materials Studio, a source
script, shell, converter or parameter generator.

PDF text parsing uses optional **pypdf 6.1.1**. Missing/wrong parser versions,
encryption, invalid pages, wrong extracted-page-text SHA-256, absent/ambiguous
snippets or missing table-row-column references remain blocked. Figure/appendix
locations are explicit alternatives. OCR status/confidence must be recorded;
this implementation does **not** run OCR. OCR-dependent rows remain unresolved
even if a supplied confidence is high. A text snippet match does not independently
prove its declared printed cell: the separate authority audit must verify that.

Archives are read in memory without filesystem extraction. Every member has
name/size/hash metadata; cited members need their own SHA-256. Duplicate names,
links, special files, unsafe paths, encryption and expansion-limit violations
are rejected. Code members are hashed, never imported or executed. Nested archives
and unrecognized member formats are not silently interpreted. Limits are 64 MiB
per input/member, 128 MiB expanded per archive, 4096 members/candidates, 256 source
artifacts and short cited excerpts. These are parser limits, not scientific gates.

## Parsing, verification and normalization

Artifacts and candidate rows retain distinct statuses:

- `source_missing`, `source_present_unparsed`, `source_parsed`;
- `literature_candidate_verified`, `literature_candidate_convention_blocked`;
- `exact_native_row`, `inferred_or_analog_value_blocked`, `end_group_scope_blocked`.

A supplied file can be parsed without acquiring source authority. Parsed claims
receive structural field-mismatch diagnostics even before authority is verified;
these diagnostics never count as authoritative convention validation. Each candidate
preserves source ID/SHA, exact PDF or line/section/record locator, raw excerpt
and checksum, ordered source types, separate elements, chemical environment,
family, original numeric strings and explicit per-value units, functional form, mixing/exclusions,
charge/torsion convention, cap/internal scope and reporting mode. Decimal
interpretation is separate from the original strings. No unit or functional-form
conversion is performed; converted values remain null. There is no 9–6-to-12–6
coercion or element-sequence-to-atom-type promotion.

Verification additionally requires a **separately reviewed and externally
hash-pinned authority audit**, supplied with `--audit` and `--audit-sha256`.
The bundle cannot set its own verification status. The audit maps each source ID
to the exact artifact metadata, captured-extraction identity, independently
checked rows, and source-bound convention assertions. `bridge()` exposes the
canonical J80-compatible audit shape. `capture_identity` is the deterministic
`digest(capture(...))`, not an unverified claim made by the file itself.
The audit's evidence kind is `independently_reviewed_primary_artifact`.

A reviewer must actually inspect the original author/publisher/release artifact,
check authority and permissions and verify every location and interpretation
before supplying that audit through the trusted administrative boundary. The
software verifies integrity and agreement with this externally trusted input;
it cannot prove authorship or honesty from a checksum. Tests manufacture such
audits only as software controls. None is installed for the delivered result.
Self-attested permission, publisher metadata, secondary summaries, reconstructed
or inferred values cannot qualify a row. GAFF, OPLS, COMPASS, CVFF, CFF91, analog
and external PDMS families cannot become native PCFF through this path.

The unchanged J80 gate implements the twelve J79/J80 convention groups: 9–6/12–6,
energy/length/angular units, sigma/epsilon/minimum meaning, mixing, exclusions and
1–4 scaling, charges and independent neutrality arithmetic, bond/angle
polynomials, torsion sign/phase/periodicity/normalization, improper/out-of-plane,
cross-term equations/order/dependencies, exact internal/cap typing and finite/
periodic electrostatics/cutoff/tail scope. Each row is also compared with that
contract. Mismatches are explicit; omitted conventions never receive defaults.
A native claim additionally needs exact source types, release/row identity,
environment, equation, phase/units and policies. An exact row never activates
runtime assignment or establishes complete native PCFF support.

A report may reach `source_parsed_ready_for_qm_authorization` only with all
supplied rows verified, complete compatible conventions and coverage of all
historical candidate IDs. That status is an intake review prerequisite only:
`qm_authorization_granted=false` in every outcome. It does not resolve the
129 J78 extension targets automatically. No pilot launcher, fitter, evaluator,
parameter acceptance, registry insertion or production export is implemented.

## Receipts and replay

The evidence includes the empty source manifest, captured inputs, provenance,
candidate normalization, convention comparison, exact remaining J78 target IDs,
result, sealed receipt/manifest, historical inventory and validation records.
`unresolved-terms.json` preserves all 129 target identities with 87 cross-term
and 69 end-group flags (48 overlap). All 365 scalar and 11 policy requirements
remain historical unresolved requirements. J80's receipt and blocked result
remain unchanged, including its 378 native unresolved instances.

Every ingestion emits a deterministic `intake-receipt.json`, including blocked
and parsed/verified status paths. The terminal release adds its separately
sealed `receipt.json` and validation records. These paths are exercised
with synthetic controls. It binds content hashes, implementation and historical
bytes. Offline replay reconstructs all reports rather than trusting re-signed
status fields. The delivered no-source snapshot has a code-pinned SHA-256.
For another bundle, keep the `snapshot_sha256` returned by ingestion outside the
output directory and provide it to replay/validation. Re-signing a snapshot or
receipt cannot replace that independent pin. Replay authenticates captured
excerpts and their pinned review; it does not pretend to re-read unavailable
copyrighted originals. Re-ingest with the original manifest pin to recheck live
source bytes. Outputs are append-only; use a fresh directory for every attempt.

All commands below use the inherited C1 process fence and J73 evidence-only
policy. If installation fails, source access stops. No scientific backend is
allowed. Optional PDF dependency installation belongs to environment setup;
normal ingestion/replay never installs packages or accesses the network.

```bash
export PYTHONPATH=.:src
J81_PY=../island-validation/phase4g1-env/bin/python

# Reconstruct the delivered blocked result, entirely offline.
"$J81_PY" -m scripts.run_phase_4j72_c1 --output /tmp/j81-replay -- \
  -m scripts.phase_4j81_ingest replay
"$J81_PY" -m scripts.run_phase_4j72_c1 --output /tmp/j81-validate -- \
  -m scripts.validate_phase_4j81_receipt

# Explicit future local intake; never automatically searches for a source.
"$J81_PY" -m scripts.run_phase_4j72_c1 --output /tmp/j81-intake-fence -- \
  -m scripts.phase_4j81_ingest ingest --directory /tmp/new-pdms-intake \
  --manifest /explicit/source-bundle.json --manifest-sha256 <reviewed-pin>
# Add --audit /explicit/independent-review.json --audit-sha256 <independent-pin>
# only after an actual independent original-artifact review.

"$J81_PY" -m scripts.run_phase_4j72_c1 --output /tmp/j81-tests -- \
  -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j81.py
```

A well-formed blocked ingestion/replay returns 0 because its software contract
succeeded; inspect `decision` for scientific status. Invalid input returns 2 and
an additive `blocked-receipt.json` when the fence and output path permit it.
No authorization is implied by a successful parser command.

## Verification and scientific boundary

Focused tests cover missing/unreadable sources, PDFs, exact lines/sections,
raw formatting and numerical transcription, SHA mutations, unit/equation/phase
mismatches, all twelve convention groups, family and atom-type separation,
end-group/cross-term blocking, independent audit/pin mutations, duplicate/conflict
identities, local-only paths, archive safety, immutable receipts, relocation and
forbidden synthetic child processes. Real QM/MD/fitting/converter/evaluator/
reserved-validation controls are explicit scope skips, never scientific passes.
J73–J80 evidence regressions, Ruff, compileall, supported dependency checks,
historical receipt validation, offline and relocated fresh-process replay,
diff checks and byte preservation are archived separately from development
attempts. Synthetic source strings and audited-ready controls do not establish
real scientific source verification.

The PDF parser dependency installation first encountered pip's blocked native
OS-discovery command. Setup was completed under the same unchanged fence with
pip's OS metadata reader configured to use files only. Failed installation,
lint and development-test attempts remain recorded with their actual outcomes.
No rejection was bypassed by broadening executable permissions.

Final scientific counts are zero verified literature/native exact rows,
20 historical blocked candidates, 87 unresolved cross terms, 69 unresolved
end-group terms, zero completed QM requests, fitted/accepted coefficients and
scientific executions. Runtime registration, automatic typing/assignment,
production export/MD and analog fallback are disabled. Production validation is
false; simulation readiness is `not_established`. J72's historical four unintended
QM attempts and 46 verified / 48 possible MD steps remain explicitly preserved
with the C1 quarantine; they are never folded into J81's zero execution counts.
