# J84: audited parameter bundle adapters and organic coverage

Decision: `audited_bundle_ingestion_ready_runtime_blocked`.
Base: `8a1c92c992951edbf34dd73c776df399443e2f06`.
Branch: `codex/phase-4j84-audited-parameter-bundle-adapters`.

J84 connects the unchanged J83 coverage registry to explicitly supplied,
hash-pinned local parameter bundles. It supplies adapters for OPLS-AA, GAFF,
GAFF2, PCFF and COMPASS, preserving each family's identity. The shared mapper
constructs J82 `TypingRule`, `TermForm`, `Convention`, `ParameterRecord` and
`EvidenceRecord` objects. J82 performs exact environment lookup, typing,
charge conservation, convention and dependency validation; J83 reports motif
and topology coverage. Neither core is modified.

No authorized scientific bundle was supplied for this implementation. The
recorded benchmark therefore reports eight organic motifs against five empty
contracts: 24 Class-I and 16 Class-II candidate reports. Missing scientific
terms remain explicit. The working ingestion path is tested with manufactured
source tables and audits; these controls establish software behavior only.
They are never scientific validation or accepted coefficients.

## Source format and trust boundary

The supported exchange format is `island.audited-parameter-bundle.v1`, with
adapter version `island.j84.canonical-json.v1`. Use `schema` below to inspect
the exact fields and `bundle_adapters.schema.validate_manifest` for strict
semantic validation. The manifest requires a bundle/profile identifier,
upstream path or URI, family and release, complete local file list, SHA-256 and
size of every file, units, complete functional forms, mixing rule, separate LJ
and Coulomb 1–4 policy, charge origin, topology scope, provenance/license
metadata, authorization and ingestion timestamp. `source_sha256` is a canonical
aggregate over sorted `(relative path, SHA-256, size)` triples, distinct from
each original-file checksum. `identity` hashes the manifest excluding its own
identity field. It includes the recorded initial ingestion timestamp.

Only files explicitly listed in that manifest are read. Paths must be relative
and remain inside the supplied bundle root; source symlinks and traversal are
rejected. Relocation changes the physical root, preserving relative names,
hashes, declared upstream identity and semantic bundle/profile identity. A new
capture timestamp is event metadata, not a new force-field release. A URI is
provenance metadata; this adapter never fetches it or bypasses access controls.

A `record_source` is a UTF-8 JSON file with fields `schema`, `family`, `version`
and `records`. Its schema is `island.parameter-source-records.v1`. Each row has
`record_id`, `kind`, exact `location`, `origin="directly_reported"` and `value`.
Row kinds are:

- `atom_type`: explicit rule ID, atom type, SMARTS, priority and overrides.
- `form`: the complete J82 form, coefficient ordering, units and equilibrium fields.
- `convention`: explicit name and value.
- `parameter`: exact ordered types, environment hash, form, original numeric
  strings, applicability, dependency keys and equilibrium references.
- `graph_scope`: exact graph identity reviewed for that bundle.

All rows retain their original JSON source bytes, row location and JSON pointer.
Internal evidence IDs and convention hashes are derived metadata; numerical
values, units and equations are passed through unchanged. No unit conversion,
9–6 to 12–6 coercion, inferred atom typing, wildcard substitution or zero-fill
is performed. Missing values fail validation. A documented non-applicability
record may omit an improper/cross term only under J82's existing evidence rules.

An original `.frc`, Amber parameter file, export, PDF or license file can be
listed as `primary_artifact`. Such attachments are hashed read-only, included
in the evidence graph and independently reviewed; this phase does not pretend
to parse undocumented native formats. An authorized, explicit record table
must supply the mapping and exact source locations. This is an evidence adapter,
not a scientific converter. Raw vendor files cannot trigger executables or
automatically supply parameter values.

A separate audit file and independently supplied SHA-256 are mandatory before
any profile is exposed to coverage. Its schema is
`island.parameter-bundle-audit.v1`. It binds the exact manifest identity,
reviewer, review time, authority basis and `policy="independent_source_review"`.
Every file requires a source review binding its path/hash, provenance and
license verification, plus the exact mapped claim hashes for every record.
The `proposals` function produces untrusted claims for review; it never approves
them. The production adapter provides no audit generation or approval function.

The caller must obtain the audit and its pin from a trusted independent review.
A locally fabricated audit with newly self-selected pins is not evidence of
independent review; software hashes cannot establish an auditor's real-world
identity or authorship. Publisher metadata, user attestations, a local copy or
a `verified_stock` label alone do not pass this gate. PolyJarvis patches,
analog origins, inferred/reconstructed rows, unknown licenses and unqualified
evidence statuses remain blocked. Original source authenticity, authorized use,
the exact mapping and conventions must be reviewed by the source auditor.

Do not publish user-supplied source assets or captured original table bytes
without their redistribution permission. The committed J84 evidence contains
no user-supplied or vendor parameter assets.

## Completeness and topology behavior

Class-I reports distinguish typing, bond, angle, torsion, improper, nonbonded,
charge, topology, mixing and 1–4 policies. Class-II additionally requires all
nine J82 cross-term families, their ordering, primary-term dependencies and
equilibrium reference values. Improper/out-of-plane equations must be supplied
explicitly. A fully audited bundle can still be incomplete for a particular
graph; ingestion is not a completeness certificate. Partial or absent bundles
keep corresponding profiles blocked, with source and missing-term diagnostics.

The benchmark includes aliphatic hydrocarbon, aromatic carbon, ether, ester,
carbonyl, amide, nitrile and common halogenated carbon. Each final explicit-H
graph is built before typing. Repeat-unit cuts, end caps, copolymer junctions,
branches and crosslinks can change graph scope, parameter environments and
site order. A changed graph invalidates the typed snapshot and requires retyping
and separately evidenced applicability. An isomorphic graph reached through a
different cut does not silently inherit the original site mapping or audit.
Equal-priority conflicting atom types remain ambiguous. Explicit priorities
and profile-ID ranking provide deterministic decisions, without guessing a type.

J83's `accepted` means audited exact-graph coverage only. Synthetic complete
controls remain `provisional`. Runtime registration, production assignment,
automatic assignment and accepted coefficients remain disabled in every case.

## Fenced commands

Run from the checkout with the supported environment. Every fence output and
new intake directory must be unused. Fence failure stops before source reading.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-schema-unique \
  -- -m scripts.phase_4j84_bundles schema

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-ingest-unique \
  -- -m scripts.phase_4j84_bundles ingest --directory /tmp/j84-intake-unique \
  --manifest /authorized/bundle/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/review.json --audit-sha256 INDEPENDENT_AUDIT_SHA256
```

The intake prints `snapshot_sha256`. Retain that external pin for supplied-source
replay; the built-in pin applies only to the committed no-source benchmark.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-custom-replay-unique \
  -- -m scripts.phase_4j84_bundles replay --directory /tmp/j84-intake-unique \
  --snapshot-sha256 CAPTURED_INPUT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-replay-unique \
  -- -m scripts.phase_4j84_bundles replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-validate-unique \
  -- -m scripts.validate_phase_4j84_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j84-tests-unique \
  -- -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j84.py
```

Offline replay reconstructs every report from the externally pinned capture,
including byte hashes, audit bindings and exact source records. It does not
re-establish the continued existence of original binary attachments; re-ingest
the supplied directory to recheck those original bytes. Relocated re-ingestion
is tested independently. J84's strict sealed receipt covers the committed
no-source benchmark and the actual fenced terminal checks. A caller's custom
intake is an auditable candidate record, not that terminal software receipt.

All J83-base tracked files, including J81/J82/J83 receipts and replay code,
remain byte-identical. Main and unrelated files remain untouched. J72's four
unintended QM attempts, 46–48 MD steps and C1 quarantine remain preserved.
J84 QM, MD, converter, physical fitting, parameter generation, energy/force
evaluation, reserved validation and production executions are zero. Explicit
scope skips are never reported as scientific passes. Readiness remains
`automatic_assignment_enabled=false`, `production_validated=false`,
`simulation_readiness="not_established"`; accepted coefficients remain zero.
