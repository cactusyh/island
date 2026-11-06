# J85: authorized bundle qualification and organic coverage

Base: `f5b9abcf459f541134d4f1162fb3af12f454a1f5`.
Branch: `codex/phase-4j85-authorized-bundle-qualification`.
Recorded decision: `source_bundle_required`.

No authorized scientific source bundle was supplied. J85 does not search the
checkout or untracked source directories for substitutes. All five profiles
remain blocked, accepted coefficients remain zero, and no evidence-backed
organic coverage is claimed. The implemented workflow can qualify explicitly
supplied bundles after independent source review. Its manufactured controls
test software behavior only.

## Qualification and trust

J85 reuses the unchanged J84 `capture` and `adapt` implementations. OPLS-AA,
GAFF, GAFF2, PCFF and COMPASS retain separate identities; no numeric overlay,
converter, fallback, parameter generation or unit conversion is introduced.
The canonical exchange schemas and independent audit contract are documented
in [J84](phase_4j84.md). Native vendor files can be hash-bound attachments;
audited canonical tables supply their explicit mappings. J85 does not claim
to parse undocumented Amber or FRC formats automatically.

The source bundle requires family/version, a complete declared file list,
each original SHA-256 and byte size, aggregate source hash, unit and
functional-form declarations, mixing rule, separate LJ/Coulomb 1–4 policies,
charge origin, license, provenance, authorization and timestamp. Original
canonical table bytes and exact row locations are preserved. Every file hash,
including opaque primary attachments, remains in the evidence graph.
Only declared files are opened; no remote source is fetched.

The caller must independently pin both manifest and audit bytes. The audit
binds source authority, licensing and every mapped claim. Hashes alone do not
prove authenticity, authorship or scientific suitability. A fabricated audit
with newly self-selected pins is not independent review. Production code
does not generate or approve audits. PolyJarvis patches, analog origins,
local patches and unqualified evidence statuses cannot qualify.

The additive `bundle_adapters.qualification` module validates all atom-type
references, parameter forms/units/order and dependency/reference keys, even
for unused source rows. It requires an additional independently audited
`convention` record named `energy_convention`. That record must describe the
energy expression and coefficient normalization, including explicit
prefactors. Existing J84 bundles without it remain valid J84 inputs but
cannot qualify under J85. No 9–6/12–6 coercion is permitted. Unit/form
compatibility is exact comparison with audited declarations; the offline
gate does not evaluate energies or establish physical dimensional correctness
of an arbitrary equation string.

Final graphs are built before typing. J82 performs exact environment lookup,
typing ambiguity detection, charge conservation, source-claim auditing,
convention checks and Class-II dependency closure. J83 supplies the motif
catalog and graph-specific coverage. J85 adds separate bond, angle, torsion,
improper, charge, nonbonded, atom typing, 1–4 and topology dimensions, plus
Class-II cross terms. The 1–4 dimension checks both LJ and Coulomb policy
evidence; exclusions/mixing and all other conventions remain independent
blocking requirements. Class-I cross terms are explicitly not required by
that contract. Missing rows are never zero-filled. Documented non-applicability
requires J82's exact audited omission records. A section header alone is
insufficient: every required final-graph instance must resolve.

Graph inputs and ranking inputs are materialized once before validation;
empty graph collections are rejected, including one-shot iterators.

`authorized_bundle_qualified_runtime_blocked` means audited coverage for the
declared exact graphs only. It does not certify chemical-family transfer,
physical accuracy, complete native PCFF, or production readiness. Unknown
motifs, missing terminal terms, unaudited graph scopes, equal-priority typing
conflicts, missing torsions/cross terms and missing energy conventions produce
`parameter_bundle_incomplete` after source authorization. Absent, unauthorized,
unaudited, unqualified-origin or synthetic-only bundles produce
`source_bundle_required`. The explicit `source_authorized` flag remains separate
from convention and graph completeness. Independent candidates are ranked by
qualification, normalized coverage, blockers and stable family/profile/bundle
identities; equal-priority conflicting atom types remain blocked.

The benchmark contains eight organic graph fixtures: aliphatic hydrocarbon,
aromatic carbon, ether, ester, carbonyl, amide, nitrile and halogenated carbon.
Ether and ester fixtures declare methyl caps, so they exercise their intended
linkages without silently including terminal alcohol/acid chemistry. The
original H-capped versions still block when the catalog cannot classify those
terminal environments. It produces 24 Class-I and 16 Class-II reports. Each supplied family is
qualified over all declared benchmark graphs; other families retain explicit
missing-source reports. Partial graph coverage remains visible even if the
overall supplied bundle is incomplete. Repeat-unit cuts, caps, copolymer
junctions, branches and crosslinks require retyping and new exact graph-scope
evidence. An isomorphic repeat cut cannot inherit the old site mapping.

## Fenced replay and intake

Every output/fence directory below must be unused. Fence installation failure
stops before source IO. The pure qualification function performs no IO or
execution. Do not publish supplied table snapshots or licensed attachments
without redistribution permission. Committed J85 evidence contains none.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-emit-unique -- \
  -m scripts.phase_4j85_qualification emit --directory /tmp/j85-evidence-unique

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-intake-fence-unique -- \
  -m scripts.phase_4j85_qualification ingest --directory /tmp/j85-intake-unique \
  --manifest /authorized/bundle/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/review.json --audit-sha256 INDEPENDENT_AUDIT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-custom-replay-unique -- \
  -m scripts.phase_4j85_qualification replay --directory /tmp/j85-intake-unique \
  --snapshot-sha256 CAPTURED_INPUT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-replay-unique -- \
  -m scripts.phase_4j85_qualification replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-validation-unique -- \
  -m scripts.validate_phase_4j85_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j85-tests-unique -- \
  -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j85.py
```

Intake prints the external snapshot pin. Each structurally valid intake,
including a source-blocked or incomplete bundle, receives
`qualification-receipt.json`; replay reconstructs
its reports and rejects re-signed semantic mutation. A byte/hash/schema
failure rejects intake with exit 2 and a structured diagnostic, before
creating any verified snapshot. The separately sealed terminal receipt binds
the committed no-source benchmark, tested implementation, actual test logs
and historical preservation. It does not certify a caller's external audit.

Relocated replay needs only captured tables/evidence and code. It does not
prove binary originals are still present; relocated source re-ingestion tests
rehash those originals separately. No parser/build pass is counted as a
scientific pass. Explicit scientific scope skips remain skips.

J85 QM, MD, converters, physical fitting, parameter generation, energy/force
evaluation, reserved validation and production execution are all zero.
Runtime registration and automatic assignment remain false, accepted
coefficients remain zero even for qualifying source candidates,
`production_validated=false`, `simulation_readiness="not_established"`.
All base-tracked files, including J81–J84 receipts and implementation, are
byte-preserved. J72's four unintended QM attempts, 46–48 MD steps and C1
quarantine remain unchanged. Main and unrelated files remain untouched.
