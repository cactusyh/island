# J85.1/J86: scoped qualification and partial-bundle safety

Base: `19fa2e76cd822987618cf03e42fca9dd9d884975`.
Branch: `codex/phase-4j86-scoped-bundle-qualification`.
Recorded decision: `source_bundle_required`.

No authorized scientific bundle was supplied. All five declared profiles remain
source-blocked. No source discovery or fetching occurs. The first authorized
profile path is implemented, but no scientific profile is qualified in this
receipt. Manufactured test bundles and audits exercise software contracts only.

## Corrected workflow

The additive `bundle_adapters.scoped_qualification` module reuses J85's
unchanged individual `qualify()` and J84's unchanged `capture()`/`adapt()`.
J81–J85 implementation, tests, source pins, evidence and receipts stay
byte-identical. Their historic aggregate/ranking APIs are retained for replay;
new qualification work must use J86's scoped API and CLI. J86 never calls the
historical cross-family ranking or aggregate logic.

Each `QualificationRequest` declares a canonical family, one exact profile ID,
and a nonempty unique set of final-graph identities. Every required slot must
have an independently inspectable per-profile result. Omitted, extra, duplicate,
wrong-family, wrong-profile and wrong-graph results are rejected. No successful
profile can hide a missing result. Missing sources retain the underlying empty
contract ID separately from the requested profile ID.

Aggregate decisions are deterministic:

- `source_bundle_required`: none of the required profiles has an authorized source.
- `partial_bundle_set`: some required profiles have authorized sources while
  others are missing or unauthorized. Incomplete authorized profiles remain
  listed alongside missing profiles.
- `parameter_bundle_incomplete`: all required sources are authorized, but at
  least one profile fails convention or exact-graph completeness.
- `family_scoped_qualification_ready_runtime_blocked`: every required profile
  passes for every graph in its declared scope.

Per-profile J85 decisions are unchanged, including
`authorized_bundle_qualified_runtime_blocked`. Aggregate readiness covers only
the explicitly declared set. A one-profile declaration cannot imply coverage
of all families. Chemical-family transfer remains `not_established`, including
when exact graph coverage passes. Class-I and Class-II reports remain separate.

`rank_qualifications(results, *, requested_family=...)` requires a family.
`opls` and `OPLS-AA` mean OPLS-AA only. GAFF, GAFF2, PCFF and COMPASS each retain
their own identity. Mixed-family input raises `cross_family_ranking_forbidden`;
there is no filtering, switching or cross-family winner. Different graph scopes
also cannot be ranked together. Reports group supplied candidates by family and exact
graph scope; missing-source declarations remain visible in the aggregate and
are never treated as rankable source candidates. Within one such group, qualification, normalized coverage,
blocker count, profile ID and bundle identity determine the stable order.

The default declaration contains one exact hydrocarbon graph and five separate
empty profile contracts, solely to expose every missing family. For the first
authorized qualification, provide one explicitly named profile and the exact
graph fixtures it covers. More extensive organic coverage requires separately
audited graphs; the default fixture does not establish transfer.

## Authorized intake and trust boundary

Sources must be explicitly supplied with external manifest and independent
audit pins. The existing J84 contract requires original source hashes, file
sizes, family/version, authorization, license/provenance, units, functional
forms, mixing, exclusions, 1–4 policy, charges, atom-type rules and exact
source locations. J85 additionally requires the independently audited
`energy_convention`, all required graph terms, exact atom-type references,
and Class-II dependency closure. Missing terms remain missing. No zero-fill,
PolyJarvis import, analog substitution, family mixing, tolerance change,
converter, scientific evaluator or numeric invention is introduced.

Hashes establish identity, not authorship. Caller-supplied independent review
and pins must come from a trusted channel; a fabricated self-signed audit is
not source authorization. Production code does not generate audits. Do not
publish captured licensed tables without redistribution permission.

A custom declaration has the same schema as
`tests/fixtures/phase_4j86/qualification-set.json`. Set `required_profiles` to
one entry with `family`, the supplied manifest's exact `profile_id`, and
`graphs` listing fixture IDs. Supply only those graph fixtures in `motifs`,
keep `bundle_captures=[]`, and independently pin the declaration SHA-256.
Ingestion receives the source separately. No bundle can silently change the
declared set. Multiple source captures can be processed through the pure
reconstruction API, without overlaying any profile.

## Offline commands

Use unused fence/output paths. The inherited C1 fence is installed before
source IO. Failure to install it stops processing. Default replay requires the
pinned committed no-source snapshot; custom replay requires its external pin.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-emit-unique -- \
  -m scripts.phase_4j86_qualification emit --directory /tmp/j86-evidence-unique

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-ingest-fence-unique -- \
  -m scripts.phase_4j86_qualification ingest --directory /tmp/j86-intake-unique \
  --qualification-set /authorized/exact-graphs.json \
  --qualification-set-sha256 DECLARATION_SHA256 \
  --manifest /authorized/profile/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/review.json --audit-sha256 INDEPENDENT_AUDIT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-custom-replay-unique -- \
  -m scripts.phase_4j86_qualification replay --directory /tmp/j86-intake-unique \
  --snapshot-sha256 CAPTURED_INPUT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-replay-unique -- \
  -m scripts.phase_4j86_qualification replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-validator-unique -- \
  -m scripts.validate_phase_4j86_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j86-tests-unique -- \
  -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j86.py
```

Every valid intake receives a snapshot-bound qualification receipt. Replay
reconstructs every aggregate, individual result, missing list, ranking, source
hash and coverage dimension; re-signing a changed report cannot bypass it.
Schema/hash failures exit with structured diagnostics before writing evidence.
The terminal strict receipt also binds implementation blobs, actual test logs,
fence events and every base-tracked file. Relocated replay requires code and
captured evidence, not external assets; a separate control re-ingests relocated
source files and rehashes them. Neither establishes scientific validity.

## Execution and readiness

All J86 QM attempts, MD steps, converter, physical-fitting,
parameter-generation, scientific energy/force, reserved-validation and
production executions are zero. `accepted_coefficients=0`,
`automatic_assignment_enabled=false`, `runtime_registration=false`,
`production_validated=false`, `simulation_readiness="not_established"`, even
for a complete declared set. No runtime selection/registry API is modified.
J72's historical four unintended QM attempts and verified 46–48 MD steps,
C1 quarantine and all historical seals remain unchanged. Main is untouched.
Explicit scientific test skips remain skips, never physical validation passes.
