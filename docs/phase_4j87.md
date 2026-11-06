# J87: first authorized OPLS profile qualification path

Base: `fe4e82de3ff1b0704d600e226e79c66cdccec2be`.
Branch: `codex/phase-4j87-first-authorized-opls-qualification`.
Recorded decision: `source_bundle_required`.

No user-supplied authorized OPLS-AA bundle accompanied this phase. No scientific
profile is qualified, no coefficients are accepted, and the eight declared
benchmark graphs have no evidenced parameter coverage. The implemented path can
qualify one explicitly supplied, independently audited profile for exact graphs.
Manufactured bundles/audits in tests exercise software contracts only.

## One profile, source-bound graphs

J87 adds `bundle_adapters.opls_qualification`; it reuses the unchanged J84
canonical source capture and independent audit adapter, J85 completeness
qualifier, J86 `QualificationRequest` and strict complete-candidate scope guard. J81–J86 code,
fixtures, receipts and evidence remain byte-identical. No runtime registry or
assignment API is changed. There is no family ranking in J87. Any request or
bundle outside OPLS-AA is rejected; source capture rejects another family before
opening its parameter files.

The pinned request uses `island.j87.opls-qualification-request.v1`:

- `family="OPLS-AA"`, exact `profile_id`, and exact `graph_identities`;
- `fixtures`: each exact graph fixture's `id`, relative `path`, SHA-256, byte
  `size`, and final `graph_identity`;
- `expected_version`: the independently reviewed OPLS-AA release;
- `expected_conventions`: every J82 Class-I convention plus J85's
  `energy_convention`, with explicit audited values;
- `expected_forms`: the six complete, audited bond, angle, torsion, improper,
  nonbonded and charge form declarations, including field order and units.

Each fixture uses `island.j87.exact-graph-fixture.v1` with `id` and `spec`.
`spec` follows J82 `PolymerSpec` through J84's existing parser. Exact file bytes
are verified before normalization; final graphs are built before typing and
must reproduce the independently declared graph identities. Missing, duplicate,
changed, extra, substituted and symlink fixture files are rejected. A topology
change requires new fixture hashes, graph identity, typing and source audit.
No repeat-unit-only typing or element-only match can qualify.

The request file itself requires an external SHA-256 pin. It cannot silently
adopt the supplied manifest's profile, graph scope, version, units or torsion
convention. All source originals are rehashed by J84 before qualification;
canonical raw table bytes are verified again during replay. Opaque primary
attachments retain hashes and audit provenance, and cannot directly supply
parsed parameter rows. No vendor-format converter is run.

The source bundle must declare a complete source file list and hashes,
OPLS-AA family/version, units, functional forms, mixing, exclusions, separate
LJ/electrostatic 1–4 policies, charge origin, atom-type rules, licensing,
authorization and provenance. The independent audit must bind the manifest,
every source hash and every mapped claim. Exact graph-scope claims must be
audited. Parser success, a local copy, a successful build or a hash alone
cannot establish scientific authority. Caller-provided audit pins must come
from an independently trusted review channel; a fabricated self-selected audit
is not authorization. Production code does not create or approve audits.

J87 compares the audited source with independently pinned expected convention
strings and complete form declarations. Unit, normalization, phase, ordering,
version, mixing and exclusion mismatches remain explicit blockers. There is no
unit conversion, equation evaluation, numeric reinterpretation or tolerance
change. Algebraic/dimensional physical validity of an arbitrary source equation
still depends on independent source review; string compatibility is not physical
validation. No 9–6/12–6 coercion, hidden patch, analog fallback, PolyJarvis numeric
import, other-family fallback or zero substitution is provided.

## Exact benchmark and decisions

The default benchmark declares eight finite graph fixtures: hydrocarbon,
aromatic carbon, ether, ester, carbonyl, amide, nitrile and one halogenated motif.
They reuse J85's graph specifications, with methyl caps for ether and ester.
The default `j82_contract_opls_aa` ID denotes the empty source-required contract,
not a new scientific profile. Its expected conventions/version are deliberately
absent. An authorized intake requires a separate explicit request naming the
actual supplied profile, reviewed release, convention contract and exact graphs.
Do not claim universal organic polymer coverage from these eight fixtures.

The output preserves atom-typing, bond, angle, torsion, improper, charge,
nonbonded, interaction_14 and topology coverage, plus Class-I cross-term
non-applicability. Missing terms, blockers, exact graph IDs, source file hashes,
source row locations/JSON pointers, claim hashes and per-instance matched
source records are explicit. No global score hides a missing term.

- `source_bundle_required`: no source or no independently audited authorized
  source; synthetic-only input also cannot qualify scientifically.
- `parameter_bundle_incomplete`: authorized source, but incomplete exact graph
  terms or source/request convention incompatibility.
- `family_scoped_qualification_ready_runtime_blocked`: every declared graph and
  required term passes source, audit, identity and convention checks.

Qualification is always `declared_exact_graphs_only`.
`chemical_family_transfer="not_established"` independently of exact graph
qualification. Blocked-source adapters retain their empty coverage-contract identity separately
from the requested source profile. That contract is an obligation report, not a
substitute parameter source. Only complete candidates enter J86's strict scope
guard. The underlying J85 result remains separately inspectable; J87
adds its request compatibility gate without rewriting historical results.

## Offline intake, replay and verification

Install the inherited C1 execution fence before source IO. Use unused fence
and output directories. If installation fails, stop. No external source is
discovered or fetched, including untracked checkout assets. Do not commit
supplied licensed tables or attachments without redistribution permission.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-emit-unique -- \
  -m scripts.phase_4j87_opls emit --directory /tmp/j87-evidence-unique

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-ingest-fence-unique -- \
  -m scripts.phase_4j87_opls ingest --directory /tmp/j87-intake-unique \
  --request /authorized/opls/request.json --request-sha256 REQUEST_SHA256 \
  --manifest /authorized/opls/manifest.json --manifest-sha256 MANIFEST_SHA256 \
  --audit /independent/opls-audit.json --audit-sha256 INDEPENDENT_AUDIT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-custom-replay-unique -- \
  -m scripts.phase_4j87_opls replay --directory /tmp/j87-intake-unique \
  --snapshot-sha256 CAPTURED_INPUT_SHA256

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-replay-unique -- \
  -m scripts.phase_4j87_opls replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-validation-unique -- \
  -m scripts.validate_phase_4j87_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j87-tests-unique -- \
  -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j87.py
```

Each valid intake has a snapshot-bound qualification receipt, including blocked
intakes. Invalid bytes/schema/identity produce exit 2 and a structured blocked
diagnostic before any evidence output is written. Replay reconstructs all
reports; re-signing a changed result does not bypass the independent snapshot
pin. The terminal strict receipt binds actual software checks, tested commit,
fence events and historical-file hashes. Relocated replay needs captured
code/evidence only; a separate test re-ingests relocated source originals,
request and fixtures, rehashing their bytes. Synthetic controls remain software
tests, not scientific validation. Explicit scientific skips remain skips.

J87 QM, MD, converter, physical fitting, parameter generation, scientific
energy/force evaluation, reserved validation and production executions are zero.
`automatic_assignment_enabled=false`, `runtime_registration=false`,
`accepted_coefficients=0`, `production_validated=false`,
`simulation_readiness="not_established"` even for accepted-looking results.
J72's historical four unintended QM attempts, verified 46–48 MD steps and C1
quarantine remain unchanged. Main and unrelated untracked files are untouched.
