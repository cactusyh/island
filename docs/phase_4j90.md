# OPLS provider registry and provisional finite-polymer parameterization (J90)

Base: `8a7933bd8088c23baeb75e9e9f3da61b4535dd67`.
Branch: `codex/opls-provider-registry-provisional-polymers`.

The offline registry accepts one explicitly selected OPLS-AA provider for one
final finite capped graph. It registers parser interfaces, not runtime force
fields. No scientific source bundle was supplied for this phase. The retained
PE, PEO/PEG, PS, PLA and PAN benchmark therefore reports `source_bundle_required`
for all 20 provider/graph combinations, with zero scientific coverage.

The benchmark declares DP2 molecules, explicit single-attachment caps, neutral
as-encoded protonation, explicit hydrogens, stereochemistry/tacticity and export
intent. PS, PLA and PAN use explicitly seeded isotactic construction. These are
graph obligations, not parameterized or scientifically validated polymers. PEG
and PEO are represented by one declared finite oligomer, not all possible PEG
end groups. Infinite-chain and chemical-family transfer remain unestablished.

## Reused code and separation

`assignment_core.PolymerSpec`, its J82 normalizer and the existing polymer builder
construct the final graph before provider assignment. The inherited graph-only
interaction inventory enumerates required bonds, angles and proper torsions.
Conservative improper enumeration requires a source-bound term or explicit
source-bound not-applicable record for every neighbor triplet. J84's pinned
intake reader, duplicate-key-safe JSON parser, path validation and canonical
identity helpers protect external files. J81–J89 and the runtime OPLS registry
are unchanged. This layer accepts site-resolved provider outputs rather than
creating another atom typer, fitting framework or force-field selection API.

`OPLSProviderRegistry` exposes `polypargen`, `ligpargen`, `oplsaa_db` and
`audited_bundle`. There is no network lookup, provider invocation, source
discovery, family ranking or provider merging. PolyParGen and LigParGen accept
self-contained GROMACS text or the canonical exchange. OPLSAA-DB means an
explicitly supplied, pinned local database snapshot and exact graph lookup;
element or substructure similarity never substitutes for an entry. An audited
bundle uses the same canonical exchange plus an independently pinned exact-graph
audit. This is a new exchange interface; it does not reinterpret legacy J84/J89
requests or relax their gates.

## Input contracts

All external IO requires the inherited J72-C1 execution fence and J73 evidence
scope. CLI requests require caller-supplied SHA-256 pins. The manifest and all
listed files, including opaque attachments, are rehashed before parameter
parsing and again during replay/export. Symlink paths, root escapes, duplicate
JSON keys and duplicate/unequal interaction candidates fail closed. Hashes
establish byte identity, not authorship or an independent reviewer's authority.
The caller must supply a trustworthy audit pin and reviewer provenance.

The manifest schema is `island.opls-provider-bundle.v1`. Its exact fields are
`MANIFEST_FIELDS` in `registry.py`. It binds provider/version, OPLS-AA family and
source version, profile identity, exact graph identity, a file manifest, license,
authorization, provenance, synthetic status and the full convention contract.
Every file has `path`, `sha256`, `size`, `role`, `format`:

| Role | Format | Meaning |
|---|---|---|
| graph | chemical_graph_json | Exact `data(finite_graph(settings))` JSON |
| parameters | gromacs_itp | Self-contained, explicit numeric GROMACS sections |
| parameters | canonical_json | Site-resolved canonical exchange |
| parameters | oplsaa_db_snapshot | Exact-graph entries in a supplied DB snapshot |
| omissions | improper_omissions_json | Explicit sourced improper omissions |
| opaque | opaque | Hash-bound original attachments; never parsed |

Exactly one graph and at most one parameter source are allowed. All files belong
to the same selected provider. License metadata requires an identifier,
`authorization="authorized_local_use"` and redistribution status. Provenance
requires citation, acquisition method and source URI. Missing authorization is
checked before opening parameter files. No source file is copied to evidence.

The request schema `island.opls-polymer-request.v1` has `provider`, `settings`,
and `graph_identity`. `settings` must supply every `SETTINGS` field in
`polymer.py`: PSMILES, positive repeat count, two explicit cap SMILES,
`protonation_state="as_encoded"`, integer charge, stereochemistry, tacticity,
integer stereo seed, and intended export engine. Protonation must be encoded in
the graph; no protonation guessing occurs. Missing/ambiguous stereochemistry,
builder limitations, radicals, multiple components and a disagreeing charge
state block construction. Any graph/settings change requires a new provider
assignment. A supplied graph snapshot must match all final bytes and atom order.

Canonical source `island.opls-provider-records.v1` contains `family`,
`graph_identity`, `forms`, `atoms`, `records`, `pairs14`, `exclusions`. Forms are
structured definitions matching the bounded `FORMS` contracts. Each atom has
site, atom type, element, source mass and source location. Each interaction has
id, family, ordered sites, form, literal string arrays, applicability, reason,
and source location. Raw values remain unchanged in the ephemeral canonical
object. Public reports retain claim hashes and source pointers, not arrays.
The source pointer combines an exact JSON pointer with the reported table/row
location, or exact section/line references for native text.

GROMACS parsing supports explicit defaults (12-6, geometric sigma/epsilon),
atomtypes, one moleculetype with nrexcl=3, atoms with explicit charge/mass,
harmonic bonds/angles, periodic/RB/OPLS torsions, harmonic impropers, explicit
pairs and exclusions. It rejects includes/macros, implicit parameter lookup,
wildcards, pair overrides, unsupported functions and repeated sections. Supply
an explicitly sourced self-contained output and original opaque attachments
when the original provider separates defaults into a different file. The parser
does not invoke a preprocessor or flatten sources. Missing improper sections do
not establish that impropers are absent; explicit omission evidence is needed.
Unsupported native formats remain blocked; parser success alone is not an audit.

The convention contract explicitly includes energy/length/charge/angular units,
12-6 sigma-zero meaning, geometric mixing, graph-distance exclusions, separate
LJ/Coulomb 1–4 scales, charge origin, typing origin, long-range/cutoff/tail rules,
torsion phase/normalization and improper ordering. Canonical ingestion currently
supports kJ/mol and nm only; incompatible units/forms are blocked. It never
converts 9-6 into 12-6. Source charges are summed as decimal literals without
renormalization or tolerance adjustment. Missing terms are never zero-filled.

The independent audit schema is `island.opls-provider-audit.v1`: reviewer,
authority basis, `independent=true`, manifest SHA, graph identity, convention
identity, and complete record/atom claim-hash arrays. It binds the exact source
scope, not a chemical-family claim. Synthetic audits can test software behavior
but never create scientific coverage or scientific evidence in this phase.

## Results and export boundaries

| Result | Meaning |
|---|---|
| source_bundle_required | Pinned authorized source/audited-provider audit absent |
| provenance_only | Verified file provenance supplied, no parameter source |
| parameter_bundle_incomplete | Missing, ambiguous, incompatible or corrupt evidence |
| provisional_parameterization | Complete exact-graph provider data, independently unaudited |
| audited_exact_graph | Complete independently reviewed source-bound exact graph |

Each report exposes typing, bond, angle, torsion, improper, charge, nonbonded,
1–4 and topology obligations separately. Complete source records can be counted
while an overall convention or charge blocker still prevents provisional status.
The benchmark aggregate preserves every blocked provider decision. No accepted
coefficient or runtime registration is created, even for audited results.

`export_parameters()` reopens and rehashes the originals, requires a complete
provisional/audited result and the explicitly intended engine, and exclusively
writes a new local artifact **outside Git worktrees**. Every artifact contains a
sealed provenance receipt, original source hashes, graph identity, raw literal
exchange and explicit non-production flags. Never commit those exports, source
tables, captured canonical payloads or proprietary attachments. Only sanitized
reports, hashes, source pointers and derived audit metadata belong in Git.

GROMACS export is an ITP parameter/topology fragment. LAMMPS export is a
coefficient fragment with explicit type-to-site mappings, not a runnable system
or coordinate file. Its exact kJ→kcal, nm→Å and harmonic half-factor conversions
are documented in the artifact; source values remain in the provenance payload.
LAMMPS RB torsion conversion is deliberately unsupported. Pair-style, cutoff and
long-range engine qualification remain separate obligations.

**OpenMM export is a data-only XML parameter exchange, explicitly marked
`loadableForceField=false`.** OpenMM's standard NonbondedForce uses arithmetic
sigma mixing, which is incompatible with OPLS geometric sigma. This phase
preserves the geometric-mixing contract, all exact terms, exclusions and 1–4
pairs in XML rather than emitting a silently incompatible ForceField XML.
Implementing and evaluating a custom OpenMM force realization requires a later
authorized phase. No OpenMM import, Context creation or engine validation occurs.

Round-trip checks are serialization checks: `read_export()` validates caller
pins, every receipt, exact graph, source hashes and the regenerated native/XML
body. They do not run GROMACS, LAMMPS, OpenMM or a converter and do not establish
physical equivalence or simulation readiness.

## Replay

Use the supported Python environment, with repository and `src` on PYTHONPATH:

```bash
export PYTHONPATH=.:src
PY=../island-validation/phase4g1-env/bin/python
$PY -m scripts.run_phase_4j72_c1 --output /tmp/j90-replay-new -- \
  -m scripts.phase_4j90_providers replay
$PY -m scripts.run_phase_4j72_c1 --output /tmp/j90-validate-new -- \
  -m scripts.validate_phase_4j90_receipt
```

For authorized local inputs (paths and hashes below are caller-supplied), use
`ingest --directory /outside/checkout/report --request REQUEST
--request-sha256 SHA --manifest MANIFEST --manifest-sha256 SHA`, optionally
`--audit AUDIT --audit-sha256 SHA`. Use the same external flags for replay and
pin the emitted `captured-input.json` with `--snapshot-sha256`. Metadata-only
replay fails clearly if originals or pins are missing. Relocation requires the
same relative source layout and bytes, not the old absolute root.

Focused tests: `pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j90.py`
inside a fresh C1 fence. The inherited portable J73–J89/graph/typing regression
inventory, strict historical receipts/replays, Ruff, compileall, dependency
check, relocated replay and byte-for-byte preservation are archived under
`docs/evidence/phase_4j90/validation`. A skipped real-source integration test is
explicit: no real provider source was supplied. Manufactured controls only test
parsing, schemas, audit gates and serialization; they are not scientific data.

All gates remain: `accepted_coefficients=0`, `runtime_registration=false`,
`automatic_assignment_enabled=false`, `production_validated=false`,
`simulation_readiness="not_established"`,
`chemical_family_transfer="not_established"`. New QM, MD, converter, physical
fitting, scientific energy/force, parameter generation, production and reserved
validation executions are zero. Historical J72 violation/C1 quarantine and all
base-tracked files remain immutable; main is untouched.
