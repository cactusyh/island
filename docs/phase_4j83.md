# Phase J83: evidence-aware coverage registry

Decision: `coverage_registry_implemented_runtime_blocked`.
Base: `c5f7c14b71b7bbcc1192c57d1f5dc02b744c1aa0`.
Branch: `codex/phase-4j83-evidence-aware-coverage-registry`.

J83 adds an offline registry over the unchanged J82 assignment core. OPLS-AA,
GAFF, GAFF2, native PCFF and COMPASS remain separate profiles. The five built-in
entries are empty scientific contracts: they contain no numeric parameters or
typing rules. No new scientific sources or PolyJarvis patches are ingested.

```python
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.coverage_registry import ForceFieldProfileRegistry

graph = normalize(PolymerSpec((("A", "*CCO*"),), 2))
registry = ForceFieldProfileRegistry.contracts()
report = registry.coverage(graph)
cell = report.capability_matrix.lookup(
    "j82_contract_opls_aa", "ether", "finite_linear"
)
```

The final explicit-H graph is constructed before motif detection or typing.
The catalog describes aliphatic C/H, aromatic, ether, ester, carbonyl, amide,
nitrile, halogenated carbon, sulfur, phosphorus, silicon, siloxane, declared
end groups and crosslink/branch anchors. All chemically unclassified atoms are
reported explicitly. Motif detection is descriptive; it establishes neither
atom types nor chemical transfer. End-group provenance comes from the existing
coordinate-free builder. Crosslink endpoints come from explicit topology edits.
Periodic metadata and unsupported stereochemistry retain the J82 blockers.

Matrix cells bind profile identity, motif and topology class to atom typing,
bond/angle, torsion, improper, nonbonded, charge, cross-term and topology
coverage. Every term touching a motif and its transitive coupling dependencies
is included. Motifs overlap. Sum neither motif term counts nor cell coverage
as unique global coverage: `profile_reports` and `unresolved_terms` retain the
authoritative per-profile site-instance inventory. Class-I cross terms are
explicitly not required by its contract; Class-II cross terms remain required
unless an exact, evidenced non-applicability record supports omission. No terms
are zero-filled or supplied by wildcard fallback.

`accepted` means independently audited coverage for this exact graph only;
it requires the complete J82 evidence, convention, charge, scope and dependency
checks. It does not establish physical validation or permit runtime assignment.
`provisional` denotes complete synthetic software controls, never scientific
evidence. All incomplete profiles are `blocked` with retained missing terms and
reasons. The original eight evidence statuses (`verified_stock`,
`verified_extension`, `local_patch`, `reported_claim`, `derived_provisional`,
`fallback`, `missing`, `blocked`) and source hashes are preserved independently
of this coverage decision. A verified label alone cannot authenticate a claim.
Catalog-unclassified cells block even when software parameter controls close.

Future authorized source bundles will be parsed by evidence adapters, such as
J81, into immutable J82 `EvidenceRecord`, `TypingRule`, `TermForm`, `Convention`
and `ParameterRecord` objects. An independent audit must bind each rule, form,
convention, exact parameter and graph scope to its source hash and location.
Supply those audited records separately to `ForceFieldProfileRegistry`.
Adding an authorized OPLS, GAFF/GAFF2 or PCFF bundle requires no resolver change.
PCFF/COMPASS Class-II profiles also require all cross-term definitions,
equilibrium references and conventions. Source authenticity, license, unit,
phase and exact typing checks remain the adapter/audit responsibility. Nothing
in this registry retrieves or licenses source files or registers runtime data.

Replay and strict validation (from the repository root, supported environment):

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j83-replay-unique \
  -- -m scripts.phase_4j83_registry replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j83-validate-unique \
  -- -m scripts.validate_phase_4j83_receipt
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j83-tests-unique \
  -- -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j83.py
```

Each output directory must be new. The inherited C1 process fence is installed
before collection; J73/J74 restrictions are layered by `scope()`. Synthetic
child controls prove forbidden executables cannot execute. Verification logs,
explicit collection, skips, source hashes, tested commit, relocated replay,
and historical byte-preservation inventory are recorded under
`docs/evidence/phase_4j83/validation/`. Relocated code/evidence replay uses the
inherited fence runtime from the original checkout and reports this explicitly.
Scientific execution tests are scope skips, never passed scientific checks.

J83 QM, MD, converter, physical fitting, parameter generation, energy/force
evaluation, reserved validation and production executions are zero. Accepted
coefficients remain zero; runtime registration and automatic assignment remain
disabled. `production_validated=false`,
`simulation_readiness="not_established"`. Historical J72's four unintended QM
attempts and 46–48 MD steps remain recorded as violations, with the C1
quarantine intact. Every tracked J82-base file, including J81/J82 replays,
validators and receipts, is preserved byte-for-byte. Main is untouched.
