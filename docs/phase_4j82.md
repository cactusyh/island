# Phase 4J82: coverage-first polymer assignment core

J82 adds an opt-in, coordinate-free **coverage analyzer**, under
`island.forcefields.assignment_core`, based on J81 commit
`279e811daece681477845c10d15169504bdbb5e8`. It adds no scientific parameters,
changes no existing runtime API or registry, and performs no scientific
execution. All five built-in contracts are intentionally empty: OPLS-AA,
GAFF, GAFF2, native PCFF and COMPASS remain separate identities.

`automatic_assignment_enabled=false`, `accepted_coefficients=0`,
`production_validated=false`, `simulation_readiness="not_established"` apply
to every decision, including complete synthetic controls. Completeness is
not chemical transfer, scientific qualification, an evaluator implementation
or production readiness. J81's blocked source-ingestion result, J76's 378
unresolved native requests and J72's four unintended QM attempts / 46–48 MD
steps remain historical facts, with their original receipts and quarantine.
J82's own scientific execution counters are zero.

## Internal contracts and topology

Frozen dataclasses define ChemicalGraph, PolymerSpec, TypedGraph,
ForceFieldProfile, ParameterRecord, EvidenceRecord, CompletenessReport and
AssignmentDecision. Auxiliary records describe rules, conventions, term
forms, requests and explicit missing terms. JSON serialization rejects
nonfinite values and hashes canonical content. Profiles bind parameter
ordering, source units, equations, evidence and complete convention identity.
Raw numeric strings are preserved; J82 performs no implicit unit conversion.

PolymerSpec represents a repeat library, explicit copolymer sequence, DP,
PSMILES connection orientation, tacticity/stereo seed, end groups, branches,
crosslinks and periodic metadata. Existing `parse_psmiles`,
`build_polymer_from_sequence`, `to_rdkit`, `RDKitSmartsAtomTypingEngine` and
`derive_interaction_inventory` are reused. The builder creates only its 2D
connectivity depiction; the analyzer discards its coordinates. There is no
3D embedding, optimization or energy/force evaluation.

The supported normalization path constructs finite, linear/copolymer chains,
adds single-dummy end groups, explicit single-dummy branches and specified
single-bond intrachain crosslinks, then completes explicit hydrogen sites.
Branch/crosslink anchors reference the pre-edit builder's one-based heavy-site
IDs. Chemically invalid attachment and duplicate bonds are rejected. PSMILES
stereo tags and builder tacticity remain part of the final graph. Connection
points must agree with the parser's head/tail orientation; they cannot silently
reorient a repeat. Periodic metadata is preserved with a blocking reason:
periodic image graphs and arbitrary multi-chain networks require an explicit
future topology adapter. They are not approximated with a finite chain.

Final graph identity binds atom order, stereo, spec, topology-edit metadata,
normalizer and RDKit version. Changing end groups, sequence, branching,
crosslinks or periodic metadata requires reconstruction and retyping. Even
an unchanged molecular formula cannot make a stale typed snapshot reusable.
Equivalent alternative input encodings are not claimed to share a graph hash;
replay is deterministic for the same declared input and dependency version.

## Typing, evidence and completeness

The existing SMARTS engine resolves explicit override graphs. The adapter
adds integer priority: a matching higher-priority rule overrides matching
lower-priority rules; equal-priority incompatible types stay ambiguous.
Override cycles and contradictions are rejected. The historical engine does
not perform chiral SMARTS matching, so stereospecific SMARTS are explicitly
blocked rather than silently matched without chirality.

Evidence statuses are `verified_stock`, `verified_extension`, `local_patch`,
`reported_claim`, `derived_provisional`, `fallback`, `missing` and `blocked`.
Only verified stock/extension records matching an independently supplied trust
set can contribute to nonsynthetic coverage. Trust entries bind the full
EvidenceRecord, including family, resolution, citation, location and source hash,
plus canonical audited-claim hashes for each rule, form, convention, parameter
and exact graph scope. Changing a row under an unchanged source hash cannot
reuse its audit. A label
or checksum alone cannot authenticate a source. Future source adapters must
hash the actual asset and obtain an independently reviewed evidence record;
this pure core does not fetch assets or grant authority. Synthetic evidence
is usable only in explicitly named synthetic-control profiles and cannot be
included in an evidence-candidate profile. PolyJarvis patches and unqualified
literature values receive no special acceptance path.

Each interaction request has a stable site-instance identity and a distinct
parameter identity based on family, ordered types and a radius-two rooted
chemical environment. Aliases are not independent parameters. No wildcard,
zero-fill, element-only fallback, implicit reversal or cross-family overlay
is performed. Source adapters must supply the exact ordered form, units and
any orientation transformation as independently evidenced records. All
conflicting candidates retain their identities; unequal rows block selection.
An explicit, evidenced not-applicable record may exempt an exact improper or
cross term without fabricating a zero coefficient.

Class-I inventories include bonds, angles, proper torsions, conservative
neighbor-triple improper candidates, per-site nonbonded and per-site charges.
Class-II inventories additionally enumerate bond–bond, bond–angle,
angle–angle, bond–bond 1–3, end/middle-bond–torsion, angle–torsion,
angle–angle–torsion and five-site torsion–torsion candidates. Every coupling
must bind its primary parameter identities and explicit equilibrium reference
values; missing fields, changed reference values and unresolved dependencies block
it. Conservative enumeration may require explicit source-supported absence
for a term. Class-I reports identify cross terms as not required by that
contract, rather than implying they passed Class-II coverage.

Convention obligations include functional forms, energy/length/angular units,
atom typing and chemical scope, charge origin and component charge conservation,
nonbonded form and sigma/epsilon definitions, mixing, 1–2/1–3 exclusions,
independent LJ/Coulomb 1–4 policy, torsion normalization/phase, improper and
cross-term ordering, electrostatics, long range, cutoff, tail and finite versus
periodic scope. A 9–6 equation cannot be coerced into a 12–6 record. No
nonbonded or charge policy is guessed. Nonsynthetic candidate scope must list
the exact graph identities authorized by its evidence adapter.

Coverage reports separate atom typing, bonded, torsion, nonbonded, charge,
cross-term and topology coverage, with missing instances, candidate source
IDs, blocking reasons, source hashes and replay identity. Profile ranking uses
completion, coverage fraction and blocker count with a stable ID tie-breaker;
it conveys no scientific preference or authorization. J82 always returns a
blocked AssignmentDecision because runtime activation is disabled.

## Extending chemical coverage

A future OPLS, GAFF/GAFF2, PCFF or polymer-family extension supplies a separate
versioned ForceFieldProfile: audited source records, atom-type SMARTS,
explicit conventions and forms, graph scope and exact environment-bound
parameter records. It can reuse the same resolver and Class-I/Class-II
validator. Family-specific source parsing remains outside the core. New
periodic or stereospecific adapters must explicitly preserve image/stereo
identity. Adding a family never patches the native PCFF registry or imports
another family's coefficients. Runtime qualification, evaluators and scientific
validation require a later authorized phase.

## Offline replay and validation

Use the supported environment and unchanged J72-C1 process fence with the
J73/J74 evidence-only layer. A missing fence stops before audit/replay. Each
command needs a fresh output directory. No real scientific fixtures are enabled.

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j82-replay \
  -- -m scripts.phase_4j82_coverage replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j82-strict \
  -- -m scripts.validate_phase_4j82_receipt

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j82-tests \
  -- -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j82.py \
  tests/test_phase_4j73.py tests/test_phase_4j74.py tests/test_phase_4j75.py \
  tests/test_phase_4j76.py tests/test_phase_4j78.py tests/test_phase_4j79.py \
  tests/test_phase_4j80.py tests/test_phase_4j81.py
```

The immutable 15-motif input covers hydrocarbon, ether, ester, aromatic,
amide, halogen, sulfur, phosphorus, silicon, end groups, branching, crosslink,
copolymer, tacticity and blocked periodic metadata. Empty contracts generate
75 incomplete candidate reports and 15 blocked runtime decisions. Synthetic
numbers in tests exist only to test software coverage, charge conservation,
conflicts and dependency gates; none enters the scientific evidence reports.

The receipt validator reconstructs every report from the pinned input, checks
an explicit implementation inventory against the tested Git commit, validates
terminal collection/JUnit and fenced events, and compares every tracked J81
base file with its original Git blob. Re-signing a changed result cannot erase
blockers or activate assignment. Relocated replay uses copied source/scripts
and evidence, no scientific source assets. Its inherited security bootstrap
still originates in the original checkout and is recorded explicitly.
Validation counts and exact tested commit are in `receipt.json` and
`validation/summary.json`. Skipped scientific controls remain explicit skips,
not passed scientific checks.

A minimal opt-in call is:

```python
from island.forcefields.assignment_core import (
    PolymerSpec, builtin_contracts, normalize, resolve,
)

graph = normalize(PolymerSpec((("A", "*CCO*"),), degree_of_polymerization=3))
decision = resolve(graph, builtin_contracts())
assert decision.status == "blocked"
# Inspect decision.reports for coverage and explicit missing terms.
```

The terminal regression inventory also explicitly includes `test_psmiles.py`,
`test_polymer_builder.py`, `test_polymer_termination.py`, `test_atom_typing.py`,
`test_atom_typing_rules.py`, `test_atom_typing_compatibility.py`,
`test_parameter_inventory.py` and `test_parameter_assignment_graph_only.py`.
No historical implementation or evidence globs are widened.
