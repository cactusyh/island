# Phase 4J20: force-field-independent crosslink networks

J20 adds an explicit crosslink planning layer on top of the J18
`FinalChemicalGraph`. It does not infer reactions and does not import a force
field backend. The caller marks reactive atoms in graph metadata, explicitly
performs any hydrogen or leaving-group removal, and then asks the planner to
select new bonds.

## Public contract

`island.crosslinking` provides:

- `ReactiveSiteRule`, an immutable rule matching one metadata key/value and a
  declared element set, bond order and site capacity;
- `ReactiveSite`, an immutable observation of the atom ID, component, current
  valence and remaining capacity;
- `CrosslinkNetworkPlan`, checksum-protected JSON with schema
  `island_crosslink_network_plan_v2` (see the corrective review below);
- `identify_reactive_sites`, `plan_crosslinks`, `apply_crosslink_plan` and
  `generate_crosslink_network`;
- `save_crosslink_plan` and `load_crosslink_plan`.

Plans bind the exact input graph identity, the complete sorted reactive-site
records, candidate-pair identities, the exact selected pair list, rule,
strategy, intercomponent policy, capacity, seed, provenance and evidence. The
only implemented strategy is `seeded_random`: eligible pairs are generated in
lexicographic atom-ID order and shuffled with a local `random.Random(seed)`;
the first feasible pairs are selected while respecting both rule capacity and
the requested maximum usage. Dictionary insertion order and coordinates do not
affect this process.

The default requires different connected components. Setting
`require_intercomponent=False` is an explicit opt-in to intramolecular pairs.
Periodic graphs are rejected because periodic-image crosslink semantics are not
defined in this phase.

For a J18 linear polyethylene chain, the caller can explicitly create the
reactive graph before planning:

```python
from island.crosslinking import ReactiveSiteRule, generate_crosslink_network
from island.graph import build_psmiles_graph, final_graph

system, original, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
carbons = [i for i, atom in system.topology.sites.items() if atom.element == "C"]
for carbon in (carbons[0], carbons[-1]):
    hydrogen = next(
        i for i in system.topology.neighbors(carbon)
        if system.topology.sites[i].element == "H"
    )
    system.topology.remove_site(hydrogen)
    system.coordinates._positions.pop(hydrogen)
    system.topology.sites[carbon].metadata["reactive"] = True
graph = final_graph(system, transformations=original.payload["transformations"])
rule = ReactiveSiteRule(
    "declared carbon", "reactive", True, ("C",), 1, 1.0,
    "caller-declared topology operation", ("fixture declaration",),
)
result = generate_crosslink_network(
    system, graph, rule, target_crosslinks=1,
    require_intercomponent=False,
    provenance="caller-selected intrachain link",
    evidence=("fixture declaration",),
)
```

This example uses a controlled software fixture. It is not an experimentally
established reaction mechanism.

## Application and identity boundaries

`apply_crosslink_plan` validates the exact graph and plan identities, copies the
input system, adds only the selected bonds, rebuilds angles and dihedrals, and
validates the resulting graph and atom valences. It does not remove atoms,
delete hydrogens, repair valence, assign atom types, assign charges, look up
parameters, or use coordinates to select a pair. The plan identity is included
in the `GraphTransformation` parameters and prior graph history is retained.

The original system remains unchanged. Coordinate-only replacements preserve
the final graph and plan identities; missing or nonfinite coordinates are
rejected when applying a plan. A changed topology, formal charge, stale graph,
duplicate/self/existing bond, unknown site, over-capacity selection, invalid
bond order, periodic input, or tampered checksum is rejected.

The saved plan checksum detects accidental or untrusted edits to its JSON.
`load_crosslink_plan(..., expected_identity=...)` additionally checks against a
trusted prior plan identity, which is needed to detect a deliberately re-signed
record. The output graph transformation also carries the plan identity.

No PCFF, OPLS-AA, GAFF or GAFF2 module is imported by the crosslink planner.
After a crosslink operation, the selected force-field preparation path must
retype and reassign the new final graph. Existing missing-interaction gates
remain authoritative.

## Acceptance fixtures and receipts

The main acceptance fixture uses two disconnected `[*:1]CC[*:2]` DP3 chains.
The caller removes one terminal hydrogen from each chain and marks those two
carbons. With seed `2026`, one bond is selected from one candidate and the
components decrease from two to one. A separate two or three fragment C=C
software fixture supplies multiple candidates and deterministic multi-bond
capacity checks. The J18 `[*:1]CC[*:2]`, DP3 fixture is also exercised
after the caller removes one terminal hydrogen from each explicitly selected
carbon; the planner itself performs no hydrogen deletion. The tests also cover
periodic rejection, stale graph/plan identities, duplicate and self pairs,
unknown IDs, invalid valence/capacity, altered seed/rule/provenance/evidence,
bundle relocation and child-process loading.

The original J20 delivery had **21 focused tests**. The correction below adds
further controls and reruns J18/J19 compatibility tests. Readiness remains
deliberately conservative:

```text
production_validated = False
simulation_readiness = "not_established"
```

This phase does not claim automatic reaction prediction, universal chemistry,
periodic-image crosslinking, packing, force-field parameter completeness,
periodic MD, or scientific validation of crosslinked polymer properties.

## Corrective review: cumulative usage, membership and chemical states

The reviewed implementation `7e1783d2d4db95d10decfbf869adca13712b41db`
contained three reproduced defects. The unmodified implementation selected
`[1,7]` with seed 2026 and then reused site 7 in `[3,7]` with seed 0;
replaced explicit `original-A`/`original-B` molecule labels by `"0"`; and
accepted a bond between two neutral NH3 nitrogen atoms, leaving each nitrogen
neutral with four single bonds. `scripts/reproduce_phase_4j20_correction.py`
runs those cases against either implementation. The original
`docs/evidence/phase_4j20.json` is retained verbatim; the separate
`docs/evidence/phase_4j20_correction.json` contains the before/after results,
commands and preservation checks.

### Cumulative capacity and record versions

The rule's `maximum_crosslinks_per_site` is a **lifetime limit**. Each v2
`ReactiveSite` records `lifetime_capacity`, `consumed_crosslinks`,
`remaining_capacity` (lifetime minus consumed), and `available_capacity`
(the smaller of remaining lifetime allowance and current chemical headroom).
The planner's `maximum_crosslinks_per_site` argument is a **per-batch limit**;
its default remains 1, and `None` means no additional batch restriction.
A seed change, a new rule name, or a new plan cannot erase prior usage.
Previously recorded lifetime limits cannot be changed for retained sites.

Accounting reconciles every crosslink provenance entry with a unique
crosslink transformation and a bond still present with the recorded order.
Missing, duplicate or conflicting entries, incorrect seeds/plan references,
missing bonds and conflicting lifetime limits reject. Existing J18 explicit
crosslink transformations and original J20 transformation evidence contribute
to usage; their serialized records are not rewritten. Usage counts every
crosslink at the site, regardless of the rule that selected it. Ordinary bonds
do not count as crosslinks. Exhausted marked sites remain observable with zero
available capacity and are excluded from candidates, so other sites can react.
An initially saturated, unused marked site is invalid input.

New plans use `island_crosslink_network_plan_v2`, and new site records use
`island_reactive_site_v2`. A plan additionally binds the exact input history
identity and the structural policy version. Original v1 plans lack these
contracts and require explicit replanning; they are rejected at the schema
gate rather than silently converted. Rule serialization remains v1. J17–J19
records, graph bundle schemas, force-field contracts and historical evidence
are unchanged.

### Original membership and current components

All three intermediate/final graph constructions in `apply_crosslink_plan`
receive the exact input `molecule_membership` mapping. These labels identify
original molecules or chains, including caller-defined labels on systems
without `chain_id` metadata. Connected components are recomputed after adding
bonds and can merge while the original molecule labels remain distinct.
Relocated graph bundles retain both mappings.

### Bounded structural policy

`island_crosslink_cho_n_structural_policy_v1` is checked against the input graph,
the proposed bonds before publishing a plan, and the copied output graph during
application. It uses element, formal charge, bond orders, aromaticity and
`radical_electrons`. It permits only integral bond orders 1, 2 or 3 and the
following nonaromatic states:

| Element | Formal charge | Bond-order valence | Reaction-site support |
| --- | --- | --- | --- |
| H | 0 | 1 | Spectator only |
| C | 0 | 2, 3, 4 | Remaining headroom below 4 |
| O | 0 | 1, 2 | Remaining headroom below 2 |
| N | 0 | 3 | Spectator only |
| N | +1 | 4, with four single bonds | Spectator only |

The incomplete C/O states support caller-declared topology fixtures. They do
not establish a reaction mechanism or an electronic state. Missing radical
metadata does not certify a closed-shell structure. Explicit nonzero radical
counts, boolean/malformed radical counts, ambiguous alternative radical keys,
aromatic atoms/bonds, other element/charge states and reactive nitrogen are
unsupported. In particular, neutral four-coordinate nitrogen rejects, and
valid NH4+ is permitted as an unchanged spectator, never as a site with spare
capacity. The operation makes no formal-charge, hydrogen or radical-state
changes.

The corrected reproduction selects `[2,4]` in round two, preserves both
custom molecule labels, and rejects the NH3 crosslink at eligibility validation.
The corrective tests cover three rounds with plan serialization and bundle
relocation between rounds, direct planning/application and the convenience API,
input preservation on success/failure, and semantically inconsistent graph or
plan records with recomputed checksums. Checksums remain integrity identifiers;
callers must retain a trusted graph/plan identity to detect wholesale replacement
of an otherwise self-consistent record.

### Executed verification

Final results: **1928 passed, 10 skipped** in the ordinary suite; **93 passed**
in the focused J20/J18/J19 suite (50 J20 and 43 J18/J19 controls). Ruff, both
active-environment pip invocations and whitespace checks passed. A read-only
comparison verified all **704** baseline files outside the two intended edits
are byte-for-byte unchanged, including the original J20 receipt.

The corrective receipt records the exact commands and results. Reproduction
against the reviewed source can be repeated without checking out or editing main:

```sh
mkdir -p /tmp/island-j20-correction/base
git archive 7e1783d2d4db95d10decfbf869adca13712b41db src | tar -x -C /tmp/island-j20-correction/base
PYTHONPATH=/tmp/island-j20-correction/base/src python scripts/reproduce_phase_4j20_correction.py --expect defects
PYTHONPATH=src python scripts/reproduce_phase_4j20_correction.py --expect corrected
```

Validation commands:

```sh
PYTHONPATH=.:src python -m pytest -q tests/test_phase_4j20_crosslinking.py tests/test_phase_4j20_correction.py tests/test_phase_4j18_final_graph.py tests/test_phase_4j19_unified.py
PYTHONPATH=.:src python -m pytest -q
ruff check .
python -m pip check
pip check
git diff --check
```

Both pip commands address the active `island` environment; they are not checks
of two different environments. Readiness remains `production_validated=False`
and `simulation_readiness="not_established"`. No force-field preparation,
implicit chemistry repair, packing or periodic MD is added by this correction.
