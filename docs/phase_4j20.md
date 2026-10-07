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
  `island_crosslink_network_plan_v1`;
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

The focused J20 suite passes **21 tests**. J18 and J19 compatibility suites
pass unchanged alongside it. Readiness remains deliberately conservative:

```text
production_validated = False
simulation_readiness = "not_established"
```

This phase does not claim automatic reaction prediction, universal chemistry,
periodic-image crosslinking, packing, force-field parameter completeness,
periodic MD, or scientific validation of crosslinked polymer properties.
