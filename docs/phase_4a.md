# Phase 4A: graph-based atom typing

Phase 4A introduces deterministic atom-type assignment while preserving ISLAND's
layer boundaries:

```text
MolecularSystem chemical graph and stable IDs
    -> AtomTypingEngine
    -> external AtomTypingResult
    -> future Phase 4B parameter assignment
```

An atom type is never written to `AtomSite` or its metadata. Coordinates,
conformation provenance, partial charges, and numerical force-field parameters are
not atom-typing inputs or outputs. Successful typing therefore does **not** mean a
system is parameterized.

## Public API

The RDKit-free interfaces are importable from `island.forcefields` or
`island.forcefields.typing`:

- `AtomTypingEngine`
- `AtomTypingRule`
- `AtomTypingRuleSet`
- `AtomTypingResult`
- `AtomTypeAssignment`
- `SiteTypingDiagnostic`
- `island_demo_v1_ruleset()`

`RDKitSmartsAtomTypingEngine` is loaded lazily when requested. It implements
ordinary RDKit SMARTS, not Foyer's complete dialect and not an established
force-field file format.

```python
from island.builders import build_linear_polymer
from island.forcefields import (
    RDKitSmartsAtomTypingEngine,
    island_demo_v1_ruleset,
)

system = build_linear_polymer("[*:1]CC[*:2]", dp=5, generate_3d=False)
result = RDKitSmartsAtomTypingEngine().type_system(
    system, island_demo_v1_ruleset()
)

for site_id, assignment in sorted(result.assignments.items()):
    print(site_id, assignment.atom_type, assignment.matched_rule_ids)
print(result.complete)
```

`type_topology(topology, ruleset)` is the coordinate-free entry point;
`type_system(system, ruleset)` is a convenience wrapper that intentionally ignores
coordinates. Both return stable-site-keyed results and leave their input unchanged.

## Graph-only RDKit conversion

`topology_to_rdkit_graph()` and `system_to_rdkit_graph()` convert only the
chemical graph and return explicit `site_id -> RDKit index` and reverse mappings.
They preserve supported atomic identity, formal charge, aromaticity, bond order,
and hydrogen representation, but do not require, generate, embed, optimize, or
attach coordinates. The existing `to_rdkit()` remains coordinate-aware and uses
the shared graph conversion internally.

The graph conversion and SMARTS backend live in the optional chemistry layer.
Importing the core or typing rule/result interfaces does not import RDKit.

## Rules and deterministic resolution

A rule has a unique ID, an output type name, a SMARTS pattern with exactly one
target atom marked `:1`, explicit `overrides`, and a description or source. Before
matching, rulesets reject invalid SMARTS targets, duplicate IDs, unknown or
self-overrides, and override cycles.

Every rule is exhaustively matched and repeated substructure matches are
deduplicated for each target site. Resolution is independent of rule declaration,
site insertion, bond insertion, and RDKit atom-index order:

1. collect all matching rule IDs for a stable site;
2. eliminate a candidate when another matching rule explicitly overrides it,
   including through a transitive override chain;
3. assign only when all surviving rules agree on one atom-type name;
4. report no survivors/matches as untyped and disagreeing survivor types as
   ambiguous.

There is no first/last-match rule, hidden specificity heuristic, priority value, or
generic fallback. Diagnostics record matched, eliminated, surviving, and selected
rules, including which matching rules eliminated each candidate.

Strict mode (the default) raises `IncompleteAtomTypingError` and attaches the
structured incomplete result as `error.result`. Use `strict=False` to inspect
untyped and ambiguous sites without fabricating assignments.

## Hydrogen policy and supported chemistry

Each ruleset declares one policy:

- `explicit_required`: every chemically required hydrogen must be an authoritative
  `AtomSite`; omitted implicit or atom-count hydrogens are rejected.
- `implicit_allowed`: SMARTS matching may use RDKit's graph-derived implicit
  hydrogen semantics without adding authoritative sites.

`explicit_required` does not require a molecule to contain hydrogen when none is
chemically required. Temporary chemistry changes and silent site addition/removal
are never performed.

Phase 4A supports atomistic `AtomSite` graphs that the graph adapter can sanitize,
with single, double, triple, and aromatic bonds and integer formal charges. It
rejects coarse-grained beads and unsupported bond/atom representations. SMARTS
chirality matching is intentionally disabled in this phase; chemical
stereochemistry remains authoritative but does not affect the demonstration type
rules.

## Results, signatures, and reuse

`AtomTypingResult` contains:

- assignments and per-site resolution diagnostics keyed by stable site ID;
- ruleset name/version and engine/dependency versions;
- untyped/ambiguous site lists and a completeness flag;
- deterministic graph, ruleset, and combined typing SHA-256 signatures.

`result.is_compatible_with(system_or_topology, ruleset)` is a reuse check, not a
cache or incremental typing mechanism. Coordinate replacement and unrelated
metadata do not invalidate it. Stable-ID, element, formal-charge, aromaticity,
hydrogen-representation, bond-order, connectivity, or rule-content changes do.
Canonical sorting makes insertion-order changes irrelevant.

## Demonstration rules and Phase 4B boundary

`island_demo_v1` is deliberately small and illustrative. Its `demo_*` labels show
aliphatic/aromatic/carbonyl carbon, several oxygen environments, neutral/positive
nitrogen, explicit hydrogen environments, and broad-rule overrides. It is not
GAFF, GAFF2, OPLS, PCFF, CVFF, or any other production force field. It provides no
partial charges, Lennard-Jones terms, bonded parameters, provenance validation
against a published force field, or chemical coverage guarantee.

Phase 4B should consume an `AtomTypingResult` only after compatibility and
completeness checks, assign independently versioned per-site and interaction
parameters into `ParameterizedSystem`, diagnose missing/ambiguous parameters, and
retain parameter/source provenance. Charge assignment, real parameter sets,
force-field file readers, engines, minimization, MD, LAMMPS export, packing, and
crosslinking remain outside Phase 4A.
