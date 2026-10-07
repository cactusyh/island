# Phase 4J17 — source-bound typed graphs

## Contract and compatibility decision (frozen before production changes)

Base: `8b1411ca338a24d0a0bede4ae2d8717faad843fb`; fetched origin main is
identical, and its tree equals reviewed `85a885d2a1a25751acac93bcddf1364c97562308`.
The requested branch already existed, clean, at that base.

Add immutable JSON-owned `PCFFTypedGraph` and `PCFFGraphCharges` records. The
external opt-in API is `bind_pcff_types`; legacy `assign_pcff_source_types` keeps
its automatic-agreement semantics. Exact source identity, final chemical graph,
stable integer IDs, provenance and evidence references are bound into identity.
Coordinates and repeat/molecule names do not select types. Source authorization
remains pinned. Structural consistency is a bounded validation result; external
chemical authority is an assertion, and automatic perception is not performed by
this contract. Unsupported structural constraints reject explicitly.

`assign_typed_pcff_charges` uses the existing increment resolver on a validated
typed graph. `provide_pcff_charges` validates a supplied vector independently of
increment availability. The latter requires elementary_charge units, exact IDs,
finite real values, provenance and evidence, and declared component totals equal
to the final graph formal charges. No redistribution or normalization is allowed.

The common class-II resolver accepts either historical automatic records or new
records without fabricated automatic/native evidence. Additive assignment/model
schemas retain old serialized identities. Public preparation receives explicit
records through distinctly named options; bundle and workflow reconstruction
must validate these records offline using only the supplied current FRC path.

## Frozen acceptance matrix

| Case | Required outcome |
|---|---|
| Historical CHO, CHN, aromatic controls | identities and results unchanged |
| PSMILES `[*:1]CC[*:2]`, DP3, seed 2026, H ends | externally established c/h labels; legacy agreement rejects; native and provided charge paths execute with existing compatibility policy |
| Same polyethylene, native increments | exact native vector; independent raw FRC increment reconstruction |
| Same polyethylene, provided charges | raw FRC reconstructed vector supplied explicitly; no native prerequisite |
| Software source with removed increments | provided path still resolves; native path incomplete |
| Saturated ether with valid labels/charges under strict policy | required missing interactions still prevent publication |
| O=O, N#N, wrong element/connectivity/valence/charge state | precise rejection |
| Unknown labels; missing/extra/bool IDs; NaN/Inf; units; totals | rejection |
| Source/topology/formal-charge/type/charge mutation | affected identities or validation invalidate |
| Rechecksummed record/bundle tampering | semantic rejection |
| Coordinate replacement; relocated bundle; separate-process continuation | supported without scientific typing imports |

Independent chemical basis for polyethylene: pinned FRC atom_types c (generic SP3
carbon), h (hydrogen bound to C/Si/H), plus explicit final-graph valence and parent
checks. This establishes bounded source-family consistency, not experimental
validation of polyethylene properties. Charge evidence must be raw source rows
read independently of production code. The converter receives those same types
and charges and does not certify their origin. No synthetic case counts as
scientific acceptance. Numerical acceptance requires normal converter output,
independent coefficient/dependency checks, LAMMPS energy/component/force agreement,
finite differences, fresh-verified minimization, relocated continuation and split
trajectory/full RNG equality at retained tolerances. Failures remain evidence.

`production_validated=False`; `simulation_readiness="not_established"`.

## Public usage

```python
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield
from island.forcefields.pcff import (
    bind_pcff_types, assign_typed_pcff_charges, provide_pcff_charges,
)

typed = bind_pcff_types(
    system, source, site_types,
    provenance="independent final-graph/source review",
    evidence_references=["source-row receipt", "chemical review receipt"],
)
# Mode B: native increments, with a versioned policy when fallback terms are used.
charges = assign_typed_pcff_charges(system, typed, resolution_policy=policy)
# Mode C instead calls provide_pcff_charges with an exact finite vector, unit,
# provenance/evidence, component_totals keyed by each component's lowest site ID,
# and total_charge. It never normalizes the vector.
prepared = prepare_forcefield(
    system,
    ForceFieldRequest("pcff", PCFFOptions(
        source_path, (0, 0, 1), (0, 0, 1), resolution_policy=policy,
        typed_graph=typed, graph_charges=charges,
    )),
)
```

`assign_pcff_source_types` remains the legacy checked-explicit API and still
requires automatic perception and agreement. `bind_pcff_types` is the opt-in
external contract. A typed graph records `origin="externally_established"`,
`automatic_perception="not_performed"`, and caller provenance separately. Native
charges record `origin="native_increments"`; provided charges record
`origin="provided"` and `native_increment_availability="not_evaluated"`.
Neither record fabricates an automatic/native wrapper. Model, prepared-facade,
checkpoint, evaluator and bundle identities include the typed graph and charge
record identities, so graph, type or charge edits invalidate dependent records.
Coordinates remain replaceable under the existing evaluator binding contract.
