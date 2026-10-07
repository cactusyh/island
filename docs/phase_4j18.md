# Phase 4J18 — unified final-graph preparation

## Base and contract

This phase is based on J17 commit `ee5ce274905704191ebf1f90fe76139cb12e3a30`.
The force-field-neutral contract is implemented in `island.graph`:

- `FinalChemicalGraph` (`island_final_chemical_graph_v1`) stores stable integer
  atom IDs, element/formal-charge/aromatic state, bond orders/connectivity,
  connected components, molecule membership, polymer/crosslink provenance,
  periodic lengths/boundaries, and separate chemical/history identities.
- Coordinates are stored in the system bundle and are deliberately excluded from
  graph identity. Coordinate changes therefore follow the existing evaluator
  binding rules.
- `GraphTransformation` (`island_final_graph_transformation_v1`) binds an
  operation, input/output graph identity, deterministic parameters and seed,
  provenance, evidence, and its own identity.
- `FinalGraphChargeRecord` (`island_force_neutral_charges_v1`) binds a graph,
  force-field family, source identity, charge method/origin, exact finite vector,
  units, provenance and evidence. It never normalizes or redistributes charges.

Public transformations are `build_psmiles_graph`, `complete_end_groups`,
`create_crosslinks`, `construct_periodic_box`, and `prepare_final_graph`.
`save_final_graph_bundle` and `load_final_graph_bundle` provide a relocatable,
checksum-verified graph bundle that can be reconstructed in a separate process
without RDKit, OpenMM, SciPy, Foyer or ParmEd.

The J17 `PCFFTypedGraph` and `PCFFGraphCharges` records are unchanged. The named
`charge_record_from_pcff` adapter checks their chemical graph against the neutral
final graph and preserves `native_increments` versus `provided` origin. Named
OPLS-AA and GAFF/GAFF2 adapters preserve their source identity and charge method;
no force-field source or charge vector is borrowed across families.

`ForceFieldRequest(final_graph=...)` validates that the requested force-field
preparation is bound to the exact final graph. Current PCFF, OPLS-AA and
GAFF/GAFF2 evaluators remain finite nonperiodic backends. A periodic final graph
therefore fails with a precise unsupported-periodic diagnostic before typing or
parameterization. A nonperiodic crosslink may proceed only when the selected
backend and pinned source contain every required interaction; missing terms still
reject publication.

## Frozen acceptance fixtures

- PSMILES `[*:1]CC[*:2]`, deterministic seed 2026, DP3: 20 sites.
- The same PSMILES, DP10: 62 sites with the same deterministic construction
  contract.
- A controlled crosslink removes one explicitly selected terminal hydrogen at
  each target, then adds one explicit crosslink bond. No implicit valence repair
  occurs. Invalid valence, duplicate bonds and existing-bond duplicates reject.
- An orthorhombic `(40, 40, 40)` periodic box with `(True, True, True)`
  boundaries is recorded with its seed and survives relocation.
- The J17 PCFF external typed/native and provided-charge path remains the
  executable nonperiodic linear slice. The periodic graph is accepted by the
  neutral graph/bundle layer and rejected by current force-field backends with
  the declared limitation.

Acceptance is software contract evidence. It does not establish scientific
crosslinked or periodic MD, full `.frc` coverage, msi2lmp parity, or universal
OPLS/GAFF parameter support. `production_validated=False` and
`simulation_readiness="not_established"` remain required.

## Validation and preserved limits

Graph validation rejects wrong IDs, missing/extra/duplicate bonds, invalid bond
orders and valence, duplicate crosslinks, changed box or boundary metadata,
changed formal charge or molecule/repeat/crosslink history, nonfinite charges,
source mismatch, stale checksums, and stale graph/bundle identities. Bundle load
never regenerates PSMILES, crosslinks or force-field typing.

The controlled normal-force numerical comparisons, fresh minimization,
relocation, workflow continuation and exact RNG checks remain the J17 evidence;
J18 adds neutral periodic/crosslink persistence only. Full-source and full-polymer
acceptance commands remain nonzero unless every declared interaction is
independently verified. No missing PCFF/OPLS/GAFF parameters are invented and no
Class-II missing term is silently set to zero.
