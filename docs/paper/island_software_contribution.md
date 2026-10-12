# ISLAND software-methods contribution

ISLAND is a provenance-preserving computational polymer construction and
molecular-simulation tool. The paper studies its reproducible software boundary,
numerical reference checks and explicitly limited physical-validation evidence.
It is not an AI prediction study.

The contribution has eight connected parts:

1. A force-field-neutral final chemical graph survives polymerization,
   crosslinking and packing transformations.
2. Unified TypedGraph and assignment records bind PCFF, OPLS-AA, GAFF and GAFF2
   preparation to one identity model.
3. Parameter and charge provenance is source-bound, checksummed and fail-closed
   when a term or source is unsupported.
4. Deterministic crosslink plans modify the final graph before typing and
   assignment, with immutable site and membership records.
5. Periodic multichain packing remaps stable atom IDs without changing graph
   identity.
6. Periodic backends export LAMMPS records and construct OpenMM Reference systems;
   J23 provides minimization, NVT, NPT, annealing and transactional recovery.
7. Checkpoint publication, relocation and separate-process resume preserve the
   accepted boundary and binary stochastic state.
8. The analysis layer validates mappings, boxes, units and source identities and
   writes checksummed observable records.

Evidence is labelled as `SOFTWARE_VERIFIED`, `NUMERICALLY_VERIFIED`,
`PHYSICALLY_VALIDATED`, `EXPERIMENTALLY_COMPARED`, `BLOCKED` or
`NOT_EVALUATED`. Software-contract tests never become physical claims. The
historical flags remain `production_validated=False` and
`simulation_readiness="not_established"`.
