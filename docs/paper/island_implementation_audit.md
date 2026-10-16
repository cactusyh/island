# ISLAND implementation audit

**Base:** `origin/main` / J23 `b303efa`  
**Branch:** `codex/phase-4j25-island-software-paper-evidence`

This audit records implementations found in the source tree and their direct
tests. The capability matrix is the machine-readable index; this document gives
the evidence boundary and limitations.

## Construction and graph

PSMILES and sequence construction are implemented in `src/island/builders/` with
linear, alternating, block and random copolymer builders, repeat-unit metadata,
stereochemical/tacticity controls and coordinate generation. Final graph records
and transformations are implemented in `src/island/graph/final.py`; identity is a
canonical digest of graph payload, transformations and membership. Tests include
`tests/test_phase_4j18_final_graph.py`, polymer builder, sequence and
stereochemistry suites.

Crosslink site identification, deterministic planning, graph application and
serialization are in `src/island/crosslinking.py`. `tests/test_phase_4j20_crosslinking.py`
and its correction companion validate membership, capacity, chemical state and
identity invalidation. Crosslinks are applied before force-field typing.

Periodic multichain packing is implemented in `src/island/packing.py` through
immutable plans, transformations, stable-ID remapping and box metadata.
`tests/test_phase_4j21_packing.py` and the correction test cover density,
placement and transactional publication.

## Assignment and backends

`src/island/graph/unified.py` provides unified source, TypedGraph, charge and
assignment records. PCFF source resolution and typed graph binding are in
`src/island/forcefields/pcff/`; OPLS-AA source, typing, charge and parameter
records are in `src/island/forcefields/oplsaa/`; GAFF/GAFF2 use the AmberTools
lineage in `src/island/forcefields/ambertools/`. `src/island/periodic.py` binds
the four periodic backend families and exports LAMMPS data/input plus OpenMM
periodic systems. Missing terms, unsupported domains and source mismatches fail
closed. Coverage is bounded by source profiles and is not universal.

## Dynamics and analysis

`src/island/periodic_md.py` implements OpenMM Reference minimization, NVT, NPT and
NVT annealing with signed manifests, accepted generation records, binary
checkpoints, relocation and child-process resume. `tests/test_phase_4j23_periodic_md.py`
contains split/uninterrupted and checkpoint-byte comparisons. The current J23
workflow intentionally retains every short frame in memory.

`src/island/analysis/mdanalysis.py` validates J23 or explicit MDAnalysis metadata,
stable maps, boxes and hashes, and writes immutable analysis records. Linear
chains are unwrapped on explicit bonded paths with orthorhombic minimum-image
displacements; branched, crosslinked and triclinic cases are rejected. Re/Rg is
calculated per chain before equal-weight aggregation. Target and measured
temperature/pressure fields are distinct. `src/island/trajectory_store.py` is an
opt-in bounded-memory v2 store with staged publication, a single-writer lock,
manifest-declared restart history and hard checkpoint intervals; valid v1 stores
remain readable. The J23 integration is an opt-in export of completed workflows,
not an online writer attached to every MD step.

J28 adds `AnalysisSelection` (`island.analysis-selection.v1`) as the graph-bound
selection contract. The verified loader retains final-graph bonds and molecule
membership, and supported explicit-H acyclic linear molecules use full declared
molecular atoms for mass-weighted Rg while endpoint-based Re/orientation use the
ordered heavy-atom path. Selection identities are stored in analysis records and
must match a trusted trajectory or expected selection at reconstruction. This is
analytical/software integration evidence only; it does not validate equilibrium
polymer properties.

J29 validates the complete retained J23 path with the real seven-frame PCFF
explicit-H fixture. The CLI produces Rg, Re, Re/Rg, gyration tensor, chain
anisotropy, orientation, nematic order, density, and connectivity records;
independent first-frame controls and same-process relocation are covered by the
acceptance tests; a fresh Python subprocess loads the relocated record and
compares the canonical payload in `tests/test_phase_4j29_real_acceptance.py`.
The reusable control script is
`scripts/validate_phase_4j29_real_workflow.py` and the acceptance tests are
`tests/test_phase_4j29_real_acceptance.py`. This validates software contracts
and numerical controls while production trajectory and physical equilibrium
gates remain blocked.

J30 adds a bounded two-chain explicit-H acceptance fixture with interleaved
stable IDs and unequal chain lengths. Per-chain bonded-path geometry is
aggregated by an explicit equal-chain or mass-weighted policy, with translation
invariance, independent first-frame controls, and relocated separate-process
record reconstruction. This remains software and numerical evidence; physical
equilibrium and production trajectory evidence are unavailable.

The J29 receipt records `validated_commit`, the exact implementation commit
whose code was executed. A separate validator checks the receipt and manifest
on disk, including schema, branch, base, SHA formatting, readiness flags, and
mutated-file rejection. This field is distinct from the corrective branch tip.

## Evidence boundary

J17-J23 receipts and identities are preserved. OpenMM Reference numerical and
workflow checks are verified within their recorded contract. The pinned LAMMPS
executable lacks KSPACE, so Ewald/PPPM runtime comparison is blocked. Short
Reference runs do not establish equilibrium or production readiness. No Tg or
mechanical claim is made.
