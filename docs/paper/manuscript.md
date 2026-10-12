# ISLAND: provenance-preserving polymer construction, force-field assignment and periodic molecular simulation

## Abstract

Polymer simulation workflows join chemical construction, force-field assignment,
periodic packing, molecular dynamics and post-processing, but their intermediate
identities and failure boundaries are often difficult to reconstruct. ISLAND is a
computational tool that preserves a force-field-neutral final chemical graph and
stable atom identities across polymerization, stereochemical construction,
crosslinking, periodic packing, source-bound parameter assignment and periodic
simulation. Unified TypedGraph, charge and assignment records represent PCFF,
OPLS-AA, GAFF and GAFF2 preparation with explicit source identities and
unsupported states. OpenMM Reference workflows provide periodic minimization,
NVT, NPT and annealing with transactional checkpoints, relocation and exact
split/uninterrupted software checks. An additive MDAnalysis-compatible layer and
chunked trajectory store bind observable records and restart boundaries to graph,
mapping, source and file hashes. Current evidence establishes software contracts
and bounded numerical behavior. Physical production validation remains open;
`production_validated=False` and `simulation_readiness="not_established"`.

## 1. Introduction

Polymer simulation combines decisions about repeat-unit representation, sequence,
stereochemistry, end groups, crosslink topology, packing, charge assignment,
parameter sources and periodic protocols. Reproducibility requires those decisions
to remain inspectable after coordinates and file formats change. ISLAND addresses
this need through immutable records and explicit failure states.

## 2. Design goals and architecture

The architecture keeps the final chemical graph independent of a force-field
family, binds every family-specific assignment to the graph and source, and treats
checkpoint publication as a transactional boundary. Coordinate changes do not
change typing or graph identities. Missing terms and unsupported domains reject.

## 3. PSMILES-to-final-graph construction

PSMILES builders support linear polymers, sequence-defined copolymers and declared
stereochemical/tacticity controls. Stable atom IDs, repeat membership, formal
charge, bond order and construction provenance enter a final graph whose digest
changes on chemical mutation.

## 4. Crosslink-network generation

Reactive-site rules identify candidates on the final graph. Deterministic seeded
planning records selected pairs, capacity and membership. Applying a plan changes
the chemical graph before typing and parameter assignment.

## 5. Periodic multichain packing

Immutable packing plans place one or more molecular units in an orthorhombic
periodic cell. Stable atom-ID remapping and packing identity are retained, and
packing changes coordinates and membership without regenerating chemistry.

## 6. Unified force-field assignment

Unified source, TypedGraph, charge and parameter-assignment records cover PCFF,
OPLS-AA, GAFF and GAFF2. Each record binds the final graph, source hash, profile,
charge origin and diagnostics. Coverage is source-profile specific.

## 7. Periodic MD and checkpoint recovery

J23 constructs periodic OpenMM Reference systems and supports bounded
minimization, NVT, NPT and NVT annealing. Signed generations retain binary OpenMM
state, box vectors and accepted boundaries. Relocation and child-process resume
reconstruct the exact numerical model. A new opt-in store separates trajectory
chunks, thermodynamic samples and durable restart checkpoints for longer runs
without changing J23 schemas.

## 8. Trajectory analysis

Analysis loads a verified J23 workflow or an external MDAnalysis trajectory with
explicit identities and stable mapping. Density uses traceable mass and the actual
frame box. Rg, end-to-end distance, Re/Rg, gyration tensor, anisotropy, RDF and
MSD have declared units and failure states. Thermostat/barostat targets are stored
separately from measured temperature and pressure.

## 9. Software and numerical validation

Graph, crosslink, packing, typing, backend, checkpoint and analysis tests establish
software contracts. Retained OpenMM Reference cases provide bounded numerical
verification, including split/uninterrupted frame and checkpoint equality. The
current LAMMPS executable lacks KSPACE, so Ewald/PPPM runtime agreement is blocked.

## 10. Representative physical validation

The first physical targets are density convergence, energy/temperature/pressure
stability, RDF, chain dimensions, anisotropy, topology preservation and replica
agreement. No converged physical cohort is currently included, so this section is
a protocol and evidence boundary rather than a property result.

## 11. Comparison with existing tools

RadonPy provides a methodological reference for automated amorphous-polymer MD.
msi2lmp, Foyer and AmberTools provide complementary conversion, typing and
parameterization capabilities. ISLAND's evaluated contribution is its graph- and
source-bound record chain and transactional recovery. Unevaluated comparisons are
reported as such.

## 12. Limitations

Universal PCFF/OPLS/GAFF coverage, complete msi2lmp parity, KSPACE-enabled LAMMPS
comparison, triclinic cells, NPT annealing and physical production readiness are
not established. Tg and mechanical properties require dedicated protocols.

## 13. Reproducibility and availability

Machine-readable capability, validation and reproducibility manifests identify
APIs, tests, sources, commands and blocked gates. Figures and tables are generated
from committed scripts and records.

## 14. Conclusions

ISLAND provides an auditable software path from polymer representation through
graph transformations, periodic parameterization, molecular dynamics, restart and
analysis. The current contribution is software-methods evidence with bounded
numerical checks; independent physical validation remains required.
