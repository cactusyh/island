# Phase 4J22: periodic force-field backends

This phase adds an additive periodic backend layer that consumes a validated
`FinalChemicalGraph`, retained native `PreparedForceField`, and explicit
periodic configuration. It performs no typing, charge generation, parameter
fitting, packing, crosslinking, or dynamics.

## Contracts

`PeriodicForceFieldConfig` uses schema `island_periodic_forcefield_config_v1`.
It binds source SHA, orthorhombic box, PME/Ewald tolerance, cutoff, mixing rule,
1-4 scales, exclusion policy, units, provenance, and evidence. The four family
adapters retain their native conventions: PCFF sixth-power 9-6 Class-II rows,
OPLS-AA geometric 12-6/RB rows, and GAFF/GAFF2 Amber Lorentz-Berthelot rows and
source 1-4 policy.

`prepare_periodic_forcefield` validates the graph, box, native source, and
native record against a box-free chemical view, then creates immutable unified
typed, charge, and assignment adapters bound to the periodic graph. Missing
terms, Wilson out-of-plane terms, source-policy mismatch, and unsupported
exclusions fail before backend construction. Historical nonperiodic
`prepare_forcefield` continues to reject periodic graphs.

LAMMPS export is deterministic (`system.data`, `in.periodic`) and records stable
atom IDs, molecule IDs, coordinates, topology, masses, charges, explicit styles,
cutoff, mixing, special bonds, long-range electrostatics, neighbor policy, and
thermo/component output. The writer never invokes `msi2lmp`. OpenMM is loaded
lazily and constructs PME/Ewald periodic systems from owned numerical records.

## Evidence boundary

The retained DP3 polyethylene fixtures under `../island-validation/phase4i2`
construct periodic OpenMM systems for PCFF, OPLS-AA, GAFF, and GAFF2. The
OpenMM Ewald receipt records finite single points, component checks, forces,
finite differences, charges, atom ordering, and box checks for all four
families. Deterministic LAMMPS data/input export, identity validation, mutation
rejection, relocation, bundle reconstruction, and separate-process loading are
separately tested.

The pinned LAMMPS executable is
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`. It lacks
the KSPACE package. Consequently, periodic LAMMPS Ewald/PPPM single-point
execution and any LAMMPS/OpenMM periodic runtime comparison are **blocked**.
The receipt status is exactly `blocked: KSPACE package unavailable`. An
auxiliary KSPACE-enabled executable was built during investigation, but it is
not the pinned executable and is not used to claim J22 LAMMPS acceptance.

A separate PME trial was incomplete/resource-limited and carries no pass claim.
No scientific production validation or readiness claim is made.

Persistence uses the existing native prepared bundle plus checksummed periodic
records and export hashes. Loading validates nested identities and never reruns
PSMILES construction, crosslinking, packing, or typing.

`production_validated=False` and `simulation_readiness="not_established"`.
