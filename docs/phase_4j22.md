# Phase 4J22: periodic force-field backends

This phase adds `island.periodic`, an additive backend layer that consumes a
validated `FinalChemicalGraph`, a retained native `PreparedForceField`, and an
explicit periodic configuration. It performs no typing, charge generation,
parameter fitting, packing, crosslinking, or dynamics.

## Contracts

`PeriodicForceFieldConfig` is schema `island_periodic_forcefield_config_v1`.
It binds source SHA, orthorhombic box, PME/Ewald tolerance, cutoff, mixing rule,
1-4 scales, exclusion policy, units, provenance and evidence. Identities are
checksummed JSON. The four families retain their native conventions: PCFF
sixth-power 9-6 Class-II rows, OPLS-AA geometric 12-6/RB rows, and GAFF/GAFF2
Amber Lorentz-Berthelot rows and source 1-4 policy.

`prepare_periodic_forcefield` validates the graph, box, native source and native
record against a box-free chemical view, then creates immutable unified typed,
charge and assignment adapters bound to the periodic graph. Missing terms,
Wilson out-of-plane terms, source-policy mismatch and unsupported exclusions
fail before backend construction. Historical nonperiodic `prepare_forcefield`
continues to reject periodic graphs.

LAMMPS export is deterministic (`system.data`, `in.periodic`) and records stable
atom IDs, molecule IDs, coordinates, topology, masses, charges, explicit styles,
cutoff, mixing, special bonds, long-range electrostatics, neighbor policy and
thermo/component output. The writer never invokes msi2lmp. OpenMM is loaded
lazily and constructs PME/Ewald periodic systems from owned numerical records;
absence of OpenMM is reported as an optional dependency error.

## Retained fixtures and verification

The retained DP3 polyethylene bundles under `../island-validation/phase4i2` were
used for PCFF, OPLS-AA, GAFF and GAFF2. All four produced finite OpenMM Reference
single-point energies and deterministic exports. The pinned LAMMPS executable
was present, but this binary lacks its KSPACE package, so the requested PPPM
single-point comparison was preserved as an unavailable dependency result rather
than converted to a pass. No production or simulation readiness claim is made.

Persistence uses the existing native prepared bundle plus checksummed periodic
records and export hashes. Loading validates all nested identities and never
reruns PSMILES construction, crosslinking, packing or typing.

`production_validated=False` and `simulation_readiness="not_established"`.
