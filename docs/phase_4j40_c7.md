# Phase 4J40-C7 — executed KSPACE reference closure

The declared LAMMPS KSPACE executable was inventoried and executed for the
existing Ewald reference controls. Its exact version, build text, SHA-256,
platform and precision are recorded in the C7 audit. The control matrix covers
neutral charges, counterion, pairwise, polymer fragment, linear, aromatic,
network, box, cutoff and reciprocal-tolerance cases.

The retained raw comparison still has an unresolved polymer Ewald component
 discrepancy. The exact input identity contract includes graph/atom ordering,
coordinates, box, masses, charges and origin, exclusions, 1–4 scaling, LJ
mixing, cutoff/switching, reciprocal tolerance, dielectric, units and boundary
conditions. Since component-level closure and force-field virial pressure do not
pass, the C7 result remains `ewald_unvalidated` and `pressure_unavailable`.

No cutoff recovery or long NPT execution is performed. Earlier C2–C6 failed
histories and hashes remain preserved. This phase makes no PPPM, equilibrium,
experimental or production-readiness claim.

`production_validated=False`; `simulation_readiness="not_established"`.
