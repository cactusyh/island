# Phase 4J40-C5 — actual backend pressure and force closure

C5 audits the retained periodic force-field records without relabeling the C2
failures. Component virial pressure is decomposed into kinetic, bonded,
nonbonded and long-range terms. Force-virial and finite isotropic volume-strain
estimators are accepted only when their residual is within a declared bound and
the long-range electrostatic virial has an independent audit.

The retained OpenMM records contain no validated full virial pressure, and the
long-range Ewald component remains discrepant against the independently executed
KSPACE control. Consequently backend pressure remains `pressure_unavailable`,
cross-engine closure remains blocked, and no recovery NPT dynamics execute.
The eleven cutoff failures remain byte-preserved. The execution gate rejects
unsafe boxes before dynamics and never changes cutoff, switching, density or
force-field semantics.

`production_validated=False`; `simulation_readiness="not_established"`.
