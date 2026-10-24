# Phase 4J40-C8 — polymer Ewald discrepancy diagnosis

C8 binds OpenMM PME and KSPACE LAMMPS frame comparisons to the exact graph,
coordinates, box, ordering, masses, charge origin/vector, exclusions, 1–4
scaling, mixing, cutoffs, switching, reciprocal tolerance, units and precision.
It compares total, Coulomb real/reciprocal, LJ, bonded, force, virial and
pressure-proxy components.

The retained C7 polymer comparisons remain outside the unchanged component bound.
The diagnosis matrix records point-charge, counterion, pairwise, fragment,
linear, aromatic, network and scaling cases as blocked pending raw paired
component records. No tolerance is widened and no charge or force-field value is
changed. Mutation of graph, box, charges, exclusions or tolerance is rejected.

Therefore `ewald_unvalidated` and `pressure_unavailable` remain explicit;
cutoff recovery and long NPT remain blocked. Earlier evidence is preserved.

`production_validated=False`; `simulation_readiness="not_established"`.
