# Phase 4J40-C10 — executed polymer Ewald closure decision

C10 executes the retained numerical comparison. The smallest retained polymer
comparison is the GAFF2 linear PE frame. OpenMM and KSPACE agree on bonded,
angle, proper, improper and LJ components to roundoff, while Coulomb differs by
`0.042948678959815245 kJ/mol`; total energy differs by
`0.04294867896004462 kJ/mol` and maximum force difference is
`0.0419290025488408 kJ/(mol nm)`. The exact graph/configuration identity and
component values are in `actual-decision.json`.

This is a real numerical failure, not a metadata-only contract result. The
retained configuration uses a 10.81388814433968 Å cubic box, 5 Å cutoff, no
switching, Ewald tolerance 1e-5, native shortest-path exclusions and GAFF2 with
exact RESP/Psi4 provided charges. No parameter sweep is selected as a repair:
there is no demonstrated one-variable implementation mismatch, so no values,
tolerances or charges are changed.

The decision is an unresolved backend/convention limitation. Ewald remains
`ewald_unvalidated`, backend pressure remains `pressure_unavailable`, cutoff
recovery and long NPT remain blocked, and all prior evidence is preserved.

`production_validated=False`; `simulation_readiness="not_established"`.
