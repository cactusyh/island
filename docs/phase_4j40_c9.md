# Phase 4J40-C9 — executed PME/KSPACE root-cause isolation

C9 records the complete ordered diagnosis matrix and controlled one-variable
experiments. The first failing system is retained from the C8 raw comparison;
its frame identity contract includes graph, coordinates, box, ordering, masses,
charges and origin, exclusions, 1–4 scaling, cutoffs, reciprocal tolerance, units
and precision. Candidate variables are charge conversion, mapping, exclusions,
1–4, box convention, cutoff, reciprocal tolerance, precision, periodic images
and network topology.

No concrete implementation mismatch was demonstrated in the retained evidence,
so no repair or tolerance change is applied. The discrepancy remains an
unsupported/non-equivalent closure case: `ewald_unvalidated` and
`pressure_unavailable`. Cutoff recovery and long NPT remain blocked. Earlier
C2–C8 evidence is preserved byte-for-byte.

`production_validated=False`; `simulation_readiness="not_established"`.
