# Phase 4J33: periodic NPT density diagnostics

J33 adds a bounded, reproducible Reference-platform NPT path on the complete
multichain graphs accepted by J32.  `island.analysis.npt.NPTConfiguration` is
immutable and identity-bound to the final graph, packing plan, source and box.
`run_npt_workflow` executes minimization, NVT pre-equilibration, an actual
`MonteCarloBarostat` NPT segment, and an optional post-NPT segment.  Every
retained frame stores coordinates, velocities, orthorhombic box vectors,
energies, kinetic temperature and opaque checkpoint bytes.

Density is calculated as total atomic mass divided by the actual retained box
volume.  Target pressure and temperature are configuration values and are
never substituted for measured values.  A backend that cannot provide a
pressure observation records `pressure_observation_status=unavailable`.
Finite-run block summaries use four declared blocks when enough frames exist;
otherwise the result is `insufficient_sampling`.  A fixed box reports
`density_relaxation_status=not_demonstrated`.  No status establishes
thermodynamic equilibrium, experimental density, material properties or
production readiness.

The on-disk manifest, frame record and analysis record each carry a SHA-256
identity and bind the backend, graph, packing, source, NPT configuration,
engine, mass and frame-content identities.  `load_npt_workflow` validates these
records without importing OpenMM, so a relocated bundle and a fresh process
can reproduce the same serialized payload.  `compare_npt_workflows` performs
an exact frame and analysis comparison for split/uninterrupted checks.

The acceptance driver is `scripts/validate_phase_4j33.py`.  Missing pinned
sources, OpenMM or AmberTools are recorded as `asset_unavailable` and leave a
required case unsuccessful.  The retained PCFF typing rejection remains an
expected blocked negative control; no PCFF coverage is inferred from it.

J33 preserves the publication boundary:

```text
production_validated=False
simulation_readiness="not_established"
```

The workflow is a software acceptance fixture and does not establish
equilibrium or any experimental or production property.
