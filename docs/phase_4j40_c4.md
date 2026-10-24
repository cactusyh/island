# Phase 4J40-C4 — validated pressure controls and safe NPT recovery

C4 adds independent controls for the virial pressure estimator and an explicit,
fail-closed recovery record. The estimator is `P=(2K+W)/(3V)` in kJ mol⁻¹ nm⁻³,
converted to bar. Controls include ideal gas, finite-volume `-dU/dV`, sign,
temperature and volume accounting. The controls do not imply that OpenMM's
retained C2 histories expose a complete virial, especially for reciprocal-space
or constrained systems.

The eleven C2 cutoff failures remain preserved. A recovery transformation may
only be a recorded repack/box transformation preserving final graph,
parameters and charge origin. The included recovery tool only writes and
preflights such a record; it intentionally executes no dynamics until pressure
and all cross-engine gates are independently validated.

Pressure remains `pressure_unavailable` for the retained backend, and scientific
acceptance remains blocked. No Ewald/PPPM claim, equilibrium claim or production
readiness claim is made.

`production_validated=False`; `simulation_readiness="not_established"`.
