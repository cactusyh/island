# Phase 4J40-C6 — long-range reference closure

The declared KSPACE executable is available at the recorded path and its
SHA-256 is preserved in the generated manifest. A pairwise Coulomb control
passes with the native kJ mol⁻¹ nm⁻¹ conversion. This toy control is not polymer
PME/Ewald validation.

The retained polymer cells remain blocked: the existing Ewald component
comparison exceeds the unchanged bound and complete long-range virial pressure
is unavailable. The closure record therefore reports `ewald_unvalidated`; no
cutoff recovery or long NPT execution starts. Original C2–C5 evidence remains
untouched.

The comparison contract binds coordinates, box, ordering, charges, exclusions,
1–4 scaling, mixing, cutoffs, reciprocal tolerance, units and precision. No
Ewald/PPPM scientific claim is made from the analytic or pairwise controls.

`production_validated=False`; `simulation_readiness="not_established"`.

The declared control matrix covers neutral point charges, explicit counterion,
pairwise reference, small polymer fragment, box scaling and cutoff/tolerance
scaling. Each comparison requires identical coordinates, box, charges, masses,
ordering, exclusions, 1–4 scaling, LJ mixing, cutoff and reciprocal tolerance;
energy, Coulomb, reciprocal, real-space, force norm, virial and finite-volume
pressure fields are compared. Missing metadata or any residual above the bound
blocks closure.
