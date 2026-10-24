# Phase 4J40-C11 — bounded GAFF2-PE discrepancy characterization

The retained first-failing GAFF2 linear PE frame has 52 atoms and a
10.81388814433968 Å cubic box. OpenMM PME total energy is 145.37771485907777
kJ/mol and Coulomb energy is -0.9164685482900288 kJ/mol. The absolute total
error is 0.04294867896004462 kJ/mol (relative 0.0002954282160899077).
The maximum force discrepancy is 0.00419290025488408 kJ/(mol nm). Charges span
-0.21333213061141376 to 0.10246676635213643 e; cutoff is 5 Å and reciprocal
tolerance is 1e-5. The complete numeric record and source identity are in
`characterization.json`.

The bounded sweep dimensions are reciprocal tolerance, cutoff, mesh,
interpolation order, Ewald parameter, FFT precision and box size. No sweep value
is selected as a repair: no independent evidence demonstrates a concrete
implementation mismatch or all-system closure. RMS/percentile force errors are
not reconstructed without an independently persisted OpenMM per-atom force
vector and remain explicitly unavailable.

The result is a bounded-difference candidate, not bit-equivalent engines.
`ewald_unvalidated` and `pressure_unavailable` remain; cutoff recovery and long
NPT remain blocked. Readiness flags remain conservative.
