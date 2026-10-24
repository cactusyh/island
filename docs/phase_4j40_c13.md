# Phase 4J40-C13: bounded numerical characterization and pressure boundary

C13 executes the corrected C12 matrix again and adds fixed, one-variable
PME/PPPM sensitivity studies and actual-polymer finite-volume pressure checks.
The outcome is **unresolved_backend_discrepancy**, **ewald_unvalidated** and
**pressure_unavailable**. No qualified numerical-equivalence policy is established.
`production_validated=False`; `simulation_readiness="not_established"`.

The comparison preserves the C12 distinction between native KSPACE LAMMPS and
LAMMPS plus an explicit analytical special-pair completion. This completion is
a single-frame composite reference, not a native LAMMPS implementation repair.
All native outputs and completion terms are retained separately. Graphs, charges,
charge origins, force-field parameters, exclusions and special-pair scaling are
unchanged. Aromatic RESP/Psi4 charges remain separately labeled. Historical
receipts and source files are unchanged.

## Executed scope

The fresh C12 rerun includes point charges, counterions, pairwise controls, a
complete retained C8H18 component, linear GAFF2 PE, OPLS PEO, OPLS-compatible
aromatic RESP/Psi4 and the packed OPLS network. Original, 1.5-times strained and
translated boxes are compared using Ewald and PME/PPPM. The unchanged C11
study adds cutoff/tolerance, box, retained replica/frame and volume-strain
controls. These are single-point evaluations, not newly equilibrated replicas.

The C13 extension declares its plan before execution: baseline mesh 32³,
LAMMPS order 5, tolerance 1e-7, alpha 8 nm^-1, cutoff 0.5 nm. Each of mesh
(24³/48³), LAMMPS order (4/6), tolerance (1e-5/1e-9), alpha (6/10 nm^-1)
and cutoff (0.4/0.48 nm) varies independently, with original and translated PE
frames. OpenMM interpolation order is fixed at 5. Fixing both mesh and alpha
makes the nominal tolerance changes numerically inactive; this is recorded,
not presented as evidence of tolerance convergence. The original C11 automatic
settings also supply independent tolerance studies. Single/mixed GPU precision
and alternative FFT builds are unavailable. The executable is the pinned
KSPACE double-FFT LAMMPS, SHA-256
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`;
OpenMM uses Reference double. Exact build, hardware and package manifests are
in the validation archive.

There are 108 characterized paired comparisons and 12 new two-sided actual-
polymer pressure derivatives. Raw outputs, XML, input decks, coordinates,
boxes, atom IDs, energy groups, forces and native LAMMPS pressure tensors are
in `artifacts.tar.xz`. `summary.json` records absolute/relative total and Coulomb
energy errors, all physical component errors, RMS/max/percentile force errors,
reference force scales, settings, mesh, atom count and box dimensions per row.
Real/reciprocal differences are descriptive: the native alpha and exception
partition differ, and completion is stored separately. Equality of those native
pieces is not a valid independent observable without aligned decomposition.

## Numerical outcome and policy

The candidate versioned record retains the existing energy bound of
1e-5 kJ/mol per total and physical component, and force bound of
1e-4 kJ/(mol nm), rtol 2e-10. It has **no pressure bound**. These are inherited
acceptance bounds, not newly justified scientific uncertainty bounds.
Thirteen of 108 comparisons pass the inherited energy/force checks. Passing a
single row does not qualify a backend or polymer workflow.

Representative unstrained direct-Ewald total-energy residuals are
1.90109e-5 kJ/mol (GAFF2 PE), 1.05628e-5 (aromatic) and 3.23410e-5 (network).
The PE fixed-mesh 24³ to 48³ sweep reduces the untranslated energy residual
from 0.00426110 to 0.0000984300 kJ/mol and RMS component force residual from
0.0119541 to 0.000300707 kJ/(mol nm). Translated 48³ residuals are
0.0000449864 kJ/mol and 0.000291721 kJ/(mol nm). The remaining residuals and
translation sensitivity do not justify promoting this to bounded numerical
equivalence. No force-field or numerical acceptance tolerance is widened.

The policy vocabulary distinguishes bitwise equality, roundoff-level equality,
bounded numerical equivalence and unresolved backend discrepancy. C13 can
report bitwise equality only on identical numerical observables; it has no
qualified roundoff or bounded-equivalence rule. All observed nonidentical rows
remain unresolved. The candidate record, summary and artifacts are hashed;
mutation of settings or numerical evidence invalidates offline reconstruction.

## Pressure and dynamics boundary

For each actual polymer, C13 executes central isotropic volume derivatives at
relative volume steps 1e-4, 1e-5 and 1e-6, scaling coordinates and box together:
`P_config = -dU/dV`. The conversion is
`1000/(N_A * 1e-27 * 1e5)` bar per kJ mol^-1 nm^-3. Velocities are zero;
these are configurational single-point pressures with zero kinetic term,
not observed thermodynamic pressures.

At step 1e-5 the OpenMM/composite finite-volume pressures (bar) are:

| system | OpenMM | composite LAMMPS | difference |
|---|---:|---:|---:|
| PE | 3303.362373 | 3303.362592 | -0.000219 |
| PEO | 4860.494266 | 4860.495280 | -0.001014 |
| aromatic RESP/Psi4 | 9379.458259 | 9379.459827 | -0.001568 |
| network | 5533.652039 | 5533.653602 | -0.001563 |

These high-pressure nonequilibrium frames do not establish pressure acceptance.
The network derivative is step-sensitive. A full independent OpenMM periodic
reciprocal force/virial tensor is not measured in this path. Therefore tensor
virial error, tensor pressure error and qualified scalar pressure error are
explicitly unavailable (`null` with diagnostic), rather than reported as zero.
The scalar finite-volume discrepancy is recorded separately. The independent
XML audit validates only the special-pair virial increment. Wrapped-position
`sum(r outer F)` is deliberately not substituted for the missing full periodic
reciprocal virial. Full kinetic, constraint and long-range pressure validation
remains required. No pressure policy is fitted from these results.

Cutoff recovery, long NVT/NPT, new replicas, equilibration diagnostics and exact
MD continuation are not executed because their prerequisite gates fail. C13
makes no equilibrium, experimental-agreement, universal force-field, PCFF,
AM1-BCC or production-readiness claim.

## Reproduction and validation

From the repository root, with `PYTHONPATH=.:src`, `OPENMM_CPU_THREADS=1` and
`OMP_NUM_THREADS=1`, run:

```bash
../island-validation/phase4e2-env/bin/python -m scripts.phase_4j40_c12_study --source /tmp/j40-validation/final/run --reference /tmp/j40-c2-validation/run/reference --lammps ../island-validation/phase4j22/lmp-kspace --output /tmp/j40-c13-reproduce/study
../island-validation/phase4e2-env/bin/python -m scripts.phase_4j40_c13_study --run /tmp/j40-c13-reproduce/study --lammps ../island-validation/phase4j22/lmp-kspace
../island-validation/phase4e2-env/bin/python -m scripts.phase_4j40_c13_audit /tmp/j40-c13-reproduce/study /tmp/j40-c13-reproduce/audit.json
../island-validation/phase4e2-env/bin/python -m scripts.validate_phase_4j40_c13_receipt docs/evidence/phase_4j40_c13
```

Both study commands exit 1 intentionally for blocked scientific acceptance;
run them separately. The portable C12 input reconstruction helper and raw
source originals are included in the archives for environments lacking the
original `/tmp` layout. The offline receipt validator exits 0 for intact
blocked evidence and reports a separate scientific acceptance exit code of 1.
Independent audits reconstruct raw residuals, force units, atom mapping,
special-pair completion and volume derivatives without using the production
comparison arithmetic. Exact commands, exit codes, suite counts, environment
manifests and historical byte hashes are recorded in the receipt and archives.
