# Phase 4J58 — executed QM pilot and Class-II identifiability

J58 executes a bounded **Psi4 1.9.1** pilot and imports real energies,
Cartesian gradients/forces and central Cartesian Hessians. It does not create a
PCFF-PDMS profile, fit a coefficient, or run PCFF MD. Parameter completeness
remains `incomplete` and scientific validation remains `not_established`.

## Base and provenance

The fetched exact `origin/main`, and actual branch base, is
`f471e06894275b5b6485461d2bc1dcd594f7dfa7`. Ancestry checks establish that main
already contains J57 `b56659d13d5e0502ebf57a7c0fc56bdc1fdb5df2` and J56
`24ad7698f64d62aeb08ed893eb8372bdbbf33819` / `1c3de73` (the full SHA is in
`plan/ancestry.json`). No prerequisite is assumed merged or cherry-picked, and
this phase does not merge main.

All changes are additive in `island.pcff_qm_pilot` and new J58 files. The entire
base tree, including historical implementation identities, receipts, bundles and
literal symlink contents, is the preservation baseline.

The 28 `author_confirmed_audited` records remain **user-attested author
confirmations**, with original author artifacts `not_yet_supplied`. They retain
`evidence_kind="user_attested_author_confirmation"` and
`authority="user_reported_author_confirmation"`. No author email, document,
author date, source-file hash or missing coefficient is invented. Confirmation
of interpretation does not supply an absent numerical coefficient.

## Predeclared pilot and engine

Python imports and executable discovery were tested in all declared environments.
Psi4 imports and runs in the existing `phase4j36-native-env`; it does not expose a
standalone command on PATH. PySCF is unavailable in the inspected environments.
No environment was modified. The other four engine adapters remain unexecuted.

The pilot selects these two J57 directed queries:

| Query | Exact fixture | Atoms | Source types | Charge / multiplicity |
| --- | --- | ---: | --- | --- |
| `angle-angle:hsi,sio,hsi,hsi` | Disiloxane | 9 | `sio`, `osi`, `hsi` | 0 / 1 |
| `bond-bond_1_3:sio,osi,sio,hsi` | Disiloxane | 9 | `sio`, `osi`, `hsi` | 0 / 1 |

Both bind graph
`7aa04b55c75e48d32abd159a51e6f6e39fe51244521118c94df49484ba1a738f`,
the original stable atom IDs/order, source typing, masses and unmodified candidate
charge vector. This is the smallest authorized exact fixture for both queries.
Its Si–H control does **not** establish transferable Si–H or methyl-PDMS support.
HMDS and DP2/DP4/DP6 are not subjected to a new QM campaign here.

The original plan, written before execution, permits 50 scan energy/gradient
jobs, two central Hessians, six sensitivity gradients and up to eight
finite-difference energies. Limits are two simultaneous jobs, two threads and
2 GiB per job, 900 seconds per job and 14,400 seconds campaign wall time. Each
Hessian has a separate internal-displacement budget. Only four external
finite-difference energy jobs are needed. The 889 historical prepared geometries
are not executed.

The explicit proposed protocol is isolated neutral singlet PBE0/def2-SVP without
dispersion, PK SCF, 1e-10 hartree energy and 1e-8 density thresholds, maximum 200
SCF cycles, integral threshold 1e-12 and an unpruned Treutler 99×590 DFT grid.
This is an exploratory backend-selected protocol, not author authority or a
validated final parameterization method. Basis files, interpreter, binary,
dependencies, LibXC version, precision and platform are hash-bound in `engine/`.
The downloaded official Psi4 v1.9.1 driver matches the installed driver
byte-for-byte; documentation/source locations for units and API calls are retained.

## Geometry and actual execution

The historical starting geometry has no clash or valence inconsistency, but its
Si–O bonds are approximately 1.70–1.72 Å and it is not an established minimum.
J58 preserves it and creates a new optimized-center scan revision. The first
optimization reaches 30 iterations: tight force criteria pass, displacement
criteria do not. Its native failure and all output remain retained. A separately
recorded, pre-execution budget amendment permits one Cartesian restart of at most
10 iterations, with **unchanged tolerances**; that restart converges.

The new center has Si–O ≈1.6554 Å. Both existing 5×5 displacement grids are used:
±0.08 radians for the angle pair and ±0.04 Å for the first/third bond pair.
`geometry-audit.json` reports every separation, bond length, achieved coordinate,
source-equilibrium comparison and new target identity. The J57 intake bridge
regenerates a target for the new exact center; it never reuses a historical target
identity. Its old generator reference label is explicitly superseded by the J58
center record. No historical geometry is edited.

All 50 scan points retain actual energy and full Cartesian gradient. Both central
Hessians use Psi4's three-point 0.005-bohr finite differences: **49 actual gradients
per Hessian**, including reference, on 24 translation-free Cartesian SALCs.
Psi4 reconstructs a full 27×27 Cartesian Hessian; rotations are not projected by
the engine. J58 projects all six rigid translation/rotation directions separately
for stability analysis. An unprojected molecular Hessian is never required to be
strictly positive definite. These two scans share the same optimized center;
their central calculations are explicit repeated calculations.

The original Hessian adapter saves both numerical matrices successfully, then
fails while calling `functional()` on Psi4's base finite-difference Wavefunction.
The original `result.json` failures remain untouched. The additive
`recovery.py` version independently reads the 49 native QCSchema blocks,
convergence markers, reference energy and completed Hessian assembly. It writes
hash-bound `metadata-recovery.json` sidecars. **No Hessian is rerun and no missing
numerical value is synthesized.** The workflow reports calculation completion
and adapter metadata-export failure separately. Two preflight launcher/validation
failures that performed no QM are also retained in validation logs.

`native-outputs.tar.xz` preserves exact inputs, native logs, indexed native
Matrix.save text, lossless no-pickle matrices, failed attempts, recovery sidecars
and execution logs. Large scratch/caches stay outside Git. `reference-datasets.tar.xz`
contains separate normalized samples, raw verification bindings and strict J57
intake results. These are numerical data, not templates or fitted parameters. Normalization v2
corrects an explanatory relative-constant estimate in v1; the actual constants and
every numerical sample remain unchanged. The v1 archive and revision/hash mapping
are preserved for audit. No QM value is rerun or replaced.

## Import integrity and numerical checks

The executed-engine importer checks native/input hashes, graph bindings, atom
order, printed and full-precision geometry, charge/spin, method/basis, SCF
convergence/settings, grid, energy, derivative shapes, Hessian symmetry/axis order
and native indexed ASCII agreement. Forces are **minus gradients**. The strict
J57 intake accepts each complete 25-point dataset, including its central Hessian,
and rejects incomplete coverage, changed graph/order/charge, missing units,
nonfinite/incomplete forces, inconsistent methods and altered files.

J57 observable normalization retains its explicit CODATA-2022/SI hartree, Bohr
and thermochemical-calorie formulas. Psi4/QCElemental geometry uses Bohr =
0.52917721067 Å, whereas the J57 normalization constant is 0.529177210544 Å.
Their relative difference, about 2.38e-10, is recorded. Native geometry is verified
with the engine constant; raw derivatives and J57 factors remain preserved. This
small representational difference is not hidden or used to tune agreement.

`independent_phase_4j58.py` imports no production importer. It independently parses
native indexed text, reconstructs all 50 normalized energies/force vectors and
both Hessians with Decimal/SI conversions, and checks a real QM energy/gradient
central difference at **0.001 and 0.0005 Å**. The force check's 0.01
kcal/mol/Å bound was declared before execution. All positive/negative energies
and sign-error comparisons are retained.

The sensitivity study executes PBE0/def2-TZVP and B3LYP/def2-SVP at the common
center, an angle corner and a bond corner. Relative deformation energies and
full-gradient differences are reported separately from large absolute electronic
energy offsets. Successful SCF convergence alone does not establish final QM
adequacy. The angle-corner deformation energy changes by 0.07995 kcal/mol with the larger
basis; the bond-corner energy changes by **0.73984 kcal/mol** (pilot 0.87629,
larger basis 1.61612 kcal/mol). B3LYP/def2-SVP changes those relative energies by
0.00095 and -0.11076 kcal/mol. Center gradient RMS differences are 3.8769
kcal/mol/Å for the larger basis and 1.1269 kcal/mol/Å for B3LYP. This supports
continued workflow/identifiability investigation, **not continued def2-SVP fitting
as an adequate final protocol**. A broader method/basis/grid study and revised
center protocol are needed before parameter inference. FD errors are 0.0002143
and 0.0002235 kcal/mol/Å: both pass, but smaller displacement does not improve
the error, so asymptotic differentiation convergence is not claimed.

The projected Hessians each have 21 internal modes, with lowest curvature about
0.00018123 hartree/bohr² and no projected eigenvalue below -1e-5. This describes
local curvature for the proposed isolated-molecule protocol, not equilibrium or
force-field validation. Full numerical outcomes are in `numerical-checks.json`.

## Identifiability and next data

The model is explicitly:

```text
QM reference = verified fixed MM contributions
             + explicitly fitted unknown MM contributions
             + model discrepancy
```

`identifiability.json` enumerates affected graph interactions, equilibrium
dependencies and nonbonded pairs. Only the verified quartic bond/angle and
adjacent bond–bond/bond–angle polynomials enter a **partial** subtraction.
Main torsions, missing site couplings and nonbonded/exclusion conventions remain
explicit other contributions. They are not set to zero. Partial energy/force
residuals are reported as unexplained contributions, not fitting errors.

The 19 directed queries are grouped by the source's exact role/reversal symmetries;
they are not assumed to represent 19 independent transferable constants. The
pilot constructs symmetry-summed numerical energy and force Jacobians, their
ranks, scaled singular values and conditioning. Source role symmetries reduce the
19 directed labels to 15 groups without proving transferability. For this control,
the four-column numerical subproblem has rank 4 and column-scaled energy-design
condition 7.5645 (energy/force design 6.2833). The partial unexplained energy RMS
is 0.18835 / 0.17644 kcal/mol and partial unexplained force RMS is 3.9369 / 4.0174
kcal/mol/Å for the angle / bond scans; these are not fitted-model errors. Even full rank of that subproblem
does not identify a complete physical residual with unbound nuisance terms.
An explicit discrepancy counterexample retains nonzero null directions.
**No uniquely isolated/transferable coefficient fit is established, and no fit is
performed.** All 19 numerical queries remain `qm_fit_required`; two now have
integrity-accepted reference data for this exact control.

J57's held-out split is one grid column and does not establish independent
validation. J58 records a separate eight-corner/edge-midpoint coverage split,
with 17 training points per scan. It is still not independent conformers.
No fitting/validation thresholds are selected after observing these points.
Before fitting, supply the complete source-bound residual/evaluator, independently
constrain other affected coefficients, introduce justified versioned joint targets,
add independent torsion/third-arm/mixed deformations and held-out conformers,
quantify QM/model uncertainty, and predeclare fit validation criteria. Methyl
HMDS/DP2 data are needed before any methyl-domain transfer claim.

## Reproduction and resume

Use the repository and `src` on `PYTHONPATH`. `--distribution` installs the pinned
source assets; it is not a profile or charge-method selection. The ordinary force
field API remains `forcefield="pcff"`.

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=../island-validation/phase4g1-env/bin/python
QM_PY=../island-validation/phase4j36-native-env/bin/python
FRC=../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files
WORK=../island-validation/phase4j58-qm-pilot

# Existing verified outputs are skipped; partial/failed jobs are never overwritten.
"$PY" -m scripts.phase_4j58 execute --work "$WORK" --python "$QM_PY"
# Explicit postprocessing repair for the preserved completed Hessians:
"$PY" -m scripts.phase_4j58 recover --work "$WORK"
"$PY" -m scripts.phase_4j58 normalize --distribution "$FRC" --work "$WORK" \
  --output /tmp/new-j58-normalized
"$PY" -m scripts.phase_4j58 diagnostics --work "$WORK" --output /tmp/j58-numerics.json
"$PY" -m scripts.phase_4j58 identifiability --distribution "$FRC" --work "$WORK" \
  --normalized /tmp/new-j58-normalized --output /tmp/j58-identifiability.json
```

For a clean execution, safely unpack the native archive, retain its original
outputs as evidence, and run the hash-bound requests into a new work directory.
The original optimization request and versioned executed adapter snapshots are
included. Run those snapshots with `PSI_SCRATCH` pointing outside Git and
`timeout 900`; run the recorded Cartesian restart separately. `prepare` binds
that successful center to the scan requests. The campaign then executes only its
62 declared jobs. A fresh execution may reproduce the visible Hessian metadata
export failure; `recover` repairs it from completed native data without QM reruns.
An actually unconverged or interrupted Hessian cannot be recovered. A new retry
requires a separate retained attempt and declared budget; no result is overwritten.

```bash
"$PY" scripts/independent_phase_4j58.py --evidence docs/evidence/phase_4j58 \
  --json /tmp/j58-independent.json
"$PY" -m island.pcff_qm_pilot.evidence --evidence docs/evidence/phase_4j58 \
  --distribution "$FRC" --json /tmp/j58-relocation.json
"$PY" scripts/validate_phase_4j58_receipt.py docs/evidence/phase_4j58 \
  --repo . --distribution "$FRC"
"$PY" -m pytest -q tests/test_phase_4j58.py
```

Relocation and fresh processes under two hash seeds reconstruct strict intake
identities and all independent statuses without a QM engine or original scratch.
The receipt pins complete ordinary-suite membership, focused J58 and J54–J57
regressions, Ruff, compileall, pip checks, diff checks and historical preservation.

QM execution, data integrity, fit readiness/completion/validation, parameter
completeness, backend execution, cross-engine equivalence and scientific validation
remain distinct. No final PCFF profile or PCFF MD is authorized. Keep
`production_validated=False`, `simulation_readiness="not_established"`.
No universal PCFF/COMPASS or complete `.frc`/PDMS coverage, production MD,
pressure validation, Ewald/PPPM closure, equilibrium or long-NPT readiness is claimed.

Completed verification: **42 focused J58 tests**, **335 J54–J57 regressions**,
and **4,786 passed / 61 existing skips** across 4,847 distinct ordinary-suite
nodes. Ruff, compileall, five pip checks, independent native/unit reconstruction,
receipt/manifest replay and historical receipt/hash checks passed. All **5,649**
base files remain byte-for-byte unchanged. No fitted parameters or complete
profile are created by those software passes.
