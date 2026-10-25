# Phase 4J40-C11 — executed GAFF2-PE convergence characterization

This continuation completes the missing numerical work in the first C11 commit.
The first C11 receipt and files are retained unchanged inside the new archive's
`previous-c11/` directory. J17–J40-C10 artifacts remain unchanged in place.
The tested implementation is `2f9c0327919d7fdf46fec6543cfac9a96b16bc89`, based on
C10 `d8a4e878d660f5ac1b1965ac077c24ad052b23e6`. Main is not modified.

## Baseline and corrections

The actual retained comparison is **OpenMM Ewald / LAMMPS Ewald**, not PME.
Separate new PME/PPPM evaluations are included below. The exact fixture is two
hydrogen-complete neutral C8H18 chains (PE DP4), 52 atoms, with independently
verified exact-graph RESP/Psi4 charges. Native AM1-BCC and native OPLS charges
are not substituted. Graph, atom IDs, numerical parameters, charge vector,
shortest-path exclusions and 1–4 factors (LJ 0.5; Coulomb 5/6) are unchanged.
The serialized source retains its inherited textual source label; its actual
provenance/configuration and charge identity remain bound to the RESP/Psi4 path.

| Quantity | Fresh native baseline |
| --- | ---: |
| Box (each orthorhombic length, nm) | 1.0813888144339681 |
| Cutoff / switching / LAMMPS skin (nm) | 0.5 / disabled / 0.2 |
| Requested reciprocal tolerance | 1e-5 |
| OpenMM total energy (kJ/mol) | 145.37771485907777 |
| OpenMM Coulomb energy (kJ/mol) | -0.91646854828999835 |
| Absolute / relative total-energy error | 0.042948678960073039 / 0.000295428216 |
| RMS force-component error (kJ/mol/nm) | 0.060332917957643922 |
| RMS atom-force-norm error (kJ/mol/nm) | 0.10449967927152397 |
| Maximum force-component error (kJ/mol/nm) | 0.41929002548727112 |
| Maximum atom-force-norm error (kJ/mol/nm) | 0.43071781778219348 |
| Maximum OpenMM atom force (kJ/mol/nm) | 2786.0723157153939 |
| Maximum LAMMPS atom force (kJ/mol/nm) | 2786.0700555446692 |
| Relative RMS force-component error | 9.85730148e-05 |
| Charge range (e) | [-0.21333213061141376, 0.10246676635213643] |

Atom-error norm percentiles (kJ/mol/nm) are {'0': 0.0004662661662496911, '50': 0.04217541739849613, '90': 0.13571717235710853, '95': 0.18546694806278857, '99': 0.3832740501527581, '100': 0.4307178177821935}.
All per-atom force vectors, differences, units and normalizations are persisted.
The baseline reproduces the retained energy to roundoff.

The initial C11 report incorrectly divided the historical per-angstrom force
error by ten. Correct conversion **multiplies by ten**. Raw LAMMPS `real`
forces convert with `4.184 × 10` from kcal/(mol Å) to kJ/(mol nm). The initial
raw maximum force was also mislabeled without this conversion. Strict native
header/atom-ID parsing and regression controls now test these conversions;
no test assertion is removed to suppress a failure. The original files remain
available for review.

## Executed bounded study

The fixed plan was frozen before the final execution. It includes 38 paired
single points (76 engine evaluations), plus 24 volume-strain evaluations:
100 new numerical evaluations in total. Each LAMMPS call executes `run 0` and
retains its input, native log, exit code and force dump. Each OpenMM call retains
its native serialized System, force vector, real/reciprocal group energies and
actual platform parameters. The requested tolerances are diagnostics, not
replacement acceptance bounds or production defaults.

| Requested tolerance | Absolute Ewald energy difference (kJ/mol) | RMS force error (kJ/mol/nm) | Max component force error (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| 0.0001 | 0.341516841 | 0.505578316 | 4.05683971 |
| 1e-05 | 0.042948679 | 0.060332918 | 0.419290025 |
| 1e-06 | 0.00391156617 | 0.00581405655 | 0.0374678578 |
| 1e-07 | 0.000102145297 | 0.000592508806 | 0.00398647793 |
| 1e-08 | 0.000102381382 | 6.3214524e-05 | 0.000501508376 |
| 1e-09 | 1.90109087e-05 | 7.04680277e-06 | 4.36303648e-05 |
| 1e-10 | 1.66771354e-05 | 4.96239259e-06 | 2.67480355e-05 |

Additional independently evaluated controls include cutoff 0.4/0.45/0.5 nm,
OpenMM PME versus LAMMPS PPPM at requested tolerances 1e-4 through 1e-7,
fixed-alpha meshes 16/24/32/48/64, LAMMPS interpolation orders 3/5/7,
Ewald alpha 6/7.3/9 nm^-1, explicit Ewald reciprocal limits 4/6/8/10,
isotropically strained box factors 1.25/1.5, and wrapped rigid translation.
OpenMM PME interpolation order is fixed at five in the pinned implementation;
the LAMMPS order study does not pretend to change OpenMM's interpolation order.
The pinned LAMMPS help records double FFT arithmetic with KISS/mpiFFT. No
single/mixed GPU platform or alternate FFT-precision executable is declared;
those precision sweeps remain unavailable. No undeclared executable is used.

OpenMM's actual Reference platform is double precision. The raw LAMMPS help,
Python/OpenMM versions, resolved native-module hashes, executable hash,
package inventories for all five environments and CPU inventory are retained.

## Energy decomposition and pressure

Real/reciprocal decomposition depends on the Ewald splitting parameter. At
baseline the analytical self contributions are about -92.2082 (OpenMM) and
-56.9896 (LAMMPS) kJ/mol, with their different native alpha choices. These
individual terms must not be mistaken for a physical total-energy error.
The net-charge background terms are below 2e-32 kJ/mol for this neutral vector;
the default three-dimensional conducting boundary has no explicit surface term.
The self formula is `-k_e alpha sum(q_i^2)/sqrt(pi)`; the uniform background is
`-k_e pi Q^2/(2 alpha^2 V)`. LAMMPS prints alpha with limited log precision,
so the analytical terms from that value are descriptive, not authoritative
native component subtraction. Actual native real/reciprocal totals are retained
separately. The pinned OpenMM alpha formula is independently checked against
`NonbondedForceImpl.cpp` at source revision b399af4 in the archive.

Six volume-strain points per engine/method estimate `-dU/dV` at relative volume
steps 1e-4, 1e-5 and 1e-6. Positions and box scale together; no dynamics is
executed. Bond, angle, proper/improper, LJ and Coulomb derivative contributions
are retained. For direct Ewald, final-step derivatives are about 3303.36236
(OpenMM) and 3303.36258 (LAMMPS) bar, compared with the LAMMPS native
configurational virial pressure about 3303.36218 bar. PME/PPPM derivatives differ
by about 0.0114 bar at the studied settings. These are high-energy packed-frame
single-point pressures with zero kinetic contribution, not observed 1-bar NPT
pressures. LAMMPS's full native configurational pressure tensor is retained.
OpenMM does not supply an independently validated full force/reciprocal virial
tensor here. No backend pressure acceptance bound is invented, and pressure
remains unavailable for scientific NPT acceptance. There are no constraints in
the serialized model; constrained pressure is outside this control.

## Interpretation and scientific boundary

The default discrepancy is strongly tolerance-sensitive. Tightening Ewald's
request from 1e-5 to 1e-9 reduces the energy difference from 0.04295 to about
1.90e-5 kJ/mol and force-component error from 0.4193 to about 4.36e-5
kJ/(mol nm). The tight residual persists around 1.7e-5 kJ/mol: the unchanged
1e-5 kJ/mol total/component energy gate still fails. The independently retained
native Coulomb prefactors differ by about 9.94e-9 relative; this small difference
is not rescaled away and is not asserted to explain the remaining residual.

Four actual retained geometries from two independently seeded replica histories
are evaluated at both 1e-5 and 1e-9. Default total errors range about
0.0212–0.1068 kJ/mol; tight errors range about 1.66e-5–2.02e-5 kJ/mol.
The shared minimized starting frame is not counted as two independent physical
samples. One favorable control passing does not qualify a production setting.

The 1.5 box-strain case has a much larger energy difference (~1.7543 kJ/mol),
including LJ (~0.40194 kJ/mol) as well as Coulomb (~1.35234 kJ/mol).
That geometry is a deliberately strained single point. It demonstrates that
this export/model combination has geometry-dependent interaction differences;
it cannot be used to infer a uniform Ewald error bound. Atom count remains 52
throughout: box scaling is not polymer system-size scaling, and growth with atom
count has not been established. No force-field, exclusion or 1–4 semantics are
changed in response. A generally valid `bounded_numerical_difference` is **not
established**.

For the intended 0.5 fs timestep, the baseline instantaneous force-difference
acceleration proxy gives a maximum one-step displacement difference ~1.89e-8 nm
and velocity difference ~7.55e-5 nm/ps. This uses the actual atom masses and
`delta r = 0.5 delta F/m dt^2`; it is neither a trajectory-error bound nor a
justification for widening minimization/MD acceptance tolerances.

The scientific classification remains `ewald_unvalidated` and
`pressure_unavailable`. No cutoff recovery, long NVT/NPT, replica equilibration,
or bit-equivalence claim is made. No equilibrium, experimental agreement, Tg,
mechanical property, universal force-field coverage or production-readiness
claim follows from these single points.

`production_validated=False`; `simulation_readiness="not_established"`.

## Reproduction

Use a fresh empty output path and the declared KSPACE binary. Inputs are included
in `artifacts.tar.xz` together with the original C11 control files. After unpacking:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python \
  -m scripts.phase_4j40_c11_study \
  --backend EXTRACTED/study/inputs/backend.json \
  --reference RETAINED-C2/reference/gaff2-pe \
  --frames-root RETAINED-J40/gaff2-pe \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output NEW-EMPTY-STUDY
```

The exact executed commands are retained in `study/invocation.json` and each
LAMMPS `command.json`. The reference and frames input files needed for a rerun
are all copied under `study/inputs/`; arranging them into the command's declared
reference and replica directory layout reconstructs those inputs without any
historical trajectory synthesis. Offline statistics and identity validation need
only the bundle, not those original paths:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j40_c11_receipt \
  docs/evidence/phase_4j40_c11/receipt.json --diagnostic-only
```

Integrity returns 0. Scientific acceptance (without `--diagnostic-only`) returns
1. A fresh-process independent audit reloads the native OpenMM XML and verifies
energy/force byte-value equality, independently recomputes force metrics from
raw real-unit dumps, and reconstructs pressure derivatives from native energies.
It does not import the production metric or characterization helpers.

[Executed tolerance-convergence plots](evidence/phase_4j40_c11/convergence.pdf)
show both engine-method pairs and the unchanged absolute acceptance bounds.

## Independent isolation of the large box-strain failure

The separate raw-parameter calculation in `exception-cutoff-audit.json` finds
**18 special pairs beyond the 0.5 nm cutoff** in the 1.5× strained box, and none
in the original or 1.25× box. Summing their source-valued 1–4 LJ terms predicts
-0.401939047921807 kJ/mol, matching the observed OpenMM-minus-LAMMPS LJ shift
-0.401939047921800 kJ/mol. The missing special-pair Coulomb correction predicts
-1.352337038039448 kJ/mol versus the actual -1.352337277938815 kJ/mol.
These residuals satisfy the unchanged 1e-5 kJ/mol diagnostic component bound.

This isolates a convention/export limitation: the OpenMM exception/bond forces
retain those pair terms beyond the nonbonded cutoff, while this LAMMPS pair path
omits them outside its neighbor cutoff. This explains the strained-box excess,
**not** the original baseline error, where all special pairs lie inside cutoff.
It is not a reason to alter 1–4 scaling, the cutoff, or force-field parameters.
No such alteration is made. The independent calculation uses native source
values and explicit minimum-image distances, not fitted correction constants.

For a fresh numerical rerun using only the bundle and declared executables,
extract `artifacts.tar.xz` to `EXTRACTED` and `validation.tar.xz` to `CHECKS`,
then execute the archived layout helper:

```bash
python CHECKS/validation/reconstruct_inputs.py EXTRACTED NEW-INPUTS
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python \
  -m scripts.phase_4j40_c11_study \
  --backend NEW-INPUTS/backend.json --reference NEW-INPUTS/reference \
  --frames-root NEW-INPUTS/frames \
  --lammps ../island-validation/phase4j22/lmp-kspace --output NEW-STUDY
```

A second, fresh-process **full numerical execution** was run after relocating and
reconstructing inputs exclusively from the archive. All 267 checks of native
energy/force/component values, derived comparisons and pressure-strain values
match exactly. `relocated-execution.json` records the commands and checks; the
second native execution is retained in the validation archive. Command paths,
runtime timings and raw log hashes naturally differ, so raw archive-byte
identity is not claimed. This is single-point reproduction, not MD checkpoint
or RNG continuation or an equilibration result.

## Final verification

The frozen implementation passed 17 focused C11 tests and all 1,079 J22–J40
regression tests. The full supported suite reported 3,047 passed and 10 skipped,
with no failures or errors. The skips require optional AmberTools execution or
archived reference manifests, retained analysis bundles, or a Foyer-environment
OPLS source configuration; they are not counted as executed integration passes.
The exact skip reasons are retained in `validation/full.xml` in the validation
archive. Ruff, compileall, all five execution-environment `pip check` commands,
and `git diff --check` passed. All 2,294 historical evidence entries through C10
remain byte-identical.

`receipt.json` and `manifest.json` bind these command records to the tested
implementation commit and the native artifact hashes. The validation archive
contains each exact argument list, environment configuration, exit code, elapsed
time and JUnit report, together with the runner used to reproduce the gates.
Two fresh receipt-validator processes returned 0 for diagnostic integrity and
1 for scientific acceptance, as expected. Their command records and receipt
SHA-256 are in `receipt-validation.json`; the acceptance failure is intentional
and does not indicate a corrupt receipt.
