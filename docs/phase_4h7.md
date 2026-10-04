# Phase 4H7 — bounded PCFF dynamics and separate-process continuation

## Repository and implementation

Branch `codex/phase-4h7-pcff-dynamics` starts at merged main
`257999e69fb0cc67158e331b94cdb3d37347e678`. Its tree exactly matches reviewed H6
`f9218fad0e2baaac00c69059be84fb4d6de1b687`:
`b2bb0da13b38172350cece7e0abc2049e67ab4f9`. No prerequisite branch was merged locally.

No production scientific engine or historical schema changed. The new acceptance
CLI and example exercise the existing `PCFFSinglePointEvaluator.open_session()`,
`initialize_velocities`, `run_nve`, `run_langevin`, `run_dynamics_segment`,
`create_dynamics_checkpoint`, save/load, and `resume_dynamics` APIs. Validation-only
helpers from G4 supply identical options, frame comparisons, segment persistence,
and Context instrumentation; they do not supply an OPLS physical model.

The H4–H6 finite model, automatic typing domain, native charges, Class II terms,
units, exact Coulomb constant, and identities are preserved. There is no new
integrator, thermostat, restart format, optimizer, force field, periodicity,
packing, crosslinking, or Amber workflow-bundle integration.

## Durable inputs and reconstruction

Acceptance reads the four retained H6 minima and original H5 numerical inputs.
Before use, it verifies every file in the committed H5/H6 evidence manifests.
Each minimum is reconstructed and integrity-checked; convergence, fresh final
verification, fmax≤0.1, and agreement with the saved optimized system are required.
The original and optimized graphs must validate against the H3 assignment and H4
model. Minimum and evaluator parameter/model fingerprints must agree, as must the
H6 summary and initial/final coordinate identities. The evaluator's parameter
fingerprint is the H4 specification identity, not a differently normalized H3 hash.

Each acceptance input directory contains data-only records:

- `base.json`, `system.json`, `minimum.json` using existing workflow encoders;
- `parameters.json` and `model.json` using existing H3/H4 signed contracts;
- local `pcff.frc`, checked against the established pin;
- `initialization.json`, preserving sampled velocities, seed, RNG/version,
  temperature, mass inventory, units, kinetic energy and 3N DOF convention;
- `inputs.json`: **pcff_acceptance_inputs_v1**, an acceptance-only manifest with
  the exact permitted relative filenames, checksums, source and model identities.

The child validates files and cross-record relationships before propagation. It
reapplies the saved minimum to a validated copy solely to compare it with the
saved optimized system. It checks initialization's system fingerprint and masses.
An additional origin check binds the checkpoint's original coordinates and
velocities to that saved setup; it does not regenerate velocities. This check was
also applied offline to all completed bundles, including the earliest child runs.
Checksums establish internal content integrity, not computational authenticity.

All operational reads use checked paths relative to the relocated directory.
Historical paths are not used to find inputs. No construction, embedding, new
charge assignment, new parameter assignment, minimization, or initialization is
run during resume. Existing record validators recompute their consistency checks;
this is not regeneration or substitution of scientific inputs. FRC files, source
artifacts and executables remain external to the repository. Nothing serializes a
Context, evaluator, or executable Python object.

Reconstruction/evaluation requires core NumPy and optional OpenMM. It does not need
RDKit, SciPy, ParmEd, Foyer, AmberTools, or LAMMPS. Each actual child reported that
none of RDKit/SciPy/ParmEd/Foyer was imported. LAMMPS is required only for the
independent acceptance oracle. Checkpoint v1's strict same-environment policy is
unchanged; same-host equality is not a cross-platform bitwise promise.

## Declaration and execution

The declaration was published before propagation under
`../island-validation/phase4h7-dynamics/declaration.json`. It binds H5/H6 file
hashes, case/site/model/coordinate identities, environment, source/executable,
settings, comparison tolerances and no outer retries.

All eight case/integrator combinations use:

- Reference, finite nonperiodic/unconstrained model;
- explicitly chosen LJ=(0,0,1), Coulomb=(0,0,1), not universal PCFF defaults;
- dt=0.1 fs, 200 steps=0.02 ps, retained interval 20;
- Maxwell–Boltzmann initialization at 300 K, seed 78123, all 3N Cartesian DOF;
- BAOAB temperature 300 K, friction 5/ps, thermostat seed 99181;
- NVE energy-deviation guard 1.0 kJ/mol relative to the original total energy;
- 202 evaluations/11 frames for each uninterrupted run;
- 102 evaluations/6 frames per 100-step segment.

No COM/rotation removal or velocity rescaling is applied. Initialization and
thermostat streams are separate. BAOAB bath energy exchange is not judged by the
NVE conservation guard.

For each combination, the CLI runs:

1. The direct NVE or Langevin API for 200 steps.
2. An uninterrupted checkpoint-capable 200-step segment.
3. A 100-step segment, durable checkpoint, directory relocation, then a genuinely
   separate Python process that reconstructs an evaluator/session and resumes
   another 100 steps.

Startup/final independent checks remain enabled. The two shared boundary frames
must have exact coordinates, velocities, step and time, compatible model
identities, and energy/force/component agreement under the existing verification
bounds. Only then is the later boundary omitted from the combined reader. Signed
segment files are not edited. Retained steps are exactly 0,20,…,200; no unstored
frames are inferred.

The split comparison uses atol=1e-10, rtol=1e-12 in canonical Å, Å/ps, kJ/mol,
kJ/(mol·Å), and ps units. Original trajectory identity and energy reference, accepted
lineage/parent checksums, cumulative calls and random-draw counts are checked.
Complete PCG64 state equality is exact. No reseeding, replay or inferred generator
advancement occurs on resume.

## Independent force-model oracle

LAMMPS checks run at steps 0, 100 and 200 for every combination. They reuse H5's
independently generated raw converter graph/parameter data. The H5/H6 AA
`ABC,ABD,CBD` equilibrium-role overrides are reconstructed from that reference's
own angle inventory; raw converter output and explicit override commands are
retained. The oracle does not consume ISLAND's compiled interaction/pair lists.

Reference revision: `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
Executable SHA256:
`ee945ead72f4997905ea42636bc987f3cc4a2a9cf89ab2af30d8e02c4c1cfebc`.
FRC SHA256:
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

The same finite all-pairs model uses nonperiodic boundaries, a cutoff beyond every
pair, no shift/tail, preserved native charges and separate explicit LJ/Coulomb
weights. The exact H4 Coulomb convention is retained. LAMMPS real energies/forces
are converted by 4.184. Totals, compatible grouped components, and every force
component must pass atol=1e-5 kJ/mol and kJ/(mol·Å), rtol=2e-10.

Split/uninterrupted agreement tests the same ISLAND integrators' continuation
consistency. LAMMPS independently tests the force model at selected coordinates;
it does not independently validate the full trajectory integrator or thermostat.

Execution used CPython 3.11.16, NumPy 2.4.6 and OpenMM 8.6.1 on Reference.
The declaration and checkpoint records retain the strict host/build identity; the
evidence manifest also records installed test/dependency versions.

## Results and accounting

All eight combinations completed 200 steps (0.02 ps). Every direct/whole and
split/whole comparison passed; all recorded split absolute differences were zero
on this host, and all four complete BAOAB PCG64 states matched exactly.

| Case | Integrator | Max energy deviation¹ (kJ/mol) | Max reference E/component error (kJ/mol) | Max reference force error (kJ/(mol·Å)) | Seconds |
|---|---|---:|---:|---:|---:|
| butane | nve | 0.014473423 | 5.51e-14 | 1.06e-12 | 59.72 |
| butane | baoab | 11.170799 | 4.97e-14 | 8.1e-13 | 61.03 |
| ethanol | nve | 0.010851061 | 1.15e-13 | 8.53e-13 | 51.45 |
| ethanol | baoab | 7.5641109 | 3.2e-14 | 1.28e-12 | 56.67 |
| PE_DP3 | nve | 0.017687313 | 3.55e-14 | 1.99e-12 | 72.52 |
| PE_DP3 | baoab | 7.3452287 | 1.87e-13 | 1.99e-12 | 70.39 |
| PEO_DP3 | nve | 0.020962325 | 1.56e-13 | 2.46e-12 | 72.50 |
| PEO_DP3 | baoab | 8.2243309 | 1.35e-13 | 2.64e-12 | 73.25 |

¹ NVE deviation is relative to the original initial total energy and remains below
1 kJ/mol. BAOAB values are bath-driven energy excursions, reported without that guard.

All 24 independent LAMMPS comparisons passed. Machine-readable identities,
measurements and file hashes are in [the evidence manifest](evidence/phase_4h7.json).

Every successful combination uses the following actual accounting:

| Operation | Evaluator calls | OpenMM Contexts |
|---|---:|---:|
| Direct 200 steps | 202 | 2: session + fresh final |
| Uninterrupted segment | 202 | 3: session + fresh startup + fresh final |
| First 100-step segment | 102 | 3 |
| Child 100-step resume | 102 | 3 |
| Split cumulative | 204 | 6 across two processes |

LAMMPS uses three external processes per combination (24 checks overall), not
OpenMM Contexts. Additional boundary calls explain the split/uninterrupted count
difference. Each continued trajectory stores 11 combined frames after validated
boundary deduplication. NVE has no thermostat RNG; equality of its absent RNG
state is not a stochastic validation claim.

Per-case elapsed times include input validation/binding, three ISLAND executions,
child reconstruction, persistence and reference operations. They are individual
wall-clock measurements, not speedup benchmarks. Commands, worker PIDs, environment,
full results, serialized initialization, minima, checkpoints, histories, logs and
artifact checksums are retained. No initial coordinates, seed, charge, policy,
budget or tolerance was altered in response to an outcome.

## Reproduction and example

```bash
PY=../island-validation/phase4e2-env/bin/python
$PY scripts/validate_pcff_dynamics.py \
  --retained ../island-validation/phase4h5-whole-system-final \
  --minima ../island-validation/phase4h6-minimization \
  --source ../island-validation/phase4h1-sources/lammps/pcff.frc \
  --lammps ../island-validation/phase4h4-upstream/build/lmp \
  --output /new/acceptance-directory
# Repeat the same command with --execute-declared to run the frozen experiment.
$PY examples/pcff_checkpoint.py \
  --inputs /completed/acceptance-directory/butane-baoab \
  --output /new/checkpoint-example
```

The first CLI invocation validates inputs and exclusively publishes a declaration;
exit zero at this stage means declaration success only. `--execute-declared` exits
zero only when all eight required combinations and input nonmutation checks pass.
Missing/changed inputs or executable are errors, never regenerated substitutes.
Failed runs retain their available signed segments and machine-readable failure
outcomes; no automatic retry follows. The example writes a new first segment,
saves, relocates, spawns a child, and reports its checkpoint at absolute step 200.

These external records are acceptance tooling, not an OPLS/PCFF extension of the
Amber high-level workflow bundle. Coordinates alone are not a restart, and reusing
a seed is not restoration of an RNG state. Only successfully saved boundaries
survive interruption; strict external input/environment compatibility still applies.

## Regressions, limitations and next steps

Synthetic software tests exercise actual PCFF/OpenMM sessions: direct and segmented
binding rejection for changed graph/element/mass/IDs/stereo; resume binding reached
through an otherwise valid analytical setup; changed real synthetic charge and
pair-policy models; startup/trial/final failures; invalidation/fresh recovery;
exact evaluator accounting; last-complete state; frame/evaluation budget stops;
checkpoint/default-application rejection for ineligible diagnostics; original
NVE reference and exact PCG64 continuation. Valid pre-trial budget stops remain
checkpoint-eligible under the existing contract.

External-input tests cover missing/altered files and wrong inventory before backend
use; valid reconstruction; rechecksummed but inconsistent initialization masses;
checkpoint origin versus a different velocity draw; nonmutation/global-RNG
preservation; and reconstruction with RDKit/SciPy/ParmEd/Foyer blocked. No generic
engine defect was found and no shared scientific code changed.

Executed verification:

- Ordinary suite: **1,281 passed, 10 skipped** (279.29 s). Skips are existing opt-in
  acceptance/dependency gates; the H7 real-source matrix ran separately.
- Focused PCFF interoperability regressions: **39 passed** (74.92 s).
- Ruff, `python -m pip check`, and whitespace checks passed.
- The PCFF checkpoint example actually resumed in a separate process after
  relocation, reaching step 200 / 0.02 ps / 204 cumulative evaluations.
- All eight declared real-source combinations and 24 LAMMPS checks executed.
- All 16 saved first/final checkpoint origins were additionally checked offline
  against saved optimized coordinates and initialization, without record changes.
- Historical hashes were rechecked after execution: 172 H5 files and
  46 H6 files remained unchanged.

No required H7 acceptance gate remains unmet in this environment. Readiness remains
limited to bounded software interoperability and force-model consistency.

This bounded interoperability milestone does not establish equilibration,
polymer-property accuracy, transferable chemistry, or scientific suitability.
`production_validated=False`, `simulation_readiness="not_established"`.
