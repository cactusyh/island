# Phase 4G4 — bounded OPLS dynamics and durable continuation

## Base and implementation scope

Branch `codex/phase-4g4-oplsaa-dynamics` starts at updated main
`81dc0da25b22506b311fd032691d9e77efcc1e70` (PR #23). Its Git tree
`fde351e3f7e00d3d11d4dde74ae472e51da4bfba` exactly matches reviewed 4G3
`b7be50e8b93a8f1b3afe1bb3c9748c4c85c05367`, including squash-merged 4G2.1.
No prerequisite merge was repeated. Main, Phase 4F5, signed historical records
and unrelated checkpoint directories were not modified.

The existing OPLS evaluator/session already satisfies the generic bound-potential
contract. **No integrator, optimizer, thermostat, checkpoint format or numerical
engine changed in this phase.** Additions are OPLS interoperability regressions,
a reproducible acceptance CLI, a separate-process example and evidence.

The physical source remains Foyer
`dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`, XML SHA256
`c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
Native charges, RB convention, geometric sigma/epsilon mixing, exclusions and
separate 0.5 LJ/Coulomb 1–4 scaling are unchanged. See [4G2](phase_4g2.md) and
[4G3](phase_4g3.md) for source scope and installation.

## API-level execution

```python
initialization = initialize_velocities(system, temperature_kelvin=300, seed=78123)
with evaluator.open_session() as session:
    segment = run_dynamics_segment(system, session, initialization.velocities, options)
checkpoint = create_dynamics_checkpoint(segment)  # enforces existing eligibility
save_dynamics_checkpoint(checkpoint, '/new/path/checkpoint.json')

# In another process: reconstruct compatible system, source, and parameter records.
# Do not repeat typing, minimization or initialization.
checkpoint = load_dynamics_checkpoint('/relocated/path/checkpoint.json')
with reconstructed_evaluator.open_session() as session:
    next_segment = resume_dynamics(checkpoint, reconstructed_system, session, budgets)
```

The unchanged direct `run_nve()` and `run_langevin()` APIs also accept OPLS
sessions. Every process owns its own session; the checkpoint contains data and
actual PCG64 state, never a Context/evaluator. Coordinates alone or a seed are
not a restart. Full-step velocities use angstrom/ps; masses use dalton; force
uses kJ/(mol*angstrom); energy uses kJ/mol. Thermal initialization and temperature
diagnostics retain all 3N unconstrained Cartesian DOFs, including translation
and rotation. A single draw is not rescaled to an exact temperature.

Changing graph, elements, stable IDs, masses or assigned stereochemistry triggers
bound-system or checkpoint incompatibility before propagation. Changed model,
parameter, backend or settings identities are rejected at the fresh startup
check before a new step. Runtime failures retain the last complete state; fresh
verification remains available through the parent evaluator after session
invalidation. Existing budget reservations, original NVE energy reference,
strict environment checks and checkpoint/application eligibility remain intact.

This is **not** support in the Amber-specific high-level workflow bundle. No
new generic workflow framework was added. The scripts use existing data-only
system/minimum encoders and parameter/checkpoint validators at the API level.

## Durable inputs and declaration

The actual inputs were the four successful 4G3 minima under
`island-validation/phase4g3-minimization`, with original systems/parameters from
`island-validation/phase4g2-energy-final`. These are retained execution locations,
not hard-coded operational API paths.

Before propagation the CLI:

1. Compares retained 4G3 outcomes against the reviewed committed evidence.
2. Verifies 4G2 source-file hashes against 4G3's original declaration.
3. Validates native parameter records against their authoritative systems and
   exact pinned XML, then checks parameter identities against the accepted minima.
4. Loads and validates each minimum, including convergence and fresh verification.
5. Applies the minimum to the original system and checks exact agreement with
   the saved optimized system, initial/final coordinate fingerprints, final energy,
   iteration/evaluation counts and parameter identity.
6. Publishes a new immutable declaration with complete consumed file hashes,
   identities, options, seeds, budgets, units, tolerances and gates.

The declaration seals generated minimum/optimized files for this experiment;
4G3's reported identities and internal minimum integrity provide the earlier
links. These are content/consistency checks, not cryptographic authenticity.
Original input hashes are rechecked at completion. Nothing is regenerated if
required retained data is missing.

### Reproducible commands

Use the existing pinned-source environment (Python 3.11, NumPy 2.4.6,
OpenMM 8.6.1 Reference, Foyer 1.2.0 pinned revision, ParmEd 4.3.1).
The full checkpoint environment record is retained in the declaration.

```sh
python scripts/validate_oplsaa_dynamics.py \
  --retained /path/to/phase4g2-energy-final \
  --minima /path/to/phase4g3-minimization \
  --xml /installed/foyer/forcefields/xml/oplsaa.xml \
  --output /new/phase4g4-run

# Inspect declaration.json, then execute exactly that plan:
python scripts/validate_oplsaa_dynamics.py \
  --retained /path/to/phase4g2-energy-final \
  --minima /path/to/phase4g3-minimization \
  --xml /installed/foyer/forcefields/xml/oplsaa.xml \
  --output /new/phase4g4-run --execute-declared

# Runnable example reusing validated external inputs and saved initialization:
python examples/oplsaa_checkpoint.py \
  --inputs /path/to/phase4g4-run/pe3-baoab --output /new/checkpoint-demo
```

The CLI refuses an existing declaration/output collision. Execution requires
exact equality to the saved declaration. Every required gate must pass for exit
0. Failures retain diagnostic result files and machine-readable outcomes with
nonzero status; there is no seed search, retry, parameter fallback or tolerance
increase. A preflight input failure is reported before creating a declaration.

### Acceptance artifacts and relocation

For each case/integrator, a new directory contains data-only external original
and starting systems, the unchanged parameter envelope, verified minimum,
initialization and local copy of the exact source XML. `inputs.json` lists these
files and SHA256 values. This small acceptance inventory is not a new public
restart or workflow schema. Checked relative paths resolve artifacts; source
paths recorded in historical metadata are not operational dependencies.

The first segment and ordinary checkpoint are published, then the directory is
moved from `*-before-relocation` to its final case/integrator name. A child Python
process loads the checksummed external inputs and checkpoint and creates a new
OPLS evaluator/session. The child records its distinct PID, output checkpoint,
evaluator count, Context count and imported optional-module status.
`whole.json`, `first.json` and `continued.json` are data-only encodings of the
existing `DynamicsSegment`; checkpoint files use the existing ordinary
checkpoint persistence API unchanged. No live Python/OpenMM objects are stored.

Resume reconstruction requires core ISLAND/NumPy and OpenMM, plus RDKit only
where assigned stereo validation needs it. For these cases the workers imported
**none of Foyer, ParmEd, SciPy or RDKit**. Those packages remain installed on the
compatible host; this observation does not override the strict checkpoint
Python/NumPy/build/CPU/OpenMM compatibility policy. Foyer/ParmEd are needed by
the independent reference oracle; RDKit/Foyer are used for the separate expected
PS rejection. Neither typing nor parameter assignment is invoked during resume.

Only the last successfully saved boundary survives process interruption. Orphaned
acceptance files are not automatically selected as newer checkpoints. MolecularSystem
alone is not a complete restart. The example performs one further demonstration
run in a new directory and does not alter its supplied artifacts.

## Predeclared experiment and accounting

All four cases ran both NVE and BAOAB, with:

- Reference, finite nonperiodic/unconstrained systems;
- 0.1 fs timestep, 200 steps = **0.02 ps**;
- retained interval 20, steps `0,20,...,200`;
- initialization at 300 K, velocity seed 78123;
- BAOAB friction 5/ps, thermostat seed 99181;
- NVE maximum absolute total-energy deviation 1.0 kJ/mol relative to the
  original trajectory's initial total energy;
- uninterrupted budget 202 evaluations and 11 frames;
- each 100-step segment budget 102 evaluations and 6 frames;
- continuation budget 100 steps / 102 evaluations / 6 frames / interval 20;
- child wall-clock timeout 300 s, no outer retries.

Every case executed the direct engine, an uninterrupted checkpoint-capable run,
and a two-segment separate-process run. The direct and uninterrupted runs each
made **202 evaluator calls**; the split run made **102 + 102 = 204**. The latter
includes the additional startup/final boundary checks. Context construction
counts per case/integrator were: direct **2**, uninterrupted segment **3**, first
segment **3**, resumed segment **3**, reference oracle **3**. Segment APIs use an
independent fresh startup and final Context in addition to their reusable Context.

The combined reader in this harness verifies the shared step-100 boundary then
keeps the earlier accepted frame, yielding 11 retained frames. It never modifies
signed segments or interpolates unstored positions. Step/time and original
trajectory identity/reference were preserved, as were cumulative counters.

Split tolerance was predeclared `atol=1e-10, rtol=1e-12`, separately in angstrom,
angstrom/ps, kJ/mol, kJ/(mol*angstrom) and ps. Complete PCG64 state comparison is
exact. Integration consistency uses the same ISLAND engines; **it is not an
independent validation of the integrator**.

## Actual numerical outcomes

All eight required runs and their separate-process continuations completed.
At common retained steps, the measured maximum split/uninterrupted differences
were **0.0** for coordinates, velocities, potential/kinetic/total energies,
forces, components and time. Direct-engine versus segment comparisons also
passed. This observed same-host equality is not a cross-platform bitwise promise.
All BAOAB PCG64 states and original trajectory identities matched exactly.

| Case | NVE max drift, kJ/mol | BAOAB max total-energy change, kJ/mol | BAOAB normal draws | NVE / BAOAB elapsed seconds |
|---|---:|---:|---:|---:|
| Butane | 0.01417585 | 4.22713656 | 8400 | 5.362 / 4.686 |
| Ethanol | 0.01176069 | 5.45919209 | 5400 | 4.046 / 3.830 |
| PE DP3 | 0.01777532 | 6.65499036 | 12000 | 5.403 / 5.253 |
| PEO DP3 | 0.02179236 | 8.83313598 | 13800 | 5.835 / 6.045 |

Whole and split maxima agreed in every row. BAOAB exchanges energy with the
bath: its energy changes are diagnostics, not NVE failures, and no conservation
guard was applied. NVE consumed no thermostat draws. Times cover all three
execution paths, process startup, persistence/reconstruction and reference
checks, not an isolated integration benchmark.

At retained steps 0, 100 and 200, the independent pinned
`Foyer.apply -> ParmEd.createSystem` geometric-mixing reference from 4G2/4G3
verified potential energies, every force and bond/angle/RB/combined-nonbonded
components. It does not use ISLAND conversion or pair-generation helpers.
Predetermined absolute tolerances stayed at 1e-5 kJ/mol for energy and
1e-5 kJ/(mol*A) for force, with relative 2e-10. Maximum measured errors across
all checks: **2.84e-14 kJ/mol energy**, **1.14e-13 kJ/(mol*A) force**.
All component checks passed. This verifies the force-model implementation
independently of the split-continuation consistency check.

PS DP3 was rejected before dynamics for the unchanged -0.23 e native-charge
residual. No chemistry was added, no source changed and no charge repaired.
There are no unmet required acceptance gates for these eight bounded runs.
The machine-readable declaration, input identities and complete per-case outcomes
are committed as `phase_4g4_acceptance.json`; detailed artifacts remain under
`island-validation/phase4g4-dynamics`.

## Verification and limits

The new OPLS-specific software tests use the explicitly synthetic H2 spring
fixture with real OpenMM sessions. They exercise direct/segmented binding
rejection for graph/element/mass/stereo changes, valid-checkpoint resume binding,
model/parameter identity mismatch, short budget exhaustion, exact split RNG,
backend failures, failed startup/final verification, count accuracy and
checkpoint/application ineligibility. Existing generic and Amber tests continue
to cover the established integrator, persistence, environment and failure rules.
No engine defect or numerical correction was required.

The pinned source has limited tested chemistry, no explicit improper section,
and the known PS native-charge failure. Neutrality and short-run completion do
not establish scientific accuracy, equilibration, polymer observables or
production readiness. No GPU, periodic cell, packing, long-chain equilibration,
new charge method, PCFF or crosslinking was added. Further OPLS chemistry or
long-time dynamics is not a prerequisite for PCFF development.

`production_validated=False`, `simulation_readiness="not_established"`.

### Final verification

- Ordinary suite: **1,108 passed, 10 expected opt-in skips**.
- New OPLS dynamics/session tests: **27 passed**, using labelled synthetic
  numerical fixtures. Existing generic/Amber regressions passed in the full run.
- Ruff and `python -m pip check`: passed; pip check ran in both ordinary and
  pinned-source environments.
- Executed examples: new `oplsaa_checkpoint.py` (200-step BAOAB, relocated child
  process, 204 cumulative calls), existing analytical `checkpoint_continuation.py`,
  `run_nve.py` and `run_langevin.py`.
- Actual retained-input acceptance: all **eight** case/integrator combinations passed,
  with expected PS rejection; all eight resume workers executed separately and
  all initial/middle/final independent-reference comparisons passed.
- No required acceptance inputs were unavailable. Historical opt-in archive/QM
  checks were not rerun and remain distinct from the executed OPLS acceptance.
