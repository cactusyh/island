# Phase 4I3: durable prepared single-chain workflows

## Base and scope

Branch `codex/phase-4i3-unified-prepared-workflow` starts from
`8d8a11aa35d1f3491c5c9455059100ecdf9fbc67`. Its tree is identical to reviewed
I2 commit `188aeddcfe133a79e30c03d47c1b631fae8b4688` (squash PR #33).
No prerequisite merge or historical workflow migration was performed.

This additive API consumes I2 bundles for GAFF, GAFF2, pinned Foyer OPLS-AA,
and the named PCFF operational model. It performs no building, typing,
parameter assignment, charge calculation, or QM. Native scientific engines,
force-field definitions, source pins, and checkpoint formats are unchanged.

`production_validated=False` and `simulation_readiness="not_established"`.
A verified local minimum and 0.02 ps trajectory are software interoperability
checks, not equilibration, independent thermostat validation, or evidence of
polymer-property accuracy. The known native chemistry/charge limitations,
including OPLS PS rejection and PCFF's limited acyclic C/H/O domain, remain.

## API

```python
from island.forcefields import PreparedForceFieldSources
from island.minimization import MinimizationOptions
from island.workflows import (
    PreparedWorkflowConfig, start_prepared_bundle_workflow,
    prepared_workflow_status, read_prepared_workflow_frames,
    resume_prepared_workflow,
)

config = PreparedWorkflowConfig(
    minimization=MinimizationOptions(max_iterations=5000, max_evaluations=10000),
    total_steps=200, segment_steps=100, recording_interval=20,
    max_evaluations_per_segment=102, max_frames_per_segment=6,
)
# Amber: sources=None. OPLS: PreparedForceFieldSources(opls_xml=...).
# PCFF: PreparedForceFieldSources(pcff_frc=...). Libraries stay external.
status = start_prepared_bundle_workflow("input-bundle", "new-run", config,
                                       sources=None, segments=1)
# Move the run directory, exit, and resolve sources explicitly in the next process.
status = resume_prepared_workflow("moved-run", sources=None, segments=1)
status = prepared_workflow_status("moved-run", sources=None)
frames = read_prepared_workflow_frames("moved-run", sources=None)
```

Configuration contains execution settings only. Family, native charge method,
source, and special-pair choices are fixed by the validated prepared bundle.
No implicit backend or charge fallback exists. Configuration validates integer
budgets/seeds, minimization options and BAOAB options before execution.
`segments` bounds each invocation; `max_segments` bounds all recorded attempts.
Evaluation/frame limits remain per segment. Artifact bytes include the copied
I2 input bundle and every manifest-listed diagnostic. Actual evaluation counts
are retained in signed segments and stage records; checkpoint counters describe
accepted lineage, while stages also retain failed-attempt counts.

Returns are owned manifest dictionaries or tuples of existing owned dynamics
frames. Preparation coordinates/provenance remain in `prepared/`; minimum and
initialized state live in a separate setup record. Compatible coordinate changes
are evaluated through native binding contracts. Results never become a generic
harmonic/LJ12-6 approximation of PCFF.

## Durable schemas and consistency

New schemas:

- `island_prepared_workflow_config_v1`: explicit numerical settings and budgets.
- `island_prepared_workflow_v1`: checksummed manifest, run identity, immutable
  configuration reference, I2 facade identity/input file inventory, stages,
  file checksums/byte counts, accepted segment references and selected checkpoint.
- `island_prepared_workflow_setup_v1`: exact configuration/facade identity,
  minimization record, applied starting system, and VelocityInitialization.

Existing tagged data-only encoders preserve integer IDs and tuples. I2 native
records are copied byte for byte, including original signatures. Original
libraries remain external and must match native source pins. No historical
absolute path is used operationally. No object deserialization or live Context
is persisted.

Start publication, status, frame reading, and resume share one validator.
It checks file coverage, contained paths, symlinks, checksums and byte budgets;
reconstructs I2 native records; binds immutable configuration, original system,
minimum, applied coordinates, initialization/masses/seeds, native parameter and
model fingerprints, and all accepted trajectory records. Each fingerprint is
compared according to its existing backend definition. Segment budgets and stage
call counts are checked against the signed diagnostic, not trusted as summaries.

The shared trajectory validator is extracted from the existing Amber workflow;
its original adapter and APIs remain compatible. It binds exact initial
coordinates and velocities, physical settings, chemical/mass inventory, original
trajectory identity, environment and lineage. Native segment/checkpoint integrity
checks bind RNG state and cumulative counters. The selected checkpoint must equal
the checkpoint derived from the last accepted segment. Completed runs undergo all
checks before returning from no-op resume.

Adjacent boundaries use the established exact synchronized-state/identity checks
and `verify_final` numerical energy/force tolerances. The combined reader retains
the earlier boundary and never changes signed segment evaluations. It returns
only actual retained frames, including irregular final endpoints.

These are internal consistency checks, not computational authenticity. They cannot
prove that arbitrary re-authored scientific evidence represents a real experiment.

## Failure and recovery

New runs refuse existing directories. `.writer.lock` provides exclusive writer
ownership using the established PID/host record. There is no automatic stale-lock
stealing: verify that the owning process has stopped before removing a stale lock.

Immutable UUID-named records and checkpoints are published before candidate
semantic validation and atomic manifest replacement. Failed publication leaves
the previous manifest/checkpoint authoritative. Unlisted orphan files are never
selected for recovery; inspect them or remove them offline with no active writer.
An interrupted initial copy may leave a directory without a usable manifest;
retain it for diagnosis and use a new directory. Storage exhaustion may prevent
saving a new diagnostic, but cannot advance the committed boundary.

Nonconverged minima retain their full diagnostic and prevent initialization.
Failed dynamics retain diagnostic segment records without adding them to accepted
history or replacing the prior checkpoint. Frame inspection excludes trial frames.
Valid setup failures remain inspectable. Resume refuses automatic retry of a
`stage_failed` run. Budget stops retain eligible verified boundaries, and never
claim completion. No automatic budget increase, reseeding, or tolerance change.

## Dependencies

Importing the workflow uses core/NumPy only. Inspection reconstructs native
records offline: Amber requires ParmEd for its copied original prmtop; OPLS/PCFF
require explicitly supplied pinned local source files. Foyer typing, AmberTools
executables and SciPy are unnecessary for inspection/resume reconstruction.
Status/frame/no-op completed resume construct no OpenMM Context and evaluate no
forces. Actual execution needs OpenMM; initial minimization additionally needs
SciPy. Existing optional chemistry imports may be available transitively but are
not required for these already prepared systems. No runtime downloads occur.

## Reproduction

```sh
PY=../island-validation/phase4e2-env/bin/python
$PY examples/prepared_chain_workflow.py \
  --bundle ../island-validation/phase4i2/gaff2/relocated \
  --output /tmp/island-prepared-example

$PY scripts/validate_prepared_workflow.py \
  --bundle ../island-validation/phase4i2/gaff/relocated \
  --output ../island-validation/phase4i3/gaff
# Repeat for gaff2. For oplsaa add:
# --xml ../island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml
# For pcff add:
# --frc ../island-validation/phase4h1-sources/lammps/pcff.frc
```

Outputs are exclusive; use new directories for a new experiment. The acceptance
CLI first checks retained I2 hashes against `docs/evidence/phase_4i2.json`, then
publishes a declaration before minimization or propagation. It retains failures
and exits nonzero on any unmet gate. The child process reconstructs from the
relocated run, and explicitly forbids minimization/velocity initialization.

Declared settings: maximum atomic force 0.1 kJ/(mol·angstrom), 5,000 iterations,
10,000 evaluations; unchanged line search/fresh verification; Reference;
BAOAB 300 K, 5/ps, 0.1 fs; velocity seed 78123, thermostat seed 99181;
200 steps, interval 20. Split budgets 102 evaluations/6 frames per 100 steps;
whole reference 202 evaluations/11 frames. The reference uses the exact saved
minimum and initialization without rerunning either. Comparisons use atol=1e-10,
rtol=1e-12 in canonical units, with complete RNG equality. No NVE criterion is
applied to BAOAB.

## Executed verification

Final measured results and file hashes are recorded in
`docs/evidence/phase_4i3.json`. GAFF/GAFF2 retained charge vectors remain labelled
synthetic provided +/-0.01 e; their live OpenMM runs do not establish chemical
accuracy. OPLS/PCFF use their existing native source charges. All acceptance
starts from retained I2 records; no preparation or QM was repeated.

The measured minima and completed verification gates follow.

| Selector | Initial → final potential (kJ/mol) | Final max / RMS force (kJ/(mol·angstrom)) | Iterations / evaluations |
|---|---:|---:|---:|
| GAFF | 66.74033303 → 11.87126967 | 0.08215756 / 0.04609018 | 125 / 130 |
| GAFF2 | 57.11223627 → 3.97594836 | 0.08935219 / 0.05891577 | 133 / 139 |
| OPLS-AA | 65.93767568 → 15.90352677 | 0.08044046 / 0.05070616 | 126 / 132 |
| PCFF | 6.84736971 → -38.68791556 | 0.07347116 / 0.04509578 | 143 / 152 |

All four minima terminated `force_converged` with fresh verification under the
predeclared budgets and unchanged 0.1 maximum-force criterion.

The acceptance wall time includes loading, copying, native semantic validation,
minimization, three dynamics executions (whole and two halves), publication,
relocation, child startup and completed-state inspection. It is not an isolated
optimizer or force-evaluation benchmark. In particular the PCFF record contract
repeats expensive deep native validation at public boundaries. This phase adds no
cache that bypasses those checks and makes no speedup claim.

Native trajectory frames do not carry a per-frame RNG snapshot. The RNG comparison
uses the complete saved PCG64 state in the final checkpoint, alongside validated
random-step/draw counters and exact origin/lineage continuity. Numerical state and
evaluation comparisons cover every common retained frame. This preserves the
existing checkpoint/trajectory schema and its limitations.

### Final gates

All four real retained-input runs passed. Each completed 200 steps (0.02 ps),
with retained steps 0,20,...,200 after boundary validation/deduplication. Every
measured maximum split/uninterrupted discrepancy was **0.0** for coordinates,
velocities, potential/kinetic/total energies, forces, components and time. Complete
final PCG64 states and trajectory origins matched exactly. This observed same-host
agreement is not a cross-platform bitwise guarantee.

Each minimization used one reusable Context plus one fresh final Context.
Each checkpoint-capable execution used three Contexts (session, fresh startup,
fresh final). Start therefore counted five; the child counted three; the whole
reference counted three. Whole dynamics used 202 evaluations; split dynamics
used 102+102=204. Completed status, frame reading and no-op resume all passed.
Original I2 input hashes were rechecked unchanged after all four runs.

Total acceptance wall times: GAFF 8.82 s, GAFF2 7.75 s, OPLS-AA 35.27 s,
PCFF 1232.19 s. Runs shared the host; these are workload timings, not backend
speed comparisons. Evidence and generated records remain external under
`../island-validation/phase4i3/{gaff,gaff2,oplsaa,pcff}/`; only the evidence
manifest is committed.

Executed environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1 Reference,
ParmEd 4.3.1, SciPy 1.17.1, RDKit 2026.3.6. Foyer was not invoked.

Final ordinary suite: **1,377 passed, 10 skipped** (470.39 s). Ruff passed.
`pip check` passed in both `phase4e2-env` and `phase4g1-env`. The 39 new focused
software tests cover all four selectors, native dispatch without preparation,
setup/config/origin contradictions, completed/paused checks, tolerated fresh
boundary offsets, failed startup/final/backend execution, budgets, publication,
writer ownership and lazy inspection. Historical Amber workflow tests also pass.
These synthetic fixtures are separate from the four retained-input numerical runs.

Executed examples: `prepared_chain_workflow.py` and
`prepared_forcefield_bundle.py` with retained GAFF2 input, plus
`checkpoint_continuation.py` (analytical NVE/BAOAB separate processes).
Existing opt-in AmberTools/QM and archived force-model matrices were not rerun;
no new independent force-field or thermostat validation is claimed. There are
no unmet gates in the declared four-selector continuation experiment. Remaining
limitations are native chemical/source scope, external pinned library resolution,
expensive PCFF inspection, checkpoint-level RNG storage, and unestablished
scientific readiness. No historical workflow or analysis schema was migrated.

Final manifest review also reproduced two rechecksummed status contradictions:
a paused segment labelled completed, and a paused manifest labelled ready.
Both tests failed before adding status derivation from signed diagnostics and
stage outcomes, and passed afterward. All four retained completed runs were
revalidated with this stricter contract without changing their artifacts.

## Corrective commit: bind rejected dynamics attempts

The reviewed I3 commit `41131cae3c5034e4113362bd153053cc9387d15b` validated
stage summaries and budgets for rejected segments, but checked setup/trajectory
relationships only for accepted segments. This left a semantic publication gap.

Reproduction used the existing explicitly synthetic OPLS fixture: pause at step 2,
independently run a step-0 segment with twice the saved initialized velocities,
inject a trial evaluation failure, and append that independently valid diagnostic
as a failed dynamics stage while retaining the step-2 checkpoint. The regression
failed before the correction (`DID NOT RAISE` at `_publish_manifest`). This is a
software-contract defect, not evidence about OPLS scientific parameters.

The shared manifest validator now processes **all recorded segment attempts in
history order**, with an authoritative boundary that advances only on accepted
segments. Every attempt is bound to the setup or that last accepted checkpoint:

- Actual system compatibility, masses, native model identities, physical settings,
  start step/time, synchronized coordinates and full-step velocities.
- First-attempt origin and saved initialization; subsequent original trajectory,
  exact parent checkpoint checksum, prior accepted lineage and environment.
- Remaining requested steps and budgets derived from the authoritative boundary,
  rather than from the diagnostic's self-reported first frame.
- Cumulative calls, random steps and normal draws relative to that boundary.
  RNG algorithm/version must match the parent. When no draws were consumed, the
  full RNG state must remain identical (or equal the initial thermostat state for
  a first attempt). For consumed draws, the existing native integrity contract
  validates the complete output PCG64 state and counts; no RNG/trajectory replay
  or additional numerical acceptance policy is introduced.

Rejected attempts never become parents, including when a diagnostic retains some
locally completed trial steps. A rejected diagnostic without a published setup is
also invalid. Existing accepted-trajectory validation remains in place.

Authoritative input validation is distinct from failed backend observations. A
first attempt may have no initial evaluation after startup failure; failed
observations need not numerically agree with the minimum or pass final
verification. For resume, native startup failure restores the accepted input
frame, while successful startup retains existing fresh-evaluation tolerances.
Final/trial failure evidence remains untouched. This does not make failed
segments eligible for checkpoints or default coordinate application.

Publication, status, frame reads and resume all use this same contract. Candidate
rejection preserves previous manifest/checkpoint bytes; explicitly rechecksummed
contradictory histories also fail public inspection before the generic
`stage_failed` resume guard. No native signatures, historical records, schemas,
force-field mathematics, scientific engines or numerical tolerances changed.

Added regressions cover the exact reproduction, a different failed trajectory
at the same resumed step, wrong parent/physical/RNG relationships, first-attempt
velocity mismatch, repeated rejected attempts retaining the same parent, public
reader rejection, legitimate first/resumed backend and startup/final verification
failures, and budget-exhausted continuation. The historical evidence file
`docs/evidence/phase_4i3.json` remains unchanged.

Verification results for this correction follow.

Executed focused command:

```sh
../island-validation/phase4e2-env/bin/python -m pytest \
  tests/test_prepared_workflow.py -k 'not pcff' -q
```

Result: **49 passed, 1 deselected**. The deselected existing PCFF software case
is included in the complete ordinary suite. Ruff and `pip check` passed in both
`phase4e2-env` and `phase4g1-env`.

Read-only retained checks used `prepared_workflow_status()` and
`read_prepared_workflow_frames()` on
`../island-validation/phase4i3/{gaff,gaff2,oplsaa,pcff}/relocated`, with the same
explicit OPLS XML/PCFF FRC paths documented above. OpenMM Context construction
was patched to fail during inspection. All four report completed step 200 and
11 retained frames. Before and after inspection, every case file checksum was
compared with the existing `docs/evidence/phase_4i3.json`: 18 GAFF, 18 GAFF2,
17 OPLS, and 18 PCFF artifacts, **71 unchanged** in total. The check's transient
summary is `/tmp/i3-correction-retained.json`; historical evidence was not rewritten.
No preparation, minimization, dynamics, force evaluation, or scientific acceptance
rerun was performed. There are no unavailable retained-workflow gates.

The earlier scientific limitations and conservative readiness flags are unchanged.
These checks establish internal input/history consistency, not computational
authenticity of arbitrary caller-supplied results.

Complete ordinary command: `../island-validation/phase4e2-env/bin/python -m pytest -q`.
Result: **1,388 passed, 10 existing opt-in skips**, 483.90 s. The skipped tests
remain the unrelated opt-in scientific matrices; execution was unchanged, so no
full scientific acceptance rerun was needed. All requested corrective gates passed.
