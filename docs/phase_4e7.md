# Phase 4E7: reproducible, resumable short-chain workflow


Current scope update: [Phase 4F1](phase_4f1.md) preserves the default 100-site
limit and adds an explicit up-to-1,000-site provided-charge budget. AM1-BCC
remains limited to 100; the original phase scope below is historical.

## Scope and base

`codex/phase-4e7-single-chain-workflow` is based on updated `origin/main`
`606d439d680527c4413f4fd633eef3b83647e4a9` (PR #14), whose tree matched reviewed
prerequisite `3a1725265e6eff42eeb0ee7ddeb67a0751ebd671`. Phase 4E6 and its
transactional-startup correction were already merged. Main was not changed.

This orchestrates existing engines for one finite capped homopolymer:
local-template construction, GAFF/GAFF2 preparation, verified local minimization,
explicit thermal initialization, and bounded BAOAB segments. It preserves the
existing supported single-bond head-to-tail inter-repeat chemistry, supported
assigned tacticity, explicit hydrogens and **100-site** AmberTools limit.
Unsupported chemistry fails through the builder/preparation rules. It does not
introduce a parameter parser, force constructor, integrator or checkpoint schema.

## API and CLI

The public module is `island.workflows`:

- `WorkflowConfig(...)`: frozen, versioned configuration with owned provided
  charges. `from_dict()` and `to_dict()` support strict JSON configuration files.
- `preflight(config, operation="start")`: checks dependencies for `inspect`,
  `start`, `prepared_start`, or `resume`. Resume never discovers AmberTools and
  does not require SciPy. All workflow imports remain lazy with respect to
  RDKit, ParmEd, OpenMM and SciPy.
- `inspect_chain(config)`: builds and returns an owned local-template chain.
  Use its actual stable IDs to supply a complete charge mapping. Inspection does
  not require AmberTools and does not assign charges.
- `start_workflow(config, segments=1)`: refuses any existing output directory;
  runs setup and at most the requested number of bounded segments.
- `resume_workflow(directory, segments=1)`: reconstructs saved inputs and continues
  only from the manifest's published boundary. It never repeats construction,
  parameterization, minimization or thermal initialization.
- `workflow_status(directory)`: checks manifest structure, checksums, progress,
  segment lineage and artifact coverage, returning an owned status dictionary.
  Phase 4E7.1 also validates the saved
  bundle/config/setup-to-trajectory relationship, including completed runs.
- `read_workflow_frames(directory)`: returns owned, validated retained frames from
  accepted segments; compatible adjacent boundaries appear once, retaining the
  earlier accepted evaluation. Conflicting duplicates or out-of-order frames are errors.
- `start_prepared_workflow(config, system, preparation, artifact_directory,
  evidence=..., segments=1)`: explicit acceptance/reuse entry for validated
  historical preparations. Evidence must be `archived_parameters` or
  `synthetic_software_test`. It still minimizes and initializes explicitly;
  it is never a fallback for failed live preparation.

Configuration exposes PSMILES, DP, optional tacticity/atactic fraction,
independent template/assembly/stereo/velocity/thermostat seeds, force-field and
charge method, exact provided charges, charge tolerance, AmberTools timeout and
installation, all `MinimizationOptions`, and BAOAB settings. Limits cover total
steps, segment steps, per-segment evaluations/frames, recording interval,
maximum segment attempts, and published artifact bytes. The byte limit applies
to immutable published artifacts, not an operating-system quota for external
tool scratch files or crash orphans. External execution retains its timeout.
All configuration remains fixed on continuation; `segments` controls only how
many bounded segments to attempt in this invocation.

Use the [runnable PE configuration](../examples/short_chain_workflow.json):

```sh
# Run from a suitable working directory, or copy/edit output_directory explicitly.
python -m island.workflows inspect examples/short_chain_workflow.json
python -m island.workflows start examples/short_chain_workflow.json --segments 1
python -m island.workflows status island-pe-run
# Exit Python/shell. In another process on the same compatible host:
python -m island.workflows resume island-pe-run --segments 1
python -m island.workflows resume island-pe-run --segments 1
```

Start needs RDKit, ParmEd, SciPy, OpenMM and a consistent AmberTools installation
on PATH/AMBERHOME. Resume needs only the reconstruction/evaluation dependencies
(RDKit, ParmEd, OpenMM), not AmberTools executables. The workflow always requests
OpenMM Reference. No default zero charges, method substitution, charge adjustment,
force-field fallback or tolerance relaxation occurs. The example explicitly
requests AM1-BCC with tolerance `1e-4 e`. Provided mode requires `provided_charges`
with canonical decimal stable-ID keys in JSON; the built inventory, not DP alone,
determines coverage. Use the same construction settings for inspection and start.

Python usage:

```python
from island.workflows import (
    WorkflowConfig, inspect_chain, start_workflow, resume_workflow,
    read_workflow_frames,
)
config = WorkflowConfig(
    psmiles="[*]CC[*]", dp=3, output_directory="my-pe-run",
    force_field="gaff2", charge_method="am1bcc",
    total_steps=60, segment_steps=20, max_segments=3,
)
chain = inspect_chain(config)
print(sorted(chain.topology.sites))
status = start_workflow(config)  # one segment; a failure is reported, not concealed
# Another Python process only needs the run directory:
status = resume_workflow("my-pe-run")
frames = read_workflow_frames("my-pe-run")
```

## Stage and success semantics

The manifest records setup outcomes (`passed`, `supplied`, `failed`,
`unavailable`) and dynamics reasons, local evaluator calls and immutable record
references. A minimization advances only if it is force converged **and** its
independent final evaluation passed. Its full result, optimizer diagnostics,
initial/final energy, forces and history are retained. Initialization records
its own temperature, seed, PCG64 identity, NumPy version, masses, units and DOF.

Workflow status distinguishes `starting`, `ready`, `paused`, `budget_exhausted`,
`stage_failed`, and `completed`. Only reaching total requested steps at a
checkpoint-eligible verified boundary produces `completed`. Budget stops are
saved when the existing checkpoint rules permit them and do not trigger an
automatic retry. Failed/unverified dynamics are retained as diagnostic segment
artifacts referenced from stage outcomes; they do not enter the accepted frame
sequence or replace the published checkpoint. An explicit later resume can retry
from the earlier published boundary, subject to the fixed attempt limit.

Segment counters retain the existing trajectory/continuation meaning. The stage
records additionally preserve calls from failed attempts, which are not part of
a later accepted checkpoint lineage. Calls from unpublished or interrupted work
cannot be inferred from the last manifest; orphan diagnostics may retain them. No step, half-step velocity, rejected
candidate, or consumed-but-unaccepted RNG state is promoted into a restart.
Fresh startup/final verification and private evaluation-session cleanup are
provided by the existing segment and session implementations.

## Durable schemas and reconstruction

- Configuration: `island_single_chain_config_v1`.
- Manifest payload: `island_single_chain_manifest_v1`, enclosed by `payload` and
  SHA-256 of canonical payload JSON.
- Prepared bundle: `island_single_chain_bundle_v1`, covered by the manifest's
  file SHA-256 and byte count.
- Segment/checkpoint records: unchanged `island_dynamics_segment_v1` and
  `island_dynamics_checkpoint_v1`, with their original content checksums.

The bundle retains original authoritative atom sites, bonds, any ordered
impropers, coordinates and metadata; original preparation record and signature;
import source/provenance/signature; logical source-file mapping; full minimum;
minimized starting system; and thermal initialization. All retained regular
parameterization artifacts are copied into the run directory, including original
prmtop, MOL2, restart, lineage, command/log files and SQM output when applicable.
Required source checksums are cross-checked against the historical preparation.
SQM output is also bound to its recorded checksum.

Tagged data containers preserve integer metadata keys, tuples and lists without
ambiguity: `{"map": [[key, value], ...]}`, `{"tuple": [...]}`, `{"list": [...]}`.
Only finite JSON scalars and these containers are accepted; repeated mapping keys
are rejected. Config charge IDs are restored through their specific schema.
There is no pickle, executable object, serialized Context or dynamic class
lookup. Fixed constructors rebuild validated ISLAND models. Checksums establish
content integrity, not computational authenticity.

On resume, original prmtop bytes are parsed using the existing Amber importer.
Reconstructed resolved content must validate against the **original** imported
signature and original signed preparation record. No changed content is re-signed
to bypass compatibility. The saved minimum must validate and reproduce the saved
starting system and provenance. Initialization must bind that system and config.
The checkpoint then applies its existing conservative environment/model/mass/
stereo/settings checks and fresh-coordinate verification before propagation.

Operational files are checked relative paths resolved inside the supplied run
directory, including symlink-escape rejection. Historical absolute paths remain
unchanged in signed records and are not used for artifact reads. Moving the
entire bundle on the same compatible host therefore works. `output_directory`
in saved config remains historical; resume uses its explicit directory argument.
Prepared input, charge, repeat-unit, stereo and historical preparation metadata
are preserved. The starting-system file intentionally remains the minimized
starting state; final dynamic coordinates and velocities live in segments and
checkpoints, not an ambiguously updated coordinates-only system.

## Publication, ownership and recovery

Publication targets POSIX filesystems supporting same-directory atomic rename,
hard links and fsync. An exclusive `.writer.lock` acquired with `O_EXCL` records
host/PID. Another writer fails with `WorkflowBusyError`; locks are never silently
stolen. After an interrupted process, an operator must confirm that the recorded
owner is no longer running before removing the stale lock. Read-only status
inspection does not need ownership.

Each segment record and eligible checkpoint is written under a unique immutable
name and fsynced before manifest replacement. The manifest is written to a
same-directory temporary file, fsynced and atomically replaced; the directory is
then fsynced. A failure before manifest publication leaves the previous boundary
authoritative. Source byte blobs use `.dat`; tagged records/segments/checkpoints
use `.json`. The logical artifact map identifies the original source filenames.

Never choose the highest numbered/newest file for recovery. Only the published
manifest selects progress. Orphan immutable files and `.publish-*` files may be
removed after verifying no writer is active and confirming that no manifest
reference names them; they are not automatically promoted or deleted. Retained
`preparation-work/` diagnostics are historical and are not resume dependencies.
If setup failed before a prepared bundle was published, resolve the failure and
start a **new** directory. If a bundle was published before the first segment,
resume uses the saved minimized system and initialization, without repeating
setup. Only the last successfully published state survives later interruption.

`WorkflowError` reports malformed configuration/bundles or I/O problems;
`WorkflowBusyError` specializes ownership conflicts. Scientific failures are
explicit stage outcomes. Dynamics publication errors are raised without advancing progress. Setup
publication problems are recorded as a failed stage if the manifest can still
be written; otherwise the I/O error is raised. Existing checkpoint eligibility and integrity checks
are retained.

## Recorded frames and restart distinction

Segment records contain only the frames the existing bounded driver retained:
initial boundary, scheduled absolute-step frames and the final accepted boundary,
without duplicates inside a segment. The combined reader deduplicates compatible
shared boundaries using exact synchronized states and existing evaluation
tolerances and never interpolates omitted steps. Diagnostic attempt frames
remain in their stage-referenced records and are not mixed into accepted history.

Frames include stable-ID coordinates (angstrom), synchronized full-step velocities
(angstrom/ps), absolute step/time (ps), U/K/total energy (kJ/mol), instantaneous
kinetic temperature (kelvin, all 3N Cartesian DOF), forces and identities. No COM
translation/rotation is removed. Segment origins/lineage bind trajectory and
segment identities. Checkpoints separately preserve the complete actual PCG64
state required for continuation. Coordinates alone, retained frames alone, or a
reused random seed are not restart state.

## Numerical acceptance

Settings were declared in `scripts/validate_workflow.py` and the example config
before execution: capped PE DP=3 (20 explicit sites), Reference, 300 K, friction
5/ps, timestep 0.1 fs, 60 steps (0.006 ps), three 20-step segments, recording every
10 steps, 22 evaluator calls/3 frames per segment. Seeds: construction 2026,
velocity 78123, thermostat 99181. Minimization: fmax 0.1 kJ/(mol*angstrom),
500 iterations/2000 calls. Uninterrupted-versus-split tolerances:
`atol=1e-10`, `rtol=1e-12` for coordinates, velocities, forces and energies.

Both actual cases completed after moving the bundle and continuing in new
Python processes:

| Evidence | Charge selection | Charge residual (e) | Minimized energy (kJ/mol) | Final minimization fmax | Minimization calls |
| --- | --- | ---: | ---: | ---: | ---: |
| Live PSMILES to AmberTools | GAFF2 / AM1-BCC, tolerance 1e-4 | -1.999967073326725e-6 | 7.506535488625179 | 0.0752321217956569 | 141 |
| Checksummed archived PE | GAFF2 / historical provided charges | 0 | 4.088844375191407 | 0.08006897945685111 | 139 |

Initial minimization energies were 60.31182752872998 and 57.183215727517144 kJ/mol.
Both independently verified minima passed the original criteria. No charge
method, charge value, parameter or acceptance tolerance was changed. The
archived provided-charge case is an execution/consistency test, not an AM1-BCC
physical substitute.

All coordinate, synchronized velocity, force and energy differences against
uninterrupted runs were **zero**; complete PCG64 states and trajectory origins
matched exactly. Each split trajectory used 66 evaluator calls versus 62 for the
uninterrupted comparison and retained seven distinct frames. Final kinetic
temperatures were 206.68574814420077 K (live) and 206.81249442263982 K (archive).
These individual instantaneous temperatures are not equilibration tests.

Actual tool environment: AmberTools 24.8, OpenMM 8.6.1 Reference, NumPy 2.4.6,
Python 3.11.16, ParmEd 4.3.1, RDKit 2026.03.6 and SciPy 1.17.1. Source signatures,
checksums and summary metrics are in [the acceptance summary](references/phase_4e7/acceptance.json).
Generated parameterization artifacts and full bundles remain outside the
repository in the user validation directory reported there.

Reproduce either path with a **new** output path:

```sh
python scripts/validate_workflow.py --output /new/live-pe --report /new/live-report.json
python scripts/validate_workflow.py --output /new/archive-pe \
  --manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" --report /new/archive-report.json
ISLAND_AMBERTOOLS_REFERENCE_MANIFEST=/path/to/references.json \
  python -m pytest -q tests/test_workflow_archive.py
```

Ordinary workflow tests use explicitly labelled synthetic preparations with real
ParmEd parsing and OpenMM evaluations. They include moved-bundle, three-process
continuation with setup functions prohibited in resumed processes, exact
uninterrupted comparison, call accounting, owned inputs, configuration/coverage/
chemistry limits, nonconvergence, malformed/missing artifacts, atomic publication
failure, concurrent writer rejection, budget stops, optional-import isolation,
and failure diagnostics. Injected failures test orchestration, not live science.

## Delivery verification

- Complete ordinary suite: **760 passed, 8 skipped** (99.79 s). Skips were the
  opt-in five-case AmberTools generator and seven archived test modules.
- Separately enabled archived single-point, minimization, NVE, evaluation-session,
  Langevin, checkpoint and new workflow acceptance: **7 passed** (191.04 s).
- New workflow tests: **26 passed**, including real synthetic OpenMM execution;
  injected stage/publication failures are labelled software tests.
- The dedicated **live PE GAFF2/AM1-BCC workflow** ran successfully, separately
  from ordinary tests. The older five-case AmberTools regeneration test was not
  enabled; no historical archived parameterization was regenerated.
- Ruff, `python -m pip check`, and checkpoint-continuation, NVE, Langevin,
  evaluation-session, minimization, local-template-builder and workflow-inspection
  examples passed.

## Remaining limits

This is thermostatted isolated-molecule execution, not bulk NVT, density, pressure,
thermal equilibration or production validation. Existing BAOAB finite-timestep
bias, nonperiodic small-system scope and checkpoint environment restrictions
remain. There is no asynchronous mid-step persistence, network/distributed lock
service, automatic stale-lock recovery, parameter-generation fallback, long-chain
charge transfer, new force-field family, PBC, packing, constraints, COM removal,
NPT or MLIP. `production_validated=False` and
`simulation_readiness="not_established"` remain in every manifest.

See [Phase 4E7.1](phase_4e7_1.md) for transactional publication validation,
setup-to-trajectory binding, boundary deduplication and inspection dependencies.
