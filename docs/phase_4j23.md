# Phase 4J23: periodic OpenMM molecular-dynamics workflows

Base: `7b5500dfff913d37b234ba723bca2e3ec36a82c6`, verified on fetched
`origin/main`. J23 adds `island.periodic_md`. All J17–J22 records, source pins,
receipts and existing nonperiodic interfaces retain their original identities.

## Inputs and public API

`PeriodicMDWorkflowConfig` selects PCFF, OPLS-AA, GAFF or GAFF2, an execution
mode, step budget, temperature, timestep, friction, velocity/integrator seeds,
and optional NPT pressure/barostat settings or annealing schedule. It is frozen;
the schedule is an owned immutable tuple. Invalid/nonfinite settings reject.

- `minimize_periodic_energy(backend, config, path)` publishes a minimized step-zero
  workflow, with an independent fresh-context force and energy check.
- `start_periodic_md_workflow(backend, config, path, steps=...)` accepts an existing
  validated J22 `PeriodicParameterizedSystem`. It optionally minimizes, initializes
  velocities, and executes NVT, NPT or annealing.
- `periodic_md_workflow_status(path)` and `read_periodic_md_frames(path)` are offline
  inspection functions. They import neither OpenMM nor RDKit/Foyer/ParmEd/SciPy.
- `resume_periodic_md_workflow(path, steps=...)` restores the saved force model and
  binary checkpoint. It performs no PSMILES generation, crosslink construction,
  packing, typing, charge assignment or parameter lookup.
- `save_periodic_md_workflow(path, destination)` validates and publishes a relocatable
  copy; `load_periodic_md_workflow(path)` validates its committed history.

Load/resume accept trusted `expected_identity`, `expected_config`,
`expected_backend_identity`, and `expected_checkpoint_identity` assertions.
Checksums detect corruption and stale nested records; a trusted external identity
is necessary to reject wholesale replacement of a valid, consistently re-signed
history. Offline validation checks all accessible JSON relationships and binary
hashes. Actual continuation additionally requires restored OpenMM state and
checkpoint bytes to equal the saved frame exactly.

`compare_periodic_md_workflows(first, second)` validates both histories and reports
exact frame/checkpoint equality. A pass requires both configured trajectories to
be complete; different generation boundaries are permitted.

## Identity and state

The manifest binds the complete J22 backend, original graph/history, initial
coordinates, packed-system metadata, ordered stable atom IDs, source/family,
typed/charge/assignment identities, periodic configuration, initial box,
integrator/barostat settings, OpenMM environment, and exact force XML hash.
The initial `FinalChemicalGraph` identity never changes. NPT box vectors live in
frames and checkpoints, separately from that initial graph and packing identity.

Each retained frame includes positions (nm), velocities (nm/ps), box vectors
(nm), step/time, energies, temperature/pressure, integrator and barostat
configuration, manifest/backend/graph/packing identities, and the complete
base64-encoded OpenMM binary checkpoint plus its SHA-256. Those bytes contain
OpenMM random and internal engine state. No JSON reconstruction of a partial
random state is presented as equivalent to a binary checkpoint.

## Supported execution contract

Execution is deliberately restricted to OpenMM **Reference** with matching
recorded version, OS/architecture and byte order. CPU/CUDA/OpenCL requests reject.
All frames are retained; this first implementation targets short software
verification runs and loads the complete history for inspection.

NVT uses `LangevinMiddleIntegrator` with explicit temperature, timestep,
friction and seed. NPT adds `MonteCarloBarostat` with explicit pressure, interval
and seed. NPT uses isotropic orthorhombic cells; a cell shrinking below twice the
nonbonded cutoff rejects the candidate segment. Triclinic cells are unsupported.
Annealing is NVT with an explicit sequence of `(exclusive_end_step, kelvin)`
boundaries. A step uses its scheduled temperature; the saved boundary records
the next step's temperature. NPT annealing is not included.

Minimization fixes the box and uses OpenMM local minimization with a finite
iteration budget. Acceptance requires a fresh-context maximum force no larger
than the requested tolerance, nonincreasing energy and finite state. The
internal RMS stopping threshold is tightened by sqrt(particle count) so the
public maximum-force criterion is checked independently. Failed candidates
retain diagnostic coordinates and error records; they are not published as
converged states. Minimization is never repeated during resume.

## Publication and relocation

The immutable root manifest is published with generation zero. Subsequent
`generation-NNNNNNNN` records append a hash-linked chain of checkpoints. Files
are written in a staging directory, flushed and fsynced; Linux
`renameat2(RENAME_NOREPLACE)` atomically publishes the directory without replacing
an existing destination. Platforms without this operation reject publication.
A writer lock serializes continuation and relocation. Failed candidates remain
in hidden diagnostic directories; the previous valid manifest and generations
remain readable. Locks left by interrupted processes require explicit inspection
and recovery; automatic lock stealing is not implemented.

Source libraries and native preparation are validated when adopting the J22
backend. Resume reconstructs OpenMM forces from the retained numerical snapshot
and verifies the resulting XML hash, avoiding family-native reconstruction or
original absolute source paths. Runtime source replacement is not an option;
use an explicit new backend/workflow for a new source or model.

## Reproducible example

```python
from island.periodic_md import (
    PeriodicMDWorkflowConfig, start_periodic_md_workflow,
    resume_periodic_md_workflow, save_periodic_md_workflow,
)

# backend is an existing validated J22 PeriodicParameterizedSystem.
config = PeriodicMDWorkflowConfig(
    family=backend.config.payload["family"], mode="npt", steps=6,
    temperature_kelvin=300, pressure_bar=1, timestep_fs=0.2,
    friction_per_ps=1, seed=2026, velocity_seed=2027,
    barostat_seed=2028, barostat_interval=1,
)
started = start_periodic_md_workflow(backend, config, "run", steps=3)
save_periodic_md_workflow("run", "relocated")
resume_periodic_md_workflow(
    "relocated", expected_identity=started["manifest_identity"],
    expected_config=config,
)
```

To reproduce the retained four-family software checks:

```sh
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 python scripts/validate_periodic_md_workflow.py \
  --run-tests --artifacts /new/empty/j23-tests --output /new/empty/j23-receipt.json
```

Fixtures use the retained J21 single-chain DP3 polyethylene packing in a 30 Å box,
seed 2026, with the existing real-source J22 assignments. NVT, NPT and annealing
run six steps per family, split after step three. Relocated child-process results
must equal every uninterrupted frame, including all checkpoint bytes. Separate
minimization cases cover all four families. These are software-contract fixtures,
not validated material-property calculations.

## Scientific and LAMMPS boundary

LAMMPS remains J22 export-only in this workflow. The pinned executable lacks
KSPACE; status remains `blocked: KSPACE package unavailable`. J23 runs no LAMMPS
Ewald/PPPM comparison. J22's resource-limited PME trial is unchanged. The retained
J23 execution fixtures use Ewald; they do not establish universal PME performance.

No scientific equilibration, validated density, Tg, mechanical properties,
universal periodic chemical coverage or production readiness is claimed.
`production_validated=False` and `simulation_readiness="not_established"`.

## Executed verification

- Focused J23: **46 passed in 73.22 s**, including 12 family/mode combinations,
  four minimizations, relocation, child-process continuation, full checkpoint/RNG
  byte equality, offline import guards and mutation/publication rejection controls.
- Ordinary pytest: **2042 passed, 10 skipped in 685.29 s**. The skips remain
  explicitly optional historical/integration checks; they are not numerical passes.
- Ruff (`ruff check src tests scripts`), compileall (`src/island scripts tests`),
  both main and isolated-Foyer-environment `pip check` commands, and
  `git diff --check` passed.
- All 12 tracked J17–J22 receipts match their base bytes; their SHA-256 values
  are recorded in `docs/evidence/phase_4j23.json`. Existing runtime/source files
  were not edited. The exact commands, retained log hashes, source identities,
  per-case initial/final boxes and checkpoint hashes are in that receipt.

The ordinary run covered the new workflow engine and existing regression suite.
The final focused run additionally exercised the public complete-trajectory
comparison helper and final checkpoint publication checks. No LAMMPS simulation
or scientific material-property acceptance was performed.
