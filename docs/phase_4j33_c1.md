# Phase 4J33-C1: NPT continuation acceptance

This additive correction starts from reviewed J33 commit
`6f34bf370746ff01fc0491490c0b3dda4a696ec0`. The original driver ran the
linear and connected-network OPLS-AA fixtures, relocated and reconstructed
finished records, but never split and resumed an NPT integration. The retained
before-gap record labels this a software acceptance gap.

The C1 driver independently minimizes and integrates each uninterrupted and
split workflow with the unchanged J33 fixture, configuration, parameters and
Reference platform. Both run eight NVT steps and 32 NPT steps, retaining every
step; the split stops after NPT step 16. Its checkpoint bundle is copied to a
new directory and resumed in a fresh Python process. No state or measurement
is imported from the uninterrupted run.

OpenMM's opaque checkpoint preserves RNG bytes but does not preserve the
MonteCarloBarostat's adaptive proposal scale and acceptance counters. Merely
loading the boundary checkpoint reproduced its state but diverged on the next
network volume move. C1 therefore also retains the **split run's own initial
NPT checkpoint**. In the child process, independent integration of that prefix
reconstructs the adaptive barostat state. The replay must reproduce the saved
boundary frame and checkpoint byte-for-byte before the boundary checkpoint is
loaded and continuation begins. This adds 16 real integration steps per resumed
case; it is explicitly recorded, not represented as checkpoint-only restoration.
A replay mismatch fails without borrowing reference values or relaxing tolerances.

Companion restart, lineage and comparison records bind the parent workflow,
configuration, backend, final graph, source, charges, typed graph and assignment;
initial and boundary checkpoints; all frame content; XML engine assets; and
path-independent identities. Loading checks these bindings and recomputes density
analysis offline. Historical J23/J32/J33 schemas and Linux publication code remain
unchanged. The new records are bounded to zero optional production steps.

The comparison checks every field of every retained frame: coordinates, boxes,
velocities, energies, kinetic temperature, volume, stage and checkpoint bytes.
Measured pressure remains unavailable, explicitly `null` in both histories.
Complete opaque checkpoint bytes are compared as RNG evidence; no separately
extracted RNG state is claimed. Density and all analysis payloads must also match.
Both absolute and relative comparison tolerances are zero.

The matrix requires both OPLS-AA cases to execute NPT, change volume, retain four
density blocks, relocate, continue in a child process and compare exactly. PCFF
must reach the frozen typing rejection on the same connected final graph with
the pinned verified source and no NPT. Runtime checks retain the exact exception
and traceback origin; receipt checks independently replay typing and verify the
exact diagnostic. GAFF2 remains `asset_unavailable` when the pinned tools cannot
be discovered; its required-case failure keeps the matrix CLI exit nonzero.
In this validation environment the pinned binary files exist, but the declared
runtime PATH does not expose them and AMBERHOME is unset. The asset-discovery
record distinguishes this environment configuration failure from absent files;
GAFF2 NPT was not started. No replacement executable was introduced.
A valid diagnostic receipt is distinct from all-family acceptance. The receipt
validator's `--require-all` option enforces the latter explicitly.

Run from the repository root with `PYTHONPATH=.:src`, `OPENMM_CPU_THREADS=1`,
`OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`:

```bash
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j33_c1 --output /fresh/empty/path
../island-validation/phase4e2-env/bin/python -m pytest tests/test_phase_4j33_c1.py -q
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j33_c1_receipt docs/evidence/phase_4j33_c1/receipt.json
```

Exact validation commands, environments, actual counts and artifact hashes are
recorded in the additive C1 evidence directory. Injected/mutated controls and the
failed checkpoint-only attempt are software diagnostics, never scientific evidence.
The retained short integrations test software continuation; they establish neither
equilibrium nor converged polymer density, experimental agreement, Tg, mechanical
properties, universal force-field coverage, KSPACE validation or production readiness.

`production_validated=False` and `simulation_readiness="not_established"`.

The tested implementation is `57a10f6394711b70aee20aefd82bc731e58fd892`.
Validation on Linux x86_64 with OpenMM 8.6.1 Reference completed as follows:

| Gate | Actual outcome |
| --- | --- |
| C1 focused | 50 passed |
| J33 | 7 passed |
| J32 | 53 passed |
| J31 including C1 | 111 passed |
| J26–J30 | 48 passed |
| J22/J23 | 91 passed |
| Full supported suite | 2,328 passed, 10 skipped |
| Ruff, compileall, diff check | passed |
| Supported/Foyer/analysis pip checks | passed |
| Unused AmberTools pip check | failed: mdtraj 1.10.3 requires absent netcdf4 |
| J22/J29/J30/J31-C1/J32/J33 retained receipt validators | passed |

Both OPLS-AA cases passed zero-tolerance continuation with 41 total frames
(32 NPT frames each), real volume changes and four density blocks. PCFF passed
its expected typing rejection with no NPT. The complete matrix exited 1 because
GAFF2 was unavailable in the declared runtime environment. The corrected receipt
validator accepts this truthful diagnostic evidence; `--require-all` exits 1.
The additive receipt records the evidence snapshot separately from the later
documentation stamp, whose branch-tip SHA is verified remotely and reported in
the delivery response. Historical evidence files were checked byte-for-byte
against the reviewed base, including the original J33 archive; its two workflows
and their original relocated copies were reconstructed read-only in fresh paths
and child processes. The ten full-suite skips remain visible in the archived
JUnit report; no required test was disabled to obtain these results.
