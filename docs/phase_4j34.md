# Phase 4J34: declared GAFF2 runtime and periodic NPT acceptance

J34 starts from `3866e5366092816df2238e96ef35a779ed2dcf5b` on
`codex/phase-4j33-c1-npt-continuation`. The corrective branch is
`codex/phase-4j34-gaff2-periodic-npt-acceptance`. Main and historical
J17–J33 evidence are preserved. This phase establishes software execution
and continuation integrity, not polymer equilibrium or material properties.

## Runtime restoration

The Linux AmberTools installation already contained `antechamber`, `parmchk2`,
`tleap` and `sqm`, but the declared workflow PATH did not expose them.
AmberTools is conda-forge 24.8, build
`cuda_None_nompi_py310h834fefc_101`. `tleap` wraps `bin/teLeap`;
`antechamber` and `parmchk2` wrap binaries in `bin/wrapped_progs`.
`netCDF4==1.7.2` and `cftime==1.6.4.post1` restore the missing dependency
reported by the existing mdtraj installation. No force-field source,
parameter, charge-conservation tolerance or historical environment was replaced.

`docs/evidence/phase_4j34/runtime.json` is the versioned runtime manifest.
Its independently pinned identity is
`cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab`.
It records executable/native executable hashes, observed banners, authoritative
package version/build, every non-bytecode AmberTools package file, conda
metadata, resolved native libraries, Python executable, installed distribution
versions/RECORD fingerprints and GAFF2 sources. Validation repeats dependency
and package checks; missing files, shadowed PATH executables, changed versions,
sources, libraries or fingerprints reject the runtime. Re-signing a replacement
manifest cannot satisfy the independently supplied identity. Python RECORD
fingerprints describe distribution metadata; they are not a hash of every file
in every Python package. The native AmberTools payload and resolved libraries
are hashed directly.

Activation is explicit and scoped to a child process. The caller and historical
validators retain their existing environment. Missing or invalid declared
runtime produces `asset_unavailable`, a stage-specific diagnostic and nonzero
CLI exit; this is never an accepted GAFF2 outcome. Wrappers retain their native
library/Python activation behavior. The archived conda explicit export and wheel
hashes describe reconstruction; an installation at another absolute prefix
requires its own reviewed manifest because executable paths and package prefix
rewrites are intentionally bound.

## Chemically complete fixture and native preparation

The fresh linear fixture is two explicitly hydrogenated PEO DP4 chains
(`[*:1]CCO[*:2]`, 60 atoms), constructed independently with seed 2032 and
unique molecule/atom IDs, then packed through production J21 packing at a
setup density of 0.3 g/cm³. The same initial conformer seed for both chains
is declared; their trajectory is integrated independently, never copied.
The connected 63-atom network uses the reviewed J32/J33 complete C/O chemistry
and crosslink provenance. Historical OPLS fixtures and acceptance bytes remain
unchanged; the J34 matrix also executes OPLS-AA on the J34 fixtures.

The existing single-molecule native Amber adapter still rejects disconnected
inputs. J34 adds an explicit multimolecule API: unwrap each complete component
with strict graph-based periodic images, run the unchanged native AM1-BCC/GAFF2
engine separately, then combine the actual resulting prmtops using ParmEd's
disjoint union. Import validates the exact packed chemical graph, stable atom
mapping, component source records and charge totals. Image shifts and all child
preparation identities are retained. The existing preparation validator only
dispatches the new additive schema; old serialized contracts are unchanged.
No missing term is manufactured or estimated. The existing independent topology,
charge, prmtop/OpenMM and force checks run before dynamics.

Earlier eligibility trials with distinct PEO/PE conformers failed native charge
conservation (approximately 0.001 e residual versus the existing 0.0001 e bound).
They remain software diagnostic controls; charges were not corrected and the
tolerance was not widened. Their artifacts are archived separately from accepted
workflows. An initial trial also exposed wrapped component geometry; strict
unwrapping addresses geometry without changing the final packed coordinates.

## Frozen acceptance and continuation contract

Exactly five cases are required: OPLS-AA linear/network, GAFF2 linear/network,
and PCFF network. Four successes must independently minimize, execute 8 NVT
steps and 32 NPT steps, retain multiple frames/four density blocks and exhibit
actual box-volume change. PCFF must satisfy the unchanged J31-C1/C33-C1 precise
source-verified typing rejection with no NPT execution. No asset failure counts
as an expected PCFF rejection or GAFF2 success.

Each successful case independently executes uninterrupted and split workflows,
splits after NPT step 16, relocates the split bundle and resumes in a fresh
Python process using the existing J33-C1 continuation API. Its declared replay
reconstructs adaptive barostat state by integrating the split workflow's own
initial checkpoint through the prefix, verifies the boundary checkpoint/state,
then continues. No uninterrupted coordinates or thermodynamic metadata feed
resume. Comparison uses the unchanged exact contract for coordinates, velocities,
boxes, volume/density, energies, measured temperature, stage labels, checkpoint
bytes and complete recorded RNG states. Measured pressure is unavailable and
retained as null, compared for matching availability; target pressure is never
substituted as a measurement. Exact native checkpoint equality and deterministic
replay provide the existing barostat-state continuation evidence, not a new
standalone serialization of hidden OpenMM barostat internals.

The matrix persists each result, diagnostic and available logs before proceeding.
A missing/unmet case makes its CLI nonzero. A signed runtime/workflow binding
connects exact runtime, graph, source, assignment, backend and NPT configuration.
The receipt validator requires all outcomes, nested bundle integrity, independent
assignment evidence, exact recomputed comparisons, source/runtime bindings and
conservative readiness. `--offline` validates archived evidence without claiming
the installed runtime is available. Negative-control fixtures are software tests,
not scientific evidence.

## Reproduction

From the repository on the declared Linux server:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
../island-validation/phase4g1-env/bin/python -m scripts.gaff2_runtime \
  --manifest docs/evidence/phase_4j34/runtime.json \
  --identity cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab --activate
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34 \
  --runtime docs/evidence/phase_4j34/runtime.json \
  --runtime-identity cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab \
  --output /tmp/island-j34-fresh-runs
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_receipt \
  docs/evidence/phase_4j34/receipt.json
```

The output path must not exist. The additive receipt archives exact gate commands,
JUnit counts, environment versions, artifact/source/executable hashes, runtime
restoration commands, accepted bundles, failed eligibility diagnostics and
historical preservation results. Its tested implementation SHA is distinct from
the later evidence commit. Full branch-tip identity is verified externally to
avoid a self-referential receipt hash.

These short runs establish neither equilibrium, converged density, experimental
agreement, Tg, mechanical properties, universal force-field coverage, msi2lmp
parity nor production readiness. Always:
`production_validated=False`, `simulation_readiness="not_established"`.

## Frozen native input and observed preparation limitation

The first frozen-code fresh matrix correctly failed GAFF2 linear preparation:
although its `input.mol2` and `sqm.in` matched the earlier successful J34 probe
byte for byte, native SQM optimization followed a different path and the resulting
rounded AM1-BCC charges summed to approximately −0.003 e for one component.
No dynamics began for that failed case. The full failed matrix, native outputs
and diagnostics are retained. This is a native preparation reproducibility
limitation, not a source-availability failure, and is not converted to success.

The accepted workflow therefore declares an explicit frozen preparation input:
`linear-native-preparation.tar.gz`, SHA-256
`6f3191258aa760c7adde1b3902809e2b79860b01c79b181d6f23faa97cd401fa`.
This contains only the genuine successful **current J34** native preparation
of the exact same graph, coordinates and runtime, including AM1-BCC outputs;
it contains no trajectory, energy history or checkpoint. The loader verifies
all root/component artifact hashes, reconstructs and validates the original
native import identity on the newly constructed exact graph, and records explicit
preparation relocation provenance. Reuse on another graph, coordinates or runtime
is rejected. The independently checked source assignment is then used for new,
independent full/split MD executions. This is a bounded reproducible prepared-input
workflow, not a claim of byte-reproducible AM1-BCC regeneration or support for
arbitrary linear systems. Direct native preparation remains available and strict;
there is no automatic fallback or charge correction.

For the frozen accepted workflow, add these arguments to the matrix reproduction
command above:

```bash
--linear-preparation docs/evidence/phase_4j34/linear-native-preparation.tar.gz \
--linear-preparation-sha256 6f3191258aa760c7adde1b3902809e2b79860b01c79b181d6f23faa97cd401fa
```

The receipt distinguishes this frozen-input acceptance from the failed fresh
preparation matrix and records the native regeneration limitation explicitly.
