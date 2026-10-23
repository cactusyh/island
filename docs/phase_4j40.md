# Phase 4J40 — periodic polymer equilibration and cross-engine gates

J40 adds a production-oriented **workflow and diagnostic contract**, exercised
with small, chemically complete retained oligomer fixtures. It does not establish
production equilibrium or physical validation. The overall acceptance command
returns nonzero for insufficient sampling or numerical-reference disagreement.

Base: fetched `origin/main` at
`2610606f8ee15b6ea64df9790ae8ef126c922b7b`, the merge of J39. Branch:
`codex/phase-4j40-polymer-equilibration-cross-engine`. No main merge is performed.
Historical source parameters, charge vectors, numerical acceptance bounds and
J17–J39 receipts are unchanged. Three unrelated untracked files were preserved.

## Frozen matrix and charge provenance

| Cell | Construction | Charge origin and parameter model |
|---|---|---|
| GAFF2 PE | Two explicit-H DP4 PE chains, C8H18 each, 52 atoms | Exact-graph RESP/Psi4, GAFF2 |
| Aromatic | Two explicit-H DP2 biphenyl chains, C12H10 each, 44 atoms | OPLS-compatible parameters with RESP/Psi4 charges |
| Multichain PEO | Two complete explicitly hydrogenated chains, 60 atoms | Native source-neutral OPLS-AA |
| Connected network | Three ancestral PEO chains, 63 atoms, three end closures | Native source-neutral OPLS-AA |

The connected fixture is an end-linked cycle (cycle rank one), **not a branched
gel**. These small systems exercise interfaces; their results do not establish
long-chain polymer behavior or transferable force-field coverage.

PSMILES, DP, end groups, seed, stereochemistry, stable atom IDs, bonds, component
and ancestry membership, final graph, packing, periodic box and setup density
are recorded in each `fixture.json`. Preparation runs on the exact packed graph.
The reviewed J38-C3 raw PE/aromatic quantum artifacts are independently audited
again before preparation. This phase imports validated graph-specific results;
it does not rerun SCF, geometry optimization or RESP fitting. Source arrays are
not copied from a different molecule. Authoritative RESP binary64 values remain
unchanged under the J36-C2 finite-precision policy; the existing GAFF2 PRMTOP
serialization bound is unchanged. Native OPLS uses its original source charges.

The native aromatic charge source remains incompatible (−0.230 e per chain).
An actual PCFF cyclic-typing negative control executes before the positive
matrix. Native AM1-BCC remains blocked by the independently established
AmberTools 24.8 precision limitation (retained neutral fixture −0.003 e); J40
hash-binds that diagnostic and does not claim a new quantum execution or success.
No normalization, redistribution, counterions, invented interactions or charge
origin relabeling is used.

## Actual staged execution contract

Each cell has two replicas with distinct integrator, initial-velocity and
barostat seeds: `(4040,4041,4042)` and `(4050,4051,4052)`. Initial packing is held
fixed, so these are stochastic replicas, not independent packing realizations.
Every uninterrupted and split trajectory starts separately from the same packed
backend, performs its own minimization, and integrates its own states.

The bounded acceptance protocol is:

1. Minimize the actual system using the existing independently force-checked policy.
2. NVT heating: 64 steps at 100 K, 64 at 200 K, 64 at 300 K.
3. NVT pre-equilibration: 256 additional steps at 300 K.
4. NPT: 1024 steps at target 300 K / 1 bar, barostat interval 25.

The timestep is 0.5 fs; NPT spans **0.512 ps**, total integration **0.736 ps**.
Reference OpenMM, LangevinMiddle thermostat (friction 1/ps), isotropic
MonteCarloBarostat and full orthorhombic periodicity are explicit. The effective
nonbonded cutoff is 0.5 nm, Ewald accuracy setting 1e−5, original family mixing
and 1–4 conventions unchanged. Frames/checkpoints are retained every four steps:
369 frames per complete replica, including minimization and 256 NPT frames.
There is no post-NPT production segment.

The additive `equilibration_run` executor retains the historical J33 executor
and serialized records. Its signed `island.heating-schedule.v1` metadata binds
the temperature stages. Historical NVT stage labels remain valid; step boundaries
and actual target temperatures identify heating versus constant-temperature NVT.
Unexpected integration failures preserve partial frames and checkpoints as an
explicit failure record rather than a completed trajectory.

The split stops at NPT step 512, relocates its complete checkpoint bundle, and
resumes in a fresh Python process. Comparisons are exact across all retained
positions, velocities, boxes, energies, measured kinetic temperatures, volumes,
densities, stage labels, checkpoint bytes and opaque RNG bytes, plus all parent,
graph, source, assignment and workflow identities. Pressure is unavailable in
both histories, not replaced with the target. The established J33-C1 independent
prefix replay reconstructs OpenMM adaptive barostat state before checkpoint
restoration; no state or thermodynamic values are borrowed from the uninterrupted
trajectory. Barostat internals are not separately exposed as a public array.

MDAnalysis runs in its declared environment in a separate process for each
replica. Its selections, PBC unwrapping, mass weighting, per-chain aggregation,
Rg/Re/Re-Rg, orientation/nematic, valid RDF, contacts, connectivity and crosslink
density definitions remain the established J38 ones. Unsupported disconnected
system moments and network end-to-end distances retain explicit blocked states.

## Equilibration is a separate, conservative gate

`island.equilibration-policy.v1` is immutable and bound into
`island.equilibration-analysis.v1`. J40 fixes burn-in in advance at 25% of NPT
frames, avoiding post-hoc window selection. Each scalar density, volume and
measured kinetic-temperature series records:

- four contiguous block means and block-mean standard deviation;
- sample standard deviation and uncertainty, taking the larger of block and
  autocorrelation standard-error estimates;
- normalized biased autocorrelation, initial-positive paired-sum integrated
  autocorrelation time, and effective sample size;
- rolling means, first/second-half difference, least-squares full-window drift;
- block independence (minimum block length at least five correlation times).

Declared software criteria require at least 512 retained post-burn-in samples,
ESS at least 50, relative drift at most 2%, temperature mean within 5% of its
target, and pairwise replica means within twice the combined standard error.
These thresholds are conservative screening criteria, not universal physical
constants or proofs of equilibration. The chosen software run has only 192
post-burn-in NPT samples and must fail the sampling floor. A constant series
cannot prove independence. Replica agreement is not accepted before individual
series pass the declared independence/sampling tests.

Density always equals total mass / actual frame-box volume. Initial target
packing density is never reused as a measured density. Instantaneous kinetic
T uses the established unconstrained 3N degrees of freedom. No instantaneous
pressure estimator is implemented in this backend; **pressure stability is
blocked**, not inferred from 1 bar. Without that gate, the combined diagnostic
cannot pass even if the remaining finite-run thresholds do.

## Independent-engine reference

The available retained J22 executable is a real KSPACE-enabled LAMMPS
10 Sep 2025 development build, source commit
`36f7c81ced8dd456430b13647e48f0c7d7bd8faa`, executable SHA-256
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
Actual help/version/package output, source hashes for Ewald/PPPM/pair kernels
and version header are recorded. The executable is reused, not rebuilt in J40.

Each actual final graph and parameter set is exported at its original packed
coordinates, preserving atom order, charges, molecule IDs, cutoff, mixing,
exclusions, 1–4 factors and box. LAMMPS executes **Ewald `run 0`**. No LAMMPS
trajectory or PPPM execution is claimed. OpenMM Reference single points include
finite-difference force sanity checks. Bonded, LJ and total electrostatic
(real + reciprocal) components, energies and all Cartesian forces are compared
under the retained J22 absolute 1e−5 kJ/mol and 1e−5 kJ/(mol Å) bounds (force
relative bound 2e−10). Target numerical Ewald accuracy and absolute comparison
bounds are different quantities; equal Ewald input tolerances do not guarantee
identical finite reciprocal sums across engines. Mismatches are retained and
keep reference acceptance nonzero. They are not hidden behind a same-engine
Foyer or Amber reader comparison, and their cause is not claimed resolved.

## Reproduction and persistence

Run from the repository root with a new, nonexistent output directory:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
mkdir /tmp/j40-reviewed
# Trusted archive hash and all internal source anchors are recorded in J38-C3.
tar -xzf docs/evidence/phase_4j38_c3/fresh-source.tar.gz -C /tmp/j40-reviewed
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j40 \
  --output /tmp/j40-fresh-run --raw /tmp/j40-reviewed/raw \
  --lammps ../island-validation/phase4j22/lmp-kspace
# Expected exit 1: finite-run equilibration/reference acceptance is unmet.
../island-validation/phase4e2-env/bin/python -m scripts.validate_phase_4j40_receipt \
  /tmp/j40-fresh-run --bundle --diagnostic-only
# Diagnostic integrity exit 0; omit --diagnostic-only for acceptance exit 1.
```

The driver selects the pinned AmberTools environment explicitly for GAFF2;
other cells use the declared Foyer/OpenMM runtime. All five execution environments,
Python/package versions and fingerprints are archived. Bundle loading validates
nested backend/source records, exact continuation, heating, replica seeds,
observables/trajectory bindings, recomputed diagnostics and acceptance flags
offline after relocation. Receipt validation additionally checks tested Git source
bytes, archived file SHA-256 values and historical byte preservation.

All execution and verification results, precise commands, implementation SHA,
JUnit counts, numerical differences, artifact hashes and remaining gates are in
`docs/evidence/phase_4j40/receipt.json` and its manifest. The evidence commit is
separate from the tested implementation; resolve its SHA using
`git log -1 -- docs/evidence/phase_4j40/receipt.json` to avoid a self-referential
commit hash.

`production_validated=False` and
`simulation_readiness="not_established"` remain mandatory. No experimental
agreement, Tg, mechanical properties, universal OPLS/GAFF2/PCFF coverage,
AM1-BCC support, universal msi2lmp parity, thermodynamic equilibrium or scientific
production readiness is established.

## Observed J40 outcome

All four cells completed both independent replicas and all eight exact
split/relocated/fresh-process comparisons. The retained initial and final NPT
volumes differ in every replica. All eight relocated MDAnalysis recomputations
reproduce the full original observable payload exactly. Scientific acceptance
remains **unmet**:

| Cell | Energy error (kJ/mol) | Maximum force error (kJ/(mol Å)) | Reference gate |
|---|---:|---:|---|
| GAFF2 PE | 0.04294867896004462 | 0.0419290025488408 | failed |
| Aromatic OPLS-compatible RESP | 0.007970133818901104 | 0.025599056887585903 | failed |
| Native OPLS PEO | 0.062370851543676054 | 0.05479818868769826 | failed |
| Native OPLS connected cycle | 0.128373825792778 | 0.03619004420235683 | failed |

Atom ordering, box export and charge values match exactly in all four reference
runs. The energy-component discrepancies are electrostatic; the network also
fails the independent packed-geometry finite-difference bound. These failures
are not waived. A supplemental, separately labeled **single-point-only** study
at the prior J22 Ewald accuracy setting 1e−10 reduces several errors but still
fails at least one unchanged bound for every cell. It neither changes the MD
profile nor supplies a passing cross-engine gate. Its code, logs, reference
artifacts and results are archived under `accuracy-study/`.

Post-burn-in density ESS ranges from about 2.94 to 10.03 and all density blocks
fail the declared independence criterion. Replica kinetic-temperature means
range from about 144.15 to 168.02 K against the 300 K target. Density uncertainty
estimates are finite-run diagnostics, not validated confidence intervals for
material density. More independent sampling, actual pressure observations and
resolution of the numerical-reference discrepancies are required before any
convergence claim.

The first prototype used one attempted barostat move per integration step; some
tiny boxes shrank below twice the cutoff. That exploratory failure is retained
and motivated the explicit 25-step coupling interval in the frozen J40 protocol.
A first focused validation also detected last-bit NumPy/BLAS dot-product and
polynomial-fit differences between execution environments (28 passed, 1 failed).
The final implementation uses deterministic Python binary64 operations with
`math.fsum` for reductions and an explicit centered least-squares formula;
a cross-environment regression protects this correction. Earlier full/regression
runs were terminated with exit −15 and are **not passing validation evidence**.
MD frames, charges, checkpoint comparisons and acceptance bounds were unchanged
by that diagnostic portability fix. Fresh dynamics and all final validation
were rerun against the final tested implementation.

The lossless evidence archive uses XZ compression and tar hard links for
byte-identical files. `artifacts.tar.xz` is the stored artifact;
`artifacts.tar.gz` is a compatibility symlink for the tested auto-detecting tar
loader. The manifest declares the actual compression. Storage deduplication
happens after independent simulations and exact comparison; it does not create
trajectory frames or replace any execution. Every extracted file is verified
against a SHA-256 inventory. This packaging keeps the artifact below GitHub's
single-blob size limit without altering scientific or workflow bytes.

The initial packing densities (roughly 0.3–0.4 g/cm³) are low-density software
acceptance starts. They are neither equilibrated condensed-phase densities nor
experimental density targets. The primary acceptance execution performs 16
minimizations and 27,648 integration-step invocations, counting independent
split-prefix barostat replay explicitly. Regression/full suites also execute
their existing short integration fixtures; those are additional test executions,
not production simulations. Development probes are separate from these counts.

## Final verification

Implementation actually tested:
`8fd8296081e42db1a891f4fea700bb1c64ef4005`.
The initial implementation was
`5bbb10ec979e0e2ec63051b47eaddf03237b1a85`; the final implementation corrects
cross-environment diagnostic reductions. The later evidence commit contains
only documentation and additive evidence.

- Focused J40: **30 passed** (49.81 s).
- J22–J39 regression selection: **929 passed** (949.41 s).
- Full supported suite: **2927 passed, 10 skipped** (1447.40 s).
- Ruff, compileall, five environment `pip check` commands, diff check: passed.
- Existing historical receipt-validator semantics and all **2222** historical
  tracked document/evidence file hashes: preserved.
- Eight fresh-process continuations and eight relocated full-observable
  recomputations: exact; eleven mutation controls: rejected with exit 2.
- Fresh extracted archive: **625** file hashes and nested offline reconstruction
  passed. Diagnostic validation exits 0; scientific acceptance exits 1.

The ten full-suite skips are existing optional reference/QM opt-ins, not disabled
J40 gates. Exact skip reasons and every command/environment/JUnit result are
archived; no interrupted run contributes to these passing counts.

For a fresh local test reproduction, set `PYTHONPATH=.:src` and point
`ISLAND_J40_RUN` to the newly generated run, then use the supported environment:

```bash
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_phase_4j40.py
../island-validation/phase4e2-env/bin/python -m pytest -q
../island-validation/phase4e2-env/bin/python -m ruff check .
../island-validation/phase4e2-env/bin/python -m compileall -q src scripts tests
git diff --check
```

The exact separately selected J22–J39 test paths, cache/basetemp controls,
source-asset environment variables, receipt commands, subprocess commands and
all five interpreter-specific pip checks are preserved in
`validation/*-command.json`; the captured helper sources use `.py.txt` suffixes.
