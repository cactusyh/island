# Phase 4J40-C2 — executed long equilibration

This phase executes the C1 protocol on the four retained molecular models. It
establishes execution and persistence evidence, **not scientific acceptance**.
One aromatic replica completes the 1 ns NPT budget; eleven replicas encounter
the unchanged periodic box/cutoff guard. Pressure remains unavailable and KSPACE
comparisons exceed the existing numerical bounds.
No temperature, pressure, force or sampling tolerance is widened.

`production_validated=False`; `simulation_readiness="not_established"`.

## Models and unchanged protocol

Base: `c7b0060594db92b7fa1c741690d9bb8c1c825323` on
`codex/phase-4j40-c1-corrective-equilibration`, fetched from origin. Main is not
modified. Historical evidence is checked against this tree, including the J40
archive symlink's link text and resolved bytes.

| Cell | Actual charge origin | Molecular scope |
| --- | --- | --- |
| `gaff2-pe` | Exact-graph RESP/Psi4 | Two hydrogen-complete C8H18 chains, PE DP4 |
| `opls-resp-aromatic` | RESP/Psi4 with OPLS-compatible parameters | Two hydrogen-complete biphenyl chains, aromatic DP2 |
| `opls-peo` | Native neutral OPLS-AA | Two complete PEO chains |
| `opls-network` | Native neutral OPLS-AA | Three ancestral PEO chains joined into a cycle; not a branched macroscopic gel |

Each of three independently seeded replicas independently prepares its exact
packed graph and parameter assignment, minimizes, heats at 100/200/300 K over
5,000/5,000/10,000 steps, and runs another 20,000 NVT steps at 300 K. Timestep
is 0.5 fs; heating plus NVT is 20 ps. The NPT request is the full C1 maximum
budget of 2,000,000 steps (1 ns), at 300 K and 1 bar, with barostat interval 25
and retention/checkpoint interval 100. This executes the maximum budget directly;
it does not claim ten completed evaluation windows. An integration exception
terminates that trajectory, while other replicas continue independently.

The seed triplets are `(4040,4041,4042)`, `(4050,4051,4052)` and
`(4060,4061,4062)`. All dynamics use OpenMM Reference, LangevinMiddle friction
1/ps, the original Ewald error request 1e-5 and cutoff 0.5 nm. The executor and
force builder are unchanged historical code. Newly prepared GAFF2 numerical
parameters and charges are checked against the independent historical source
control; regenerated source receipts retain their actual distinct identities.
Those identities are never overwritten to satisfy the historical replica
comparator. Their descriptive numerical replica comparisons are recorded, while
identity comparability remains explicitly blocked.

## Executed outcomes

See the machine-readable receipt and matrix for per-replica retained frame
counts, exact failure diagnostics, last retained integration step, NPT duration,
actual observables and statistics. The final retained frame precedes an exception;
it is not necessarily the last attempted integration step. The independent
continuation attempt records its exact attempted NPT step when available.

The OpenMM failure is: `The periodic box size has decreased to less than twice
the nonbonded cutoff.` The box is not expanded and the cutoff is not shortened
to conceal this boundary. Partial trajectories are persisted as failures, never
published as completed workflow bundles.

The split run starts independently from the same packed coordinates, box,
configuration and seeds, checkpoints at NPT step 100, relocates the bundle and
continues in a fresh Python interpreter. An additional failure-preserving child
runner reconstructs adaptive barostat state through independent prefix replay,
restores the boundary checkpoint and retains the genuinely integrated suffix.
It borrows no coordinates or thermodynamic values from the uninterrupted run.
Equality includes every retained frame field and opaque checkpoint/RNG bytes.
Equal partial histories do not establish completion of the requested workflow.
The historical continuation child is also attempted. Eleven integration failures
are retained separately. The redundant legacy child for the 1 ns aromatic run
is intentionally interrupted while its independently prepared parallel C2 child
continues; the interruption, command and exit status are retained explicitly.
The parallel split restart is byte/identity equal to the independently generated
primary split. Relocating that child's own completed artifact never borrows
coordinates from the uninterrupted history.

## Observations and scientific gates

Actual temperature is computed from retained kinetic energy, never from the
thermostat target. Density uses actual mass and actual frame box volume.
The sidecar observations contain the volume, density, temperature, potential and
kinetic energies, stage, time and box traces. Component energies are recomputed
in independent Reference contexts on each actual retained geometry and box;
the sum is compared with recorded potential energy at the original 1e-5 kJ/mol
bound. Component aggregation uses `math.fsum` so canonical JSON key sorting
cannot change the recorded comparison error. A preliminary ordinary-sum
bookkeeping failure and the original component records are retained in the
validation archive; individual component energies are unchanged. These evaluations do not advance or change trajectory RNG states.

The unchanged J40 policy discards the first quarter of each NPT history and
requires four contiguous blocks, at least 512 retained samples, ESS at least 50,
block length at least five autocorrelation times, relative drift at most 2%,
temperature within 5% of target and replica differences within two combined
standard errors. The records preserve block means, uncertainties, biased ACF,
initial-positive-sequence integrated autocorrelation, ESS, rolling means and
drift. Short or correlated histories fail these criteria explicitly.

Measured instantaneous pressure is `None`, with observation source `unavailable`.
An OpenMM barostat target is not a pressure measurement. No incomplete position
force virial is substituted for the full periodic, reciprocal-space virial.
NPT scientific acceptance therefore remains blocked even if a temperature or
sampling subcriterion improves. No equilibrium claim follows from run duration.

The completed 20,401-frame aromatic history also exposes a historical v1 NumPy
convolution portability limitation. Strict v1 continuation loading succeeds in
the declared execution environment, while a few rolling means differ by one
last bit in the supported test environment (for example density differences
2.78e-17 g/cm³). Both outcomes and an independent value comparison are retained.
The original analysis bytes and v1 loader remain unchanged. The additive C2
completed-history loader validates the signed workflow, all restart asset
hashes, frame/graph/source/assignment/checkpoint/RNG bindings and actual frame
content; C2 diagnostics independently recompute deterministic fsum statistics.
It does not label the cross-environment v1 summary reconstruction as passing,
and it introduces no serialization or scientific tolerance.

## Independent numerical audits

The retained KSPACE-enabled LAMMPS executable actually runs four fresh Ewald
single points. Graph, coordinates, box, ordering, masses, charges, parameter
terms, units, mixing, exclusions, 1–4 scaling, cutoff and electrostatics settings
are audited from the actual exported input files. Input-file hashes match the
historical exact-model control. The original absolute energy/component and
force bounds remain unchanged. Electrostatic discrepancies remain unresolved;
matching requested relative Ewald accuracy is not proof of absolute component
agreement. Unit and force-sign checks, periodic mapping, neighbor/cutoff
boundaries and reciprocal electrostatics are recorded separately. PPPM is not
claimed.

A separately executed mathematical two-charge unit control measures native
Coulomb prefactors of 138.93545764438198 (OpenMM) and 138.935456264
(LAMMPS), in kJ nm/(mol e²), a relative difference of approximately 9.94e-9.
Installed OpenMM headers, the declared LAMMPS source object, executable hashes,
actual inputs and outputs are archived. No prefactor or charge is rescaled.
This unit control is not molecular dynamics or periodic KSPACE acceptance.
The larger periodic discrepancies remain unresolved. The diagnosis record
separates input units, native constants, images, neighbor/cutoff boundaries,
force signs and reciprocal electrostatics rather than assigning an unproved
single cause.

The minimized network audit is re-executed on the **fresh C2 minimization frame**.
Every atom and Cartesian axis is tested with central differences over
1e-4 through 3e-9 nm. A component needs two adjacent passing step sizes at
1e-5 kJ/(mol Å). Full-model bond, angle, torsion and nonbonded contributions are
checked. Ordinary and crosslink-associated bonded/angle/torsion terms are also
partitioned by the actual three crosslink edges and tested independently.
The source assignment contains no improper terms; none are fabricated.
Constraints, periodic-image margins, cutoff margins and angle conditioning are
recorded. The minimized audits pass; the separately retained packed-geometry
failure is not relabeled. This finite-difference result does not repair the
cross-engine Ewald discrepancy.

## Reproduction and validation

All commands, environment overrides, command exit codes, elapsed times, JUnit
files, executable/package inventories and validation helpers are archived in
`docs/evidence/phase_4j40_c2/validation.tar.xz`. Use a fresh empty output directory.
A representative replica invocation is:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  -m scripts.validate_phase_4j40_c2 \
  --case opls-network --replica 0 --output NEW/opls-network/replica-0 \
  --raw /path/to/independently-validated-J38-C3-raw
```

Repeat for all four case names and replicas 0–2. For GAFF2, use the declared
AmberTools environment returned by `declared_environment` from the retained
J34 runtime manifest; do not substitute PATH tools. Run the additive finalizer
only after all replica attempts are available:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  -m scripts.phase_4j40_c2_finalize --output NEW
```

The validation archive contains the exact independently executed KSPACE and FD
commands and helpers needed to add those audit records. Scientific execution
and finalizer commands return nonzero while acceptance is blocked.

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j40_c2_receipt \
  docs/evidence/phase_4j40_c2/receipt.json --diagnostic-only
# Expected integrity exit 0; omit --diagnostic-only for acceptance exit 1.
```

The validator checks nested source/graph/configuration identities, actual frame
and checkpoint content, statistics, partial continuation equality, component
sums, minimized FD intervals, archived file inventories and historical bytes.
The declared historical J40 archive supplies independently audited raw source
assets and comparison controls, never new dynamics. Fresh-process validation
after relocation uses path-independent record identities.

Native aromatic OPLS remains source-incompatible; unsupported PCFF remains an
expected typing rejection; native AM1-BCC remains blocked by its documented
precision limitation. None is reclassified as RESP/Psi4. This phase establishes
neither universal force-field coverage nor experimental agreement, Tg,
mechanical properties, production equilibrium or scientific readiness.

## Final measured run table

The durations below end at the last retained frame; failed integration attempts
may occur up to one retention interval later. Every row includes 20 ps of prior
heating/NVT. Temperature means and density ESS use the declared 25% NPT burn-in.

| Cell | Replica | Retained NPT (ps) | Frames (all stages) | Mean T (K) | Density ESS | NPT request completed |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `gaff2-pe` | 0 | 92.05 | 2242 | 296.763 | 6.610 | no; box/cutoff failure |
| `gaff2-pe` | 1 | 55.80 | 1517 | 300.472 | 4.672 | no; box/cutoff failure |
| `gaff2-pe` | 2 | 254.65 | 5494 | 295.479 | 8.010 | no; box/cutoff failure |
| `opls-network` | 0 | 4.00 | 481 | 295.574 | 3.706 | no; box/cutoff failure |
| `opls-network` | 1 | 2.55 | 452 | 314.807 | 4.625 | no; box/cutoff failure |
| `opls-network` | 2 | 2.35 | 448 | 320.674 | 2.934 | no; box/cutoff failure |
| `opls-peo` | 0 | 42.20 | 1245 | 305.991 | 3.104 | no; box/cutoff failure |
| `opls-peo` | 1 | 19.80 | 797 | 305.343 | 3.239 | no; box/cutoff failure |
| `opls-peo` | 2 | 19.15 | 784 | 312.783 | 3.541 | no; box/cutoff failure |
| `opls-resp-aromatic` | 0 | 103.15 | 2464 | 300.174 | 4.622 | no; box/cutoff failure |
| `opls-resp-aromatic` | 1 | 35.15 | 1104 | 291.310 | 2.994 | no; box/cutoff failure |
| `opls-resp-aromatic` | 2 | 1000.00 | 20401 | 301.357 | 7.069 | yes |

All density/volume histories remain insufficient under the declared criteria.
Full-budget execution is distinct from equilibration and scientific acceptance.

The trajectory archive may be stored as ordered `artifacts.tar.xz.part-NNN`
files with `artifact-parts.json`. Each part and the combined stream has a SHA-256
identity; concatenation reconstructs the original archive byte-for-byte. This
is lossless storage of independently executed histories, not frame synthesis.

[Density, kinetic-temperature and volume traces (PDF)](evidence/phase_4j40_c2/npt-traces.pdf)
are generated from the retained observations with Matplotlib 3.11.2 in the declared
analysis environment. The dashed temperature line is the requested target;
it is not a measured trace. Logarithmic time axes expose the shorter failed
histories alongside the full-budget aromatic run. Input and output hashes and
the plotting source identity are in `plot-provenance.json`.

## Final verification record

The implementation actually tested is
`9b2bf120cc24db2d42f125d564896c5f7a960173`. The documentation/evidence commit
is separate; resolve it with the receipt's declared `git log` command.

| Gate | Actual result |
| --- | --- |
| Focused J40-C2 | 29 passed |
| J22–J40-C2 regressions | 1,025 passed |
| Full supported suite | 2,993 passed, 10 skipped |
| Ruff, compileall, diff check | Passed |
| Five execution-environment pip checks | Passed |
| Historical byte preservation | 2,267 entries unchanged |
| Independent trajectory-archive inventory | 773 files hash-verified |

The ten full-suite skips are the declared opt-in real AmberTools integration
(`ISLAND_RUN_AMBERTOOLS_INTEGRATION`), retained Phase 4E7.1 analysis
(`ISLAND_ANALYSIS_REFERENCE_ROOT`), the real Foyer source permutation
(`ISLAND_OPLSAA_XML`), and seven archived AmberTools reference tests
(`ISLAND_AMBERTOOLS_REFERENCE_MANIFEST`: checkpoint, NVE, evaluation sessions,
Langevin, minimization, singlepoint and workflow). Their exact names and reasons
remain in the archived JUnit file. These skips are not scientific passes.

Final fresh-process receipt validation passes integrity (exit 0 with
`--diagnostic-only`) and correctly rejects scientific acceptance (exit 1 without
that option). Both report `integrity="passed", scientific_acceptance=false`.
The exact commands, elapsed times, tested commit and receipt input hash are in
`docs/evidence/phase_4j40_c2/receipt-validation.json`; their output logs are
committed alongside the receipt. These post-build validation results are bound
by the evidence commit, avoiding a self-referential receipt hash.

The final staged diff check also covers generated publication assets. Matplotlib's
SVG trailing whitespace was stripped without changing coordinates or plotted
values; the generator, plot hashes and receipt were regenerated. The original
SVG/provenance and first receipt-validation results remain in the validation
archive under `preliminary/publication-first`.
