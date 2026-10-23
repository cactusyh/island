# Phase 4J36-C1: lossless native RESP precision audit

Base branch: `codex/phase-4j36-executable-psi4-resp-gaff2`.
Exact fetched and remotely verified base:
`a53dcd4bc1222dbb3ded245aa52df2f134726f7c`.
Corrective branch: `codex/phase-4j36-c1-precision-preserving-resp`.

**RESP/Psi4 GAFF2 acceptance remains blocked.** The original two fixture failures
were reproduced before production changes. Charge import loses neither atom
ordering nor native float64 bits. The pinned native RESP constrained solver
returns small nonzero exact binary charge totals. Expanding a float into its full
exact decimal does not change that float, increase solver precision, or repair
its residual. No charge values, charge tolerances, force-field parameters or
historical records are changed.

## Audited path and real native ESP experiment

The original path is retained intact:

1. Complete graph components and stable IDs generate explicitly ordered quantum
   input, formal charge and multiplicity. Graph unwrapping is deterministic.
2. J36 optimized the isolated components at HF/STO-3G, then evaluated HF/6-31G*
   ESP with Psi4 1.9.1. The pinned CDS RESP conda package is 1.0, while its module
   reports 0.8. Both version fields retain their actual meanings.
3. The original native RESP driver writes grid coordinates and loads Psi4 ESP
   through ten-decimal text. Its `results.out` charge table has eight decimals.
4. The authoritative **charge** vector is the native fitter's float64 ndarray,
   already saved losslessly as `.npy`. JSON float/repr rows round-trip to the same
   IEEE-754 bits. ISLAND never uses the rounded table as charge authority.
5. The previous provided-charge gate sums decimal representations exactly for
   these fixtures and rejects the residual. C1 additionally sums the exact
   binary source values using `Fraction.from_float`, independently of Decimal
   context, and verifies total molecular and component charges.

An earlier intermediate precision limit was identified: ESP text quantization.
The pinned Psi4 native API
`ESPPropCalc.compute_esp_over_grid_in_memory(Matrix) -> Vector` exposes the ESP
without this rendering. C1 **actually used this official API** for fresh
HF/6-31G* SCF/ESP calculations and fresh native `resp.espfit.fit` calls on both
retained fixtures. The optimized geometries are the retained native J36 arrays;
geometry optimization was not rerun. Native VDW grid generation, single-stage
hyperbolic RESP options, component constraints and SCF settings remain those of
J36. No Amber fit or external OPLS-AA vector is substituted.

The new independent worker imports no ISLAND typing/charge implementation.
Native ESP and grid arrays, the exact ESP array entering the fit, fitted charges,
geometry, source and executable fingerprints, input/settings, stdout and raw
results are archived. A comparison to Psi4's official GRID_ESP export verifies
API units/semantics: maximum differences are approximately 5e-11 atomic units,
as expected for ten-decimal rendering. This declared check (5.1e-11 au) is an
**ESP export check**, not a charge-conservation tolerance. The rounded export
never supplies the fit or charge vector.

| Fixture/component | Original exact binary total (e, abbreviated) | Fresh native ESP/RESP total (e, abbreviated) |
| --- | --- | --- |
| Linear component 0 | 7.8062556418956319e-17 | -4.4235448637408581e-17 |
| Linear component 1 | -6.0715321659188248e-17 | -8.0664641632921530e-17 |
| Network component 0 | 6.9280518821823733e-17 | 9.9855020085914958e-17 |

All formal component totals are zero. Exact fractions and full terminating
binary-decimal expansions are recorded, together with the original shortest-
decimal totals (1.06e-16, -2.65e-17, 3.307e-17 e). These are different exact
representations of the same individually round-trippable float values; none is
a repaired charge vector. Both component and system gates fail. Eliminating
ESP text truncation therefore does not resolve exact conservation.

Source inspection finds `np.zeros` float64 matrices and `np.linalg.solve` in the
pinned native fitter. It returns float64 charges with 53 significand bits. No
higher-precision **charge** output was identified for this inspected pinned
path. AmberTools 24.8's standalone RESP executable has a six-decimal `(8f10.6)`
charge output; its help and binary fingerprint are recorded. This is not an
improvement over the native ndarray and is not used as a substitute. The
inspected upstream Amber source excerpt is identified as rc5 source used by the
conda recipe before updates, rather than asserted to be an exact 24.8 source
tree. No universal claim about other RESP implementations or versions is made.

## Loader, provenance and failure semantics

`island.resp-native-precision.v1` and
`island.resp-native-memory-precision.v1` bind the trusted original execution,
exact graph, stable atom IDs/elements, initial and optimized geometry, runtime,
worker, source files and artifact hashes. They preserve unchanged native values,
full decimal expansions and hexadecimal binary representations. Loaders
reconstruct records rather than trusting self-signed labels. Mutations of
ordering, elements, coordinates, formal charge, sources, method/basis, ESP,
finite values, charge files or audit identities are rejected.

Charge conservation uses **zero tolerance**, exact rational arithmetic. When
calling historical Decimal-based validators, a local arithmetic context large
enough to represent all finite binary64 values prevents ambient Decimal
rounding; it neither alters values nor relaxes any gate. Historical code and
serialized identities remain unchanged.

`PrecisionRESPCharges.load` is the strict accepted-path bridge. It requires an
unchanged native audit, exact native component/system totals, a live pinned
runtime hash/probe check, and the existing independently reviewed J35-C1
RESP/Psi4 result. Its reviewed rows must be the same native values with full
binary-decimal expansions. No failed vector is published through this bridge.
`provided()` reloads and revalidates before using the existing GAFF2 provided
mode. Actual retained fixtures are rejected at
`native_component_and_system_total_validation`, classified
`authoritative_native_precision_insufficient`.

Consequently GAFF2 parameter assignment, minimization, NVT/NPT, MD relocation,
fresh-process MD continuation, split/uninterrupted histories, and byte-identical
MD checkpoint/RNG comparison are **not_executed** for this correction. Quantum
record relocation/reconstruction is explicitly distinct from MD continuation.
The native AM1-BCC blocked path and external OPLS-AA provided-charge control
remain separate historical origins, not RESP evidence. A software adapter error
in the first new API attempt was corrected before the successful native ESP
runs; its original diagnostics are archived as software failures, not physical
or unsupported-domain evidence.

## Reproduction and validation

Use fresh empty output directories and the retained pinned Psi4 environment
manifest/explicit lock under `docs/evidence/phase_4j36`. From repository root:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j36_c1 --output /NEW/c1-audit
# Expected exit 1: both original and native-memory charge gates are unmet.
```

This command reconstructs genuine archived calculations and rechecks live
runtime hashes; it does not claim to regenerate them. `--offline` explicitly
limits it to artifact reconstruction. To reproduce fresh native ESP and fits,
extract J36's `artifacts/native-quantum.tar.gz` into `/NEW/native`, then for each
`KIND` (`linear`, `network`) execute in a new scratch working directory:

```bash
/DECLARED/PSI4/PREFIX/bin/python /ABS/REPOSITORY/scripts/phase_4j36_c1_native_esp.py \
  --native /NEW/native/quantum-KIND --output /NEW/native-esp-KIND
# Expected exit 1 for nonzero native charge residuals. Preserve raw outputs.
```

Use the scoped environment from `psi4_runtime.environment` and validate the
retained manifest/identity first, as the archived `run_native.py.txt` does.
Exact commands, thread settings, versions, executable and file hashes are
archived. `phase_4j36_c1_source_audit --output /NEW/source-audit.json` reproduces
source/runtime inspection. New fitting results are separate from the frozen
J36 control, and no historical artifact is overwritten.

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j36_c1_receipt \
  docs/evidence/phase_4j36_c1/receipt.json --diagnostic-only
# Exit 0 means intact blocked evidence only.
# Omit --diagnostic-only: acceptance returns 1.
```

The receipt records the implementation commit actually tested and a resolver for
the later evidence commit (`git log -1 --format=%H --
docs/evidence/phase_4j36_c1/receipt.json`), avoiding a self-referential Git hash.
JUnit counts, exact verification commands, environments, known skips and
historical byte-preservation checks are included with the receipt.

## Recorded verification

Tested implementation: `8e12858c0a4e41200f87981415ef85ed477bb4bf`.
Evidence and this verification summary are committed separately.

| Gate | Actual result |
| --- | --- |
| C1 focused | 27 passed |
| J22–J36 analysis, compatibility and corrective regressions | 667 passed |
| Full supported suite | 2,635 passed, 10 skipped |
| Ruff, compileall, all five execution-environment pip checks | Passed |
| Historical byte comparison | 1,108 evidence/manifest files unchanged |
| Historical receipt and retained-bundle checks | Passed, including expected nonzero blocked-acceptance exits |
| C1 diagnostic receipt integrity | Exit 0 with `--diagnostic-only` |
| C1 scientific/software workflow acceptance | Exit 1; native component totals remain nonzero |

The full-suite skip records exactly match the retained J36 skip records. They
cover separately configured older archived integrations and the opt-in real
AmberTools reference test; no new skip or disabling flag was introduced.
See `artifacts/skipped-tests.json` for exact reasons. Existing regression tests
run their established short integration fixtures; they do not establish fresh
RESP-based MD acceptance. All five environment versions/package inventories,
command lines, JUnit files, source/runtime fingerprints, re-signed software
negative controls, and read-only relocation receipts are retained under
[`docs/evidence/phase_4j36_c1`](evidence/phase_4j36_c1/receipt.json).

The recorded runner uses `PYTHONPATH=.:src`, `OPENMM_CPU_THREADS=1`,
`OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and the established
`ISLAND_J33_C1_RUN`, `ISLAND_J34_RUN`, `ISLAND_J35_RUN` fixture paths.
It invokes pytest with `-p no:cacheprovider -q --tb=short`, a fresh `--basetemp`
and `--junitxml` for each suite. The regression selection is every
`tests/test_phase_4j*.py` file whose phase number is 22 through 36 inclusive;
the full suite has no path restriction. See `artifacts/gates.py.txt` and each
`*-command.json` for executable reproduction commands and actual exit codes.

`production_validated=False`; `simulation_readiness="not_established"`.
This establishes no equilibrium, experimental validity, AM1-BCC support or
scientific production readiness. Resolving the source solver's exact residuals
would require a separately reviewed source/representation policy, not charge
normalization, redistribution or hidden tolerance relaxation in ISLAND.
