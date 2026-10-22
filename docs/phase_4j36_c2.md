# Phase 4J36-C2: finite-precision RESP/Psi4 conservation

Base: `codex/phase-4j36-c1-precision-preserving-resp`, remotely verified at
`2a62eeb0b26e8225a0d46db00ba02f44f1f386ba`. Corrective branch:
`codex/phase-4j36-c2-finite-precision-resp`. Main is not modified or merged.

## Numerical finding and declared policy

The C1 exact-equality diagnostics remain intact. Every original J36 and native
in-memory C1 RESP fixture is evaluated, per component and for the whole system.
The authoritative source is the native binary64 NPY array (53 significand bits),
not the eight-decimal RESP table or ten-decimal ESP export. The new records retain
all individual source values, full binary-decimal expansions, hexadecimal values,
array hashes, graph/optimized-geometry mapping and source/executable fingerprints.
There is no normalization, redistribution, or correction of any individual charge.

The nonzero sums are **not solely naive summation errors**: they remain in an
exact rational sum of the source floats. `math.fsum` and a 2,200-digit local
`Decimal.from_float` reconstruction independently agree with the exact binary
sum. There is no authoritative decimal source for these fits, so Decimal is
explicitly a reconstruction of binary values, not a claim about decimal output.
An independent standard-library NPY parser uses header parsing and `struct`
without NumPy or production charge imports, and cross-checks JSON value bits,
ordered stable atom IDs, elements, components and source array fingerprints.
C1's trusted source audits independently validate the geometry mapping.

For binary64 round-to-nearest, let `u=2^-53`, `n` be the atom count, and
`L1=sum(abs(q_i))`. The immutable versioned policy
`island.resp-binary64-conservation-policy.v1` defines:

- Absolute allowance in e: `atol=sum(0.5*ulp(q_i))`, the source representation
  envelope. This uses the exact rational half-ULPs, including subnormal values.
- Dimensionless relative allowance: `rtol=gamma_n=n*u/(1-n*u)`.
- Declared finite-precision threshold: `atol + rtol*L1`.
- Decision: exact rational `abs(sum(q_i)-formal_charge)` does not exceed that
  threshold, for **each component and the entire system**.

`gamma_n` is the standard forward-error bound for an n-value floating-point
reduction (Higham, *Accuracy and Stability of Numerical Algorithms*, second
edition, section 3.1). This is a numerical consistency budget for the native
constrained float64 calculation and representation, not a new force-field
accuracy criterion or a rigorous forward-error bound on the complete RESP solve.
The threshold is derived from precision, atom count and charge magnitudes,
independently of the observed residual, with no residual-fitted constant.
In fact every retained component and system passes even the smaller half-ULP
source envelope alone. These residuals are consistent with source rounding;
they do not establish chemically meaningful charge loss or prove exact arithmetic
inside the solver. A 0.003 e deficit and a 1e-8 e deficit are rejected by controls.

Relative residuals are reported against L1. Relative residual against formal
charge is undefined (`null`) for these neutral components; division by zero is
never used to manufacture a favorable relative result. Full exact values and
rational bounds are in each policy record. Abbreviated absolute residuals:

| Source / fixture | Component sizes | Component residual magnitudes (e) | System residual magnitude (e) |
| --- | --- | --- | --- |
| J36 original linear | 30, 30 | 7.8063e-17, 6.0715e-17 | 1.7347e-17 |
| J36 original network | 63 | 6.9281e-17 | 6.9281e-17 |
| C1 native ESP linear | 30, 30 | 4.4235e-17, 8.0665e-17 | 1.2490e-16 |
| C1 native ESP network | 63 | 9.9855e-17 | 9.9855e-17 |

For the selected native ESP linear vector, component thresholds are
1.3852e-14 and 1.4265e-14 e; system threshold is 5.5633e-14 e.
The native ESP network threshold is 5.0851e-14 e. L1-relative residuals range
approximately 2.1e-18 to 1.9e-17 across the retained cases. Source and output
precision have not increased: accepting an appropriate finite-precision residual
does not make the individual charges more accurate.

A second standalone, standard-library numerical review independently recomputes
both the half-ULP and `gamma_n` bounds from raw NPY bytes, without importing
production policy code. All four sources pass that review. Its source and signed
`independent-policy.json` are archived; the final fresh acceptance command records
the review identity and runs only after that review completed. This is additional
numerical evidence, not a substitute charge provider.

## Persistence, origins and fresh execution

`island.resp-finite-precision-charge.v1` binds the policy, native audit, independent
raw parser result, stable mapping, formal/component totals, unchanged rows and
pinned runtime (Psi4 1.9.1; CDS RESP conda package 1.0 / module 0.8). Loading
requires the caller's independently established record identity and recomputes
all numerical policy fields and nested identities. Changing the vector, source,
executable, tolerance, mapping, totals or origin invalidates the record. A small
schema-specific dispatch in the existing provided-charge wrapper supports this
new record; the historical exact-decimal provided and RESP loaders retain their
original behavior and serialized records.

All four charge records must pass independent parsing before any dynamics starts.
The selected vectors for the fresh GAFF2 cases are the **C1 native in-memory ESP
RESP fits**, not external OPLS-AA charges or the frozen J35 vector. C2 loads those
retained authoritative quantum results; it does not claim a new quantum fit or
geometry optimization. It then performs fresh GAFF2 native typing/parameter
preparation, charge import, independent Amber energy/force verification,
minimization, NVT, NPT, split/uninterrupted execution, checkpoint relocation and
fresh-process resume for the linear and connected network fixtures.

The first linear preparation exposed an independent precision-loss bug: ParmEd
AmberParm component addition reconstructed atom charges from rounded bridge
`parm_data`, discarding edits made before combination. The fix assigns the
unchanged authoritative vector **after** native structure combination/conversion,
before existing `3E24.16` native CHARGE serialization. No tolerance is changed.
The accepted source vector stays bit/value-identical to the native fit; native
PRMTOP scaling/unscaling has the existing <=1e-14 e serialization bound and its
actual maximum error is separately reported. Both complete real workflows are
regression tests for this path; the linear case fails before that fix. Initial
failed preparation artifacts are retained separately from corrected execution.

The inherited J33-C1 continuation comparison uses zero numerical tolerances and
byte-identical checkpoints, including complete opaque checkpoint RNG state.
It compares retained coordinates, velocities, box vectors, volume, density,
energies, measured temperature, stage and all workflow/source/assignment bindings.
Pressure observation remains unavailable (`null`), not copied from target
pressure. OpenMM barostat proposal state uses the established independent prefix
replay before boundary-checkpoint restoration; this is not a claim that checkpoint
bytes alone contain every adaptive barostat implementation detail. Both NPT runs
must change box volume and retain four density blocks. They are short software
acceptance fixtures, not equilibrium or production trajectories.

AM1-BCC remains blocked by its documented native precision failure. External
OPLS-AA electrostatics remains an explicitly separate historical/control origin.
This finite-binary64 policy is specific to authoritative RESP/Psi4 records; it
cannot authorize either alternative origin.

## Reproduction and validation

Tested implementation: `af1a0e77b1bab90fadf1cbf4633bee6046b16ab3`.
An initial receipt check found a software canonicalization defect: NPT outcome
identities include a terminal newline, while charge/receipt identities do not.
The validator now uses each schema's established convention; a dedicated test
rejects cross-domain and mutated identities. Initial passing suites and the
failed receipt check are preserved in a separate archive. All final verification
gates and fresh acceptance runs were repeated against the corrected commit.
The final evidence commit is separate and resolvable using
`git log -1 --format=%H -- docs/evidence/phase_4j36_c2/receipt.json`.

The frozen protocol uses 8 NVT steps and 32 NPT steps, a split at NPT step 16,
0.5 fs timestep, 300 K, 1 bar and the OpenMM Reference platform. There are four
NPT density blocks and no production segment. These very short segments test
software execution and continuation; they do not establish physical equilibration.
The actual native PRMTOP charge roundtrip errors are 1.7347e-18 e (linear) and
5.5511e-17 e (network), separately from the unchanged authoritative source vector.


Use fresh empty output directories. Exact commands, versions, package manifests,
source and executable SHA-256 hashes, JUnit files, numerical policies and workflow
artifacts are retained with `docs/evidence/phase_4j36_c2/receipt.json`.

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  -m scripts.validate_phase_4j36_c2 --output /NEW/j36-c2-run

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.validate_phase_4j36_c2_receipt \
  docs/evidence/phase_4j36_c2/receipt.json
```

`--audit-only` produces policy diagnostics but returns nonzero for incomplete
workflow acceptance. Any preparation, execution, comparison or asset failure
preserves diagnostics and cannot satisfy the matrix. Offline reconstructed
policy records attest artifact integrity, not present executable availability;
the real acceptance driver additionally verifies the pinned quantum and Amber
runtimes before execution.

Historical J17–J36-C1 evidence is preserved byte-for-byte. Since J35's historical
receipt validator also freezes its implementation source files, that original
validator is run from a detached read-only historical checkout when validating
its original implementation. The additive C2 validator binds the current changed
provided-charge loader and assembly fix; no historical receipt is re-signed.

## Final recorded verification

| Gate | Result |
| --- | --- |
| C2 focused suite | 29 passed |
| J22–J36 regression suite | 696 passed |
| Full supported suite | 2,664 passed, 10 skipped |
| Ruff, compileall, all five environment pip checks, diff checks | Passed |
| Historical evidence/manifest byte comparison | 1,232 unchanged files |
| Fresh GAFF2 linear and network workflow cases | Both accepted |
| Relocation, fresh-process continuation and reconstruction | Passed for both |
| Split/uninterrupted frames, checkpoints and checkpoint RNG hashes | Exact comparison passed for both |

The ten skips exactly match J36-C1's retained optional/archive test configuration;
`artifacts/skipped-tests.json` lists every reason. No new skip or disabling flag
was introduced. The focused suite includes two real short GAFF2 workflow cases,
which execute again in the regression and full-suite gates. They are acceptance
integrations, not new quantum fitting or production MD.

`artifacts/gates.py.txt`, each `*-command.json`, and the JUnit files record exact
commands and counts. The established configuration uses `PYTHONPATH=.:src`,
`OPENMM_CPU_THREADS=1`, `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1`, and the
retained `ISLAND_J33_C1_RUN`, `ISLAND_J34_RUN`, `ISLAND_J35_RUN` paths. Pytest uses
`-p no:cacheprovider -q --tb=short`, a fresh basetemp and a separate JUnit output
per gate. The regression selection includes every `test_phase_4j*.py` with phase
number 22 through 36 inclusive; the full suite has no path restriction.
Package inventories and executable hashes for the five environments are in
`artifacts/execution-environments.json`.

`production_validated=False`; `simulation_readiness="not_established"`.
No AM1-BCC support, universal GAFF2 coverage, thermodynamic equilibrium,
experimental agreement, Tg, mechanical-property or production-readiness claim.
