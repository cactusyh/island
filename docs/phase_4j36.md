# Phase 4J36: executable Psi4/RESP, failed exact charge acceptance

Base: `5bd1498745bab65a57ca1d2b752c7c1638aa5a86`, from
`codex/phase-4j35-c1-resp-psi4-provenance`. Work uses the exact base;
main is not modified. J17–J35 evidence remains byte-for-byte unchanged.

A dedicated Linux environment now actually executes Psi4 geometry optimization,
ESP generation and the native CDS RESP fitter for both retained GAFF2 fixtures.
**J36 workflow acceptance fails at exact component-charge conservation.** No
RESP-based GAFF2 parameterization or MD is claimed or executed in this phase.
The generator and evidence loader are usable, but a validated J35-C1 provider
publication and downstream workflow remain gated and unverified. The acceptance
CLI always returns nonzero until these downstream requirements are implemented
and demonstrated as well as exact charge conservation.

## Restored runtime and real calculations

The dedicated prefix is recorded in `evidence/phase_4j36/environment.json`.
It contains Python 3.11.14, Psi4 1.9.1 and conda RESP package 1.0
(`pyh31278eb_0`, Psi4 channel). Its Python module reports version **0.8**;
these are different version fields, both recorded, not relabeled.
LibXC 7.1.2 from an initial solver result could not import the pinned Psi4:
`XC_GGA_XC_TH_FL` was absent. Pinning the documented working LibXC 7.0.0,
Libint 2.9.0 and optimization dependencies restored execution. The separate
RadonPy installation was audited but neither modified nor used for final runs.
The explicit conda lock, package metadata, executable/module/library/basis
fingerprints, command logs and failed restoration attempts are archived.

Fresh independent worker calculations used each graph component, stable atom
IDs, deterministic graph-based PBC unwrapping, isolated-component HF/STO-3G
optimization (`GAU_LOOSE`, at most 100 iterations), and HF/6-31G* ESP.
The **single-stage** native RESP fit used four VDW shells, unrestrained H atoms,
A=0.0005, B=0.1, iteration tolerance 1e-8 and at most 200 iterations.
This is not a claim of a conventional two-stage RESP protocol. The exact settings
are in each input and native fitting output. Geometry optimization completed for
all three components; the network optimizer's intermediate recovery diagnostics
are retained. Nothing was replaced with an external OPLS-AA charge vector.

The authoritative output is the unchanged float64 ndarray returned by
`resp.resp(...)[-1]`, saved as `.npy` and lossless JSON float/repr rows.
Printed `results.out` totals are rounded and are **not** charge authority.

| Fixture | Atoms / components | Native exact decimal component totals (e) | Expected |
| --- | --- | --- | --- |
| Linear, two complete PEO DP4 chains | 60 / 2 | 1.06E-16; -2.65E-17 | 0; 0 |
| Connected network | 63 / 1 | 3.307E-17 | 0 |

These small floating-point residuals fail the existing exact Decimal provided-
charge contract. They also have nonzero sums when each binary float is expanded
exactly. This is an acceptance/serialization limitation, not a demonstrated
physical scientific failure. No tolerance, parameter, charge or historical
record is changed. No normalization, redistribution, rounding or correction is
performed. Required downstream minimization, NVT, NPT, split/resume and MD
comparisons are **not_executed**, because the charge prerequisite fails.
A separately reviewed numerical/representation decision would be needed to
resolve this conflict; this phase does not make that decision.

An independent read-only normal-equation reconstruction from archived ESP grids
and optimized geometry reproduces the native vector within about 1.5e-11 e per
atom. Grid serialization has ten decimal places; this comparison is diagnostic
and does not replace the native vector or accept a nonzero total. Its code does
not import the native fitter or production charge extractor.

## Evidence and validation semantics

`island.psi4-resp-runtime.v1` binds an explicit prefix, exact versions, conda
metadata and asset hashes. Runtime execution never falls back to system PATH.
`island.psi4-resp-execution.v1` binds graph, atom mapping, input/optimized
geometry, method/basis, worker/runtime identities, raw arrays and native
intermediates. Loading reconstructs the record from artifacts. Reordered,
rounded, truncated, stale or mutated inputs are rejected. Failed charge records
are diagnostic records, never published as validated provided charges.

Native output relocation and fresh-process reconstruction validate identical
serialized records offline. They are **quantum evidence reconstruction**, not
MD continuation. The frozen J35 external OPLS-AA (`oplsaa_provided`) control
remains distinct. Native AmberTools AM1-BCC retains its historical -0.003 e
blocked diagnostic. No J17–J35 evidence is regenerated or overwritten.

The additive receipt validator returns 1 for unmet acceptance even when its
failed-evidence integrity checks pass. `--diagnostic-only` returns 0 solely for
valid diagnostic evidence. It must never be described as scientific acceptance.
Tests use actual failed native outputs and clearly separate software mutation
controls; no all-blocked synthetic fixture counts as positive acceptance.

## Reproduction

Use a new empty prefix and output directory. Commands are from repository root:

```bash
../miniforge3/bin/conda create --yes --prefix /NEW/psi4-runtime \
  --file docs/evidence/phase_4j36/environment-explicit.txt
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python - <<'PY'
import json
from island.forcefields.ambertools.psi4_runtime import capture_runtime
r = capture_runtime('/NEW/psi4-runtime')
open('/NEW/runtime.json', 'w').write(json.dumps(r, indent=2) + '\n')
print(r['identity'])
PY
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j36 --runtime /NEW/runtime.json \
  --runtime-identity FULL_IDENTITY_PRINTED_ABOVE --output /NEW/quantum-run
```

The final command intentionally returns 1 for the retained failed gate. Preserve
its result, native files and logs. For artifact-only reconstruction, extract
`artifacts/native-quantum.tar.gz` to a new directory and add
`--reconstruct-native /EXTRACTED` with the retained declared runtime/identity.
This mode still checks the declared runtime; the receipt's offline artifact
loader does not require an installed quantum runtime.

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j36_receipt docs/evidence/phase_4j36/receipt.json
# Expected exit 1: acceptance failed.
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j36_receipt docs/evidence/phase_4j36/receipt.json \
  --diagnostic-only
# Expected exit 0 only when the complete failed-evidence record is intact.
```

Verification on implementation `ac6e3387f89c8c5b829411eaee5cda16105b695a`:
47 focused tests passed; 593 J22–J35 regressions passed; the full supported suite
passed 2,608 tests with 10 existing skips. Ruff, compileall, five execution-
environment pip checks and diff checks passed. Historical receipt validators
passed their integrity checks; historical blocked acceptance remains nonzero.
All 992 historical evidence/manifest files remained byte-identical.
Exact executed commands, environments, JUnit results and skipped-test reasons
are archived in the receipt artifacts. Implementation SHA and later evidence
commit remain distinct; resolve the evidence commit using
`git log -1 --format=%H -- docs/evidence/phase_4j36/receipt.json` to avoid a
self-referential receipt hash.

`production_validated=False`; `simulation_readiness="not_established"`.
No equilibrium, experimental agreement, Tg, mechanical, universal force-field
coverage or production-readiness conclusion follows from this work.
