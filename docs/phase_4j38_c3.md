# Phase 4J38-C3: exact-graph charge completion

C3 completes the two missing **GAFF2 software workflow** cells with newly
calculated, independently audited RESP/Psi4 charges. Native aromatic OPLS-AA
remains `source_incompatible`: its authoritative source gives −0.230 e for
each neutral C12H10 chain. No source charge was repaired or relabeled.
The overall benchmark acceptance command still exits **1**. Short-trajectory
software execution is distinct from thermodynamic equilibration and scientific
validation.

Base: `codex/phase-4j38-c2-source-complete-benchmark-paths`, fetched and verified
at `3d298cccd7f430340c82d8e069f387a352e4f19b`.
Tested implementation: `36466826c831866586b5e8385f1fff4bf923b6db`.
Later evidence commits change documentation/evidence only. No merge into main.

## Exact chemical systems

| Cell | Construction | Complete components | Atoms | Final graph SHA-256 identity |
|---|---|---|---:|---|
| GAFF2 PE | `[*:1]CC[*:2]`, DP4, explicit H, two chains | C8H18 × 2 | 52 | `1baea92607ac7d99dd85d04d55c22808c6bae9d6501824f7828088a2fca0b66b` |
| GAFF2 aromatic | `[*:1]c1ccc([*:2])cc1`, DP2, explicit H, two chains | C12H10 × 2 | 44 | `f7730225baafb0a2164cdcf3db2f8813d83d3b8d145923053b0fc4ffd193f052` |

The frozen J38 builder, seeds 2038/2039, end-group completion, molecule ancestry,
production packing, coordinates, box and atom IDs are unchanged. Both systems
have formal charge zero in each component. All atoms of each component enter
its own quantum calculation. No DP3 PE, benzene, retained PEO or different-graph
charge vector is transferred. Each component's actual conformation is optimized
and fitted independently; ISLAND does not average symmetry-related atoms or
chains.

## Fresh quantum evidence and independent audit

The pinned J36 runtime passed live executable/source/library hashes and an
import/version probe. It is Psi4 **1.9.1**, CDS RESP conda package **1.0** with
module version **0.8**, Python **3.11.14**. Its full environment and executable
fingerprints are in `raw/runtime.json` inside the artifact archive, identical
to the historical J36 environment manifest. AmberTools is the pinned **24.8**
runtime; all five execution environments are recorded separately in
`execution-environments.json` with executable hashes and package inventories.

Fresh calculations executed before the implementation freeze. The external
workers import Psi4/RESP, not ISLAND. Exact commands, logs, return codes, source
copies and SHA-256 hashes are archived. `phase_4j36_quantum_reference.py` runs
HF/STO-3G optimization (`GAU_LOOSE`, at most 100 iterations), followed by an
HF/6-31G* ESP/RESP fit. Its legacy exact-decimal gate returns 1; this is retained,
not reclassified as successful CLI termination. C3 then recomputes HF/6-31G*
ESP at those new optimized geometries through `ESPPropCalc`'s native in-memory
API and calls the pinned CDS `espfit.fit` on the native float64 arrays.
The C3 native worker records the correct fresh-geometry provenance explicitly.

The C3 raw native worker also retains the legacy exact-zero diagnostic/exit 1.
Those original failed zero-tolerance records remain unchanged; the separate
versioned J36-C2 policy supplies the numerical acceptance decision. An external
CDS RESP refit independently reconstructs inverse-distance inputs and reproduces
all four component charge arrays bit-for-bit, without rerunning SCF or using the
production importer. Its commands and results are archived separately.

The fitting settings retain J36's one-stage RESP protocol: four vdW shells
(1.4, 1.6, 1.8, 2.0), point density 1.0, RESP_A=0.0005, RESP_B=0.1,
IHFREE=True, tolerance 1e-8, at most 200 iterations. This is the declared
software acceptance protocol, not a claim of universally transferable or
experimentally validated GAFF2 charges. The official rounded ESP export is
checked for API units/semantics but never used as charge authority.

`phase_4j38_c3_reference.py` performs a separate standard-library raw audit,
independent of the production charge importer and conservation implementation:

- reconstruct graph identity and connected components from raw bonds;
- compare stable IDs, elements, formal charges and membership;
- compare original centered coordinates, optimized NPY coordinates, Psi4 input
  and printed coordinate witnesses;
- inspect actual Psi4 version, HF/basis/convergence log witnesses;
- parse native NPY binary64 bytes without NumPy and verify JSON charge bits;
- reconstruct exact rational and Decimal residuals plus compensated sums;
- independently reconstruct the unchanged J36-C2 conservation bound;
- check executable/environment, raw source, intermediate, geometry, charge and
  worker hashes against newly reviewed C3 anchors.

Trust roots are `scripts/phase_4j38_c3_anchors.json` in the tested implementation,
not self-signed charge labels. The archive is a newly generated C3 source asset,
not a copy of an earlier charge receipt. File hashes prove integrity and bind
execution provenance; they do not establish experimental validity.

## Unchanged finite-precision policy

The authoritative vector is the native `charges-native.npy` array. JSON import
is bit-exact. No normalization, redistribution, truncation, residual subtraction
or tolerance relaxation occurs. The J36-C2 bound remains
`sum(0.5*ulp(q)) + gamma_n*sum(abs(q))`, with `u=2^-53` and
`gamma_n=n*u/(1-n*u)`. Exact rational arithmetic makes the decision; Decimal
reconstruction and `math.fsum` are independent witnesses. The residual exists
in the native fitted vector, rather than being solely a naive-summation error.

| Vector | Absolute residual (e) | Allowed residual (e) |
|---|---:|---:|
| PE component 0 | 4.16334e-17 | 2.35323e-15 |
| PE component 1 | 2.16840e-18 | 3.59205e-15 |
| PE system | 3.94650e-17 | 1.17403e-14 |
| Aromatic component 0 | 2.25514e-16 | 6.65688e-15 |
| Aromatic component 1 | 1.00614e-16 | 6.58138e-15 |
| Aromatic system | 3.26128e-16 | 2.60574e-14 |

The exact values, signs, fractions, relative residuals, source precision, policy
identities and summation witnesses are serialized in the charge records and
independent audits. These are the final retained C3 fits; exploratory fits are
not substituted into this evidence. PRMTOP serialization remains subject to
its existing 1e-14 e bound and is checked against the unchanged source vector.

## Actual periodic execution

After independent source and charge checks passed, each cell underwent fresh
provided-charge GAFF2 preparation and independent assignment checks. Each then
executed minimization, 64 NVT steps and 256 NPT steps at 300 K/1 bar,
0.5 fs timestep on OpenMM 8.6.1 Reference. LangevinMiddle and isotropic
MonteCarloBarostat use the existing deterministic seeds 2033/2034/2035.
There are **321 retained frames** and **four 64-frame NPT density blocks** per
completed history; the split boundary is NPT step **128**.

Uninterrupted and independently executed split histories use the same actual
packed inputs/configuration. The split bundle is relocated, then resumed in a
fresh interpreter. Existing J33-C1 continuation independently integrates the
128-step NPT prefix to reconstruct OpenMM's adaptive barostat proposal state
before restoring the boundary checkpoint. This extra integration is real work;
no trajectory coordinates or thermodynamic values are copied from the full run.
Every retained coordinate, velocity, box, stage, energy and measured temperature,
volume/density, identity, opaque RNG/checkpoint byte sequence compares exactly.
The pressure observation is **unavailable (`null`)** in the existing backend;
1 bar is a target, never reported as measured pressure.

Both cases exhibit actual box-volume changes. Existing workflow config records
retain `cutoff_nm=0.9`; as documented for earlier phases, execution takes its
**effective 0.5 nm cutoff** from the bound periodic backend, not this unused
configuration field. C3 changes neither historical schema nor physical settings.

Independent `AmberPrmtopFile` Ewald energy/force/component construction agrees
at the first and last actual frame under the retained 1e-6 kJ/mol energy and
1e-5 kJ/mol/angstrom force bounds. Both paths use OpenMM Reference: this is an
independent construction check, **not an independent engine or LAMMPS KSPACE
validation**.

MDAnalysis reconstruction, portable trajectory export and existing selection,
PBC, mass weighting and aggregation policies are preserved. Density comes from
mass/actual frame volume. Short-run block means are 0.3039026 g/cm³ (PE) and
0.2933646 g/cm³ (aromatic), block standard errors 0.0024543 and 0.0008569 g/cm³.
These are finite-run diagnostics, not equilibrium densities or experimental
comparisons. Rg, Re/Re-to-Rg, RDF, orientation/nematic order, contacts and
connectivity diagnostics are in each `observables.json`, with selection and
sampling limitations.

## Matrix and compatibility

`matrix.json` updates only the two GAFF2 cells. Graph construction, packing,
charge completeness, typing, parameter completeness, executable workflow,
numerical reference, equilibration and scientific validation remain separate.
The five other executable J38 cells remain read-only historical evidence;
native aromatic OPLS retains its previous row unchanged. All seven executable
cells still have insufficient equilibrium sampling and blocked scientific
validation. PCFF rejection and the blocked native AM1-BCC control remain
historical controls; RESP is never labeled AM1-BCC or external OPLS-AA.

C3 adds scripts/records only; no production force-field parameters, charges of
historical records, source profiles, tolerances, J23/J32 identities, or Linux
atomic-publication behavior are changed. All **1961 tracked historical docs,
evidence files and `.gitattributes`** at the base are byte-identical.

## Reproduction and validation

From the repository root, use the declared environments and set:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export ISLAND_J33_C1_RUN=/tmp/j33-c1-validation/runs
export ISLAND_J34_RUN=/tmp/j34-validation/runs
export ISLAND_J35_RUN=/tmp/j35-validation/runs
```

Extract `docs/evidence/phase_4j38_c3/fresh-source.tar.gz` into a new directory,
then replay from its `raw` subdirectory into a fresh empty output namespace:

```bash
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c3 \
  --output /tmp/new-c3-run --from-reviewed /tmp/extracted-c3/raw
# exit 1: corrected GAFF2 software cells pass, native aromatic OPLS still blocked
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c3_receipt \
  --bundle /tmp/new-c3-run --diagnostic-only
# exit 0: strict evidence integrity; this is not overall acceptance
```

Omitting `--from-reviewed` executes fresh graph-specific quantum generation.
New raw outputs require an independently reviewed anchor update before the
worker allows MD; the driver cannot silently trust a newly self-signed vector.
The exact original quantum/RESP command lines and logs are in the source archive;
repeat them into new directories rather than overwriting retained evidence.
Fresh electronic-structure runs may differ in low-order floating bits; this
phase's exact continuation assertion concerns MD from the same validated bundle.

`validation/*-command.json` records exact test/check argv, environments, tested
SHA, exit codes and durations; JUnit files record actual counts. The runner
sources are archived as inert `.py.txt` files. Tests exercise actual subprocess
exit codes, raw-source mutations, independently re-signed provenance/policy
mutations, wrong graph/order/coordinates/formal charges, missing runtime,
relocation, changed MD frames/checkpoints, and failed child continuation.
Injected and mutated records are software negative controls, never scientific
evidence. The supported suite includes its existing short integration runs;
none are disabled. Detailed final counts are in the additive receipt. An initial full-suite run
was deliberately terminated to correct missing-reviewed-asset classification;
its nonzero termination and earlier passing checks are archived separately in
`pre-classification-fix-validation.tar.gz`. Required final gates were rerun
against the final tested implementation; the interruption is not a passing gate.

Receipt acceptance exits 1; `--diagnostic-only` exits 0 after full integrity
validation. An invalid/mutated receipt exits 2. A structurally readable failed
source record is never enough to prove corrected-cell software acceptance.
Historical receipt validators retain their original blocked/nonzero semantics.

## Final verification results

| Gate | Actual result |
|---|---|
| C3 focused | 49 passed |
| J22–J38-C2 regressions | 837 passed |
| Supported full suite | 2854 passed, 10 skipped; zero failures/errors |
| Ruff / compileall / diff check | passed |
| Five declared environment pip checks | passed |
| Historical receipt semantics / 1961 historical byte hashes | passed |
| Independent raw audit / relocated fresh-process bundle reconstruction | passed |

Six consistently re-signed receipt/manifest software negative controls were
rejected with exit 2; their actual diagnostics are in
`software-negative-controls.json`. They are not scientific evidence.

The 10 full-suite skips are unchanged optional legacy gates in the established
supported configuration. Exact reasons are preserved in JUnit: one real
AmberTools/QM opt-in (`ISLAND_RUN_AMBERTOOLS_INTEGRATION`); seven archived
AmberTools reference cases needing `ISLAND_AMBERTOOLS_REFERENCE_MANIFEST`;
one retained PE analysis needing `ISLAND_ANALYSIS_REFERENCE_ROOT`; and one
noncontiguous OPLS source test needing `ISLAND_OPLSAA_XML` in its prescribed
environment. Those optional gates were **not executed** in this full-suite run;
they are not promoted to success. No marker or configuration was disabled for
C3. The focused acceptance executions and historical validators above provide
the separately recorded current source/workflow checks.

## Remaining limits

The native aromatic OPLS source still needs an independently established,
chemically neutral source for this exact final graph. C3 introduces no external
OPLS charge substitution. Native AM1-BCC remains blocked by its documented
precision limitation. No universal OPLS-AA/GAFF2/PCFF or msi2lmp coverage,
experimental agreement, Tg, mechanical properties, thermodynamic equilibrium,
LAMMPS Ewald/PPPM, or production readiness is established.

`production_validated=False`

`simulation_readiness="not_established"`
