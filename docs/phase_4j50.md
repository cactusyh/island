# Phase 4J50: source conventions and explicit LAMMPS execution authority

J50 separates the source model from converter equations, backend numerical
choices and scientific verification. The resulting identity includes **four
separately checksummed, immutable records**, the canonical assignment, exact
graph/order/coordinates/box and source-bound charges. Integrity validation
rebuilds the records from the pinned local archive and rejects re-signed edits.

**Research LAMMPS MD remains blocked.** Explicit opt-in cannot supply missing
source evidence. All four fixtures execute diagnostic static/component/force/FD
comparisons, and relocate/reconstruct in fresh processes. No J50 minimization,
NVT, NPT, checkpoint or trajectory is produced. J49's bounded diagnostic MD
evidence remains historical evidence; it is not promoted to J50 research MD.

`production_validated=False`, `simulation_readiness="not_established"`.
Pressure/virial qualification, polymer Ewald/PME, cutoff recovery, sampling,
long NPT and scientific equilibrium remain unestablished.

## Exact base and preserved model

Fetched `origin/main` and exact branch base:
`b207a3741f0b6c2e1a6c632a283d82fefdd5e0d9`, merged J49.
Branch: `codex/phase-4j50-opls-pdms-lammps-convention-reconciliation`.
Every tracked base file is byte-preserved, including J17–J49 receipts, bundles,
serialized identities, parameters, native sources and implementation files.
J50 adds a policy facade without modifying historical canonical assignments.

The public family is **`forcefield="opls"`**. The resolver still selects
`island_external_opls_pdms_xiang2025_dft_local_v1` automatically, with its exact
charge origin, source hashes, internal version and bounded TMS-capped siloxane
domain. It is external literature-derived OPLS-compatible local research data.
It supplies neither native universal OPLS-AA nor PCFF/COMPASS PDMS support.

Source: Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density,
Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*,
JPCB 129 (2025), 1864–1873, DOI `10.1021/acs.jpcb.4c08471`.
The serialized author **DFT/NPA** branch is retained; the paper's Qeq results
are not scientific validation of this branch. No charge, coefficient, acceptance
bound or physical profile is changed. No parameter row is added or zero-filled.
Authorization remains user-requested local research; redistribution permission
is unestablished. Source ZIP/PDF bytes are not redistributed.

## Four authority layers

`source_authoritative_settings` records literal source commands, documented
interpretations of those commands, their exact file hashes/line/raw text and
absence evidence. `converter_settings` records the equation-preserving
real-unit/OpenMM conversions and pinned implementation/engine definitions.
`backend_execution_settings` explicitly records each selected numerical value
and its owner. `scientific_verification_settings` records unchanged acceptance
bounds and the still-open scientific gates. Changing any layer changes the
combined identity; merely recomputing checksums cannot authorize execution.
Prior diagnostic success is tied to the exact J49 receipt and its four fixtures;
preparation alone never marks a new graph numerically verified.

Each setting has exactly one of the six requested classifications:
source-stated, converter-stated, family-inherited, documented derivation,
backend-selected, or unresolved. No family inheritance is needed in this audit.
The source record never promotes backend-selected values or unresolved fields
to source authority. Units and style equations are documented interpretations
of the author's commands, not independent author quotations.

The pinned OPLS include chain explicitly specifies `units real`, harmonic
bonds/angles, OPLS dihedrals, harmonic improper style, CHARMM LJ/Coulomb-long,
geometric mixing, `special_bonds lj/coul 0 0 0.5`, 9/11 Å cutoffs and PPPM
`1e-4`. Bounded topology contains no impropers; no fictitious improper row is
introduced. Sigma denotes the LJ zero crossing, not r-star. The OPLS torsion
retains its four half-factored terms.

The archive has 32 files. Other example scripts explicitly write neighbor
settings and pressure output, but include generic
`system.in.init/system.in.charges`; the pinned OPLS DFT path instead includes
`opls_system.in.init`. Those examples do **not** establish the conventions of
this DFT/NPA branch. Exact occurrences are retained with an unproven-scope
diagnostic. Independently, the paper's Methods on page 3 explicitly states
periodic boundary conditions in x, y and z for all simulations. That establishes
`p p p` directly, with a hash/page/raw-text citation. No COMPASS/Qeq settings,
charges or thermophysical-validation claims are borrowed.

Ten source declarations remain unresolved: dielectric, Coulomb prefactor,
periodic-image policy, precision, neighbor-list policy,
pressure convention, virial convention, pair-table bits, erfc kernel and LJ tail
correction. Each diagnostic identifies the missing OPLS DFT-scoped declaration
and the searched include-chain hashes. Backend values, including dielectric 1,
periodic conducting images and a configurational pressure proxy, are explicit
diagnostic policy choices. They do not close these source gaps.

The machine-readable source/converter/backend table is
`docs/evidence/phase_4j50/audit/source_replay_v1-convention-table.json`;
the exact-erfc table records its numerical overrides separately. Each row also
references the converter and scientific-record identities. All 37 required
settings appear in these tables.

## Execution modes and public use

| Mode | Permission | Current result |
| --- | --- | --- |
| `diagnostic` | Static/component comparisons after every selected backend setting resolves | Available; no MD authority |
| `research_opt_in` | Short minimization/NVT only after explicit caller opt-in, complete source/backend settings and independent qualification | Blocked by source conventions |
| `production_candidate` | Existing Ewald, pressure, cutoff, sampling and cross-engine gates must all pass | Unavailable |

The generator gates before creating files. Execution gates before inspecting or
starting an engine. It emits explicit dielectric, mixing, special bonds,
cutoffs, torsion style, solver/tolerance, precision-bound output and pressure
compute commands. An unresolved selected backend value cannot reach LAMMPS.
The diagnostic pressure compute is `pressure NULL virial`: no kinetic term,
six ordered tensor components, in atm, with an explicit conversion to a virial
proxy. This is not validated thermodynamic pressure.

The additive J50 facade preserves ordinary canonical selection:

```python
from island.polymer import OPLSAssets, use_opls_assets
from island.opls_execution import prepare_polymer_forcefield, run_static, LAMMPSOptions

assets = OPLSAssets(
    "island_ffref/pcff/PDMS_benchmark-main.zip",
    "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4",
)
with use_opls_assets(assets):
    prepared = prepare_polymer_forcefield(system, forcefield="opls")
report = prepared.payload
# selected_profile, charge_origin, domain, four convention records,
# unresolved_settings, execution_mode and scientific readiness are visible.
# prepared.assignment is the unchanged canonical J48 assignment.

result = run_static(system, prepared, "new-diagnostic-directory", engine=LAMMPSOptions(
    "../island-validation/phase4j22/lmp-kspace",
    "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191",
))
```

No profile or charge-method argument is required. Assets are installed explicitly
once; no network access, file discovery or family fallback occurs. The old
`island.polymer` and `island.opls_lammps` identities/boundaries remain unchanged.

Requesting `execution_mode="research_opt_in"` and passing
`research_opt_in=True` to `run_md` still raises `source_conventions_unresolved`,
listing every missing declaration and available evidence. Without that Boolean
it raises `explicit_research_opt_in_required`. Diagnostic mode rejects MD;
production-candidate mode rejects every execution. Opt-in is necessary and never
sufficient. No NPT adapter is authorized.

## Reconciliation and measured effect

HMDS and TMS-capped DP2/DP4/DP6 use identical final graphs, atom order,
coordinates, masses, charge literals, 40 Å boxes, exclusions, 0.5 LJ/Coulomb
1–4 scaling, geometric mixing, CHARMM 9–11 Å LJ switching and unit conversions.
The solver/kernel differences are explicit experimental variables, so these
runs are not represented as identical-numerical-contract research acceptance.

| Fixture | Source-command replay Coulomb error, kJ/mol | Explicit exact-erfc backend error, kJ/mol | Exact backend max force error, kJ/(mol Å) |
| --- | ---: | ---: | ---: |
| HMDS | 0.0699658923 | 0.00000978363 | 0.00000190698 |
| DP2 | 0.1436606402 | 0.0000100741 | 0.00000262941 |
| DP4 | 0.4018122427 | 0.00000949844 | 0.00000219366 |
| DP6 | 0.5085865935 | 0.0000147126 | 0.00000388548 |

Original source-command replay retains PPPM `1e-4` and the pinned kernel/table.
The explicit backend control retains J49's local exact-erfc executable,
`table 0`, Ewald `1e-7`, matched splitting alpha and `kmax 20 20 20`.
These are **representational numerical changes**, not source-authoritative
settings or new physical parameters. Coulomb real/reciprocal, total energy,
forces and FD are affected. Missing source scope/pressure authority remains
**unresolved**. Charges and parameters are unchanged; no output correction is
applied. The unchanged energy/force bounds are `3e-4` and FD bound is `1e-4`.

The Coulomb prefactors remain LAMMPS `332.06371` kcal Å/(mol e²) and OpenMM
`1389.35457644382` kJ Å/(mol e²); their small precision difference is retained.
No source-authority claim is made for either engine-owned constant.
The combined pair style has no standalone per-atom Coulomb-force observable:
the report explicitly labels total-minus-independent-bonded/LJ decomposition.

All eight workflows have deterministic data/input, actual pinned LAMMPS runs,
raw component energies/forces/virial proxies, sampled three-step-size FD,
relocation and two fresh-process reconstructions under different hash seeds.
The independent audit reconstructs archive and PDF citations, engine-definition
citations, 288 FD samples and raw force dumps outside production parsers.
Only the explicit local PDF check uses `pypdf` (6.20.0 in this environment);
the numerical log reconstruction uses the standard library.
Research minimization/NVT/continuation/checkpoint/RNG acceptance is **not run**
because its prerequisite source gate is open. Diagnostic success is separate:
diagnostic gate 0; research-opt-in and production-candidate gates 1.

## Reproduction and validation

Use explicit local assets and new output directories:

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J50_PY=../island-validation/phase4g1-env/bin/python

"$J50_PY" scripts/phase_4j50.py audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --exact-kernel /tmp/island-j49-exact/build/lmp \
  --exact-kernel-sha256 3a989d8e36345bd4508d932832a78e0eca3313b134d5c414c280065d615029ea \
  --output /tmp/j50-audit
# Expected exit 1: unresolved source/research authority, even if diagnostics pass.

python scripts/audit_phase_4j50_rows.py --directory /tmp/j50-audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --lammps-source ../lammps/lammps --paper island_ffref/pcff/jp4c08471.pdf \
  --output /tmp/j50-independent.json

"$J50_PY" -m pytest -q tests/test_phase_4j50.py
"$J50_PY" -m pytest -q tests/test_phase_4j49.py
"$J50_PY" -m pytest -q
ruff check src tests scripts
"$J50_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$J50_PY" scripts/validate_phase_4j50_receipt.py docs/evidence/phase_4j50 --repo .
"$J50_PY" scripts/validate_phase_4j49_receipt.py docs/evidence/phase_4j49 --repo . \
  --sha256 8bcff997fc807c7fc45cadb9aca3a5357638d844bebfb79de3c8aa954fb5c4e9
```

The exact-erfc executable is the explicitly modified J49 diagnostic build;
`scripts/build_phase_4j49_kernel.py` and its checksummed build receipt document
reconstruction. Compiler/path changes may change executable bytes; this J50
policy pins the measured build, and never discovers or substitutes another one.

Focused tests cover all four fixtures and all requested convention, family,
archive, graph, atom-order, box and charge mutations; re-signing layer identities
cannot bypass integrity. Verification totals and all full-suite node identities
are recorded in the receipt, JUnit files and reproducibility manifest.

Final verification: **61 J50 focused tests passed**, **52 J49 regressions
passed**; the complete collected suite has **4168 tests: 4107 passed and 61
existing skips**, with exact node accounting across nonoverlapping file shards.
Ruff, compileall, five pip checks, source/PDF/FD reconstruction and all 4424
historical-file checks pass. Research and production acceptance remain nonzero.
