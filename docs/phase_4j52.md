# Phase 4J52: final PDMS convention audit for the pinned supplied corpus

The nine unresolved conventions inherited from J51, in their original order, are:

1. `dielectric`;
2. `coulomb_prefactor`;
3. `periodic_image_policy`;
4. `precision`;
5. `neighbor_list_policy`;
6. `pressure_convention`;
7. `virial_convention`;
8. `pair_table_bits`;
9. `erfc_kernel`.

**All nine are `unresolved_permanently` for the exact supplied, hash-bound
corpus.** This records exhaustion of that corpus. New author evidence can
supersede this result through a new audited version. It does not mean that
future evidence is impossible. Each diagnostic backend implementation remains
separately `backend_resolved`; that status does not establish author authority.
Research minimization/NVT and production execution remain blocked.

Fetched base: `78a516cdf9f373b905ad8fd729bc82453ac8a886`, the exact
`origin/main` containing merged J51 (#75). Branch:
`codex/phase-4j52-opls-pdms-final-convention-resolution`.
All tracked base files, including J17–J51 receipts, source identities, bundles
and implementation files, remain byte-identical. This phase is additive.

## Final nine-convention table

| Convention | Evidence and scope | Exact missing artifact | Final author result |
| --- | --- | --- | --- |
| dielectric | OPLS input has no declaration; `dielectric 1.0` belongs to Dreiding; present engine definition is retained separately | DFT/NPA-scoped declaration or author executable/restart proving the value | `unresolved_permanently` |
| Coulomb prefactor | `units real` is explicit; the current engine uses 332.06371 kcal Å/(mol e²), divided by dielectric | Author engine revision/build and its electrostatic constant | `unresolved_permanently` |
| periodic images | Paper p3 explicitly states PBC in x/y/z; pinned backend documents its conducting sum | DFT/NPA surface/slab/image-sum policy | `unresolved_permanently` |
| precision | Supplied data contain binary64-formatted charge literals; current diagnostic is CPU double precision | Author binary hash, build precision, compiler and accelerator manifest | `unresolved_permanently` |
| neighbor list | Generic inputs include `check yes`; viscosity example has `neighbor 2.0 bin`; selected OPLS entry point does not include those drivers | Complete driver binding skin, rebuild and sorting policy to the DFT/NPA branch | `unresolved_permanently` |
| pressure | Generic scripts request `press`/`pxy`; current diagnostic explicitly computes a configurational proxy | DFT/NPA estimator including kinetic/constraint terms | `unresolved_permanently` |
| virial | Current engine tensor/conversion is auditable; no author tensor/image/kinetic contract is supplied | DFT/NPA tensor, contribution and conversion definition | `unresolved_permanently` |
| pair-table bits | No OPLS author interpolation setting; table changes have a measured numerical effect | Author `pair_modify table` declaration or exact executable | `unresolved_permanently` |
| erfc kernel | Present original polynomial and explicitly modified J49 exact-erfc implementation are pinned | Author binary/revision establishing energy/force kernel | `unresolved_permanently` |

The complete machine table is
`docs/evidence/phase_4j52/audit/final-nine-conventions.json`. Each row contains
all searched file hashes, exact page/line/raw matches, DOI/citation, interpreted
value (null where absent), units, conversion status/version, independent
reconstruction, checksum, missing artifact, affected modes/observables and the
reason a default would be unjustified. The permanent-blocker report is
`audit/permanent-blockers.json`. Lack of author evidence has no calculable
unique numerical effect; no hypothetical author value is invented to measure it.

J51 already resolved units, LJ form/parameter conversion, torsion order and
half-factor, special pairs/exclusions, geometric LJ mixing, switching/cutoffs,
PPPM tolerance and x/y/z boundary conditions. Their complete 37-setting table
and four authority layers remain unchanged in each workflow. J51's documented
CHARMM no-tail restriction remains `converter_resolved`; no zero parameter is
inserted. The nine-item audit does not reopen those resolved declarations.

## Evidence searched

The independent offline reader searches all bytes/lines of **64** explicitly
supplied or pinned file objects: all **32** author ZIP members, **3** PDFs
(37 pages), **10** converter/export/build files and **19** backend documents
and source files. It records every query, every match, empty search results,
file size, SHA-256, page/line count, extracted-text hash and scope. All **90**
relevant J47–J51 JSON comparison/result/receipt records are read and hashed.
No other force-field paper, network discovery or downloaded default is used.

Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat
Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*,
JPCB 129 (2025), 1864–1873, DOI `10.1021/acs.jpcb.4c08471`:

| Supplied asset | SHA-256 |
| --- | --- |
| `pcff/PDMS_benchmark-main.zip` | `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4` |
| `pcff/jp4c08471.pdf` | `04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc` |
| `pcff/jp4c08471_si_001.pdf` | `6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f` |
| `pcff/ja11225.pdf` (original OPLS supporting tables) | `451cf754daab3fcf6e357a3cdfb64e5d8ab8eb800ce8afdad53897c0bbedcc5a` |

The last asset identifies Jorgensen/Maxwell/Tirado-Rives, JACS 118 (1996),
11225–11236, DOI `10.1021/ja9621760`. Its general OPLS tables cannot establish
the later PDMS DFT/NPA execution conventions. It supplies no new parameter.
The author's original Moltemplate version/source is not supplied. The pinned
ISLAND conversion implementation is audited without claiming that it is the
author's converter version. Source PDFs/ZIPs are not redistributed.

SI p6 Table S4 says `pair_style lj/cut/coul/long 11 11`; the selected author
OPLS include, line 8, says `pair_style lj/charmm/coul/long 9.0 11.0`.
The source does not reconcile those scopes for the serialized DFT/NPA branch.
No precedence, combination or source-model substitution is invented. Same-scope
competing claims remain conflicting under the unchanged J51 reconciler.
The paper's Qeq thermophysical validation is not transferred to DFT/NPA.

## Public API and execution policy

The profile remains
`island_external_opls_pdms_xiang2025_dft_local_v1`: external literature-derived
OPLS-compatible, limited to its declared TMS-capped PDMS graph domain. It is
not native universal OPLS-AA, PCFF or COMPASS PDMS. No charge method, physical
parameter, tolerance or historical public facade is changed.

Install the explicitly supplied local assets once; ordinary selection is only
`forcefield="opls"`:

```python
from island.polymer import OPLSAssets, use_opls_assets
from island.opls_final_resolution import ConventionSources, use_convention_sources
from island.opls_final_execution import prepare_polymer_forcefield

with use_opls_assets(OPLSAssets(
    "island_ffref/pcff/PDMS_benchmark-main.zip",
    "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4",
)), use_convention_sources(ConventionSources(
    "island_ffref/pcff/jp4c08471.pdf",
    "island_ffref/pcff/jp4c08471_si_001.pdf",
)):
    prepared = prepare_polymer_forcefield(system, forcefield="opls")
report = prepared.payload
```

No profile or charge-method argument is required. The immutable result exposes
family, internal profile, original DFT/NPA charge origin, exact domain,
source/converter/backend/scientific settings, final unresolved rows, execution
mode and readiness. It binds the unchanged canonical assignment and J51
identity to the complete-corpus audit identity. Runtime integrity checks rehash
the original ZIP and both Xiang PDFs. The ordinary runtime needs no PDF parser.

| Mode | Current permission |
| --- | --- |
| `diagnostic` | Explicit pinned static/component/FD comparisons; no MD authority |
| `research_opt_in` | Disabled, including explicit opt-in; exact nine blockers returned before engine access or output creation |
| `production_candidate` | Disabled; source, Ewald, pressure, sampling and cross-engine gates remain open |

`unresolved_permanently` cannot be re-signed into authority or cleared through
a backend policy selection. Neither diagnostic policy authorizes minimization,
NVT or NPT. No new checkpoints, trajectories or RNG continuation are claimed.
Historical J49 diagnostic continuation remains preserved without promotion.

## Repeated numerical evidence

HMDS and capped DP2/DP4/DP6 each have two actual static/FD workflows:
`source_replay_v1` and `exact_erfc_diagnostic_v1`. Every workflow has repeated
deterministic input generation, total/component energies, all-atom forces,
central differences at 2e-4/1e-4/5e-5 Å, independent special-pair partitions,
affine virial records and native configurational pressure proxy. ZIP/PDF and
workflow relocation is reconstructed in fresh processes at hash seeds 1 and
99397. The independent reader reconstructs all **288** raw FD samples, the
17 source rows, graph/type/mass/charge mapping, interaction enumeration,
eight special-pair/pressure partitions and both exact-kernel source edits.

Physical inputs and deterministic exports are byte-identical to J51: graph,
coordinates, atom order, masses, charge literals, 40 Å periodic box,
exclusions, 0.5 special-pair scaling, geometric LJ mixing, CHARMM 9–11 Å
switching, Coulomb 11 Å cutoff and unit conversions. Pressure/virial remain
proxies without an inherited scientific acceptance bound. Native combined
LAMMPS pair style does not expose standalone 1–4 or Coulomb force components;
physical equation partitions and force subtraction are labeled explicitly.
No output correction is applied.

The source-command replay uses PPPM `1e-4`, unlike OpenMM Reference Ewald
`1e-7`; this is an explicitly unequal numerical contract. The separate J49
control uses explicit table=0, Ewald `1e-7`, fixed alpha/reciprocal extent and
the locally modified exact-erfc kernel. This diagnostic control cannot become
source authority. Both engines retain their own Coulomb constants unchanged.

The existing bounds remain energy 3e-4 kJ/mol, force 3e-4 kJ/(mol Å), and FD
1e-4 kJ/(mol Å). Source replay still fails the Coulomb/FD bounds; the explicit
exact-erfc control passes. Diagnostic acceptance is zero, research/production
acceptance remains nonzero. There is no Ewald/PPPM closure or pressure claim.

Pinned original LAMMPS executable:
`../island-validation/phase4j22/lmp-kspace`, SHA-256
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
Pinned exact-erfc diagnostic executable: `/tmp/island-j49-exact/build/lmp`,
SHA-256 `3a989d8e36345bd4508d932832a78e0eca3313b134d5c414c280065d615029ea`.
The recorded backend is LAMMPS 10 Sep 2025 development, revision
`36f7c81ced8dd456430b13647e48f0c7d7bd8faa`, serial CPU double, KSPACE/MOLECULE.
Version/help/build information and hashes are retained per workflow.

## Reproduction

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J52_PY=../island-validation/phase4g1-env/bin/python

"$J52_PY" scripts/phase_4j52.py audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --supplement island_ffref/pcff/jp4c08471_si_001.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --exact-kernel /tmp/island-j49-exact/build/lmp \
  --exact-kernel-sha256 3a989d8e36345bd4508d932832a78e0eca3313b134d5c414c280065d615029ea \
  --output /tmp/j52-audit
# Expected exit 1: nine final source blockers, no research MD authorization.

# Explicit offline pypdf 6.20.0 reader; no network or runtime discovery.
python -m scripts.reconstruct_phase_4j52_evidence \
  --directory /tmp/j52-audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --references island_ffref/pcff --backend ../lammps/lammps \
  --exact-source /tmp/island-j49-exact/source \
  --output /tmp/j52-independent.json

"$J52_PY" -m pytest -q tests/test_phase_4j52.py
"$J52_PY" -m pytest -q tests/test_phase_4j49.py tests/test_phase_4j50.py tests/test_phase_4j51.py
"$J52_PY" -m pytest -q
ruff check src tests scripts
"$J52_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$J52_PY" scripts/validate_phase_4j52_receipt.py docs/evidence/phase_4j52 --repo .
"$J52_PY" scripts/validate_phase_4j51_receipt.py docs/evidence/phase_4j51 --repo . \
  --sha256 2e7f4838ad01ef2994b395a22c354cc807fb6b30b64d2a8573ebe22adfdacf8a
"$J52_PY" scripts/check_phase_4j52_receipt.py docs/evidence/phase_4j52 --repo . \
  --output /tmp/j52-receipt-controls.json
```

The exhaustive source-only search can be reproduced with
`python scripts/audit_phase_4j52_sources.py --archive
island_ffref/pcff/PDMS_benchmark-main.zip --references island_ffref/pcff
--backend ../lammps/lammps --output /tmp/j52-corpus.json --verify
docs/evidence/phase_4j52/audit/supplied-corpus.json`.
The receipt/manifest can be deterministically sealed from the completed audit,
independent reconstruction, collection and validation artifacts with
`python scripts/seal_phase_4j52_receipt.py docs/evidence/phase_4j52 --repo .`.
Exact commands and disjoint full-suite shard membership are retained in evidence.

Verification: **83 focused tests passed**, **237 J49–J51 regressions passed**,
and **4,314 ordinary-suite tests passed with 61 existing skips** (4,375 nodes,
no failures/errors). Ruff, compileall, five pip checks, historical receipt
validators and all 4,897 historical file hashes passed. Independent reconstruction
verified 2,102 raw convention matches, 304 backend citation lines, all 288 FD
samples and 1,017 J47–J51 raw files. Workflow reconstruction used 16 relocated
fresh processes. Re-signed convention/FD/readiness tampering and relocated
receipt validation are retained separately in the validation evidence.

Readiness stays `production_validated=False` and
`simulation_readiness="not_established"`. No native universal OPLS-AA, universal
PDMS/PCFF coverage, production MD, pressure validation, thermodynamic equilibrium
or long-NPT readiness is established.
