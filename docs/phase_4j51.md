# Phase 4J51: scoped external PDMS source-convention closure

The ten unresolved **source** declarations inherited from J50 are, in order:

1. dielectric;
2. Coulomb prefactor;
3. periodic-image policy;
4. precision;
5. neighbor-list policy;
6. pressure convention;
7. virial convention;
8. Coulomb pair-table bits;
9. erfc kernel;
10. LJ tail correction.

**One closes by a documented converter rule; nine remain source-unresolved.**
The selected CHARMM switched LJ style explicitly cannot use tail corrections.
This is a functional-form restriction, not inference from a default flag and
not a new parameter or inserted zero term. The remaining nine cannot be proven
from the supplied DFT/NPA-scoped source. Their pinned diagnostic backend
implementations are documented independently. Research LAMMPS MD remains
blocked, including explicit opt-in. No J51 minimization, NVT, NPT, MD trajectory,
checkpoint or RNG-continuation claim is made.

Fetched exact `origin/main`: `72dbb07a0d61a15436a04720357d76f63e41d927`,
merged J50 (#74). Branch:
`codex/phase-4j51-opls-pdms-source-convention-closure`.
All tracked base files remain byte-identical, including implementations and
every J17–J50 receipt, bundle and serialized identity. J51 is additive.

## Ten-item resolution table

| Setting | Evidence checked | J51 source/converter result | Explicit diagnostic backend result |
| --- | --- | --- | --- |
| dielectric | Pinned OPLS includes contain no declaration; Dreiding's `dielectric 1.0` has another scope | `unresolved` | `backend_resolved`: explicit 1; `force.cpp` is backend evidence only |
| Coulomb prefactor | `units real` is explicit; author's exact engine implementation/version is not | `unresolved` | `backend_resolved`: 332.06371 kcal Å/(mol e²), divided by dielectric; OpenMM's constant is retained separately |
| periodic images | Paper page 3 states x/y/z PBC; no DFT-scoped surface/slab/image-policy declaration | `unresolved` | `backend_resolved`: 3-D conducting sum, no slab correction; pinned KSPACE code |
| precision | Archive does not establish author binary/build precision | `unresolved` | `backend_resolved`: pinned serial CPU double build, `%.17g` output; build/help captured per run |
| neighbor list | Generic example scripts say `check yes`; no binding to serialized OPLS DFT include chain | `unresolved` | `backend_resolved`: explicit 2 Å bin skin, rebuild every step, `check no`, no sorting |
| pressure | Generic examples request pressure; no DFT-scoped pressure definition | `unresolved` | `backend_resolved`: `compute pressure NULL virial`; configurational proxy only |
| virial | No DFT-scoped tensor/kinetic/image convention | `unresolved` | `backend_resolved`: xx/yy/zz/xy/xz/yz, pinned atm-to-virial conversion; scientific gate open |
| pair table | No author interpolation/bit declaration | `unresolved` | `backend_resolved`: explicit 12 for source-command replay, 0 for exact-erfc control |
| erfc kernel | No author binary/kernel/version | `unresolved` | `backend_resolved`: pinned original polynomial versus explicitly modified J49 `std::erfc` diagnostic build |
| LJ tail | Exact `lj/charmm/coul/long` command plus pinned style documentation forbidding tail correction | **`converter_resolved`**: no extra LRC; not source-authoritative | `backend_resolved`: explicit `tail no` |

Both policy-specific complete tables are under
`docs/evidence/phase_4j51/audit/*-ten-conventions.json`. All **37** conventions
are retained, with exact citations, raw text, hashes, units, interpreted values,
conversion formula/version, reconstruction method and evidence checksum.
J50's four records remain unchanged inside the new result; its original ten-gap
list is historical evidence. J51's new scoped table records the one closure and
nine remaining gaps without rewriting that list.

For an unresolved author value, its unique numerical effect cannot be calculated.
The audit does not assign a hypothetical author value to measure it. Measured
effects of table/kernel/solver differences are retained component by component.
Pressure/virial differences remain scientifically unresolved.

## Source and scope review

The existing profile is unchanged:
`island_external_opls_pdms_xiang2025_dft_local_v1`.
It remains external literature-derived OPLS-compatible local research data,
limited to its declared TMS-capped PDMS domain. It supplies no native universal
OPLS-AA or PCFF/COMPASS PDMS parameters. Canonical selection remains
`forcefield="opls"`, with automatic original DFT/NPA charge selection.

Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat
Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*,
JPCB 129 (2025), 1864–1873, DOI `10.1021/acs.jpcb.4c08471`:

| Asset | SHA-256 |
| --- | --- |
| `pcff/PDMS_benchmark-main.zip` | `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4` |
| `pcff/jp4c08471.pdf` | `04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc` |
| `pcff/jp4c08471_si_001.pdf` | `6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f` |

The archive's OPLS input specifies `real`, harmonic bonds/angles, four-term
half-factored OPLS torsions in coefficient order K1/K2/K3/K4, geometric sigma
and epsilon mixing, `special_bonds lj/coul 0 0 0.5`, CHARMM LJ 9–11 Å switching,
and PPPM `1e-4`. Explicit 1–2/1–3 exclusions and 1–4 scales follow those commands,
not customary OPLS defaults. Sigma is the LJ zero crossing, with
`r_min=2^(1/6)*sigma`. Real-to-OpenMM distance/energy conversions and harmonic
half-factor conversions remain unchanged. No improper interaction is declared
for these bounded graphs; no fictitious coefficient is added.

SI **page 6, Table S4**, instead says `pair_style lj/cut/coul/long 11 11`.
This differs from `param/OPLS/opls_system.in.init` line 8,
`pair_style lj/charmm/coul/long 9.0 11.0`. SI's OPLS-AA table is not explicitly
mapped to the archive's serialized DFT/NPA branch. J51 records both statements
and the missing branch reconciliation; it imports neither the SI charge table
nor its pair-style instruction into the registered profile. Treating both as
claims for one execution contract produces `conflicting` and rejects execution.
No source-precedence rule is invented. The paper's Qeq property results do not
scientifically validate this DFT/NPA branch.

The independent reader re-reads the relevant Xiang/Sun originals, including
Sun's polysilane paper (`ma00107a006.pdf`, DOI `10.1021/ma00107a006`), Sun's
COMPASS paper and aromatic-polyester CFF93 paper, and Smith/Borodin/Smith's
separate PDMS model (`jp047434r.pdf`). Their chemical domains/functional forms
cannot provide missing OPLS DFT execution conventions. Six source reviews and
hashes, exact excerpts and J47–J50 comparison-record hashes are retained in
`independent-reconstruction.json`. No cross-family row is imported.
Only supplied local files are read. No discovery or download occurs. Local
research authorization does not establish redistribution rights; no source
PDF/ZIP is committed.

## Resolution and public execution boundary

`island.opls_resolution` adds immutable, replay-validated records with statuses
`source_resolved`, `converter_resolved`, `backend_resolved`, `derived_resolved`,
`unresolved`, and `conflicting`. Evidence classifications are source statement,
source equation/table, source input/code, converter-defined, backend-defined,
independently derived, and unresolved. Equal claims retain all citations;
different same-scope values retain both candidates and reject execution.
The reconciler is an audit primitive, not a parameter-registration mechanism.

The new identity includes the unchanged canonical assignment, exact graph,
coordinates/order/box, charges and all four J50 layers, plus the source PDFs,
resolution table, execution mode and scope review. Re-signing altered fields
cannot supply authority. Both PDFs and the original ZIP are rehashed during
integrity validation. Source installation is explicit and can be done once:

```python
from island.polymer import OPLSAssets, use_opls_assets
from island.opls_resolution import ConventionSources, use_convention_sources
from island.opls_resolved_execution import prepare_polymer_forcefield

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

No profile or charge-method argument is needed. The result exposes selected
profile, charge origin, domain, source/converter/backend/scientific records,
every unresolved convention, mode and readiness. Historical facades retain
their existing byte-identical behavior.

`diagnostic` permits pinned static comparisons only. `research_opt_in` requires
explicit opt-in plus resolved source/converter and backend conventions and
independent qualification; it currently rejects before engine access/files,
listing all nine remaining source gaps and the SI/archive scope issue.
`production_candidate` remains unavailable. Neither mode permits NPT.

## Numerical evidence and conditional acceptance

HMDS and capped DP2/DP4/DP6 each have two new actual pinned LAMMPS diagnostic
runs, independent repeat data/input export, all-atom force comparison,
component FD at 2e-4/1e-4/5e-5 Å, virial/proxy records and two relocated fresh
process reconstructions (hash seeds 1 and 99397), with ZIP and PDFs relocated.
The eight workflows retain raw input/log/force dumps and all 288 FD samples.

Physical inputs match: final graph, coordinates, atom order, masses, charges,
40 Å periodic box, exclusions, 0.5 special-pair scales, geometric LJ mixing,
CHARMM switch, 11 Å Coulomb cutoff, unit conversions and double precision.
Original PPPM `1e-4` differs from OpenMM Reference Ewald `1e-7`; it is explicitly
an unequal numerical contract, not an identical-solver acceptance claim. The
separate exact-erfc/Ewald control records each J49 numerical override and the
two engines' retained Coulomb constants. It cannot replace source authority.

The unchanged bounds are 3e-4 kJ/mol energy, 3e-4 kJ/(mol Å) force, and 1e-4
kJ/(mol Å) FD. Source-command replay still fails Coulomb/FD; the explicit J49
exact-erfc control passes those numerical bounds. Full bonded/angle/torsion/LJ,
real/reciprocal Coulomb and total components are retained. 1–4 Coulomb/LJ
physical partitions are independently reconstructed and clearly labeled:
the combined LAMMPS style has no native standalone 1–4 or Coulomb-force dump.
Coulomb force subtraction is identified explicitly. No output correction is
used for acceptance. Virial is compared to affine OpenMM energy derivatives
at two strain sizes, with raw configurational pressure proxy; no inherited
virial/pressure acceptance bound is invented.

Research minimization/NVT/continuation and exact checkpoint/RNG equality are
**not run**, because source prerequisites remain open. Historical J49 bounded
diagnostic continuation is preserved without promotion. Acceptance stays
nonzero for research/production even when software and diagnostic checks pass.

## Reproduction

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J51_PY=../island-validation/phase4g1-env/bin/python

"$J51_PY" scripts/phase_4j51.py audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --supplement island_ffref/pcff/jp4c08471_si_001.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --exact-kernel /tmp/island-j49-exact/build/lmp \
  --exact-kernel-sha256 3a989d8e36345bd4508d932832a78e0eca3313b134d5c414c280065d615029ea \
  --output /tmp/j51-audit
# Expected exit 1: nine unresolved source conventions; no research MD authority.

# Explicit offline pypdf 6.20.0 reader; ordinary preparation does not need pypdf.
python scripts/audit_phase_4j51_sources.py --directory /tmp/j51-audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --reference-directory island_ffref/pcff --lammps-source ../lammps/lammps \
  --exact-kernel-source /tmp/island-j49-exact/source \
  --output /tmp/j51-independent.json

"$J51_PY" -m pytest -q tests/test_phase_4j51.py
"$J51_PY" -m pytest -q tests/test_phase_4j49.py tests/test_phase_4j50.py
"$J51_PY" -m pytest -q
ruff check src tests scripts
"$J51_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$J51_PY" scripts/validate_phase_4j51_receipt.py docs/evidence/phase_4j51 --repo .
"$J51_PY" scripts/validate_phase_4j50_receipt.py docs/evidence/phase_4j50 --repo . \
  --sha256 f61eda5b5b423b2501d29d665a196da9be6cc1374ab4678458b0160e8d31d7b5
```

The J51 receipt can be deterministically regenerated from the completed
`audit/`, `independent-reconstruction.json`, `collection.json` and
`validation/` artifacts, using
`python scripts/seal_phase_4j51_receipt.py docs/evidence/phase_4j51 --repo .`.
The sealer rejects failed/incomplete tests and changed historical bytes.
The full suite is accounted once across eight disjoint file shards plus the
new focused file; exact shard commands are retained in validation evidence.
The separate J49/J50 regression run is an additional check, excluded from
that union to avoid duplicate counts.

Verification: **124 focused tests passed**, **113 J49/J50 regressions passed**,
and **4,231 full-suite tests passed with 61 existing skips** (4,292 collected,
no failures/errors). Ruff, compileall and all five pip checks passed. Independent
reconstruction checked 792 convention citations, 352 backend citations, all 288
raw FD samples, eight 1–4/pressure partitions and both exact-kernel source edits.
All 4,652 historical files retain their hashes. Receipt/manifest validation and
relocated separate-process validation pass; research acceptance remains 1.

The exact diagnostic executable can be rebuilt with the preserved J49 build
script/receipt; compiler/path differences may change its bytes. No engine is
discovered or substituted. Raw validation logs, exact full-suite node accounting,
all-source reproducibility manifest and historical hashes accompany the receipt.

`production_validated=False`, `simulation_readiness="not_established"`.
Polymer Ewald/PME/PPPM closure, backend pressure validation, cutoff recovery,
sampling, equilibrium and long NPT remain unestablished. No universal OPLS/PDMS,
native PCFF PDMS, production readiness or experimental-agreement claim is added.
