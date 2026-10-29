# Phase 4J48: canonical OPLS selection

OPLS selection now has a canonical public entry point: users specify
`forcefield="opls"`; the resolver selects the internal profile and its charge
method from the validated final graph. The bounded methyl-PDMS cases select
`island_external_opls_pdms_xiang2025_dft_local_v1` and `author_dft_npa`.
These remain **external, literature-derived, OPLS-compatible research records**.
They are not native universal OPLS-AA, PCFF or COMPASS coverage.

The software selection gate passes. **LAMMPS end-to-end acceptance remains 1**:
automatic selection reproduces the J47 static electrostatic energy, force and
finite-difference failures. No new LAMMPS minimization, NVT, NPT, trajectory,
checkpoint or MD continuation is performed. The historical OpenMM minimization,
342 DP4 finite differences and 1000-step research NVT records are preserved; they
are not substituted for LAMMPS acceptance.

The branch starts at fetched `origin/main`
`fdf117aadf613912786fdfbebe653ad58695bea3`, the merge of J47 (#71).
It is `codex/phase-4j48-canonical-opls-registry`; main is not merged into it.
All implementation is additive. The manifest checks every tracked base file,
including J17–J47 receipts, serialized identities, sources and bundles.

`production_validated=False`; `simulation_readiness="not_established"`.

## Public API and explicit installation

The canonical public facade is `island.polymer`. The older expert facades remain
available with their historical explicit-selection semantics and unchanged bytes.
Local source installation is separate from selecting a force field. The resolver
never scans directories or downloads files. Configure authorized local assets
once, with an independently supplied SHA-256, then use only the public family:

```python
from island import SimulationBox
from island.polymer import OPLSAssets, prepare_polymer, use_opls_assets

assets = OPLSAssets(
    archive_path="island_ffref/pcff/PDMS_benchmark-main.zip",
    archive_sha256="8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4",
)
with use_opls_assets(assets):
    prepared = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="opls",
        dp=4,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        box=SimulationBox(40, 40, 40),
        # coordinates=coordinates,  # optional exact atom-order angstrom array
    )
print(prepared.payload["selected_profile"])
print(prepared.payload["charge_origin"])
# prepared.system is the explicit-H final system;
# prepared.assignment is the immutable parameter-assignment record.
```

Degree of polymerization, terminal groups, box and coordinates define the
actual system. They do not select charge or parameter models. An already built
linear, multichain or crosslinked system uses the same graph-first boundary:

```python
from island.polymer import (
    ForceFieldRequest, capture_final_graph, prepare_forcefield,
    prepare_polymer_forcefield,
)

with use_opls_assets(assets):
    prepared = prepare_polymer_forcefield(system, forcefield="opls")
    request = ForceFieldRequest(
        family="opls",
        final_graph=capture_final_graph(system),
        typed_graph=prepared.payload["typed_graph"],
        graph_charges=prepared.payload["graph_charges"],
    )
    checked = prepare_forcefield(system, request)
```

Graph/type/charge records constrain the source result. They cannot replace its
charges or parameters. Records returned by the canonical API can be used directly.
Historical graph objects are accepted only if their own integrity/element contract
validates the system; no historical schema is expanded to bypass Si restrictions.
Unsupported multichain end mappings or crosslinked PDMS graphs fail closed.

`Preparation.validate_integrity(system)` re-resolves the graph, rehashes source
bytes and reconstructs the assignment. A re-signed changed record cannot confer
execution authority. Profile choice and the chemical graph identity exclude names, coordinates and
box; the full assignment also carries execution bindings. Chemical identity includes stable atom order, bonds, components, molecule IDs
and the chemical builder end mapping. The separate execution identity binds the
actual coordinates and box, so stale execution input is rejected.

## Registry and charge policy

`island_opls_registry_v1` is distributed as immutable code-backed data and an
exact JSON export in `docs/evidence/phase_4j48/canonical-registry.json`. This
avoids a dependence on source-tree JSON installation or filesystem order.
Every entry records source hashes, kind, environments, graph/end restrictions,
charge method/origin, bonded forms, nonbonded model, mixing, 1–4/exclusions,
long-range conventions, backend restrictions, verification state and priority.

| Public family | Internal profile | Source kind | Default charge method |
|---|---|---|---|
| opls | island_oplsaa_foyer_v1 | native, bounded pinned Foyer source | native source atom charges |
| opls | island_external_opls_pdms_xiang2025_dft_local_v1 | external_literature | author_dft_npa |

Data-only J43–J45 candidates and arbitrary user manifests are not authorized
execution profiles and are not automatically promoted. No additional complete
literature extension has qualified. Explicit user-provided adapters remain on
the historical expert path until independently reviewed for registry inclusion.

Resolution first checks exact chemistry/domain, then checks the selected
profile's complete source typing, charges, rows and conventions. Specificity and
priority determine compatible candidate order. Equal-priority conflicting
parameter, charge, convention or source declarations produce `profile_conflict`.
Unapproved registry changes fail before backend construction. Reordering the
canonical entries does not change the result. Failure of an exact matched domain
never triggers broader incomplete coverage or another force-field family.

The PDMS recognizer requires neutral, explicit-H, nonisotopic, linear methyl
siloxanes with two trimethylsilyl ends. HMDS has a unique oxygen. Longer chains
use the validated builder's graph-bound head oxygen and chemical end-group
provenance for the existing asymmetric E/M assignment. Arbitrary symmetric chains
without that mapping return `typing_ambiguous`; no end is selected by atom name,
coordinates or an arbitrary minimum atom ID. Si–H, Si–OH, Si–O–C, cycles,
networks and P extensions remain unsupported by this profile.

The exact charge vector preserves atom order, component membership/totals,
source identity, profile, units, origin, generation method and evidence. Exact
source literals have tiny residual totals for some fixtures; these are retained
as rational sums. There is no charge normalization or redistribution. Missing
source charges or mismatched supplied charge records fail closed. There is no
zero filling, interpolation, cross-family borrowing or new physical parameter.

## Source identity and numerical outcome

The external source is Xiang et al., DOI `10.1021/acs.jpcb.4c08471`, JPCB 2025,
129, 1864–1873. The active branch is the author's serialized **DFT/NPA** model;
the main paper Results use Qeq. Source access authorizes this local research
exercise; redistribution permission for the source archive is not established.
No paper or source ZIP is committed.

| Asset | SHA-256 |
|---|---|
| PDMS_benchmark-main.zip | 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 |
| jp4c08471.pdf | 04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc |
| pinned LAMMPS | 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 |
| native Foyer oplsaa.xml | c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4 |

The registry also preserves all five author companion file hashes. Runtime reads
and hashes the explicitly installed ZIP/companions; the PDF is an independently
pinned audit reference, not a runtime discovery requirement. Seventeen profile
coefficient rows are independently reconstructed out of 3108 source coefficient
candidates. Original source zeros remain independently evidenced original rows.

`island.opls_lammps.run_static` validates the canonical assignment, then uses the
unchanged J47 deterministic exporter and pinned engine. It seals the canonical
assignment together with original source preparation, system, data/input,
engine/build record, raw output, forces, energies and virial diagnostic.
Relocation requires explicit source/executable paths and independent hashes.
Fresh processes re-resolve chemistry, regenerate original inputs, verify file
membership/checksums and reparse raw outputs. Two process hash seeds are checked
for every fixture. Split/uninterrupted MD and exact checkpoint/RNG comparison
remain explicitly blocked behind failed static gates; they are not claimed.

| Fixture | Atoms | Total energy error kJ/mol | Max force error kJ/(mol Å) | Max LAMMPS FD error kJ/(mol Å) |
|---|---:|---:|---:|---:|
| HMDS | 27 | 0.0699658923 | 0.0839530434 | 0.2647347082 |
| capped DP2 | 37 | 0.1436606402 | 0.2561358095 | 0.4239769054 |
| capped DP4 | 57 | 0.4018122426 | 0.4421014032 | 0.3366947118 |
| capped DP6 | 77 | 0.5085865935 | 0.3382671888 | 0.5528278625 |

Energy and force bounds remain `3e-4`; finite-difference bound remains `1e-4`.
Each fixture has 24 OpenMM and 24 LAMMPS FD comparisons. The OpenMM FD checks
pass; the source-settings LAMMPS FD checks fail. Bond, angle, proper torsion and
LJ energy components pass. Improper absence is source-declared for this domain.
Special-pair identities and separate 1–4 reconstructions are retained.

The reproduced discrepancy classification is **long-range solver**: source
LAMMPS uses PPPM `1e-4` and its default pair table; OpenMM uses Ewald `1e-7`.
Identical numerical splitting/mesh treatment is not specified by the profile.
Real/reciprocal partitions therefore differ and their Coulomb sum is compared.
J47's separate exact-equation control also identifies **precision** (LAMMPS erfc
polynomial) and **unit conversion** (Coulomb constant 1389.35456264 versus
OpenMM 1389.35457644382 kJ mol⁻¹ Å e⁻²). Its raw failures and source hashes are
preserved. Explained or corrected residuals are not used as acceptance values.

Closing this requires an explicitly reviewed, matched numerical solver contract
and an engine implementation that passes the unchanged force/FD gates. J48 does
not assume authorization to replace the source's reciprocal/tables convention.
Virial has no inherited acceptance threshold and stays unresolved; pressure is
unavailable. No NPT, density, equilibrium, Ewald/PPPM closure, universal coverage,
production PDMS LAMMPS MD or scientific readiness claim is made.

## CLI and reproduction

An explicit installation can also be a local JSON configuration with keys
`archive_path`, `archive_sha256`, and optionally `native_source_path` and
`native_source_sha256`. Point `ISLAND_OPLS_ASSETS` at that file and supply its
independent SHA-256 in `ISLAND_OPLS_ASSETS_SHA256`. No default asset search occurs.
Ordinary inspection does not need profile/charge-method flags:

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J48_PY=../island-validation/phase4g1-env/bin/python

# After explicitly configuring the local asset JSON and independent hash:
"$J48_PY" -m island.inspect_forcefield \
  --psmiles '[*:1]O[Si](C)(C)[*:2]' --dp 4 \
  --head-end-group '[*][Si](C)(C)C' --tail-end-group '[*]C' \
  --forcefield opls --json /tmp/j48-inspection.json
"$J48_PY" -m island.inspect_forcefield \
  --forcefield opls --registry --json /tmp/j48-registry.json
```

The inspector uses the declared 40 Å cubic audit box for a siloxane export
inspection, with deterministic builder coordinates; it never runs MD or infers
chemistry from those coordinates. Inspection exit 0 means parameter assignment
completeness, not cross-engine or scientific acceptance. Expert-only explicit
`--profile`, `--charge-method`, `--source` invocations delegate to the unchanged
J46 module. The additive package entry point preserves that historical file's
bytes and regression behavior.

```bash
"$J48_PY" scripts/phase_4j48.py \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --output /tmp/j48-reproduction
# Expected exit 1, with a complete classified static matrix and no MD.

PYTHONPATH=. python scripts/audit_phase_4j48_inputs.py \
  --directory /tmp/j48-reproduction \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --output /tmp/j48-independent-inputs.json

"$J48_PY" -m pytest -q tests/test_phase_4j48.py
"$J48_PY" -m pytest -q tests/test_phase_4j47.py
"$J48_PY" -m pytest -q
ruff check src tests scripts
"$J48_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$J48_PY" scripts/validate_phase_4j48_receipt.py docs/evidence/phase_4j48 --repo .
"$J48_PY" scripts/validate_phase_4j47_receipt.py docs/evidence/phase_4j47 --repo . \
  --sha256 b74d0caca18859cc5072989418ecf3478ff729d27c206e4a5229ced7eefa08c7
```

The new validator verifies the exact collected test-node inventory across full
suite shards plus focused tests, rather than inferring coverage from a total.
Receipt, manifest and relocated fresh-process reconstruction checks are software
integrity checks. Numerical acceptance remains a separate, nonzero result.
