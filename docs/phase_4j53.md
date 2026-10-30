# Phase 4J53: authoritative OPLS-PDMS source intake

The supplied directory contains **62 files whose bytes already have hashes in
committed evidence**. There is no new authoritative source package. The source
status is `not_supplied`. This phase implements intake, precise diagnostics and
permanent-blocker management; it creates no profile and runs no new LAMMPS
static calculation, minimization, NVT or NPT. Existing component comparisons
are referenced by their original hashes, without repeating the exhausted audit.

Fetched exact base: `c5a06496a9c27560ef627baf21bb95aeaef40080`, merged J52
(#76). Branch: `codex/phase-4j53-authoritative-opls-pdms-source-intake`.
All tracked base files, including J17–J52 receipts, bundles, identities and
implementation code, remain byte-identical. The implementation is additive.

## Nine permanent blockers

These are the exact J52 settings, in their original order. Permanence applies
to that supplied corpus. A new source can support a new reviewed version;
it cannot rewrite the old conclusion.

| Setting | Artifact required to reopen |
| --- | --- |
| `dielectric` | DFT/NPA-scoped declaration or author executable/restart proving the value |
| `coulomb_prefactor` | Original engine revision/build and real-unit electrostatic constant |
| `periodic_image_policy` | DFT/NPA surface/slab/image-sum policy; xyz PBC alone is insufficient |
| `precision` | Author binary hash and precision/compiler/accelerator manifest |
| `neighbor_list_policy` | Complete DFT/NPA driver binding skin, rebuild and sorting policy |
| `pressure_convention` | Actual DFT/NPA estimator, including kinetic and constraint contributions |
| `virial_convention` | Tensor, images, constraints, kinetic terms and conversion definition |
| `pair_table_bits` | Author `pair_modify table` declaration or exact executable |
| `erfc_kernel` | Author binary/revision establishing the energy/force kernel |

`audit/blocker-reopen.json` contains every original permanent-blocker record
verbatim, its exact missing artifact, affected modes/observables and the new
intake identity. No blocker is reopened. J51's converter-resolved CHARMM LJ
no-tail restriction remains unchanged. Backend-resolved diagnostic settings
remain separate from source authority.

## Source package contract

`python -m island.intake_opls_source schema` emits the versioned JSON Schema.
A package is an explicitly installed directory containing `source-intake.json`
and an independently supplied SHA-256 pin for those manifest bytes. The
manifest inventories **every other physical file**, including hidden files,
empty files, source text extractions, reconstruction code and permission records.
The manifest itself is pinned externally, avoiding a self-hash cycle.

Every entry records relative path, byte count, SHA-256 and role: paper,
supplement, author parameters, input/data, converter, source version,
permission, extracted evidence, independent reconstruction or exact system.
ZIP/TAR members are additionally inventoried and hashed, including nested
archives. Links, unsafe paths, duplicate members/records/JSON keys, corrupt
archives and undeclared files are rejected. Archives are inspected without
extraction or execution, with explicit size/nesting limits. Source code from
the package is never executed by intake.

Required metadata includes DOI/citation, source version and its evidence,
external OPLS identity, local-research authorization, rights holder and separate
redistribution status. Local use does not establish redistribution permission.
The article license does not automatically authorize a separate parameter ZIP.
Authorization must explicitly cover the paths and hashes of every non-permission
file in the package. Missing archive coverage remains a conflict even after
the manifest is re-pinned.

All 37 existing execution settings are required, together with explicit
parameter functional form, charge method and exclusions. This covers units,
dielectric/prefactor, LJ form and sigma/epsilon conversion, torsions, special
pairs, cutoffs/switching, long-range solver/tolerance, periodic images,
precision, neighbors, virial and pressure. Only `source_resolved` or
`converter_resolved` declarations may enter the candidate contract; a
converter declaration requires inventoried converter evidence. Backend defaults
and family conventions cannot clear a missing source setting.

Each convention/interaction record retains its family, labels, functional form,
original/converted value, units, conversion formula and implementation/version,
source path/hash, DOI/citation, page/table/equation/line locator, exact raw UTF-8
lines and evidence checksum. Binary PDFs need an inventoried text extraction.
Reconstruction code/result locations, version, value and units are required;
the machine-readable result must agree exactly. This checks provided records,
not the authenticity of an unexecuted reconstruction. Qualification must
independently reproduce the extraction and conversion before selection.

For each declared exact graph, the contract includes a checksummed system,
typing identity, atom types in exact graph order, the complete charge vector,
origin/method/units and source-bound rows for every required interaction.
Interaction IDs and dependencies are generated independently from the final
graph. Missing/conflicting rows, type/charge/mass inconsistency and unresolved
units are rejected. Charges are neither normalized nor redistributed.

## Status and selection policy

| Computed status | Meaning |
| --- | --- |
| `not_supplied` | No explicit new package installed |
| `supplied_unreadable` | Directory, manifest, extraction or archive cannot be read |
| `supplied_incomplete` | Required metadata, permission, conventions, rows or reconstruction are absent |
| `supplied_conflicting` | File/hash/evidence, graph, charge, type or duplicate record conflict |
| `source_complete_not_executed` | Structurally complete evidence contract; no execution or registration |
| `source_complete_executed` | Reserved for a pinned source-specific adapter's actual qualification runs |
| `source_complete_independently_verified` | Also requires independent numerical and workflow reconstruction |

J53 has **no new qualified adapter or registry entry**. Intake can reach
`source_complete_not_executed`, but cannot reach either executed status.
Package-supplied status/qualification flags are forbidden. A complete synthetic
contract in tests exercises the validator; it supplies no scientific PDMS data.

Automatic selection still requires authorized registry membership, exact graph
domain and typing, complete source-bound charges/interactions, source/converter
conventions and independent source-specific qualification. A candidate package
cannot add itself to the registry, replace the existing profile or trigger a
fallback. No new source is automatically selected in this phase.

The unchanged internal prerequisite is
`island_external_opls_pdms_xiang2025_dft_local_v1`, external literature-derived
OPLS-compatible, bounded to TMS-capped methyl PDMS. The DFT/NPA branch is not
the paper's Qeq Results model. The original SI/author-input pair-style scope
conflict is preserved. The profile is not native universal OPLS-AA, PCFF or
COMPASS. No charge, parameter or acceptance tolerance is changed.

## Public API and modes

Configure the existing local assets once. Ordinary callers still select only
`forcefield="opls"`:

```python
from island.polymer import OPLSAssets, use_opls_assets
from island.opls_final_resolution import ConventionSources, use_convention_sources
from island.opls_source_intake import prepare_polymer_forcefield

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

No profile or charge-method argument is required. The immutable report retains
the selected historical profile, source/typing/charge/domain identities, four
convention layers, J52 resolution identity, source intake, blocker-reopen
report, execution mode and scientific readiness. Paths do not enter identities.
Installing a new package uses `SourcePackage` and `use_source_package` at the
administrator boundary. It does not grant execution rights.

`island.opls_source_execution` provides guarded diagnostic export/static
wrappers and a report-bundle save/reconstruct path. With no new package, an
explicit historical diagnostic remains subject to J52's unchanged policy.
An installed unqualified candidate cannot drive even static execution.
`research_opt_in` minimization/NVT and `production_candidate` remain blocked
before engine access or output creation. NPT is blocked. Saved preparations
are reports, not executable checkpoints.

HMDS and capped DP2/DP4/DP6 each have a family-only preparation, unchanged
physical fixture bytes, three pre-engine blocked-action checks and relocated
reconstruction in separate processes at hash seeds 1 and 99397. These tests
qualify the intake workflow only. No new energy/force/FD, MD, checkpoint or
RNG-continuation result is claimed. J49–J52 component/FD/virial records remain
referenced by their original hashes in `audit/prior-comparison-records.json`.

## Reproduction

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J53_PY=../island-validation/phase4g1-env/bin/python

"$J53_PY" -m island.intake_opls_source schema --json /tmp/intake-schema.json
"$J53_PY" -m island.intake_opls_source inspect \
  --json /tmp/intake.json --blocker-reopen-json /tmp/blockers.json
# Expected exit 1: not_supplied; no source discovery or engine execution.

# For a genuinely new directory, prepare the manifest using the emitted schema:
"$J53_PY" -m island.intake_opls_source inventory --package /path/to/new-package
"$J53_PY" -m island.intake_opls_source inspect --package /path/to/new-package \
  --manifest-sha256 <independently-reviewed-manifest-sha256> --json /tmp/new-intake.json
# Exit 0 means structurally complete, not authorized or scientifically verified.

"$J53_PY" -m scripts.phase_4j53 audit \
  --sources island_ffref --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --supplement island_ffref/pcff/jp4c08471_si_001.pdf --output /tmp/j53/audit
# Expected exit 1: no new source; four blocked reports and eight fresh processes.

"$J53_PY" -m pytest -q tests/test_phase_4j53.py
"$J53_PY" -m pytest -q tests/test_phase_4j49.py tests/test_phase_4j50.py \
  tests/test_phase_4j51.py tests/test_phase_4j52.py
"$J53_PY" -m pytest -q
ruff check src tests scripts
"$J53_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
git diff --cached --check
"$J53_PY" scripts/validate_phase_4j53_receipt.py docs/evidence/phase_4j53 --repo .
"$J53_PY" scripts/check_phase_4j53_receipt.py docs/evidence/phase_4j53 --repo . \
  --output /tmp/j53-receipt-controls.json
```

All source hashes, profile identities, historical bytes and implementation
hashes are in `manifest.json`. Exact test commands and disjoint full-suite
membership are retained in validation artifacts. The seal command is
`python scripts/seal_phase_4j53_receipt.py docs/evidence/phase_4j53 --repo .`.
Historical J49–J52 validators are run with the same Python environment.

Validation completed: **95 focused tests passed**, **320 J49–J52 regressions
passed**, and **4,409 ordinary-suite tests passed with 61 existing skips**
(4,470 distinct nodes, zero failures/errors). The ordinary-suite union contains
all 4,375 baseline nodes and all 95 J53 nodes exactly once; the separate
regression run is excluded from that count. Ruff, compileall, five pip checks,
diff checks, J49–J52 receipt validators and all **5,165** historical file hashes
passed. Eight relocated fresh processes reconstructed the four preparations.
Twelve re-signed receipt/evidence mutations were rejected. The retained control
record identifies the preliminary sealed receipt it tested; the final receipt
is then re-sealed to include those results and validated again after relocation.

Readiness remains `production_validated=False` and
`simulation_readiness="not_established"`. There is no universal OPLS-AA/PDMS
coverage, native PCFF PDMS, Ewald/PPPM closure, pressure validation, production
equilibrium or long-NPT readiness claim.
