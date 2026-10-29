# Phase 4J46: explicit unified force-field selection

J46 adds a common preparation and inspection boundary for PCFF, OPLS,
GAFF and GAFF2. It preserves the selected family, final atom order,
components, chemical construction metadata, source hashes, typing, charges,
parameter assignments and backend identity. Missing interactions remain visible
and prevent execution. No new physical parameter, native-source override,
Class-II zero coupling or production PDMS simulation is introduced.

`production_validated=False` and `simulation_readiness="not_established"`.
Polymer Ewald/PME scientific closure remains unvalidated; backend pressure is
unavailable; cutoff recovery and long NPT remain blocked. No equilibrium,
experimental agreement or universal source/backend coverage is claimed.

## Exact base and preserved history

After `git fetch origin`, `origin/main` was
`fdf8e8058173faa9c013703aa4910b3065820e7c`. Main has not merged the later
J43–J45 and PDMS work. The requested branch starts at that exact main.
Prerequisite commit `127f4d1` carries the 319 added prerequisite files unchanged
from published snapshot `8257a2ef25042494d0e4541b8da117cd8e97072a`.
The J46 implementation starts in `61856c7`; the tested implementation is
`2fda1cc`.
There is no merge or change to main.

The J43, J44 and J45 reports, receipts, source catalogs and gap reports were
read before implementation. Their findings are preserved: J43 had no supplied
assets; J44 found 52 local files and 136 numerical candidates without a complete
profile; J45 retained missing conventions for the Xiang Table S4 candidate.
The later source-bound author DFT/NPA research adapter has separate provenance
and does not rewrite those outcomes. J46 reuses these contracts and does not
repeat the literature inventory or distribute the user papers/author archive.

The reproducibility manifest protects all **3,567 tracked entries** in the
published prerequisite snapshot, including every J17–J45 historical receipt,
bundle, source and identity. Symlink identities are checked against their Git
link-target bytes. The older PDMS validator separately preserves its 3,476-entry
historical inventory, including local assets.

## Public interface

```python
from island.builders import build_linear_polymer
from island.forcefields.preparation import OPLSOptions
from island.forcefields.selection import (
    SourceSelection, prepare_polymer_forcefield, create_evaluator,
)

system = build_linear_polymer(
    "[*:1]CC[*:2]", dp=2,
    head_end_group="[*]C", tail_end_group="[*]C",
    add_hydrogens=True, random_seed=46,
)
result = prepare_polymer_forcefield(
    system, forcefield="opls",
    source_profile=SourceSelection(
        "island_oplsaa_foyer_v1", OPLSOptions("/explicit/path/oplsaa.xml"),
    ),
    charge_method="native", resolution_policy="strict_source",
)
report = result.payload
# Inspect report["stages"], report["requirements"] and report["diagnostics"].
# This raises UnifiedForceFieldError with the complete result when blocked.
backend = result.require_executable(system)
# Construct an evaluator only after that same final graph passes integrity checks:
evaluator = create_evaluator(system, result)
```

Existing `ForceFieldRequest`, `prepare_forcefield`, scientific gates and numerical
implementations remain unchanged. The new package delegates to those adapters
and to the existing opt-in research/external adapters. A returned native backend
can be passed to its existing workflow. A PDMS backend retains the bounded
research workflow described in `docs/pdms_user_workflow.md`; J46 does not authorize
additional MD modes or LAMMPS execution.

There is no implicit source discovery. Options must name the existing explicit
source asset. Native sources retain their established source pins; external
manifests and the PDMS archive require independent supplied hashes. Every
external companion XML/license asset is also checked. Serialized JSON reports
are inspection evidence, never reconstructed execution authority. After
relocation, repeat preparation with explicit assets and compare identities.

The same API accepts multichain and crosslinked **final graphs**. It derives
interactions from final bonds, rather than carrying repeat-unit atom types
through topology edits. Changed atom IDs/order, component metadata, masses,
chemical graph, source bytes or resolver implementation invalidate old results.
Coordinates and documented coordinate-only metadata may change independently.

## Capability and assignment reports

`capability_registry()` returns owned, deterministic declarations for each
selection profile and the historical PCFF registrations. Each entry includes
elements/environments, bonded and nonbonded forms, all Class-II couplings,
charges, mixing, exclusions, 1–4 conventions, electrostatics and backend limits.
An empty element list means source-dependent coverage, never all-element support.
Historical operational profiles requiring legacy fallback policy are explicitly
marked `selection_adapter=False`; the new strict boundary does not enable them.
Their existing APIs and finite-silane evidence remain unchanged. Older broad
PCFF graph rules cannot bypass the J41 phosphorus oxidation/connectivity
boundary: P-containing inputs remain blocked with the audited chemical-state
diagnostics, including for graph profiles v1–v6.

Reports independently expose `typing_complete`, `charge_complete`,
`parameter_complete`, `executable`, `numerically_verified` and
`scientifically_validated`. Registry declarations alone do not establish any
stage for an arbitrary graph. In particular, `prepared` establishes an existing
backend software capability, not numerical or scientific qualification of a new
domain. The last two flags remain false in these inspection reports.

Every required interaction has a stable identifier, site IDs, status and
source/equilibrium dependencies. Successful assignment retains the underlying
source-row/equivalence evidence. Unresolved requirements and diagnostics remain
available together in deterministic order. Atom types, exact source charges,
masses, policies and graph-distance special pairs are inventoried separately.
Native PCFF records retain ordinary-equivalence/source selection and all
bond/angle/torsion/improper couplings. Explicit source zeros remain distinct from
missing rows. OPLS records retain RB versus OPLS torsion forms, geometric mixing,
0.5/0.5 scaling, exclusions and the selected electrostatics policy.

Blocked statuses distinguish `asset_unavailable`, `typing_ambiguous`,
`charge_incomplete`, `parameter_incomplete` and `source_incompatible`.
`create_evaluator` raises `UnifiedForceFieldError` with `.status`,
`.diagnostics` and the complete `.result` before constructing a numerical engine.
Backend charge exceptions are classified through their original exception chain.

GAFF/GAFF2 keep AmberTools AM1-BCC and provided-charge paths. Provided vectors
require an explicit origin. `resp_psi4` additionally requires `RESPSource` with
an independently reviewed identity and every existing Psi4/ESP/RESP asset;
the exact verified vector must match the options. A provided vector named RESP
cannot bypass that loader. External charges remain explicitly external.
AmberTools binaries still require the caller's explicit PATH/environment; the
selector does not change tool discovery or charge tolerances.

## Missing-term policies

| Policy | Accepted evidence | Execution boundary |
| --- | --- | --- |
| `strict_source` | Complete existing source rows | Every missing required row blocks; no completion evidence accepted |
| `derived_from_declared_formula` | `derive_geometric_pair`: explicitly declared geometric sigma/epsilon mixing, two exact cited self-pair rows, units, hashes, equation and independent Decimal/binary64 reconstruction | Data-only `derived` record; cannot replace a native coupling or create an execution profile |
| `fit_from_external_reference` | `fit_harmonic_bond`: explicit hash-bound external isolated harmonic-bond scan, labels, provenance, objective/weights/constraints, residual covariance and independent least-squares/normal-equation reconstruction | Data-only `fitted_external` candidate; not an original source row or MD authority |

There is no defensible general analytical formula for the missing PCFF
Class-II coupling rows. J46 does not derive them. The fitting adapter rejects
whole-molecule total energies presented as an isolated bond contribution,
ambiguous labels, unsupported forms, unresolved units, invalid weights and
rank-deficient data. It reports residual-based statistical uncertainty explicitly;
this does not measure physical model uncertainty. Analytical test datasets in
the tests are software fixtures, never supplied PDMS physical parameters.

Completion records must match the explicit policy and selected profile, survive
source/raw-row/unit/evidence reconstruction, and remain separate from native
assignments. A complete existing native assignment may carry a reconstructed
candidate record for review; the candidate does not change that assignment.
Attaching completion evidence cannot make an incomplete native model executable.
Any eventual physical fit still needs a complete separately authorized external
bundle and independent numerical qualification under the existing contracts.

## Bounded source inspection

The reproduced audit contains 42 cases. Its source assets are:

| Asset | SHA-256 |
| --- | --- |
| pinned `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| pinned Foyer `oplsaa.xml` | `c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4` |
| user `PDMS_benchmark-main.zip` | `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4` |

Examples of the deliberately different outcomes:

- Native PCFF methane is prepared; the strict aliphatic/multichain/PEO cases
  retain missing coupling rows. Ethane has nine missing `bond-bond_1_3` paths
  with exact source keys and bond-equilibrium dependencies.
- Native OPLS alkane, benzene, PEO, two-component and bounded bridged-network
  controls are prepared. The methyl-capped aromatic PSMILES chain is
  `charge_incomplete` (native component residual `-0.23000000000000004 e`),
  without repair or another family's charges.
- Native OPLS Si/P cases remain `typing_ambiguous`. The existing separate
  Xiang author DFT/NPA external profile prepares only its supported neutral
  explicit-H trimethylsilyl-capped siloxanes (bounded DP0, DP2 and DP4 here).
  This is not the article's Qeq main model, native OPLS/PCFF completion,
  universal Si support or universal msi2lmp parity.
- GAFF/GAFF2 cases without explicit options are `asset_unavailable`; source
  discovery and speculative Si/P support are not inferred from family names.

No production PDMS MD or new profile qualification is performed by this audit.
The focused tests compare unified and original finite PCFF/OPLS evaluators
exactly. Existing historical workflow tests may exercise their previously bounded
research simulations; they do not change any scientific acceptance flag.

## Reproduction

Use the existing validation environment; the paths below are explicit local
assets, not runtime defaults:

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=../island-validation/phase4g1-env/bin/python
PCFF=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/tools/msi2lmp/frc_files/pcff.frc
OPLS=../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml

$PY -m island.inspect_forcefield --registry --json /tmp/j46-registry.json
$PY -m island.inspect_forcefield --smiles CC --forcefield opls \
  --charge-method native --source "$OPLS" --json /tmp/j46-ethane.json
$PY -m island.inspect_forcefield --smiles CC --forcefield pcff \
  --charge-method native --source "$PCFF" --lj 0 0 1 --coulomb 0 0 1 \
  --json /tmp/j46-pcff-blocked.json
# The preceding strict PCFF command must exit 1: missing bond-bond_1_3 rows.

$PY -m island.inspect_forcefield --psmiles '[*:1]CC[*:2]' --dp 2 \
  --head-end-group '[*]C' --tail-end-group '[*]C' --forcefield opls \
  --charge-method native --source "$OPLS" --json /tmp/j46-polyethylene.json

$PY scripts/phase_4j46.py --output /tmp/j46-audit \
  --pcff-source "$PCFF" --opls-source "$OPLS" \
  --pdms-archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --pdms-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4

$PY -m pytest tests/test_phase_4j46.py -q
$PY -m pytest -q
ruff check src tests scripts
$PY -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
$PY scripts/validate_phase_4j46_receipt.py docs/evidence/phase_4j46 --repo .
```

The audit command exits zero when it successfully records both prepared and
blocked inspection cases. Individual inspection CLI commands exit zero only
for an executable assignment and one for a classified blocked result. Neither
is a full-source PCFF acceptance or scientific-readiness command. Existing
native full-source acceptance remains blocked.

`receipt.json` and `manifest.json` record the exact tested implementation,
source identities, commands, test counts and historical byte hashes.
`summary.json` retains every case's blocked keys; `audit.tar.xz` contains the
complete assignment/source-row evidence without the source publications or ZIP.
The complete ordinary suite is partitioned by test file into eight independent
processes, plus the focused J46 file. The collection and JUnit files account
for every collected node exactly once; the validator rejects omissions,
duplicates, failures, changed artifacts and promoted readiness.

Historical J42's validator pins its pre-C1 implementation. It is reproduced
against its original tested detached snapshot `2302d5275c924aefd840756b1ac4256fce95f038`;
J42-C1 and later validators run against the preserved current tree. This avoids
rewriting old receipts to match their later intentional corrections.

## Verification results

- Focused J46: **80 passed**.
- Complete collected suite: **3,867 passed, 61 skipped, zero failures/errors**
  across the historical file shards and focused J46 file (**3,928 nodes**).
- Ruff, compileall, all five pip checks, diff checks, J46 receipt/manifest and
  historical hash validation pass.
- J42 in its original tested snapshot, J42-C1, J43, J44, J45 and the existing
  externally pinned PDMS research receipt validate. The deliberately nonzero
  scientific and native full-source acceptance gates remain unchanged.
- Fresh-process CLI relocation reproduces exactly the same native OPLS report;
  changed sources, final graphs, atom order, components, charges, evidence and
  profile identities are rejected. Unified/native finite evaluator results agree
  exactly in the bounded PCFF and OPLS controls.
