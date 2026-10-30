# Phase 4J54: separate external PCFF/COMPASS Class-II PDMS audit

J54 starts at exact fetched main `b296a06488abca891ae46511aba3a55b075beebf`
(merged J55); the implementation branch is
`codex/phase-4j54-external-pcff-pdms-class2`. The receipt records the full base
identity. Historical J17–J55 files, scientific claims, serialized identities and
native force-field implementations remain byte-for-byte unchanged.

**No complete external Class-II PDMS domain qualifies from the supplied assets.**
The proposed identity `island_external_pcff_pdms_class2_v1` is a candidate label,
not a registered or executable profile. J54 adds a source-bound audit and exact
missing-row diagnostics. It does not provide PCFF PDMS MD. The existing external
OPLS PDMS profile is neither used nor modified.

## Separate source candidates

The earlier source inventories, Sun numerical catalog, PDMS FRC completion study
and `pcff/main.pdf` review are reused. Only known relevant literature hashes are
checked; there is no blind repeat of the literature inventory. The two original
Sun PDFs are independently reread: 105 numerical records reconstructed, with all
8 printed angle blanks retained. Their file hashes, citations, page/row locations,
OCR corrections and evidence checksums remain linked to the historical catalog.
No paper/FRC coefficient overlay is performed.

1. **One pinned PCFF distribution.** `pcff.frc` SHA-256
   `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`,
   plus independently pinned RLB/templates, supplies audit-only typing, masses,
   increments and available Class-II rows. All numerical candidates in this
   assessment come from this single version. No FRC DOI is invented or replaced
   with the DOI of a numerically different paper. Each used row retains its
   internal reference, source revision, original line, labels, decimals, original
   and converted units/values, conversion formula/code hash and evidence checksum.
2. **The author's separate COMPASS input example.** The supplied
   `PDMS_benchmark-main.zip`, SHA-256
   `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4`,
   contains `param/COMPASS` as well as other models. J54 inventories every archive
   member recursively and audits only the COMPASS inputs. No OPLS rows, charges,
   torsions or conventions enter this assessment. Its explicit `units real`,
   Class-II styles, LJ9–6, sixth-power mixing, `special_bonds lj/coul 0 0 1
   dihedral yes`, 10 Å cutoff and PPPM tolerance 0.0001 belong to that example;
   they cannot authorize missing conventions for the other FRC version.

The COMPASS example contains **53 active coefficient commands**, including
literal zero coupling coefficients. These are preserved as author-exported input
cells; they are not zeros inserted by ISLAND. Source atom-type paths and
command-specific dimensions accompany the cells. Their generating converter and
transfer rules are not supplied. Its 32,448-atom topology has improper/AA terms
centered on carbon but none centered on silicon: **12,800 Si-centered neighbor
triples** lack an author-topology improper/AA mapping or explicit absence rule.
This is a graph-site count, not a count of independent fitted constants. An
inactive coefficient with only a numeric type ID does not establish a transferable
Si-centered rule. New HMDS/DP2/DP4/DP6 typing and terminal-charge mapping, exact
converter provenance and archive authorization remain incomplete. The bulk charge
vector is preserved as source evidence; it is not transferred or normalized.
The archive has no license/permission file; the article's license is not promoted
to an archive redistribution grant.

Sun 1995 silane parameters, Sun–Rigby 1997 siloxane parameters, the FRC version and
the author COMPASS example remain separate. Smith 2004's Buckingham/sine-squared
model cannot supply this LJ9–6 Class-II domain. The aromatic polysiloxane
application `pcff/main.pdf` provides no new transferable coefficient rows. The
known Table S3 coupling-equilibrium discrepancies remain linked to the J44
findings, even though some corresponding author-archive entries differ. A common
family name is not authority to merge these versions.

## Exact graph and dependency results

| Fixture | Atoms | FRC typing / increments | Missing FRC coupling requests |
| --- | ---: | --- | ---: |
| Disiloxane | 9 | Complete, separate Si–H diagnostic only | 54 |
| HMDS | 27 | Complete candidate reconstruction | 378 |
| TMS-capped DP2 | 37 | Complete candidate reconstruction | 528 |
| TMS-capped DP4 | 57 | Complete candidate reconstruction | 828 |
| TMS-capped DP6 | 77 | Complete candidate reconstruction | 1,128 |
| `[*:1]O[Si](C)(C)[*:2]` | Open descriptor | No final typing or charges | Explicit DP and caps required |

There are **2,916 missing FRC site requests**, grouped into **47 directed source
queries**, not 2,916 independent coefficients. The graph audit independently
reconstructs all selected source rows, source ordinary equivalences, graph paths
and Decimal charge increments. No wildcard, lower-order fallback, interpolation,
charge normalization or coefficient zero filling is allowed. Source-printed
symmetric rows may repeat their explicit coefficient block under the recorded
converter rule; a blank coefficient is never synthesized.

For HMDS the 49 bond–bond, 49 bond–angle and all 60 proper main terms have source
rows. Missing requests include 78 angle–angle entries and 60 each of bond–bond
1–3, end-bond–torsion, middle-bond–torsion, angle–torsion and
angle–angle–torsion. DP6 extends the exact inventory without expanding chemical
typing rules. Equilibrium dependencies are retained for every coupling and missing
bond/angle dependencies propagate to dependent rows. The empty source-declared
torsion–torsion section and domain-specific improper applicability remain explicit
model gaps. The source's linear-Phi angle–angle–torsion comment is not replaced
with a backend cos(Phi) model. Units or source exclusions/1–4/cutoff/long-range
rules are never inferred from general PCFF defaults.

The ordinary candidate domain is one finite, linear, explicit-H,
trimethylsilyl-capped methyl-siloxane chain. Si–H is not included in ordinary
selection: the exact disiloxane control has its own explicit diagnostic switch
and source evidence. Si–OH, Si–O–C, cyclic siloxanes, crosslinks, phosphorus,
charged/radical/isotopic states and open attachment points reject domain typing.
These methyl-only fixtures contain no C–C bond; C–C coverage is not extrapolated.
There is no periodic execution convention or terminal-group inference.

## Public inspection and independent J55 statuses

J54 is an additive audit facade. The canonical registries and historical J55
implementation stay unchanged because no new profile qualifies for registration.
Administrators install the pinned source directory once; ordinary inspection
requires only the public family, without profile or charge-method arguments:

```python
from island.pcff_pdms import (
    PCFFPDMSAssets, use_pcff_pdms_assets, prepare_polymer_forcefield,
)

with use_pcff_pdms_assets(PCFFPDMSAssets("/local/original/frc_files")):
    result = prepare_polymer_forcefield(system, forcefield="pcff")
report = result.payload
# family="pcff", candidate_profile=..., selected_profile=None
# parameter_completeness="incomplete", typing_status="complete",
# charge_status="complete", cross_engine_equivalence="not_claimed"
```

`ISLAND_PCFF_PDMS_SOURCES` provides the equivalent administrator directory setting.
`prepare_polymer` accepts PSMILES, explicit DP and explicit end groups; a repeat
string alone does not become a final molecule. Results expose source identities,
typed graph, exact charges, per-site source-to-row mapping, all missing requests
and dependencies, required conventions and the independent J55 statuses.

```bash
python -m island.pcff_pdms --smiles 'C[Si](C)(C)O[Si](C)(C)C' \
  --forcefield pcff --distribution /local/original/frc_files --json report.json
# Explicit independent claims return nonzero for this incomplete candidate:
# --require-parameter-complete
# --require-cross-engine-equivalence
```

Backend export/build/runtime are `blocked` by missing source rows/conventions.
Cross-engine equivalence is `not_claimed`; its failure/unavailability is never a
parameter-completeness gate. Scientific validation and production readiness are
`not_established`. Immutable records replay live sources, graph, order, charges,
units, rows, source/code hashes and authority; re-signing mutated JSON is
insufficient. Saved audit bundles contain no prepared MD system or engine input.

## Evidence and reproduction

`docs/evidence/phase_4j54/audit` contains the Class-II dependency matrix,
source-to-row mapping, exact missing-row report, separate author COMPASS audit,
known-literature context and per-fixture independent reconstructions/bundles.
The manifest protects all source hashes and historical tracked files. Five
relocated bundles reconstruct identically in ten fresh processes with distinct
hash seeds. The original PDFs/archive are not redistributed.

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J54_PY=../island-validation/phase4g1-env/bin/python
"$J54_PY" -m scripts.phase_4j54 audit \
  --distribution ../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --references island_ffref --pdf-python ../miniforge3/bin/python \
  --output /tmp/j54-audit
# Add --require-parameter-complete for expected exit 1 after writing diagnostics.
"$J54_PY" -m pytest -q tests/test_phase_4j54.py
"$J54_PY" -m pytest -q tests/test_phase_4j55.py
"$J54_PY" -m pytest -q
ruff check src tests scripts
"$J54_PY" -m compileall -q src tests scripts
# Five pip checks and exact ordinary-suite membership are retained in validation.
"$J54_PY" scripts/validate_phase_4j54_receipt.py docs/evidence/phase_4j54 \
  --repo . --distribution ../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip --references island_ffref
git diff --check
```

No backend static, energy/force, finite-difference, minimization, NVT or NPT check
is run for the incomplete Class-II candidates. This is a source-completeness
block, not a requirement for trajectory agreement. No new complete profile is
fabricated. Keep `production_validated=False` and
`simulation_readiness="not_established"`. No universal PCFF/COMPASS or PDMS
coverage, complete `.frc` coverage, production MD, pressure validation,
Ewald/PPPM closure, equilibrium or long-NPT readiness is claimed.

Validation completed: **55 focused tests**, **74 J55 regression tests**, and an
ordinary-suite union of **4,538 passed / 61 documented skips** (4,599 distinct
collected nodes). The initial ordinary run had one historical J31 GAFF2
crosslinked-network AmberTools charge-consistency failure; that exact test passed
on standalone retry. Both raw runs are retained. `retry-selection.json` binds the
failure, exact-node retry and unchanged passed-node excerpt; the union does not
hide the initial failure or count the retried node twice. No historical code,
charge or tolerance was changed to obtain that retry.

Ruff, compileall, five pip checks, staged diff checks, J48–J53/J55 receipt
validators and all **5,397 historical file hashes** passed. Nine re-signed audit
mutations were rejected by raw-source replay. Receipt/manifest validation and
relocated receipt validation also passed. The profile remains unregistered and
parameter-incomplete regardless of those software checks.
