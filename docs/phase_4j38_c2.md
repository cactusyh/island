# Phase 4J38-C2: aromatic native-source incompatibility and exact-graph charge paths

This phase clarifies the physical model and charge-source boundaries of the
existing eight benchmark cells. It does not make the three failed J38-C1 cells
executable. The acceptance driver and normal receipt validator return **1**.
Diagnostic integrity is separately checked with `--diagnostic-only` (exit 0).
No charges, parameters, numerical tolerances or historical outcomes changed.

Verified base: `4e9ce54af161669577bdfda9360e64d6e1ef5361`, fetched and compared
with `origin/codex/phase-4j38-c1-benchmark-charge-provenance` and `git ls-remote`.
Tested implementation: `0795329aed2c3bf75b7121d30aecfcd68bbc6bb0`.
The final evidence commit contains documentation/evidence only. Main was not
modified or merged. Unrelated pre-existing untracked files were preserved.

## Chemical interpretation of the native aromatic charge failure

The unchanged C1 driver was rerun before C2 implementation. It reproduced:
`OPLSAssignmentError: Native component neutrality failed: -0.23 e`.
The reproduction output (gzip-compressed without changing its bytes) and fresh
source audit are retained. These are software
source/failure diagnostics, not scientific validation evidence.

The exact declared aromatic final graph is
`f7730225baafb0a2164cdcf3db2f8813d83d3b8d145923053b0fc4ffd193f052`:
PSMILES `[*:1]c1ccc([*:2])cc1`, DP2, two independently built explicit-H finite
chains. Each chain is biphenyl, C12H10, not benzene or an infinite repeat cell.
The independently reconstructed graph matches a separate RDKit biphenyl
reference. Each repeat contributes C6H5; terminal valence is completed with real
hydrogens. Atom ordering, repeat/end-group metadata, packing remaps, molecule and
component membership are checked against the exact graph. Every atom's formal
charge is zero. There are no ionic repeat units, charged end groups or counterions.

The repository's declared original native source is the pinned Foyer OPLS XML:
SHA-256 `c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
C2 reuses the independent raw-graph/Foyer/XML audit unchanged, verifies upstream
Foyer implementation hashes, then reconstructs the decimal/component totals
again from the audited atom rows. Native stable-ID assignments are separately
compared with the external typing result. Source definitions and provenance cite
DOI `10.1021/ja9621760`; the separate fused-ring rule cites
`10.1021/ct200731v`. These citations do not establish an unreviewed replacement
charge set for biphenyl.

The selected source rows give all twelve aromatic carbons `opls_145`, −0.115 e,
and ten hydrogens `opls_146`, +0.115 e:

`12 × (−0.115) + 10 × (+0.115) = −0.230 e` per component/chain.

Two inter-ring carbons per chain have no H: IDs 4/7 and 26/29. The exact decimal
system total is **−0.460 e**, compared with formal charge **0**. Compensated sums
agree. The residual is not rounding, missing end groups, omitted/duplicated atoms,
misordered charges or incorrect component assignment. It is not an intended
ionic repeat representation of the declared C12H10 chains.

The pinned selected charge set is therefore **source-incompatible with this
neutral condensed-phase model**. A periodic engine can define a net-charged-cell
background convention, but that would be a different physical model; no such
background model is declared or accepted here. Successful engine construction
could not turn this residual into evidence for the declared neutral benchmark.
This is not a conclusion about every OPLS-AA implementation or aromatic polymer.
No independently reviewed charge-complete replacement native source has been
established. In particular, the zero-charge `opls_147` naphthalene ring-fusion
row is not a supported type substitution for nonfused biphenyl ipso carbons.

C2 adds an immutable `source_incompatible` decision and an explicit
`SourceIncompatibleError` before any new MD. The original native OPLS rejection
and all C1 records remain unchanged. All seven fresh C1 diagnostic files are
byte-identical to the retained C1 files; `c1-preservation.json` records both hashes. Individual charges were never normalized,
redistributed, subtracted, rounded or otherwise repaired.

## Exact-graph GAFF2 RESP/Psi4 asset requirement

The PE graph is
`1baea92607ac7d99dd85d04d55c22808c6bae9d6501824f7828088a2fca0b66b`;
the aromatic graph is the one above. Neither is identical to the reviewed PEO
linear graph
`1480f6da41f766852310f979bb8f820de11fe1832061b61157acd8577e52172c`
or the PEO network graph
`01ec38d6b0177fe2f3f40e747f6c53315034be9b7902fd168b21e69a83e47407`.
Even a similarly named polymer or different degree of polymerization is not an
acceptable charge source for these exact graphs.

The current declared reviewed inventory is re-extracted into a fresh temporary
location and independently audited using the J37 raw-source implementation.
Archive/input/geometry/ESP/native-charge/source/executable hashes, atom IDs and
ordering, elements, coordinates, formal and component totals, binary64 vector
identity and the finite-precision policy are checked. The inventory and reviews
are persisted in `run/c1-recheck/resp-inventory.json` and bound to the C2 matrix.
No new reviewed PE/aromatic graph-specific raw RESP artifacts were found in this
declared inventory. This is not a claim of an exhaustive search of all external
filesystems or of unavailable Psi4 software.

Both cells remain **asset_unavailable**, before parameter preparation and MD.
The required external assets are independently reviewed raw Psi4/RESP geometry,
ESP and native charge results, stable atom-ID/component/coordinate mapping,
formal totals, source hashes and executable provenance bound to **each exact
requested graph**. There is no transfer of benzene, PEO, another PE oligomer,
external OPLS or historical frozen-bundle charges. Installation of Psi4 alone
would not satisfy this graph-specific evidence requirement.

Because neither source gate passed, C2 performed **zero new benchmark MD
workflows**: no preparation, minimization, NVT/NPT or continuation for these
cells. If new exact-graph assets become available, they require an independently
reviewed source record and actual workflow/continuation evidence before the
frozen blocked outcomes can change. This phase does not introduce a bypass or a
new generic unreviewed RESP-file acceptance mechanism.

## Explicit benchmark source paths

| Cell | Selected source | C2 source/workflow outcome |
| --- | --- | --- |
| OPLS-AA PE | Native pinned OPLS | Retained J38 software evidence, read-only |
| OPLS-AA PEO | Native pinned OPLS | Retained J38 software evidence, read-only |
| OPLS-AA aromatic | Native pinned OPLS | Source-incompatible; no new MD |
| OPLS-AA network | Native pinned OPLS | Retained J38 software evidence, read-only |
| GAFF2 PE | Exact-graph RESP/Psi4 requested | Asset unavailable; no parameter preparation or MD |
| GAFF2 PEO | Exact-graph RESP/Psi4 | Retained independently reviewed J38 evidence, read-only |
| GAFF2 aromatic | Exact-graph RESP/Psi4 requested | Asset unavailable; no parameter preparation or MD |
| GAFF2 network | Exact-graph RESP/Psi4 | Retained independently reviewed J38 evidence, read-only |

The matrix remains eight cells and retains all nine C1 stage outcomes separately:
graph construction, packing, typing, charge completeness, parameter completeness,
executable backend, numerical reference, equilibration and scientific validation.
Each row now additionally distinguishes:

- `native_oplsaa`: selected only for native OPLS cases; aromatic is incompatible.
- `externally_provided`: explicitly not selected and not enabled in this tranche.
  No new external vector or adapter was introduced. An external vector could
  never be labeled native OPLS-AA, RESP/Psi4 or AM1-BCC.
- `resp_psi4_exact_graph`: retained independently verified exact-graph PEO
  sources, or unavailable PE/aromatic assets; no implied charge generation.
- Source-artifact availability, the requested graph, review identity where
  available, parameter completeness and workflow status.

Five previously successful software cells remain `retained_j38_readonly`, with
C1/J38 graph, packing, row and review bindings. They are not new C2 simulations.
The historical validators revalidate their retained workflow/analysis artifacts,
including relocated separate-process reconstruction, split/uninterrupted frame,
checkpoint and complete RNG equality. Source/trajectory hash checks retain the
existing numerical reference strategies; no new independent engine or LAMMPS
KSPACE execution is claimed.

The diagnostic bundle itself also survives relocation and fresh-process
reconstruction. Loading recomputes the raw audits, exact fixtures, chemical
source decision and matrix. Mutating a vector, origin, graph, source hash,
component total, review claim, stage or readiness flag invalidates it even when
JSON identities are consistently re-signed. Six fresh-process controls also re-sign the full receipt/manifest after
mutating source compatibility, external/native origin, RESP graph identity,
parameter completeness, acceptance or charge values. All are rejected for
the semantic mismatch. The archived mutated records are **software negative
controls, never scientific evidence**. Diagnostic validation is not an
acceptance success. The unchanged PCFF rejection and native AM1-BCC blocked
controls remain bound through the retained J38 matrix.

## Reproduction and validation

On the declared Linux server, use a new nonexistent output directory:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c2 /tmp/j38-c2-new-output
# Expected exit 1: source/asset gates unmet, diagnostics retained.
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c2_receipt /tmp/j38-c2-new-output --bundle --diagnostic-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c2_receipt docs/evidence/phase_4j38_c2/receipt.json --diagnostic-only
# Omit --diagnostic-only for the nonzero acceptance gate.
```

The exact command ledger, JUnit files, environment package versions, executable
hashes and reproduction helpers are archived with the receipt. Tests retain the
established `ISLAND_J33_C1_RUN`, `ISLAND_J34_RUN` and `ISLAND_J35_RUN` fixture
configuration. The J22–J38-C1 regression group excludes only the new C2 module
because that module is run separately as the focused suite and included again
in the full supported suite. No required test or marker was disabled.
The suites execute established short integration MD runs; those are testing,
not a new scientific campaign or new success for the blocked cells.

Historical standalone validators retain their original diagnostic/acceptance
semantics. J35's original source-hash validator uses its retained original
implementation checkout. All **1,831 tracked documentation/evidence/`.gitattributes` files**
present at the base are checked against Git blobs and current bytes. Historical JSON readability
and identities are preserved; no historical receipt or scientific test count is
rewritten. All checks used the frozen tested implementation:

| Gate | Actual result |
| --- | --- |
| Focused J38-C2 tests | 34 passed |
| J22–J38-C1 regressions | 803 passed |
| Full supported suite | 2,805 passed, 10 established skips |
| Ruff / compileall / diff check | Passed |
| All five environment pip checks | Passed |
| Historical diagnostic/receipt/hash checks | Passed with original blocked acceptance semantics |
| C2 acceptance | Exit 1: three source/charge gates remain unmet |

Exact commands, JUnit counts and skip reasons are recorded in the additive receipt.
No new test exclusions or skips were introduced.

Unmet gates remain the incompatible aromatic native source and missing exact
PE/aromatic RESP assets. Short retained trajectories do not establish equilibrium,
experimental agreement, Tg, mechanical properties, production readiness,
universal OPLS-AA/GAFF2 or PCFF coverage, AM1-BCC support, universal msi2lmp parity
or LAMMPS Ewald/PPPM validation.

`production_validated=False`

`simulation_readiness="not_established"`
