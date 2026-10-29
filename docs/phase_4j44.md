# Phase 4J44: local literature parameter audit

The supplied collection contains numerical parameter data, but this audit does
not establish a complete, independently reconstructed and authorized executable
profile for any requested fixture. Acceptance exits **1** (`parameter_incomplete`).
No literature profile, prepared numerical system, MD trajectory or checkpoint is
created. Candidate transcriptions are evidence, not usable force-field rows.

## Base and historical preservation

`origin/main` was fetched and verified as
`fdf8e8058173faa9c013703aa4910b3065820e7c` before J44 edits. Main contains
J42-C1 but does not contain J43. The nine published J43 contract/evidence files
were carried unchanged from `be8d94ec6e526d8cf8251ceb7c1baef86d609ee5` in
prerequisite commit `e737567`. This preserves the requested main base while
reusing J43's provenance and decimal reconstruction contracts. The manifest
compares every historical file under `docs`, `src`, `scripts` and `tests` with
that J43 snapshot, including all earlier receipts and serialized identities.
No native PCFF/OPLS source, C5 schema or existing workflow was edited.

## Local collection

All **52 files** were inventoried recursively, including the hidden notebook
checkpoint and duplicate PDFs. There are **47 PDFs, four ZIPs and one XLSX**.
Each filename has its own size and SHA-256. The checked-in `source-audit.json`
records title, DOI, authors, year, journal, access/license evidence, category,
parameter-row presence, PDF page/table locations and archive member hashes.
Unknown bibliographic fields remain null; a filename, download stamp or open
access to a paper does not establish redistribution rights for another file.
No original PDFs, spreadsheets or parameter archives are republished here.

Sources are grouped by actual contents into Class-II PCFF/COMPASS, OPLS-AA,
PDMS/siloxane, other extensions, explanatory sources and incompatible/restricted
sources. Source family is retained separately from eligibility. For example,
`pcff/ja11225.pdf` is the original OPLS-AA supplement, with an explicit
personal-use/no-redistribution notice. The Rigby PEO paper's title page gives
1997 despite its filename saying 1999. The XLSX contains published cross-engine
energies, not force-field coefficients. The Murzyn ZIP contains PDB/GRO
coordinates and benchmark states, not parameters.

## Numerical evidence and exact blockers

The PDMS benchmark supplement (Xiang et al., DOI
`10.1021/acs.jpcb.4c08471`) supplies **136 coefficient-line candidates** across
separate UA, HUANG, COMPASS, OPLS and Dreiding tables. Every candidate preserves
source hash, DOI, page/table, extracted one-based line, exact raw text, labels,
numeric columns, original decimal strings and evidence checksum. No numerical
conversion was attempted where units or policy are unresolved. Explicit source
zeros remain in the evidence; no missing coefficient is replaced by zero.
J43/C5's independent decimal readers reconstruct these tokens, and a separate
command re-extracts their original PDF pages and compares them directly.

The source-specific findings are retained in `dependency_findings`:

- COMPASS Table S3 gives `Si-O = 1.64`, but the C-Si-O bond-bond and bond-angle
  dependencies both use `1.8995, 1.8995`. C-Si-C coupling instead uses unlike
  equilibrium lengths `1.8995, 1.64`. These inconsistencies are preserved.
- The HUANG table lacks the mixed C-Si-C, C-Si-O and H-C-Si angle coupling rows.
  Its OPLS torsions cannot supply omitted Class-II torsion cross terms.
- The COMPASS table declares `improper_style class2`, then supplies
  `dihedral_coeff` lines without a complete improper/angle-angle inventory.
- The external OPLS-PDMS table does not establish explicit units, special-pair
  scaling, mixing, exact terminal atom/charge mapping and improper applicability.
  The OPLS name alone cannot supply those missing declarations.
- These methyl-PDMS tables do not define Si-H types/charges, Si-H bonds and
  associated angles/torsions/couplings for disiloxane or Si-H termination.
- Smith et al. (`10.1021/jp047434r`) use Buckingham exp-6 with pair-specific
  parameters. Jorge et al. (`10.1021/acsphyschemau.1c00014`) use PolCA united-atom
  alkyl sites and polarization corrections. Neither is an all-atom Class-II
  bundle, and neither is simplified or relabelled to make it executable.
- Phosphazene validation/ESP charge tables do not supply a full state-specific
  Class-II library. Murzyn's dimethyl phosphate anion and neutral DMPH remain
  distinct; neither transfers parameters to other P states.
- The supplied PCFF archive contains exactly the historical FRC bytes:
  `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
  It adds no source coverage. Its README's missing-improper-to-zero instruction
  is rejected. No license file is present in that archive.
- The OPLS-based API library explicitly declares geometric mixing and 0.5
  LJ/Coulomb 1–4 scaling, but its README replaces charges for different models.
  It supplies no Si/P labels or graph-bound ethane/benzene assignment. It is
  isolated from native OPLS and from the other papers.

The graph-only cases are disiloxane, hexamethyldisiloxane, capped explicit-H PDMS,
PSMILES PDMS DP2/DP4, methylsilane, ethane, benzene, ethanol and the explicitly
anionic dimethyl phosphate graph. The existing C5 builder validates PDMS
termination, valence and hydrogen completion. Requirements enumerate each
atom/type/mass/charge, graph interaction, special pair and Class-II equilibrium
dependency. Nothing is marked covered merely because a table contains numbers.
Native historical controls remain covered by the regression suite; they are not
relabelled as literature-profile successes.

## Data contracts and boundaries

`LiteratureAudit`, `LiteratureSourceReview` and `LiteratureBundleAudit` are
immutable, checksum-bound data records. Candidate bundles require an explicit
path and independent SHA-256. Row kinds cover typing, charge, bond, angle,
torsion, improper, cross term, nonbonded, special 1–4, exclusion and mixing
records. Accepted candidate tokens require source citations, functional form,
units, conversion identity/version and declared conventions. They still confer
no chemical or execution authority. Conflicting duplicates, wildcards,
family/paper mixing, unsupported units/conversions, stale graphs, changed
component/molecule memberships, source mutations and evidence mutations reject. Review reproduced an additional label-evidence gap:
changing an atom label while preserving its cited number initially passed. The
corrected validator requires the exact ordered labels in the primary raw row,
outside the numerical columns; a dedicated regression preserves this rejection.

Acceptance distinguishes inventory, citation, typing, charge, parameter,
backend, independent numerical and scientific stages. Only the inventory is
complete. Bibliographic review was performed, but unresolved citations/licenses
and source dependencies prevent a complete citation/qualification gate.
Executable and numerical gates remain false. `require_executable()` raises a
classified `LiteratureBlockedError` preserving the complete result. There is no
new numerical adapter or registry entry. Therefore energy/force, minimization,
NVT/NPT and continuation qualification are not run for incomplete literature
models. Relocation and fresh-process reconstruction test the audit records.

## Reproduction

Use the repository's supported environment for the audit and tests; the offline
PDF reconstruction uses the locally installed `pypdf==6.20.0`. No command
performs network discovery or downloads.

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
REFS=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/island_ffref
CATALOG=docs/evidence/phase_4j44/source-audit.json
python scripts/phase_4j44.py --references "$REFS" --catalog "$CATALOG" --output /tmp/j44-study
# Expected exit 1; writes classified blocked outcomes and exact missing dependencies.
python scripts/phase_4j44.py --references "$REFS" --reconstruct /tmp/j44-study
# Run with the local Python containing pypdf; verifies all original files and archive members.
python scripts/audit_phase_4j44_sources.py --references "$REFS" --catalog "$CATALOG"
python -m pytest -q tests/test_phase_4j44.py
python scripts/validate_phase_4j44_receipt.py docs/evidence/phase_4j44
```

The receipt and manifest contain exact commands, tested implementation commit,
focused/full results and the J22–J43 regression subset, five pip checks, Ruff,
compileall, diff checks, independent source reconstruction, relocated/fresh-process
reconstruction and historical byte comparisons. Receipt-only reconstruction
checks archived evidence without requiring redistribution of the original
assets; the separate raw-source command requires the supplied local originals.

`production_validated=False`, `simulation_readiness="not_established"`.
Polymer Ewald/PME remains unvalidated; backend pressure is unavailable; cutoff
recovery and long NPT remain blocked. No universal coverage, converter parity,
experimental agreement, equilibrium or scientific production claim is made.
