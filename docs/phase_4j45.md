# Phase 4J45: external literature domain gap closure

No complete executable domain can be established from the supplied local assets.
The selected target is **Xiang et al.'s external OPLS-compatible methyl-PDMS
model, Table S4**, DOI `10.1021/acs.jpcb.4c08471`. Its 17 coefficient rows and
four two-column charge rows are independently reconstructable, but do not
establish complete units, typing, charges, pair conventions or execution
provenance. Acceptance exits **1**. No profile, parameter assignment, converted
physical row, charge vector, prepared system, numerical MD, checkpoint or
trajectory is created.

## Exact base and preserved prerequisites

After fetching origin, `origin/main` is
`fdf8e8058173faa9c013703aa4910b3065820e7c`, through J42-C1. J43 and J44 have not
merged. Branch `codex/phase-4j45-complete-external-domain` starts at that exact
main; prerequisite commit `2c36c62` carries all 20 committed J43/J44 files
byte-for-byte from J44 snapshot `95531e5509b262265583e43c078afc34daa5c54a`.
Main is not modified or merged. The reproducibility manifest protects every
tracked file in that J44 snapshot, including all J17–J44 receipts, bundles,
serialized identities, sources, tests and validators. Native PCFF, native OPLS,
GAFF/GAFF2, C5 and all existing workflows remain unchanged.

## Selection from the committed audit

J44's recursive inventory and gap reports are reused: 52 files, 136 coefficient
candidates, ten graph-only fixtures. J45 checks those known source hashes without
rediscovering files or repeating their bibliography. The J44 catalog and new
selection record are independently checksum-bound. The selection record retains
the original 17 Table S4 candidates unchanged, plus 34 exact one-based PDF
page/line citations from the original files. Source metadata and all 52 original
hashes remain available through the pinned J44 catalog and J45 manifest.

Table S4 is the strongest bounded **all-atom methyl-PDMS** target because its
four LJ pairs, three harmonic bonds, six harmonic angles and four OPLS torsions
cover all literal bond/angle/proper paths of hexamethyldisiloxane. Literal element
matches are explicitly candidate associations, never assigned types or accepted
physical coefficients. No global numerical coverage percentage or optimality
claim is inferred from this count.

The alternatives are considered separately in `selection.json`. COMPASS S3 has
contradictory bond-equilibrium dependencies and unresolved improper/cross terms;
HUANG S2 omits mixed angle couplings. Smith 2004 explicitly uses Buckingham
exp-6, quartic bonded terms and a sine-squared bend–torsion coupling; it requires
its own functional form and cannot supply LJ or ordinary Class-II rows. Its
Table 4 lacks Si–H valence rows and leaves higher-order cells blank; those cells
are not converted to zero. PolCA uses a united-atom/polarization representation.
The original OPLS supplement has restricted redistribution and incompletely
reconstructed graph-specific conventions. The supplied API library requires
model-specific charge replacement and off-atom mapping. The PCFF archive adds
no coverage beyond the historical pinned FRC. DMP refinements and coordinate
files do not establish a complete state-specific P model. None supplies a
complete alternative domain without unresolved dependencies or cross-source
borrowing. Their original J44 blocked findings remain unchanged.

## Exact remaining gaps

The six final-graph cases are disiloxane, hexamethyldisiloxane, explicit-H capped
PDMS repeat, PSMILES PDMS DP2 and DP4, and methoxytrimethylsilane as a mixed
Si–O–C negative control. The repeat is `[*:1]O[Si](C)(C)[*:2]`; the existing
builder uses explicit methyl and hydroxyl terminal groups, then completes H.
Those actual terminal environments are audited; the repeat string is not a
force-field typing rule.

Each fixture's `gap-report.json` lists every atom/type/mass/charge request,
bond, angle, proper path, improper applicability and special pair, with graph
identity, exact atom order, molecule/component membership, element labels,
matching candidate evidence checksums and missing dependencies. All assignments
remain missing even when a literal coefficient token is found.

- Si–H types, charges, bonds, angles and torsions are absent for disiloxane.
- O–H and O–C terminal environments in the capped repeat/DP2/DP4 lack exact
  source typing, charge and bonded terms. Methoxytrimethylsilane exposes the
  same O–C boundary. These cases return `typing_ambiguous`.
- Hexamethyldisiloxane returns `parameter_incomplete`: it fits the table's
  elemental methyl/Si/bridging-O motifs, but that does not prove typing.
- The interior/end charge columns lack exact atom-by-atom terminal rules.
  Both raw columns are preserved, unassigned. No charge is changed, balanced
  or normalized; a plausible neutral vector is not independent authority.
- Table S4 does not explicitly bind energy/length/force-constant units, exact
  masses, unlike-pair mixing, 1–2/1–3 exclusions, 1–4 LJ/Coulomb scaling or
  improper/inversion applicability. These cannot be borrowed from native
  OPLS, another article or a converter default.
- `dihedral_style opls` is explicit; independent equation/conversion/reference
  reconstruction is still unqualified. There is no silent RB conversion.
- The main article supplies the OPLS cutoff **1.1 nm** and PPPM context. This
  partial evidence is preserved; it does not establish a complete units and
  backend-input policy for all boundaries. The article's CC-BY notice does not
  automatically establish the supplement's redistribution authorization.

The immutable `DomainSelection` and `DomainClosure` records reuse J43 pinned
reads and J44 graph/provenance contracts. Every result distinguishes
source-complete, typing-complete, charge-complete, parameter-complete,
executable, numerically verified, workflow verified and scientifically validated;
all eight gates remain false. Missing source files return `asset_unavailable`;
changed sources return `source_incompatible`. `require_executable()` always
raises a classified `DomainClosureBlockedError` preserving the full result and
missing dependencies. This data layer registers no engine or external profile.

## Reproduction and verification

```bash
export PYTHONPATH=.:src PYTHONDONTWRITEBYTECODE=1
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
REFS=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/island_ffref
python scripts/phase_4j45.py --references "$REFS" \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --selection docs/evidence/phase_4j45/selection.json --output /tmp/j45-study
# Expected exit 1: writes the bundle-generation refusal and exact missing rows.
python scripts/phase_4j45.py --reconstruct /tmp/j45-study --references "$REFS"
# Offline independent source reader; local pypdf==6.20.0 is required.
python scripts/audit_phase_4j45_sources.py --references "$REFS" \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --selection docs/evidence/phase_4j45/selection.json
python -m pytest -q tests/test_phase_4j45.py
python scripts/validate_phase_4j45_receipt.py docs/evidence/phase_4j45
```

The receipt/manifest record the tested implementation commit, exact commands,
focused/full tests and J22–J44 regression subset, Ruff, compileall, five pip
checks, diff checks, prior receipt validators, source reconstruction and every
historical hash. Relocation and fresh-process reconstruction validate blocked
records and source bytes. Mutation controls reject changed graphs, atom order,
component membership, rows, units, charges, profile identity, source hashes,
evidence, stages and readiness. No numerical qualification is claimed for an
incomplete model; energy/forces, minimization, MD and checkpoint/RNG continuation
remain unrun behind the completeness gate.

`production_validated=False`, `simulation_readiness="not_established"`.
Polymer Ewald/PME remains unvalidated, backend pressure unavailable, cutoff
recovery and long NPT blocked. This phase establishes no universal coverage,
converter parity, experimental agreement, equilibrium or production readiness.
