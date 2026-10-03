# Phase 4G1 — pinned Foyer OPLS-AA typing and native charges

## Base and identity

Branch `codex/phase-4g1-oplsaa-typing-charges` starts from `origin/main`
`3fe9de0132c497746961240b97eb7cb1f15ea24e` (the Phase 4F6 squash merge). Its Git tree is identical to reviewed
`8eb0e75d4f9151cd407490315cd55f8e4687f23f`. Phase 4F5 remains on its separate
research branch. No prerequisite branch was merged again; main and unrelated
checkpoint directories were not changed.

The selected implementation is **Foyer-distributed OPLS-AA**, not “latest OPLS”,
OPLS3/4 or LigParGen/CM1A-LBCC:

- Repository: <https://github.com/mosdef-hub/foyer>
- Pinned revision: `dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`.
- Foyer implementation declares `1.2.0`; revision dated 2026-09-26.
- File: `foyer/forcefields/xml/oplsaa.xml`; no includes/dependent XML files.
- XML SHA256: `c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
- XML metadata: name `OPLS-AA`, version `0.1.0`, combining rule `geometric`.
- [Pinned XML](https://github.com/mosdef-hub/foyer/blob/dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd/foyer/forcefields/xml/oplsaa.xml)
  records per-type DOIs, including Jorgensen et al., `10.1021/ja9621760`, and
  provenance from GROMACS revision `aa66efd3179a671aace62b827c54aadc3d89c1ee`.

The pin manifest in `src/island/forcefields/oplsaa/pin.json` records the XML,
license and implementation-file hashes. Runtime checks the installed Foyer
Python implementation against these hashes **before import**. No download,
source substitution or update occurs during assignment.

### License and acquisition

Foyer's [license](https://github.com/mosdef-hub/foyer/blob/dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd/LICENSE.rst)
is MIT (copyright Vanderbilt University); its README warns that sub-portions
may have independent terms. The XML identifies imported GROMACS parameter data
but does not supply an independent complete data-license statement in the file.
ISLAND therefore ships hashes and acquisition instructions, **not the XML or a
copied parameter library**. Users supply/install the exact local upstream file
under applicable upstream terms. Code was independently written; no Foyer
implementation was copied into ISLAND.

Clone upstream, checkout the exact revision, and use a separate environment:

```sh
git clone https://github.com/mosdef-hub/foyer.git /path/to/foyer
git -C /path/to/foyer checkout dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd
# In an isolated conda environment, install the pinned source's dependencies:
conda create -p /path/to/env -c conda-forge python=3.11 'numpy>=2,<2.6' \
  'gmso>=0.16.2' 'openmm>=8.4' parmed networkx lark lxml ele requests rdkit pip
/path/to/env/bin/python -m pip install --no-deps /path/to/foyer -e /path/to/island
```

RDKit is needed by the example's builder, not the typing adapter. Upstream
imports OpenMM/ParmEd/GMSO even for typing, but the adapter constructs no OpenMM
System or Context. Core ISLAND and offline record validation import none of these
optional packages. Installation does not downgrade the existing ISLAND environment.

## Public contract

```python
from island.forcefields.oplsaa import (
    load_oplsaa_source, type_atoms, assign_native_charges,
    save_opls_result, load_opls_result,
)
source = load_oplsaa_source('/path/to/foyer/foyer/forcefields/xml/oplsaa.xml')
typing = type_atoms(system, source)
charges = assign_native_charges(system, typing, source)
charges.validate_integrity(system, source)
save_opls_result(charges, system, source, 'new-charges.json')
```

The source identity binds family, variant, revision, XML hash, declared version,
typing implementation and native-charge policy. Matching only `"opls"` or type
names is insufficient. Both typing and native charges are implemented; complete
bonded assignment, ISLAND energy evaluation and dynamics are explicitly false.
Availability of the algorithm, installed optional dependencies, actual chemical
coverage and scientific suitability are separate claims.

`OPLSTypingResult` retains stable IDs, explicit external-index mapping, graph
identity, actual Foyer white/blacklists and resolved types, and dependency
versions. Foyer interprets its own SMARTS dialect: recursive type labels `%type`,
connectivity/ring predicates and override precedence are not translated into
RDKit SMARTS. The adapter uses `AtomTypingRulesProvider` and `find_atomtypes`,
with explicit bond orders and no residue-name reuse. The independent acceptance
path uses `Forcefield.run_atomtyping(..., use_residue_map=False)` on a separately
constructed ParmEd graph.

This separate result avoids pretending Foyer rules satisfy the existing
RDKit-specific `AtomTypingResult` contract. `OPLSChargeResult` binds the complete
typing result identity to each type's exact `NonbondedForce/Atom[@type=...]` charge,
unit `elementary_charge`, per-component totals and source identity. It cannot be
converted to `ImportedAmberResult` or a fully parameterized system.

The initial graph adapter conservatively requires explicit H, nonperiodic
atomistic systems, neutral uncharged sites, closed-shell valences, supported
H/C/N/O/F/S/Cl/Br/I elements and finite positive masses. Charged/isotopic/radical,
implicit-H and unsupported-valence input is rejected. This is narrower than all
825 XML types. Existing assigned atom-stereochemical metadata is preserved and
bound to compatibility; the selected source does not add stereospecific rules.
Bond semantics are preserved, not inferred from coordinates. No input topology,
metadata or coordinates are changed, and typing requires no coordinates.

Errors are `OPLSAssignmentError`; missing/incompatible installation raises
`OPLSDependencyError`. `OPLSTypingError.diagnostics` preserves the upstream error
and explicit index-to-stable-ID map for untyped/ambiguous graphs. There are no
fallback types, charges or method substitutions.

### Charge tolerance and integrity

Native charges are decimal library entries with at most four decimal places
and maximum absolute value 4 e. They are not noisy QM results. The declared
component and molecular tolerance is **N × 1e-12 e**: a conservative binary
conversion/summation allowance, far below the source's smallest decimal quantum
for the tested inventories. It is not permission to redistribute charge or round
away an actual source imbalance. Charges are copied unchanged, and every
connected component is checked independently.

Schemas `island_foyer_typing_v1` and `island_opls_native_charges_v1` use the existing
strict tagged-JSON/checksum format. Records are owned immutable strings with
owned payload access. Validation checks actual coverage, source identity,
graph/stereo/mass identity, override evidence, native charge lookup, totals and
readiness restrictions. Publishing validates first and refuses existing files.
Loading requires the compatible authoritative system and local pinned source,
but does not require Foyer, RDKit, OpenMM, ParmEd or SciPy. Coordinate changes
are permitted; graph/type/source/charge contradictions are not.

Checksums establish content consistency, not computational authenticity.
Offline validation does not rerun SMARTS or prove that arbitrary caller-supplied
match evidence was genuinely computed. Actual adapter fidelity is tested by
independent upstream invocation.

This differs from GAFF fragment AM1-BCC: OPLS charges here are source-native
per-type constants selected from the **final molecular graph**. There is no
fragment cap collapse, conservation projection, AM1-BCC/RESP fallback, or PE/PEO
fitted table. No historical GAFF contract or signed record is changed.

## Declared acceptance and execution

```sh
python -m scripts.validate_oplsaa --xml /path/to/oplsaa.xml --output /path/to/NEW
python -m scripts.validate_oplsaa --xml /path/to/oplsaa.xml --output /path/to/NEW \
  --execute-declared
python examples/oplsaa_typing_charges.py --xml /path/to/oplsaa.xml \
  --psmiles '[*:1]CC[*:2]' --dp 3
```

The declaration precedes execution and pins source, dependency versions, the nine
cases (ethane, butane, ethanol, dimethyl ether, benzene, PE DP3/DP50, PEO DP3,
polystyrene-type DP3), numerical tolerance and checks. It compares complete types
and charges exactly against upstream, checks charge entries independently using
XML Decimal parsing, and tests noncontiguous IDs/reversed atom and bond insertion,
coordinate independence and nonmutation. Alkane terminal/interior types and
alcohol/ether oxygen distinctions are explicitly checked. No embedding/QM or
AmberTools is invoked. Failed cases remain in a nonzero-exit report.

### Actual execution

The matrix ran in `island-validation/phase4g1-acceptance/` after publication of
its declaration. All nine graphs typed identically to direct upstream Foyer.
Eight native-charge assignments passed; the ninth failed neutrality and remains
failed. The aggregate CLI exited **1**, honestly preserving incomplete coverage.

| Case | Explicit sites | Native total (e) | Result |
| --- | ---: | ---: | --- |
| Ethane | 8 | 0 | Passed |
| Butane | 14 | 0 | Passed |
| Ethanol | 9 | -8.33e-17 | Passed |
| Dimethyl ether | 9 | -2.78e-17 | Passed |
| Benzene | 12 | 0 | Passed |
| PE DP3 | 20 | 0 | Passed |
| PE DP50 | 302 | 0 | Passed |
| PEO DP3 | 23 | -8.33e-17 | Passed |
| Polystyrene-type DP3 | 50 | -0.23 | **Charge assignment rejected** |

Successful cases have exactly zero maximum charge discrepancy against upstream.
Terminal alkane carbons select `opls_135` (-0.18 e), interiors `opls_136`
(-0.12 e); alcohol O selects `opls_154` (-0.683 e), ether O `opls_180` (-0.4 e).
The direct XML decimal checks confirm these constants. All cases also reached
the insertion-order/noncontiguous-ID/coordinate comparison before the PS charge
gate. The retained PS typing record remains valid; no successful PS charge
record or normalization was fabricated.

The PS source arithmetic is informative: 18 ring carbons receive `opls_145`
(-0.115 e) and 15 ring hydrogens `opls_146` (+0.115 e). One terminal benzylic CH2
receives `opls_149` (-0.005 e); two backbone benzylic methines receive generic
`opls_137` (-0.06 e). Together with the remaining alkane sites these sum to
-0.23 e, including with independent decimal summation. The pinned source has a
special CH2 benzyl rule; the observed generic CH methine typing does not supply
analogous compensation. This diagnoses the selected source's limitation for this
graph, not permission to invent a new charge/type. Complete PS coverage requires
separately reviewed source/typing evidence; 4G2 must not bypass this failed gate.

[Committed evidence](phase_4g1_acceptance.json) records the complete declaration,
identities, dependency versions, output checksums, continuous charge totals and
failed case. Retained outputs were revalidated offline with Foyer, GMSO, RDKit,
ParmEd, OpenMM and SciPy blocked: nine typing and eight charge records passed.

The isolated environment used Python 3.11, Foyer 1.2.0 at the pinned revision,
GMSO 0.17.0, NumPy 2.4.6, NetworkX 3.6.1, Lark 1.3.1, lxml 6.1.3, ele 0.2.0,
ParmEd 4.3.1 and OpenMM 8.6.1. An explicit conda package listing is retained with
the outcomes. The ordinary environment remained unchanged.

Verification: **1,013 ordinary tests passed, 9 existing opt-in tests skipped**
(83.28 s). The 22 new tests use labelled synthetic sources and mocked failures
for contracts, rather than claiming upstream numerical evidence. Ruff and pip
check passed (pip check in both environments). The actual OPLS PE DP3 example
and existing cached-fragment example passed. No AmberTools, embedding, QM,
force-field application, OpenMM Context, energy calculation or dynamics ran.

## Source forms and bounded Phase 4G2 proposal

Inspection found 825 atom-type/nonbonded records, 307 harmonic bond records,
964 harmonic angle records and 1,089 Ryckaert–Bellemans records, **all tagged
Proper**. No explicit improper, periodic-torsion, virtual-site or constraint
sections occur in this selected XML. Only 228 atom types have a SMARTS definition;
listing a type does not establish typing coverage.

[OpenMM's XML specification](https://docs.openmm.org/7.7.0/userguide/application/05_creating_ffs.html)
and the pinned Foyer loader define the relevant conventions:

- Harmonic bonds: length in nm, stiffness in kJ/(mol nm²), harmonic one-half.
- Harmonic angles: equilibrium angle in radians, stiffness kJ/(mol rad²), one-half.
- RB torsions: six coefficients c0–c5 in kJ/mol; preserve ordering and explicitly
  validate the RB cosine-angle convention during conversion.
- LJ 12-6: sigma in nm, epsilon in kJ/mol; geometric sigma **and** epsilon mixing.
- Fixed Coulomb: elementary-charge units; explicit 1–4 Coulomb and LJ scaling
  both 0.5. Shorter bonded exclusions and ring/multiple-path exception semantics
  require independent validation in the next phase.

A bounded 4G2 should resolve this exact source's bond/angle/RB/LJ/charge records,
prove complete coverage and unit conventions, then compare NoCutoff energies and
forces against independent upstream construction with correct geometric mixing.
It must test RB sign/phase, finite differences, zero-LJ charged sites and ring
exclusion/1–4 counting. Do not invent missing impropers or assume Amber's mixing
rules. Unsupported interaction families remain errors. None of that conversion
or energy/dynamics functionality is implemented in 4G1.

Upstream agreement establishes adapter fidelity only. It does not validate every
polymer or establish scientific suitability. All records retain
`production_validated=False` and `simulation_readiness="not_established"`.
