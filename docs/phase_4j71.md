# Phase 4J71 — PCFF literature reading and parameter extraction

The decision is **source_conflict**. The supplied literature was read, and 467 source rows containing 1,470 nonblank numerical cells were extracted. No coherent source contract resolves the 378 J70 interaction instances. The reference bundle remains **incomplete**, with zero accepted coefficients, no execution profile and no automatic assignment.

The actual base is fetched origin/main, **4c9066e4501591211049c44401e82b3110118f58**, containing required J70 **9040f4bf6b9948c68dabcd85d7dd25a40e5a6536**. The branch is codex/phase-4j71-pcff-literature-parameter-extraction. Main was neither merged into nor pushed by this phase. Historical J54–J70 evidence, original observations/producers, source rows and charges remain unchanged.

## What was actually read

The pre-reading [inventory](evidence/phase_4j71/source-inventory.json) froze 25 files and their original hashes. Its initially unknown citation fields remain preserved. The [enriched inventory](evidence/phase_4j71/literature-file-inventory.json), [reading ledger](evidence/phase_4j71/reading-status.json) and [document notes](evidence/phase_4j71/document-reading-notes.json) provide established metadata and page coverage.

All **15 primary PDF documents, 287 pages**, were read completely: 12 main papers and three supporting documents. Six checkpoint PDFs have byte-identical content to read originals. One checkpoint is truncated at 2,097,152 bytes and could not be read; the complete Sun1995 original was read. The adjacent nitrate-esters PDF has a different hash but duplicate paper content. These copies do not increase unique document coverage.

| Document | Pages | Parameter-bearing material and consequence |
|---|---:|---|
| Sun1995, polysilanes, Macromolecules28,701–712, DOI10.1021/ma00107a006 | 12 | SectionIII equations, Table5 and Appendix p11; 77 rows. “PDMS” here means poly(dimethylsilane), a Si–Si backbone. |
| Sun/Rigby1997, siloxanes, Spectrochim.ActaA53,1301–1323 | 23 | Eq1a–b and Appendix pp21–22; 28 rows. General cross-term equations do not supply omitted coupling coefficients. |
| Sun1994, aromatic polyesters, J.Comput.Chem.15,752–768 | 17 | CFF93 ester/aromatic environments; table/equation index retained, no siloxane substitution. |
| Sun/Mumby/Maple/Hagler1994, polycarbonates, JACS116,2978–2987 | 10 | CFF93; Table10 p9 and carbonate equations do not establish siloxane authority. |
| Sun1998, COMPASS overview, JPCB102,7338–7364 | 27 | Explicitly distinguishes COMPASS from PCFF; AppendixB pp24–25 is alkane/benzene data. |
| Sun/Ren/Fried1998, phosphazenes, CTPS8,229–246 | 18 | COMPASS P/N/organic parameters, Appendix pp16–18; not Si/O rows. |
| Rigby/Sun/Eichinger1997, PEO, PolymerInternational44,311–330 | 20 | Ether/alcohol appendix pp17–20. Types o2e/o2h and c4o do not become siloxane osi. Filename year1999 is not the printed publication year. |
| Bunte/Sun2000, nitrate esters, JPCB104,2477–2489 | 13 | COMPASS nitrate tables and AppendixA; nitrate O is not bridging Si–O–Si O. |
| Patil/Mark2000, aromatic polysiloxanes, CTPS10,189–195 | 7 | pcff2 means COMPASS1.0; crystal/RIS results are not new Class-II coefficients. |
| Smith/Borodin/Smith2004, PDMS, JPCB108,20340–20350 | 11 | Table3 HMDS charges, Table4 constants, Eq6 torsion/bend modulation; another model. |
| Xiang et al.2025, PDMS benchmark, JPCB129,1864–1873 | 10 | Main paper and archive association, ref22. Charge plots were not digitized. |
| Xiang et al. supporting information | 14 | TablesS1–S5, pp2–7: UA, HUANG, COMPASS, OPLS and Dreiding; commands remain separated by model. |
| Jorge/Milne/Barrera/Gomes2021, PolCA, ACSPhysChemAu1,54–69 | 16 | Tables2–6, Eqs9–15; united-atom, rigid bonds, half-prefactor harmonic angles and RB torsions. |
| PolCA supplementary material | 76 | Entire SI read, including rotated pages; TableS72 p65 repeats16 RB rows. Property tables are not PCFF coefficients. |
| Jorgensen/Maxwell/Tirado-Rives1996 OPLS SI, DOI10.1021/ja9621760 | 13 | ja11225.pdf is SI, not the main paper. Every page was rendered and read manually without OCR. Table6 compares AMBER/CHARMM22. No Si parameters. |

The associated ZIP archives were inventoried by member hash. The benchmark README, small inputs, style declarations and applicable coefficient rows were read. Complete static topology parsing selected used coefficient IDs: **124 archive numerical rows** were independently checked against exact member lines. Large unused generic parameter tables, production coordinates and charge-command files are not claimed as manually read papers. No archived simulation input was executed. The engineering README's suggested zero fallback is excluded.

No full papers, full extracted texts, page images or large archives are committed. The benchmark main paper explicitly states CC-BY4.0. An initially unsupported CC-BY label for PolCA was corrected: its supplied text establishes authors' copyright/open access, but the exact redistribution grant was not located. Originals remain external. Contract v2 records that correction, the exact PolCA SI title and complete SI text versus the absent OPLS main paper. No numerical source cell changed.

## Extracted values and source authority

The [record catalog](evidence/phase_4j71/parameter-records.json) binds every candidate row to a source SHA, citation, page/table/row/field coordinates or archive member/line. It separates raw decimals/units from normalized values and conversion factors. Scalar potentials and their analytic Cartesian chain-rule forces share one equation. Equilibria, phases, multiplicities, orientation, missing fields and transfer concerns remain explicit.

The 105 historical Sun rows were compared after fresh full reading and manual Appendix checks. Their historical review producer is preserved; J71 invents no native/QM execution identity. Examples confirmed from the actual documents:

- Sun1995 Si–C: r0=1.899 Å; K2=189.65, K3=-279.42, K4=307.51, with energy ΣK_n(r-r0)^n.
- Sun1997 Si–O: 1.640 Å,350.1232,-517.3424,673.7067; Si–O–Si:159.0°,8.5000,-13.4188,-4.1785, with angular displacements in radians.
- Sun1997 distinct Si4 and Si4c LJ9-6 rows: (rmin,epsilon)=(4.405 Å,0.198 kcal/mol) and (4.290 Å,0.131 kcal/mol). Different typed sites are not values to average.
- Smith2004: Si–O r0=1.651 Å; Si–O–Si theta0=137.63°,10.305,-18.101,10.100; exp6 nonbonded terms, half-strength1–4 and negative-half-prefactor torsions. Its sine-squared bend/torsion modulation is not a Class-II angle-torsion coefficient.
- PolCA: Si–OB0.1640 nm; Si–OB–Si149.5° and k=61.3 kJ mol^-1 rad^-2 in k/2 Δθ²; UA RB torsions and LJ12-6. The printed OH–H length **0.945 nm** remains extraction_unresolved, not silently changed to0.0945.
- Sun1995 AAT uses cosine; it does not resolve native FRC literal-Phi/AAT authority. Sun1997's phase-containing general equation and coefficient-only Appendix rows are separately qualified.

Blank cells remain null. Printed zeros remain source evidence, without authority to zero-fill another source or graph. Some locations use logical table row/column coordinates because text extraction has no unambiguous glyph box. Numeric transcriptions are not verbatim prose quotations.

The [human table](evidence/phase_4j71/source-evidence-table.md) groups all extracted rows into the requested bond, angle, coupling, torsion, charge and nonbonded families. The [machine table](evidence/phase_4j71/source-evidence-table.json) adds target IDs, compatibility, missing conventions, author confirmations and fitting restrictions. Non-siloxane parameter tables are indexed in reading notes rather than transferred into HMDS.

## Conflicts and the 378-instance mapping

Two independently checked conflicts determine the decision:

| Declared environment | Benchmark SI TableS3 | Associated archive, checked topology |
|---|---|---|
| C–Si–O angle | 113.1855,36.2069,-20.3939,20.0172 | type102, O–Si–C:114.9060,23.0218,-31.3993,24.9814 |
| C–Si–C angle | 114.906,23.0218,-31.3993,24.9814 | type29, C–Si–C:113.1855,36.2069,-20.3939,20.0172 |

Archive atom types6/13/39/44 were checked against masses and the complete atom/bond/angle/dihedral topology, rather than guessed from IDs. Source records, archive/member hashes and line numbers are retained in the [conflict report](evidence/phase_4j71/source-conflicts.json). No value is averaged, selected for a better fit or relabeled to conceal disagreement.

Other obstructions include reversed Si–C–H reference lengths, the SI's dihedral-coefficient commands after its improper-style declaration, and its UA hybrid description versus the archive's LJ12-6 style. These are source/converter convention questions, not proof that native PCFF is physically wrong. Sixteen equal PolCA main/SI RB duplicates remain separate records.

The [mapping report](evidence/phase_4j71/mapping-report.json) preserves all **378 exact instance IDs**, stable sites, ordered types, graph identity and equilibrium dependencies:78 AA sites and60 sites in each of five torsional coupling families, including BB13. Sparse candidate lists contain 654 compatible_but_non_native comparisons and54 ambiguous comparisons. Every absent pair is explicitly defined as not_applicable. Individual angle-angle mappings preserve both the central and shared atom; only I/L exchange is allowed. LAMMPS AA triplets remain aggregate non-native candidates. These are conservative screening results, not exact assignments. **Zero instances are exactly resolved.** UA, OPLS, GAFF, COMPASS and chemically different oxygen environments are not promoted to native PCFF.

The Sun/Rigby primary appendix is the strongest supplied authority for its printed siloxane scalar terms, but supplies no complete missing-coupling contract. Published benchmark code remains external COMPASS evidence, not original PCFF author-primary material. The 28 confirmations remain user-attested; original author artifacts remain not_yet_supplied.

## Bundle, API and bounded decision

The [partial bundle](evidence/phase_4j71/partial-literature-bundle.json) is a non-executable reference catalog. Its executable native, author-primary and external inventories are empty; all 378 remaining IDs are listed. No complete PCFF/PDMS profile is created.

The [converter report](evidence/phase_4j71/strict-converter-report.json) and [static closure report](evidence/phase_4j71/static-closure-report.json) record **nonactivation**. J71 made zero converter attempts because no complete coherent bundle exists. J70's strict returncode17 failure at endbonddihedral hc c3 sio c3 remains historical evidence, not a new J71 run. No permissive conversion or manufactured zero row was used. Historical available-subset evaluator closure is not complete literature-bundle closure.

The ordinary prepare_polymer(..., forcefield="pcff") interface remains unchanged. Actual HMDS preparation still reports 378 missing interactions and refuses an execution-ready profile. No manual profile or charge-method requirement is introduced. The [capability report](evidence/phase_4j71/canonical-capability.json) lists exact blockers.

The bounded next dependency is **source clarification**: resolve TablesS3/archive labels and reference ordering with the benchmark authors, and obtain exact authoritative Si/O coupling rows with a native Phi/AAT and equilibrium-dependency contract. No messages were sent. Source absence is not automatically a QM-fit requirement. [Target references](evidence/phase_4j71/qm-targets.json) preserve J70's exact unexecuted templates; no fitting is authorized. No generic sampling, QM or MD allocation follows.

J68's bound-free incompatibility and J67/J68 failures remain independent, unchanged results. Parameter completeness alone would not establish reproduction of those QM targets.

## Verification and replay

All 93 J71 tests pass. They cover re-signed source/citation/unit/equation/location/convention mutations, family/resolution separation, ambiguity, zero/wildcard rejection, duplicate/conflict preservation, converter gating, actual fail-closed preparation, receipt/manifest semantics, untouched validation and fresh-process replay without PDFs. These controls are software verification, not new QM evidence.

The first focused run retained 27 failures from treating canonical JSON strings as bytes and a brittle decimal-format assertion. A corrected overlapping run retained one formatting-assertion failure. The terminal J71 run passes93/93 after the shared-neighbor regression was added. Failed XML/log records and initial Ruff findings are retained. Two failed orchestration syntax calls did not execute shell commands. Optional PDF tooling was isolated; an unused incompatible installation was excluded. The ambient shell Python's pip check also reports existing Torch2.5.1 requirements (missing fsspec and incompatible sympy1.14.0). That environment is not dependency-clean and was not repaired. The supported project environments and a clean isolated PDF-reader environment pass pip checks; the latter independently reproduces all readable source page hashes. The ambient failure is retained as an advisory, not counted as a passing check.

The unique terminal total is **5,660 tests: 5,247 passed, 413 skipped, zero failures/errors** (5,567 historical-suite nodes plus93 J71 nodes). Historical byte checks preserve8,722 files. Focused reruns add zero unique coverage.

The supported full suite, J54–J70 receipt regressions, Ruff, compileall, pip checks, diff checks and historical preservation checks are recorded in validation/. Unique terminal coverage is the union of the historical-suite collection and disjoint J71 nodes; focused reruns add no coverage. Native/MD integration tests excluded by the phase guard are skips, not new campaign evidence. Older validators with broad implementation globs use pinned historical snapshots; J71 has an explicit implementation-file inventory.

Offline replay needs repository source and retained parent evidence, not local PDFs, network access, converter or optional solvers:

~~~bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/replay_phase_4j71.py
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/validate_phase_4j71_receipt.py docs/evidence/phase_4j71
~~~

Use the supported interpreter (../island-validation/phase4g1-env/bin/python here) if shell Python lacks dependencies. Add --deep-parent to replay for strict J70 receipt validation. Optional verification of mounted originals and exact archive lines:

~~~bash
PYTHONPATH=.:src python scripts/read_phase_4j71_sources.py
~~~

extract_phase_4j71.py rebuilds reports from pinned reviewed cells; it does not automatically read a new paper or infer a parameter. A matching digest alone never establishes scientific declarations or source authority: validation reconstructs expected tables and typed declarations.

Software verification passes the recorded checks. Source completeness remains incomplete; converter completeness and complete static closure are blocked. No parameter is accepted. Exact PCFF equivalence, chemical transfer, scientific validation and production readiness are not established. New QM attempts, QM fits, MD steps and reserved-validation executions are zero. **production_validated=False; simulation_readiness="not_established".**
