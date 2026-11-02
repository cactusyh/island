# Phase 4J69 — authoritative Si–O source and Class-II dependency audit

The decision is **`source_incomplete`**. The available pinned PCFF Si–O bond, Si–O–Si angle and fixed cross terms close numerically against the retained research evaluator. No newly demonstrated parsing, equilibrium, units or prefactor error justifies changing that evaluator. The pinned distribution lacks native coefficients for 378 applicable coupling interactions, and original author artifacts remain unavailable. Closure of the available subset cannot establish native PCFF equivalence or a unique physical explanation for the QM residuals.

The fetched main commit `8784955934ef6b750d54b3cfea49d824a1218c76` contains required J68 `c957d9f327a025e57052531baed184c110e313aa`. The branch `codex/phase-4j69-authoritative-sio-source-audit` starts from that fetched main. No main merge or push is performed. Unrelated and untracked files are preserved.

## Frozen inputs and reconstruction

`plan.json` was frozen before the new numerical audit, SHA256 `a72ab77e0e4ebc0ff52a9b7287cfef01569b33aa0394a5dfc8ab2f92c96cf36a`. It binds the source archives and exact J54/J67/J68 inputs, 27-atom neutral-singlet explicit-H graph, stable order, all 26 observations and their original producer mapping, fixed evaluator, 58 physical columns, one energy-only intercept, original bounds, protected22/full26 limits, J68 diagnostic columns and original failures. Strict J68 receipt reconstruction passed before execution. Historical observations are reused with their original identities; no execution identity is fabricated.

The retained target remains PBE0/def2-TZVP, DF/def2-universal-jkfit, unpruned 250×1454 grid, SCF energy/density thresholds 1e-12/1e-10, integral threshold 1e-12, no dispersion, Psi4 1.9.1. J69 executes **zero new QM attempts, zero MD steps, zero validation calculations and zero research optimization calls**. It reconstructs existing energies, gradients and arrays. The three reserved validation configurations remain untouched. Si–H and capped-DP2 investigations remain inactive.

The J68 bounded rho≈1.196374166 and bound-free rho≈1.196230308, their independent SCS checks, checked lower bounds above one, bridge-angle and Si–O failed diagnostics, and complete evaluator comparisons are preserved. J66 original22 feasibility, J67 expanded-data incompatibility and the disclosed J67 timing deviation are unchanged. Removing the bounds did not make the existing research span feasible. This audit neither re-solves nor reclassifies those fits.

## Source authority and dependency graph

`source-assets.json` inventories the exact external paths, uncompressed and archived SHA256 hashes, pinned revision and provenance. All required available assets are archived additively under `sources/` for offline replay. The source is the LAMMPS-distributed **pcff.frc**, SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`, revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`. Its `cff91` namespace is part of this exact PCFF distribution. It is not a substitution from another force-field family. Internal references, including the silicon row reference4, are retained verbatim. No paper DOI is substituted for distribution identity. All 27 inventoried assets are available. `inherited-source-authority.json` also rechecks the J60 authorities already pinned by the frozen model, including real-unit constants, dielectric, boundary and mixing sources. The RLB, typing templates, complete converter C/header sources, relevant LAMMPS Class-II implementations and official equation documentation are also pinned. Original author artifacts remain `not_yet_supplied`; the 28 confirmations remain user-attested.

`source-row-inventory.json` preserves the actual J54 row IDs, source locations, raw text, units, conversions, citations and parser identities. The independent Decimal ASCII oracle is checked against the actual ISLAND source parser and exact/ordinary-equivalence resolution. No permissive wildcard or lower-order fallback is added. `source-dependency-graph.json` binds each graph interaction to source rows, stable sites/types, scalar energy and force equation, orientation, multiplicity, units, transformed values, equilibrium dependencies, generated LAMMPS representation and evaluator implementation. Metadata, converted fixed terms, absent native coefficients and research unknowns remain distinct. `pair-source-dependencies.json` binds all351 atom pairs to their source LJ and charge records, graph exclusions and derived mixing. LJ96 uses minimum-energy rmin, sixth-power mixing, dielectric1 and Coulomb332.06371; graph distances1/2 are excluded and distances≥3 carry full weight. These are unchanged declared isolated-research conventions; the FRC alone does not establish their native exclusion authority. Missing terms are never silently zeroed.

The 378 absent interaction rows comprise 78 angle-angle sites and 60 sites in each of five torsional coupling families. Type symmetry reduces these to the historical 56 unknown coupling coefficients; two further quadratic Si–C/C–H research columns yield58. These are explicit research variables, not discovered native rows. Wilson terms use the declared three-connected-center convention: HMDS has no such center. Degree-four angle-angle slots remain explicit. The empty torsion-torsion section supplies no numeric zero authority. Literal FRC linear-Phi AAT and the selected research cosine AAT remain explicitly different.

## Executed numerical closure and converter result

The audit independently reconstructs available source terms directly from raw rows, uses scalar coordinate equations, differentiates their conservative potentials, and compares every fixed family with the unchanged complete evaluator. Central differences at 1e-5 Å check all 81 Cartesian components at all 26 geometries. Reports retain family and interaction extrema, responsible stable atoms/components, configuration IDs and units. Analytic sums cover Si–O and other quartic bonds, Si–O–Si and other quartic angles, bond-bond, bond-angle, fixed angle-angle, source torsions, LJ96 and Coulomb. Full MM array comparisons explicitly supply all 58 research coefficients and the single energy intercept.

The maximum available-source energy difference is 1.78e-15 kcal/mol; maximum force-component difference is 4.76e-14 kcal/mol/Å and maximum vector-force difference is 4.89e-14 kcal/mol/Å. The largest family finite-difference force-component discrepancy is 8.70e-8 kcal/mol/Å, within the frozen 2e-5 derivative tolerance. This is numerical closure of the declared source subset and research conventions, not validation against physical experiment.

One **non-QM converter attempt** compiled and ran unmodified msi2lmp v3.9.11 / 6 Sep 2024 from pinned sources. The exact retained center, supplied source types and unchanged charges were exported through the existing CAR/MDF interface. Strict Class-II conversion, without `-ignore`, terminated with code17: `Unable to find endbonddihedral data for hc c3 sio c3`. Commands, input files, compiler output, executable hash, stdout, stderr and runtime are retained. No complete msi2lmp-generated coefficient set exists for this graph. The 173 graph interactions represented by the printed bond/angle/BB/BA tables match the independently generated representation to 3.56e-15 in their respective coefficient units (`converter-partial-closure.json`). These printed values provide a partial conversion check, not authority for unprinted sections. The generated per-interaction LAMMPS coefficient representation in the dependency graph is therefore labelled a research representation; it is not presented as a successful converter export. No LAMMPS simulation or MD is run.

## Equilibrium, units, phases and propagation

The Si–O source bond row is `quartic_bond:cff91:1775`:

```
r0 = 1.6562 angstrom
E = 306.1232 dr^2 - 517.3424 dr^3 + 673.7067 dr^4  kcal/mol
```

It applies to (Si2,O5) and (O5,Si6). The Si–O–Si row is `quartic_angle:cff91:2365`:

```
theta0 = 157.026 degree = 2.740620711236616 radian
E = 9.0740 dt^2 - 19.5576 dt^3 + 8.5000 dt^4  kcal/mol
```

There is no extra1/2. Angular displacements and angular coefficient dimensions use radians; only equilibrium angles and phases convert from degrees. Internal energy-bearing coefficients multiply by 4.184 to kJ/mol; execution divides them back. Equilibrium lengths and angles are not divided by 4.184. The source references propagate unchanged into each recorded bond-bond, bond-angle, angle-angle and torsional coupling dependency. Bond/angle ordering and end-block exchanges are tested; angle-angle swaps only the outer arms at fixed center/shared arm. Tests also check conservative forces under proper rotation and translation.

O–Si–O is not an HMDS interaction: each silicon in this exact graph has only one oxygen neighbor. The available O–Si–O source row is `quartic_angle:cff91:2417`, with equilibrium 110.693°; it is inventoried without assigning a nonexistent HMDS interaction. Any applicable source query outside this graph is recorded separately; no O–Si–O evaluator coverage or PDMS interior coverage is inferred.

The FRC prints plus-cos main torsions while LAMMPS uses minus-cos. The converter copies source phases unchanged. Under the *same declared Phi*, literal converter phases differ from the printed source scalar potential: the largest atom-vector force contrast over the retained 26 is about6.69 kcal/mol/Å. This historical convention issue was already disclosed by J60/J63, whose research representation adds180° to preserve the source plus-cos equation exactly, including constants. Source geometry authority is insufficient to turn this algebraic contrast into a newly established native-PCFF implementation error. FRC AAT linear-Phi versus research cosine remains unresolved as native source authority, and no applicable numerical AAT row is invented.

J68's shared angular and Si–O polynomials use `a sum(q-q0)+b sum(q-q0)^2`, with forces from that scalar potential. They leave source equilibria and all other dependencies unchanged. They are **not globally propagated equilibrium revisions**. Their reduction of rho does not prove that row1775 or2365 is physically wrong. No new independently motivated source issue was found that would authorize a J69 hypothesis, so **zero hypotheses and zero new minimax solves** are activated. No diagnostic coefficients are accepted.

## Residual localization and decision

`residual-family-localization.json` retains vectors and residual-direction projections for O5, Si2, Si6, C3 and H10 at every one of the26 configurations, including joint002/032/034 and retained joint020/028/040. It separates total fixed force, the explicit58-column contribution, Si–O bonds, bridge angle, other bonds/angles, BB/BA/AA, torsions, LJ96 and Coulomb, plus projections along both J68 diagnostic feature pairs. Original bounded, bound-free and both diagnostic vectors are reported separately. Sum-of-family norms, norm of their sum and cancellation differences prevent magnitude-only causal interpretation.

For the retained bound-free vector at joint-002/O5, the local residual is 2.392461 kcal/mol/Å. The Si–O bond force has norm 19.5810 and projection −19.4686 along the residual direction; the bridge-angle force has norm 5.0006 and projection +0.4369. The sum of separate family norms is 46.217 while the net MM force norm is 18.893 kcal/mol/Å. This cancellation and opposite projection show why the largest family norm alone cannot identify a source error. At joint-020/Si6, the family norms sum to 58.864 while their net is 0.250 kcal/mol/Å.

At C3 and H10, both Si–O diagnostic feature-force vectors are exactly zero across all 26 configurations. Their residual reduction under the J68 Si–O amendment therefore comes through refitting the other58 columns, not a direct force from that local polynomial. A fully propagated equilibrium revision could affect neighboring atoms through cross terms, but that is a different potential and has not been justified or tested here.

The available bond/angle baseline has no demonstrated source-to-evaluator discrepancy. The existing research Hamiltonian remains incompatible with its frozen training limits. Native source incompleteness and unresolved original conventions prevent attributing the limiting O/Si/C/H residuals to a uniquely wrong source parameter or asserting closure of a complete native PCFF model. No implementation correction is justified by this audit.

Supplying native numerical values for the already free56 coupling columns alone cannot remove the checked bound-free incompatibility. A future source finding would need to justify a baseline, convention or dependency change outside that unchanged affine problem. Source incompleteness limits attribution; it is not presented as proof that merely filling those coefficient cells will fix the fit.

The next dependency is specific: obtain authoritative, hash-bound missing Si/O Class-II numerical rows and native Phi/AAT/pair-convention interpretation, with original author artifacts, and reconcile them with this dependency graph. Then freeze any genuinely source-supported conservative model amendment separately and retest the *existing* protected22/full26 constraints. Additional QM cannot remove the demonstrated incompatibility of those existing constraints. No new QM allocation, generic sampling, bound widening, extra parameter family or automatic combination of J68 diagnostics follows from this result.

## Verification and offline replay

The implementation is additive and research-only. Explicit implementation inventories prevent future phase modules from changing historical receipt scope. Strict validation checks scientific declarations independently of outer digests, source/archive identities, exact inventories, original producer mappings, historical bytes and unique terminal test nodes. Tests mutate source SHA/rows, parser conversions, phases, equilibria, ordering, prefactors, missing terms and re-signed declarations. They check derivative closure, reversal/rotation, duplicate/null diagnostic columns, historical bounded/free classification, reserved-validation nonactivation and fresh-process relocation without importing Psi4 or solvers.

Preliminary failures are retained in `validation/`: metadata rows were initially sent to the energetic-row resolver; initial focused tests were started before reports completed; reversal exposed an audit-only family-label orientation bug; hard links across filesystems failed during relocation; a later metadata-edge schema inconsistency was caught by fresh-process replay. The resolver, audit label, metadata-edge records and relocation copy were corrected. The partial converter-table intake initially omitted reverse-ordered sites; it now checks all 173 applicable graph interactions. An early report represented pair fields in source order; the final generated LAMMPS representation correctly uses epsilon,rmin. Neither error changed historical source values or the MM evaluator. Preliminary reports and failed test logs remain separate from terminal checks. The missing Ruff executable at the initially assumed environment path was corrected to the available Ruff command.

The supported full collection is run under the explicit J69 zero-QM/zero-MD policy. Retained verified native fixtures are reused. Tests that require fresh native calculations or dynamics steps, including manufactured integrator execution, are scope skips rather than passes. Static energy/force and synthetic conic software controls remain enabled. Terminal unique coverage is 5,482 cases: 5,069 passed and 413 skipped, including 68 J69 cases and 944 J54–J68 cases. 352 cases are explicit J69 execution-scope skips. All 8,499 base files are byte-preserved. Original full-worker J69 results remain in the raw archive; the final 68 J69 cases replace their overlap exactly once, while every historical case is retained byte-for-byte. Terminal unique coverage and all scope skips are recorded in `collection.json` and `validation/terminal-test-counts.json`; focused reruns add no unique coverage. Applicable J54–J68 regressions, strict parent receipts, Ruff, compileall, environment pip checks, diff checks and byte preservation are recorded in the receipt.

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j69.py

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j69.py --deep-parent

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/validate_phase_4j69_receipt.py \
  docs/evidence/phase_4j69
```

Software verification, source completeness, source-to-evaluator closure, diagnostic correction, parameter acceptance, chemical transfer and scientific validation remain separate. This exact-HMDS protocol-specific audit establishes no complete-basis claim, native PCFF equivalence, capped-DP2 coverage, methyl-PDMS transfer, condensed-phase accuracy or production readiness. Ordinary `forcefield="pcff"` behavior is unchanged; no source row is overwritten and no research candidate is registered.

```
accepted_coefficients = 0
automatic_assignment = False
production_validated = False
simulation_readiness = "not_established"
```
