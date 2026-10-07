# Phase 4J14 — source-family-aware aliphatic amines

## Scope and baseline

Branch: `codex/phase-4j14-pcff-amine-typing`. Base:
`eaf292ec9aae8ea6dafffbb388d4e4b031713e24` (merged J13 correction, PR46).
After fetching origin, its complete tree compared equal to reviewed
`ed210c5981c83e4e65c57daedf9de544d0af4655`. Main was not modified or merged.

New, **opt-in typing profile**: `island_pcff_source_graph_v5`, with the bounded
`bounded_ref1_aliphatic_amines_v1` convention. It fixes an incompatible hydrogen
family selection for ordinary neutral and protonated aliphatic amines. This is
an explicit choice of a coherent source family, not a claim that every commercial
PCFF implementation or every nitrogen environment has a unique equivalent label.
Historical profiles v1–v4 remain unchanged and readable. Typing profile v5 is
independent of resolver policy v3/v4; neither resolver policy was modified.

Seven newly charge-complete cases have independently checked executable models
under the existing **opt-in converter-compatibility v4** policy. Every declared
case remains incomplete under strict FRC-source resolver v3. Full-source PCFF
coverage remains **false**. `production_validated=False` and
`simulation_readiness="not_established"` remain unchanged.

## Chemical authority and deliberate disagreements

The sole runtime source remains [LAMMPS pcff.frc](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc),
revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
No library or external implementation is bundled. The compact
[chemical evidence](evidence/phase_4j14_chemical_evidence.json) records exact atom,
ordinary-equivalence and positional auto-equivalence rows, versions, references
and external file hashes.

| Label | Source atom line; interpretation used here |
|---|---|
| `hn` | 105, v1/ref1: descriptive N-bound H, ordinarily equivalent to `h*` in every family. Selected on bounded neutral aliphatic `na` N. |
| `h*` | 99, v1/ref1: broad N/O-bound H family; the resolver reaches it through `hn`, not by a new alias. |
| `hn2` | 106, v2.1/ref8: another amino-H family, with its **own** ordinary and increment equivalences. It is not interchangeable with `hn`. Ref7 increment line816 pairs it with urethane `n_2`. This does not prove all possible `hn2` environments are urethanes. |
| `h+` | 100, v1/ref1: cationic hydrogen family, selected only on the explicitly formal +1 tetrahedral center described below. The H atom's formal charge stays zero. |
| `na` | 130, v1/ref1: saturated amine N; no carbonyl, aromatic or delocalized-cation interpretation is added. |
| `n4`, `n+` | 123/118, v1/ref1: four-connected ammonium N family. `n4` ordinarily maps to `n+`; no new alias is introduced. Quaternary ammonium is separately exercised and already worked historically without N–H bonds. |
| `n_2` | 129, v1/ref1 atom declaration, with later ref7 refined urethane increments. Its resonance/carbonyl environment is outside the new rule. |

The source's ref1 entry identifies the Biosym 25-Dec-1991 CFF91 family; ref7
identifies the 1993 urethane refinement, and ref8 the 1993 estimated-parameter
addition. Versions choose competing numerical rows **within** a pattern; they
are not a rule to replace one chemical type by another.

Independent implementation evidence was pinned and inspected:

* [pysimm PCFF typer](https://github.com/polysimtools/pysimm/blob/fb33814128189a99d3c8d7c4eddac2dfe44262fe/pysimm/forcefield/pcff.py),
  revision `fb33814128189a99d3c8d7c4eddac2dfe44262fe`, SHA256
  `4ed9372aef3446537e1388c1902922d953fc9d0b32bdbe3ca1c5d2b726a459c9`:
  its H-on-N branch selects `h*`. This independently corroborates the neutral
  hydrogen family, **not** every N rule or protonated assignment.
* [LUNAR PCFF typer](https://github.com/CMMRLab/LUNAR/blob/67dabeda9e6bd3cc8f968aa0c6a88730886138ef/src/atom_typing/typing/PCFF.py),
  revision `67dabeda9e6bd3cc8f968aa0c6a88730886138ef`: its sp3-amine hydrogen
  branch precedes generic N-bound H and selects `hn2`; protonated N selects
  `n4`, with attached `hn`. **J14 deliberately differs.** No LUNAR agreement is
  claimed. Its urethane recognition also uses whole-formula conditions, which
  are not adopted as final-polymer graph rules. GPL implementation code was not
  copied. The downloaded pysimm license and LUNAR hashes remain external.

The charged convention is grounded in the source's cationic atom descriptions,
formal +1 four-coordinate local graph, ordinary/automatic equivalences and the
ref1 N–C/N–H increment system. Parameter availability and neutrality alone did
not select the types. Neither independent typer establishes universal charged
PCFF typing. Aromatic amines, amides/carbamates, pyridinium and guanidinium remain
separate environments; this phase does not resolve their disagreements.

### Exact graph predicates

`amine_domains.py` starts from validated v4 decisions and changes only verified
N/H sites in already resolved components:

* N is nonaromatic, acyclic, standard mass, zero isotope/radicals, with the
  correct atomic number. Formal charge 0 requires three single bonds; formal +1
  requires four. At least one neighbor is carbon.
* Every carbon neighbor is neutral, nonaromatic, standard-state and saturated
  with four single bonds. Every other neighbor is an explicit neutral H with
  one single bond. No N–N or other heteroatom neighbor is admitted by this rule.
* Neutral N/H becomes `na`/`hn`; cationic N/H becomes `n4`/`h+`. Four carbon
  neighbors are the distinct quaternary case. NH3/NH4 without a carbon neighbor
  are not newly authorized.
* Assigned stereochemistry, masses, bonds, coordinates and provenance are not
  modified. The rule does not use names, site order, repeat index or coordinates.
  Explicit supplied labels undergo the same chemical checks.

Rule IDs, neighbor stable IDs, charge state, source rows, profile evidence and
all downstream oriented increment/parameter paths are owned and signed. The new
module participates in the existing validation-cache rule/function invalidation.
Old cache and graph-binding guarantees remain intact.

## Native charge evidence

Every authoritative bond is processed once by the unchanged resolver; no repair,
normalization, base-charge change or tolerance relaxation was made.

| Oriented source pair | FRC line | Endpoint increments (e) |
|---|---:|---:|
| c → na | 513 | +0.0827, −0.0827 |
| h* → na | 804 | +0.2487, −0.2487 |
| c → n+ | 508 | +0.4071, −0.1571 |
| h+ → n+ | 813 | +0.2800, −0.0300 |

The neutral pairs each sum to zero. The cationic pairs each sum to +0.25, so the
four N bonds contribute exactly +1.0000 in independent Decimal summation;
ordinary hydrocarbon bonds cancel componentwise. Both endpoint values and their
orientation are retained. No hydrogen acquires formal +1 merely from its label.
Native float sums meet the existing **1e-12 e** component/molecular tolerance.

The former `na–hn2` and `n+–h*` searches were absent from the exact source. Under
v5, declared aliphatic N–H environments reach the appropriate ref1 rows. This is
a versioned chemical typing change; it is not an alias or resolver fallback.

## Declaration, before/after outcomes and remaining failures

The [15-case declaration](evidence/phase_4j14_declaration.json) preceded production
edits and acceptance. All coordinate seeds are 2026. Original v4 results were
retained before the update, and both unchanged resolver policies were inspected.
[Outcomes](evidence/phase_4j14_outcomes.json) contain endpoint contributions,
component sums and exact missing interaction families; full per-interaction
searches remain in the external strict-inspection and acceptance directories.

| Case | SMILES or PSMILES | v4 charge | v5 charge | v5 model, v4 compatibility |
|---|---|---|---|---|
| Primary | `CN` | incomplete | complete, 0 | complete |
| Secondary | `CNC` | incomplete | complete, 0 | complete |
| Tertiary | `CN(C)C` | complete | complete, 0 | complete (existing chemistry) |
| Primary protonated | `C[NH3+]` | incomplete | complete, +1 | complete |
| Secondary protonated | `C[NH2+]C` | incomplete | complete, +1 | complete |
| Tertiary protonated | `C[NH+](C)C` | incomplete | complete, +1 | complete |
| Quaternary | `C[N+](C)(C)C` | complete | complete, +1 | complete (existing chemistry) |
| Amine oligomer | `[*:1]CCN[*:2]`, DP2 | incomplete | complete, 0 | complete |
| Closed-shell sidechain | `[*:1]C(C[NH3+])C[*:2]`, DP1 | incomplete | complete, +1 | complete |

The rebuilt declared oligomer and sidechain files were **byte-identical** to the
retained J10 and J13 inputs respectively:
`32e981fd9e21fed5f3613b146bea6b925001bdb61a3565d6e7c9bb601c50d4d3`
and `4ba78264aa8a60f6a50afd85a18f266dc09a8292bc34ad061021bdd2cbfa43d0`.
The oligomer contains both terminal primary and internal secondary nitrogen,
actual hydrogen chain ends and inter-repeat bonds. Assignment uses the final graph.

Six declared negative near-misses remain unchanged:

* Acetamide `CC(=O)N`: typing/charges complete; BB/BA/MBT missing and proper-torsion
  candidates conflict. Public rejection is at model assembly.
* Carbamate `COC(=O)N`: existing `n_2/hn` family mismatch, two missing increments.
  No new `hn2` carbamate override was inferred here.
* Aniline `Nc1ccccc1`: existing `nn/hn2`, two missing increments; no neutral
  aliphatic rule may consume aromatic carbon neighbors.
* Pyridinium: two `nh+–cp` and one `nh+–hn` increments remain missing; no `nh`
  replacement.
* Guanidinium: printed source contributions remain **0.9999 e**, residual
  **−0.0001 e** against formal +1; complete row coverage is not charge validity.
* `C[NH]`: radical typing rejection, before native charge/model work.

The new acceptance CLI requires all 15 names exactly once, a bound declaration,
correct stage and chemistry diagnostic, and independent positive labels/charges.
Unexpected success, unrelated exceptions, wrong reason or missing cases fail.
The original J13 negative contracts still run with **historical v4 typing**:
peroxide fails model assembly, sulfoxide/charged-terminal radical fail typing,
and the old closed-shell sidechain still fails charges. Their original receipts
are unchanged; new J14 observations do not rewrite them.

Counts (same 15 fixtures): typing **14/15 → 14/15**; native charges **3/15 → 10/15**;
compatibility models **2/15 → 9/15**; strict-v3 models **0/15**. Seven new cases
have independent numerical/workflow evidence. The v5 ledger contains bounded
predicates for **87/133** labels (v4: 86/133), not 87 globally verified labels.
The [additive ledger](evidence/phase_4j14_coverage.json) keeps all 133 labels,
family statuses and exact unresolved next actions. The full source remains
incomplete, including charged/aromatic families, context-dependent metals and
surfaces, missing cross terms, nonzero Wilson equilibrium interpretation and
empty torsion–torsion scope. v4's non-cp BB13 zeros remain explicitly
**converter-derived**, never source rows or strict coverage.

## Independent numerical and durable-workflow verification

[Declarations](evidence/phase_4j14_vertical_declarations.json) and
[numerical receipts](evidence/phase_4j14_numerical.json) retain all seven new cases,
identities, source rows, commands, artifact hashes and the failed first reference.
External root: `../island-validation/phase4j14-declared`.

`pcff_j14_reference.py` supplies independent bounded C/H/N labels. The independent
Decimal raw-FRC reader rebuilds interaction inventories, source selections,
coefficients, equilibrium dependencies and charges without importing the native
PCFF typer/resolver. Those independently supplied types/charges enter optional
CAR/MDF for the compiled msi2lmp reference; no native coefficient list is exported.
msi2lmp does **not** validate automatic typing or increments and is never required
by runtime preparation, evaluation, persistence or resume. No `-ignore` is used.

Pinned executable SHA256 values:

* msi2lmp: `a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`.
* LAMMPS: `db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`
  (30 Sep 2026 build of the pinned revision).

Reference settings: finite, unconstrained OpenMM Reference / LAMMPS all-pairs;
explicit LJ=(0,0,1), Coulomb=(0,0,1). Source constants and units unchanged.
Initial and `0.015*sin(arange(3*N)+0.3)` Å perturbed coordinates, then minimum
and steps 2/4 are compared. Energy and grouped-component atol **1e-5 kJ/mol**;
force-component atol **1e-5 kJ/(mol·Å)**; rtol **2e-10**. Central differences use
**1e-4, 1e-5, 1e-6 Å** and converge to the analytical forces. No tolerance,
seed, charge, parameter or starting coordinate was adjusted after failure.

### Retained converter failure and separate correction

The first protonated-tertiary comparison failed by **0.028352771945 kJ/mol** in
the angle–angle/total grouping. `tertiary_protonated-vertical/` remains failed.
The converter assigned the `c,n+,c,c` triple of couplings as
`[-1.5155,-4.2781,-4.2781]`; center-preserving raw-source matching requires
three instances of **line4143**, `-1.5155`. Generic reversal moved the N apex.
This is the already documented converter central-role problem, now exercised
by an asymmetric charged amine.

Evidence: pinned [GetParameters.c lines1007–1053](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L1007)
invokes generic `find_match` for AA patterns; its generic reverse can move B.
Pinned [improper_class2.cpp](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/improper_class2.cpp#L621)
uses K1 ABC×CBD, K2 ABC×ABD, K3 ABD×CBD. The separate
`raw_source_center_preserving_j14_v1` reference reconstructs those three lookups
and theta dependencies from raw source and the converter's own inventories.
Raw output, each source match and explicit override are retained. The H5
ABC/ABD/CBD equilibrium-role correction is also retained. A focused synthetic
asymmetric regression rejects central-atom movement. Native formulas and model
identities did not change. The separately declared corrected-reference run
passed all gates; it does not erase the failed converter comparison.

| Newly complete case | Min iterations/evaluations | Final max force | Max E/component difference | Max force-component difference |
|---|---:|---:|---:|---:|
| Primary | 34 / 40 | 0.03731 | 2.49e-14 | 1.54e-12 |
| Secondary | 60 / 65 | 0.08804 | 8.79e-14 | 1.12e-12 |
| Primary protonated | 51 / 60 | 0.09864 | 7.64e-14 | 1.22e-12 |
| Secondary protonated | 60 / 66 | 0.07086 | 8.53e-14 | 1.53e-12 |
| Tertiary protonated, corrected oracle | 57 / 62 | 0.08605 | 5.68e-14 | 1.37e-12 |
| Amine oligomer | 116 / 123 | 0.09748 | 1.07e-13 | 2.63e-12 |
| Protonated sidechain | 99 / 102 | 0.08824 | 7.11e-14 | 1.49e-12 |

Energy units are kJ/mol; force units kJ/(mol·Å). Worst final-displacement FD error
is **1.17e-7 kJ/(mol·Å)**. All seven minima met the unchanged 0.1 maximum atomic
force criterion and independent fresh verification within 5,000 iterations and
10,000 evaluations. Initial/final energies, RMS force and all counts are in the
receipt. A local minimum is not equilibration.

Every new case ran four BAOAB steps (0.1 fs, 300 K, 5/ps; velocity seed78123,
thermostat seed99181), whole versus 2+2 split. Per segment: 4 evaluations/3 frames;
whole: 6/5; split cumulative: 8 evaluations/5 deduplicated frames. Bundles and
runs were relocated; continuation ran in separate child processes. All common
frame differences were zero on this host with exact complete RNG equality
(atol1e-10, rtol1e-12). Context counts were 5 for setup/minimum/first segment,
3 for whole reference and 3 for child continuation. This is orchestration
consistency, not independent thermostat validation or a bitwise portability claim.
Completed status/frame/no-op resume worked with scientific imports blocked.
Rechecksummed model-source and policy tampering were rejected.

## Public use and reproduction

The native API accepts constructed `MolecularSystem` objects. Use
`type_pcff_atoms(system, source, profile="island_pcff_source_graph_v5")` followed
by `assign_automatic_pcff_charges(..., resolution_policy=...)` for inspection.
For preparation select the same typing profile and **explicit** resolver policy
in `PCFFOptions`. Prepared bundles retain the profile and validated native
identities; reconstruction does not rerun scientific preparation. OpenMM/SciPy
remain lazy; optional CAR/MDF and executables are reference-only.

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
ROOT=../island-validation/phase4j14-declared

$PY examples/pcff_amine_typing.py --source "$FRC"
$PY examples/pcff_amine_typing.py --source "$FRC" --strict-source  # expected 1
$PY scripts/validate_pcff_amines.py --source "$FRC" --output NEW/acceptance
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" --declaration "$ROOT/strict-declaration.json" --output NEW/strict

# Each case declaration (also stored in phase_4j14_vertical_declarations.json)
# binds the already checksummed input coordinates. Use a new output directory.
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/primary-declaration.json" --output NEW/primary
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/secondary-declaration.json" --output NEW/secondary
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/primary_protonated-declaration.json" --output NEW/primary-protonated
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/secondary_protonated-declaration.json" --output NEW/secondary-protonated
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/tertiary_protonated-declaration.json" --output NEW/original-tertiary  # retained reference failure
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/tertiary_protonated-corrected-reference-declaration.json" --output NEW/corrected-tertiary
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/amine_oligomer-declaration.json" --output NEW/oligomer
$PY scripts/validate_pcff_profile.py --source "$FRC" --declaration "$ROOT/protonated_sidechain-declaration.json" --output NEW/sidechain

$PY scripts/audit_pcff_completion.py --source "$FRC" --declaration docs/evidence/phase_4j11_declaration.json --retained-matrix ../island-validation/phase4j10-declared/parity-matrix --output NEW/full-source.json --require-full-source  # expected 1
$PY scripts/inspect_pcff_tranche_controls.py --source "$FRC" --declaration docs/evidence/phase_4j13_declaration.json --output NEW/j13-main
$PY scripts/inspect_pcff_tranche_controls.py --source "$FRC" --declaration docs/evidence/phase_4j13_charge_declaration.json --output NEW/j13-sidechain
$PY scripts/revalidate_pcff_tranche.py --source "$FRC" --root ../island-validation/phase4j8-declared/acceptance --root ../island-validation/phase4j9-declared/vertical-slice --root ../island-validation/phase4j10-declared/fused-vertical --root ../island-validation/phase4j12-declared/vertical-slice --root ../island-validation/phase4j13-declared/peo-vertical --root ../island-validation/phase4j13-declared/thioether-vertical --output NEW/historical.json

$PY -m pytest -q tests/test_pcff*.py tests/test_psmiles*.py tests/test_prepared_workflow.py tests/test_prepared_bundles.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

## Verification and preservation

Final executed counts and log hashes are recorded in
[verification](evidence/phase_4j14_verification.json). Real-source receipts are
separate from synthetic parser/typing/gate tests. No QM was run. The four-step
workflows and numerical checks were confined to the seven newly complete cases.
Six historical J8/J9/J10/J12/J13 workflows/bundles were revalidated read-only,
with OpenMM/RDKit/SciPy/Foyer/ParmEd imports blocked. Historical source/evidence,
failed diagnostics, bundles and checkpoints were checksummed before and after.
Unrelated directories were not changed. The failed original J14 reference is
preserved alongside its separately corrected experiment.

The remaining gates are full-source chemical/interaction coverage and scientific
suitability. In particular, v3 still requires source rows for missing BB13
couplings; v4's converter compatibility does not satisfy that strict gate.

Executed final results: **1,735 passed / 10 skipped** ordinary tests (333.12 s),
**593 passed** affected regressions (235.92 s), and **30 passed** new focused
regressions. Ruff and both pip checks passed. The new example passed; its strict
variant exited 1 with the exact missing BB13 list. The historical PEO example
retained model identity `b0057359aa7a07700259412f98ce5c0dfb98246eb3f43006005f7c95c93d793a`
and energy 82.36346573972186 kJ/mol. All **3,543 historical path entries / 2,739
unique files** matched their original hashes. The full-source command exited 1;
its unchanged historical-v4 audit still reports 47 unresolved labels, while the
additive v5 bounded-predicate ledger has 46 globally unresolved labels. Neither
count is a count of globally validated chemical domains.
