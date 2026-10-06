# Phase 4J2 — explicit executable source fallbacks

## Baseline and completion contract

Branch `codex/phase-4j2-pcff-executable-source-fallbacks` starts at fetched
`758449bb125ed706d2e982b34b4a2c2d8f27d56c`. Its tree
`e0fa331c321e61bd106ff6f55ba7064c03ea4204` exactly matches reviewed correction
`f142e606705b588b6bc66629461829e5a7b60df4`. Corrected I3 `787deb0` is an ancestor
of main. There is no unmerged prerequisite. Main and unrelated notebook/checkpoint
work are preserved.

Full-source PCFF coverage remains **partial**. This phase adds executable
supplementation, not a new claim that all parsed rows have chemical rules or
complete models. Source interpretation, automatic typing, native charges,
parameter coverage, executable models and independent verification have separate
gates. Numerical agreement is not scientific suitability.

`production_validated=False`; `simulation_readiness="not_established"`.

## Exact source and reference audit

The external source remains LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
The local checkout and complete file bytes were verified. Independent raw scanning
again finds 133 distinct atom types, 134 ordinary equivalences, 108 automatic
equivalences and 564 bond increments. J1's full numerical inventory and all 133
labels remain in the separate J2 ledger. No source file or external code is bundled.

Pinned implementation evidence:

- [README](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README),
  LIMITATIONS: no automatic-equivalence supplementation or bond increments.
  Atom types are supplied by the input. Successful conversion cannot attest typing,
  charge assignment or fallback semantics. `-ignore` is never used.
- [InitializeItems.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/InitializeItems.c),
  lines 33–75 select quadratic versus quartic and torsion_1 versus torsion_3
  by force-field class; this is not mixed-model automatic supplementation.
- [ReadFrcFile.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/ReadFrcFile.c),
  `ReadFrcFile` validates the Class II r-eps nonbonded representation and loads
  explicitly named families. [MakeLists.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/MakeLists.c)
  builds interaction inventories from connectivity, not coordinates.
- [GetParameters.c](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c):
  bonds 119–128, angles 175–185 and torsions 299–311 search direct then ordinary
  family equivalence; cross-term assignments 220–233 and 396–515 preserve
  directional blocks and equilibrium dependencies. Non-cp BB13 initialization/
  lookup remains the explicitly named historical converter policy.
  `find_match`/`match_types`, lines 1055–1175, try exact then leading-asterisk
  wildcard matching, including numbered wildcards. The converter chooses by file
  order. J2 deliberately rejects unresolved competing coefficients instead.
  `find_angleangle_data`, lines 1007–1052, still calls generic reversal that can
  move the central atom: J1's original failed comparison and separate central-role
  correction remain unchanged. ISLAND does not adopt that defect.
- FRC automatic-equivalence header at lines 323–330 distinguishes bond-increment,
  bond, angle-end/apex, torsion-end/center and OOP-end/center columns. It does not
  authorize a universal equivalence map.
- FRC equations at lines 1010, 1781, 2431, 3234 specify quadratic bond/angle,
  plus-cos single torsion and harmonic Wilson forms. For executable checks,
  [bond_harmonic.cpp](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/MOLECULE/bond_harmonic.cpp),
  [angle_harmonic.cpp](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/MOLECULE/angle_harmonic.cpp),
  and [dihedral_fourier.cpp](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/EXTRA-MOLECULE/dihedral_fourier.cpp)
  preserve the no-half prefactor and plus-cos convention (Fourier computation,
  lines 156–180). Wilson retains the H4/J1 signed mean-of-three-angle equation
  from pinned `src/CLASS2/improper_class2.cpp`.

The pinned LUNAR `ff_functions.py` (67dabeda9e6bd3cc8f968aa0c6a88730886138ef)
was also inspected: `charge_from_bond_increments` lines 217–285 has the previously
reported unreachable automatic `elif`; it is not used as a charge oracle.
Its explicit positional roles and limited wildcard tests inform the audit, not
an imported implementation. No external implementation code is copied.

## Opt-in interpretation and exact preservation of historical behavior

`resolution_policy="island_pcff_positional_fallbacks_v1"` is required to enable
supplementation. Omitting it preserves historical CHO/J1 selection, schemas and
fingerprints. The public `resolve_pcff_source_record` accepts the same explicit
policy for inspectable queries; query results alone are not complete models.

For each ordinary family the new policy tries direct labels, then ordinary
family equivalence. Within each tier exact patterns precede wildcard patterns;
more constrained patterns take precedence, with conflicting equally specific
oriented coefficients reported as ambiguous. Versions are selected per ordered
source key, never across unrelated patterns. A numbered wildcard is recognized
only at the beginning (`*`, `*1`, ...); `h*` remains a literal type label.
This is an explicit conservative ISLAND policy, not a claim to recover every
commercial heuristic or the converter's file-order behavior.

Only **missing** ordinary base records trigger the automatic family:

| Requested base | Supplemental namespace/family | Positional columns |
|---|---|---|
| Quartic bond | cff91_auto / quadratic_bond | bond / bond |
| Quartic angle | cff91_auto / quadratic_angle | end / apex / end |
| torsion_3 | cff91_auto / torsion_1 | end / center / center / end |
| Wilson | cff91_auto / wilson_out_of_plane | end / center / end / end |

Direct automatic labels precede positional automatic equivalence. Ambiguity does
not trigger a lower-priority guess. Every result retains stable IDs, original and
transformed labels, equivalence row evidence, candidate rows, versions, legal
permutations, selected coefficients, namespace and preceding failed lookup.
Angle reversal preserves the apex. Proper reversal swaps directional cross-term
blocks. Wilson permits six peripheral permutations with center fixed. Angle-angle
permits outer-arm exchange with center/shared-arm fixed, never full reversal.

Lower-order equilibria feed the existing explicit dependencies, including source
row identities. **All existing cross-term requirements remain active**. A missing
coupling does not become zero because its base became quadratic. Only the existing
named non-cp BB13 policy supplies its narrowly defined zero. Empty torsion–torsion
remains an operational scope exclusion; nonzero Wilson equilibrium remains
unsupported. These unresolved cases cannot be prepared as complete models.

Quadratic terms remain explicit two-coefficient forms, without pretending that
absent cubic/quartic coefficients were source zeros. `torsion_1` uses
`K*(1+cos(n*phi-phase))`, cis=0, preserving the constant even for n=0. Existing
Class II `torsion_3` keeps its historical minus-cos interpretation. Kcal converts
by 4.184; lengths remain angstrom; equilibrium angles/phases convert to radians,
angular force constants do not receive degree factors. OpenMM expressions use
10×distance(nm), with returned forces converted to kJ/(mol*angstrom).

Native charge supplementation uses the dedicated automatic **bond_increment**
column only after missing direct and ordinary bond-equivalence rows. Endpoint
values and reversal orientation are retained. Zero base, component/formal totals
and 1e-12 e tolerance are unchanged. No wildcard charge rule, normalization,
external charges or hn2→hn alias is introduced. The actual na/hn2 automatic labels
remain na/hn2 (source lines 395/373); the only hn2 increment is hn2–n_2 at line 816.
Both declared amine cases therefore remain incomplete.

## Records and public path

Additive v2 typing/charge/assignment/model records preserve v1 contracts. Parameter
and model validators recompute the chosen policy and reject rechecksummed changes.
`island_pcff_source_graph_v2` retains J1 predicates and adds neutral closed-shell
homonuclear singly bonded halogens, using valence-one source atom labels and
source-matched increments/parameters. This is graph recognition, not a molecule
name dispatcher. Dichlorine's independently audited rows are atom type line 88,
increment line 693 (explicit 0/0), auto-equivalence line 356 and quadratic bond
line 1287 (1.988 Å, 236.5339 kcal/(mol Å²)). Other environments remain separately
bounded; no new isotope/radical or unqualified generic fallback is enabled.

```python
request = ForceFieldRequest("pcff", PCFFOptions(
    source_path=frc_path, lj=(0, 0, 1), coulomb=(0, 0, 1),
    typing_profile="island_pcff_source_graph_v2",
    resolution_policy="island_pcff_positional_fallbacks_v1",
))
prepared = prepare_forcefield(system, request)
```

The same keyword is available on `assign_automatic_pcff_charges` and
`assign_pcff_parameters`. `assign_pcff_source_types(..., profile=...)` shares the
validated graph contract for explicit labels. It does not bypass chemical checks.
Source parameters remain separate from provided charges.

Facade adoption, native persistence, I2 bundles and I3 workflows retain their
schemas. They delegate to v2 native validators. The shared offline
`pcff_evaluation_identity()` derives the new implementation/profile identity;
evaluator and workflow do not duplicate it. Historical formulas are unchanged.
Offline inspection/loading uses core Python/NumPy and explicit hash-verified
external source files, with no OpenMM Context, Foyer, AmberTools or QM.

## Acceptance declaration and commands

The committed `phase_4j2_declaration.json` was published before implementation and
numerical execution. It revisits all 14 J1 fixtures and declares water, silane,
dichlorine and formaldehyde; construction seed 2026, original and asymmetric
0.015 Å sinusoidal perturbation, unchanged source charges, and explicit LJ/Coulomb
(0,0,1). Cross-term failures are retained. No QM, alternate seeds or force-field
source modifications are involved.

Canonical comparison tolerances: energy/forces atol 1e-5, rtol 2e-10. Kernel finite
differences use 1e-4, 1e-5, 1e-6 Å. The workflow uses the first newly complete case,
0.1 force criterion, 5000 iterations/10000 evaluations, and unchanged line-search/
verification defaults; BAOAB four steps, two segments of two, 300 K, 5/ps, 0.1 fs,
seeds 78123/99181, segment budgets four evaluations/three frames and uninterrupted
six evaluations/five frames. Split tolerances remain 1e-10/1e-12 with exact RNG.

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY examples/pcff_source_fallbacks.py --source "$FRC"
$PY scripts/validate_pcff_fallbacks.py --source "$FRC" \
  --lammps /external/pinned/build/lmp --output /new/j2
$PY scripts/validate_pcff_fallback_workflow.py --source "$FRC" \
  --bundle /new/j2/dichlorine/bundle --evidence /new/j2/dichlorine/outcome.json \
  --output /new/j2-workflow
```

The assignment/term CLI's default gate is explicitly bounded; the separate workflow
command must also pass for J2 operational acceptance. `--require-full-source` and
the historical full-source CLI remain nonzero while required coverage is missing.

## Measurements and remaining gates

Measured receipts and test counts follow below.
Original J1 numerical comparisons, failures, corrected reference and workflow
receipts remain unchanged. J2 does not establish full-source completion, universal
PCFF chemistry, periodic condensed-phase MD or crosslinked-network readiness.

### Executed source/model results

| Measure | Reviewed J1 (14 declared cases) | J2 on the same 14 | J2 including 4 additional cases (18) |
|---|---:|---:|---:|
| Automatic typing complete | 14/14 | 14/14 | 18/18 |
| Native charges complete | 12/14 | 12/14 | 16/18 |
| Complete operational parameter models | 8/14 | 8/14 | 12/18 |
| Models with all required executable term forms | 8/14 | 8/14 | 12/18 |
| Independently checked whole systems | 7 original passes + 1 separately corrected reference | Same 8 prior results, not rerun | Prior 8 + 1 new dichlorine check; 3 additional models not numerically attested |

The 133-label recognition ledger still contains **66 implemented labels and 67
unimplemented/ambiguous labels**. The new graph predicate adds an elemental-halogen
environment, not another label. Kernel/selection checks exercise representatives
of all four supplemental families, not all 1,186 automatic numeric rows
(628 bonds, 330 angles, 216 torsions, 12 Wilson rows). Parsing those rows is not
1,186 independent physical validations.

Exact unresolved model interactions after the new searches:

- Acetal: **16**, down from 19. `quadratic_angle:3,4,6` selects source line 1899;
  proper sites `[2,3,4,6]` and `[3,4,6,7]` select automatic torsion line 2544.
  Independent raw-row checks confirm the angle 109.5°/70 kcal and torsion
  0.13 kcal, n=3, phase=0. Required cross terms remain missing. No absent-coupling
  policy was invented to complete it.
- Acetone: **6** missing bond–angle interactions, all type pattern
  `c_0–c3–hc` (ordinary `c_0–c–h`).
- Chloromethane: **15** missing bond–bond, bond–angle and angle–angle couplings,
  with Cl/C/H roles retained. The existence of the C–Cl base parameter is
  insufficient.
- Ester: **39** unresolved interactions/dependencies, retained individually in
  the J2 gap ledger.
- Methylamine and amino-alcohol: respectively two and one na–hn2 bonds with
  no increment after direct, ordinary bond and automatic bond-increment searches.
  No successful parameter model is fabricated from their partial charge vectors.

Water, silane and formaldehyde have complete models with ordinary records; they
are not counted as recoveries caused by automatic fallback. Dichlorine is newly
complete through its automatic quadratic bond and the explicitly documented
valence-one graph extension. This small halogen case does not establish missing
organic cross-term coverage.

All eight previously complete J1 models have exactly unchanged numerical
coefficient, ordered-site, equilibrium, charge, nonbonded and pair-policy content
under the new policy. Their new interpretation/provenance identities deliberately
differ; their original records and historical identities were not replaced.

### Numerical evidence

A separate LAMMPS executable was built from the verified, clean pinned checkout:

```sh
cmake -S ../island-validation/phase4h4-upstream/lammps/cmake \
  -B ../island-validation/phase4j2-tools/build -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=off -D BUILD_OMP=off -D PKG_CLASS2=on \
  -D PKG_MOLECULE=on -D PKG_EXTRA-MOLECULE=on
cmake --build ../island-validation/phase4j2-tools/build -j 8
```

It reports `30 Sep 2026 - Development - patch_30Sep2026-65-ge891a3e10`.
Build logs, help/package listing, commands, executable checksum, reference inputs,
outputs and force dumps are retained separately. The historical executable was
not rebuilt. Its initial missing-Fourier failure remains in the first receipt.

| Check | Maximum energy error (kJ/mol) | Maximum force error (kJ/(mol Å)) |
|---|---:|---:|
| Dichlorine whole system, two geometries | 8.71e-14 | 9.10e-13 |
| Automatic quadratic bond, asymmetric term | 0 | 2.28e-13 |
| Automatic quadratic angle | 0 | 2.85e-14 |
| Automatic torsion_1 | 1.39e-16 | 6.67e-16 |
| Automatic zero-equilibrium Wilson | 6.40e-14 | 1.43e-13 |
| Synthetic torsion with nonzero 37° phase | 3.56e-15 | 5.33e-15 |

Raw source keys/rows are independently specified for the term reference; no ISLAND
selected coefficient list builds the LAMMPS oracle. The selected source IDs and
normalized coefficients are checked separately against that raw scan. For the
whole dichlorine graph, source atom/equivalence/increment rows and exact site/bond
coverage are independently checked. No msi2lmp auto-assignment claim is made.
The tiny static LAMMPS reference uses unit masses, which do not affect potential
energies/forces; ISLAND binding and workflow retain actual chlorine masses. No
LAMMPS trajectory or kinetic-energy comparison is claimed.

Finite-difference errors decrease toward roundoff as displacement shrinks from
1e-4 to 1e-6 Å; at the smallest displacement maxima are below 1.87e-7 in canonical
force units. The intermediate 1e-5 displacement is sometimes more accurate because
of subtraction roundoff. Production forces remain analytically differentiated by
OpenMM. Fresh and session results/identities agree for both whole-system frames.

### Remaining full-source work

The gap ledger identifies every still-blocked interaction and all 67 remaining
labels with source-row descriptions and next actions. Required work includes:

1. Missing source cross terms for the four organic/halogen fixtures above, or a
   separately justified and versioned physical model. A lower-order base alone
   cannot justify dropping those terms.
2. Independent chemical predicates and explicit charge-state conventions for
   metal/ionic labels, guanidinium/arginine/imidazole nitrogen variants, aromatic
   amines, phosphorus, sulfur oxidation states, zeolite environments and isotopes.
   Source comments are not complete graph-typing rules.
3. Resolving generic/specialized aliases (`h`, `h*`, `o`, `n`, `n+`, `ho2`, `s`,
   `o=`, `oe`) with independently typed examples. Neither generic catch-alls nor
   the old ho/ho2 discrepancy are silently converted into coverage.
4. Nonzero Wilson equilibrium signed-angle/permutation semantics. Empty
   torsion–torsion remains a named model exclusion, not a recovered parameter.
5. Wider whole-system evidence for mixed quadratic/quartic models and automatic
   torsion/Wilson use. Isolated equation checks establish a bounded numerical
   capability; they do not supply missing coupling data or universal applicability.

The low-level assignment's inherited `conventions` describes the original H3
source-record tier; v2's explicit `resolution_policy` and the model's full policy
record govern supplementation. Historical v1 payloads retain the original tier
without supplementation. Native consistency checks establish internal agreement,
not computational authenticity of arbitrary supplied records.

### Durable workflow result and verification record

The checked dichlorine bundle passed public preparation, save/load and the I3
workflow. Minimum: 8.17929899845 → 0 kJ/mol, **1 iteration, 4 evaluations**, maximum
and RMS force both 0, `force_converged`, independently fresh-verified. The exact
saved initialization was used for the uninterrupted comparison; it was not redrawn.
Four BAOAB steps retained `[0,1,2,3,4]` at 0.0004 ps. Split/uninterrupted differences
in coordinates, synchronized velocities, potential/kinetic/total energies,
components, forces and time were all 0 on this host; complete PCG64 state matched
exactly. This is not a cross-platform bitwise guarantee or an independent
validation of the thermostat.

The first process used five Contexts (two for minimization, three for the first
segment); the uninterrupted comparison used three, and the relocated child used
three. Evaluations were six uninterrupted versus eight over the two segments,
including the existing boundary verification calls. Completed status, frame reads
and completed no-op resume passed. Bundle hashes remained unchanged. Total harness
wall time was 472.00 s, dominated by reconstruction/semantic validation for the
source-bearing records, not four-step propagation. No speedup is claimed.

Identities:

- Prepared facade: `b1e3aeff559df6063fdfbcf112ee8cafeedc3611dfc29d99ecd5be519591767b`.
- Native parameter: `16cbb326b133b46e1a2de9320b8b481fbd110b789b53fa64b09da8d23418e8e3`.
- Evaluator model: `4bbf673154c4d5aa0971989c2ea415c2b8263e862a039e7a5a940a880946da7f`.

All new artifacts are under `../island-validation/phase4j2-declared/`:

- `acceptance/`: all 18 declared source/coverage results, including the initial
  harness attribute error and unavailable-Fourier executable failure.
- `corrected-checks/`: only affected dichlorine and isolated-term checks rerun
  with the corrected result attribute and separately compiled executable.
- `workflow/`: first attempt's valid paused checkpoint/minimum plus its wrapper
  failure (required evidence argument omitted from child invocation). Kept intact.
- `workflow-checked/`: identical declared inputs/settings in a new run; corrected
  child command, relocation and successful continuation. No tolerance, seed,
  force criterion or physical parameter change.
- `historical-readonly.json`, `j1-committed-hashes.json`: 12 historical J1 model
  records, two historical/J1 bundles and two completed historical workflows
  validated with OpenMM/RDKit/ParmEd/Foyer/SciPy imports blocked. All 252 files
  checked around reading remained unchanged; all 243 files listed in the original
  committed J1 artifact manifest matched its original hashes.

A separate J2 receipt, coverage ledger and exact interaction/type gap ledger are
committed in `docs/evidence/phase_4j2*.json`. No generated full parameter catalog,
FRC library, binary or trajectory is committed. The original J1 evidence files
are unchanged.

Executed verification:

- Ordinary suite: **1,436 passed, 10 skipped**, 1030.77 s. Existing optional
  archive/backend skips remain; no synthetic test substitutes for real evidence.
- Final focused fallback regressions: **9 passed** after adding the strict
  charge/parameter-policy mismatch guard. They cover selected-row precedence,
  wildcard conflicts, legal roles, equilibrium dependencies, missing coupling
  rejection, finite-difference forces, charge orientation, tampering, owned
  persistence and optional-import isolation. Existing corrected I3 rejection
  tests remain in the passing ordinary suite.
- Earlier historical PCFF expanded/model/Class II subset: **55 passed**.
- Ruff and `git diff --check` passed. `pip check` passed in the main validation
  environment and retained Foyer environment.
- New fallback example and historical expanded example passed. The historical
  full-source CLI returned **1**, as required for incomplete coverage.
- Actual LAMMPS checks, fresh/session comparisons and separate-process workflow
  ran as reported above. No new QM, unrelated force-field matrix or long trajectory.

Primary evaluation environment: Python 3.11.16, NumPy 2.4.6, RDKit 2026.3.6,
OpenMM 8.6.1 and SciPy 1.17.1; all force evaluations used Reference. Foyer is not a
PCFF runtime/reconstruction dependency. The independently built LAMMPS executable
is needed only for validation.

```sh
$PY scripts/summarize_pcff_fallbacks.py \
  --root ../island-validation/phase4j2-declared --source "$FRC" \
  --output /new/j2-summary
$PY -m pytest -q tests/test_pcff_fallbacks.py
$PY -m pytest -q
ruff check .
$PY -m pip check
```

**Bounded J2 operational acceptance passed. Full-source completion remains unmet.**
The next useful work is source-supported coupling coverage and independently
justified remaining chemical/charge-state rules; conservation or converter success
cannot replace that evidence.

The separate final offline inspection also passed for the **new J2 completed
workflow**, with OpenMM, SciPy, RDKit, ParmEd and Foyer imports blocked: status,
steps `[0,1,2,3,4]`, frame reading, completed no-op resume and all file hashes were
unchanged. The evidence summarizer returned 0 for the bounded J2 gate and **1**
with `--require-full-source`, without rerunning scientific calculations.

Source ingestion remains 4,359/4,359 numerical rows, plus the four semantic
catalogs, before and after J2. All 22 data section/namespace entries (including the
empty torsion–torsion section) are retained. The change is executable lookup and
assembly for four supplemental family/namespace paths, not an assertion that
every source row now has a justified chemical rule or complete physical model.

A bounded **actual msi2lmp control** also ran on checksum-verified retained ethene
`.car`/`.mdf` inputs, with a declaration published first. Command:

```sh
../island-validation/phase4h5-tools/msi2lmp reference -class II \
  -frc /explicit/pinned/pcff.frc -p 3 -nocenter
```

All 19 mass, coefficient, atom and interaction sections exactly matched the
retained raw converter result. The separately saved command/log/output/checksum
receipt is `msi2lmp-control/`. This checks the converter's ordinary Class II path;
it is not automatic typing, charge assignment or automatic fallback evidence.
`SearchAndFill.c`, pinned lines 143–217, also confirms version replacement per
ordered source key and replication of omitted symmetric coefficient halves.
ISLAND keeps equal-version coefficient conflicts explicit instead of file-order
selection. The original J1 angle-angle failure and correction are untouched.

Automatic nonbonded-equivalence supplementation beyond existing ordinary matches
was not added here; all declared cases that reach parameter assignment have
ordinary nonbonded coverage. Such a future missing nonbonded case must be diagnosed
as an unimplemented lookup path until its automatic mapping is checked, rather
than immediately called a missing source parameter. No claim of full positional
fallback coverage outside the explicitly implemented families is made.

## Review correction: declared typing scope and fail-closed acceptance

This correction is on the existing J2 branch, based on reviewed commit
`ddbdb6c4e814b281c59bcb0334055f8f641f5a51`; the merged base remains
`758449bb125ed706d2e982b34b4a2c2d8f27d56c`. Main and the remote feature tip were
verified before editing. The source pin, physical model, identities for valid
inputs, numerical tolerances, seeds, budgets and historical evidence are unchanged.

### Reproductions

The pinned real source and construction seed 2026 reproduced both unintended
assignments. `O=O` and `N#N` were incomplete under source-graph v1 but complete
under v2, with `o_1/o_1` and `nt/nt` respectively. Public preparation and fresh
Reference evaluation also succeeded (energies 2.0563544223281984 and
1.2116917048074658 kJ/mol). These are evidence of the software defect, not valid
chemical models. Their incorrect typing/model bytes are retained in the separate
`phase4j2-review-correction/` directory and must now fail native validation.

The second reproduction ran actual water preparation with `--cases water` and
replaced **only** `terms()` with the five retained successful term records from
committed J2 evidence. No dichlorine or new LAMMPS calculation ran. The old CLI
reported `assignment_and_term_gate=True` and returned 0. This is an explicitly
labelled software control, not additional scientific evidence.

### Corrections and contracts

`elemental_halogens` now affects only the declared neutral homonuclear single-bond
halogen branch. Oxygen double-bond and nitrogen triple-bond recognition again
require a carbon parent. O=O and N#N remain unresolved in both automatic and
explicit-type paths; supplied `o_1`/`nt` labels cannot bypass the chemical checks.
Carbonyl, nitrile, water and ClCl controls retain their behavior. No new profile
or signature is introduced: this enforces the already-declared v2 scope. Valid
records retain their identities; previously accepted invalid records are rejected
without repair or re-signing.

The acceptance and summary CLIs share offline checks in
`scripts/pcff_fallback_gates.py`. The assignment/term gate requires:

- Exactly one named dichlorine case and no duplicate case names.
- Successful typing, native charges and complete model, without execution errors,
  including errors after numerical comparison but before bundle publication.
- Both declared numerical comparisons and distinct retained reference geometries.
- A published bundle with the expected file set and checksums, validated native
  reconstruction, matching facade/model/assignment identities and a Cl–Cl graph.
- All five distinct declared `(family, fixture kind)` term checks, finite numerical
  errors, real-row selection evidence and the declared finite-difference checks.
- The unchanged-source result and an explicit local source matching its pin and
  experiment declaration.

Legacy comparison records were written only after successful original allclose
assertions. Their recorded absolute errors independently satisfy the declared
1e-5 absolute bounds; the reader requires that sufficient condition rather than
inventing missing per-component relative-error evidence. The original numerical
comparison still uses atol=1e-5 and rtol=2e-10. These checks establish internal
consistency of retained evidence, not computational authentication.

The summary resolves the corrected case by name, verifies its original case,
SMILES, source inventory and model identity, and recomputes the assignment gate.
Array position and a saved success boolean cannot establish acceptance. Workflow
acceptance remains separate and must match the validated parameter/model
identities, unchanged-input and RNG checks. Original harness failures remain
identified separately from successful corrections. Full-source completion remains
false; exploratory partial runs retain diagnostics and return nonzero for the
bounded gate if mandatory evidence is absent.

Summary inspection now takes an explicit source path (no historical-path search):

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY scripts/summarize_pcff_fallbacks.py \
  --root ../island-validation/phase4j2-declared --source "$FRC" \
  --output /new/review-summary
$PY -m pytest -q tests/test_pcff_review_typing.py tests/test_pcff_fallback_gates.py
```

The separate corrective receipt is `docs/evidence/phase_4j2_review_correction.json`.
It records executed checks, read-only artifact hashes and identities, post-fix
real-source boundary checks, and original versus corrected CLI outcomes. Original
J1/J2 receipts, failures, bundles and checkpoints are not overwritten. This
correction does not complete the remaining 67 unresolved/unimplemented source
labels or the documented charge/cross-term gaps, and does not start J3.

Executed corrective software checks include **33 new focused tests** and **98
historical PCFF/expanded-model/prepared-workflow regressions**, including the
corrected I3 rejected-attempt tests. The focused typing tests also ran against the
reviewed pre-fix module in an isolated process: both O/N negatives failed, while
all four positive controls passed. The post-fix water-only CLI control returns 1
and records `Missing mandatory dichlorine case`; it still retains the actual water
preparation and the labelled reused term-check records.

Read-only native inspection passed for the retained J2 dichlorine bundle and
completed workflow (5 frames), the J1 corrected workflow (5 frames), and the
historical I3 PCFF workflow (11 frames), including completed no-op resume. OpenMM,
SciPy, RDKit, ParmEd and Foyer imports were blocked during inspection. Both the
historical CHO and J1 bundles retain their native/facade identities. All **228
files** in the original J2 experiment tree and all **8 original committed J1/J2
evidence files** remain byte-identical.

The real-source nitrile control `CC#N` retains complete typing and accepts the
same checked explicit labels, but does **not** have complete parameter coverage.
Public preparation reports the existing ambiguous angle candidates and
missing cross terms. This is an unchanged full-source coverage limitation, not a
newly successful nitrile model. The initial corrective verification script's
assumption of full nitrile preparation was wrong; its failure log is preserved,
and the continuation checks this outcome against the reviewed implementation.
No types, parameters or tolerances are changed to make it pass.

Bounded real-source post-fix checks completed: formaldehyde, water and dichlorine
public preparation and fresh Reference evaluation succeed. All three model
identities equal their original J2 identities; dichlorine also preserves its
original facade identity `b1e3aeff559df6063fdfbcf112ee8cafeedc3611dfc29d99ecd5be519591767b`.
The nitrile typing identity and exact public incomplete-model diagnostic match
the reviewed code. O=O/N#N public preparation and checked explicit typing reject
as unresolved, and their saved invalid typing records reject as contradictory.
No new QM, LAMMPS comparison or dynamics was executed in this correction.

The complete ordinary suite ran once after the corrections: **1,469 passed,
10 skipped** in 1012.72 seconds. Ruff, `git diff --check`, and pip checks in the
main validation and retained Foyer environments passed. The corrected summary
returns 0 for the bounded J2 gate and 1 with `--require-full-source`. All requested
corrective gates are satisfied; the explicitly reported full-source/nitrile
parameter-coverage limitations remain. Readiness flags remain
`production_validated=False` and `simulation_readiness="not_established"`.
