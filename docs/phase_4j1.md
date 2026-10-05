# Phase 4J1 — PCFF full-source coverage expansion (partial)

## Base and honest completion boundary

Branch `codex/phase-4j1-pcff-full-source-coverage` starts at reviewed I3 correction
`787deb0b7c93e3c1ab08e29bc5a7325111ec863c`. At fetch, `origin/main` was
`8d8a11aa35d1f3491c5c9455059100ecdf9fbc67` (I2). Corrected I3 is an **unmerged
prerequisite**. No prerequisite merge, main edit, force push, or historical
scientific-data rewrite was performed.

**The complete FRC target is not finished.** This is executable progress beyond
an inventory: a separately versioned graph profile, native charged-component
bookkeeping, expanded parameter/model records, Wilson energy/analytic forces,
and real whole-system reference comparisons. The full-source gate remains false.
The old saturated acyclic CHO profile remains readable and numerically unchanged;
it is not the permanent coverage target.

The machine-readable [coverage ledger](evidence/phase_4j1_coverage.json) contains
all 133 source types and every parameter namespace, not just successful fixtures.
Its conservative full-source completion counters are separate from **bounded
fixture coverage**. A passing fixture does not validate every use of its types.

`production_validated=False`, `simulation_readiness="not_established"` throughout.

## Source inventory and interpretation

Source remains external:

- [LAMMPS](https://github.com/lammps/lammps), revision
  `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- `tools/msi2lmp/frc_files/pcff.frc`, SHA256
  `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
- Installed here: `../island-validation/phase4h1-sources/lammps/pcff.frc`.
- Author notices and data redistribution handling remain as documented in H1.
  No full FRC or external implementation code is bundled.

An independent line scanner confirmed **133 distinct atom types, 134 ordinary
equivalence rows, 108 automatic-equivalence rows, and 564 bond-increment rows**.
The numerical inventory is:

| Namespace | Family | Rows |
|---|---|---:|
| cff91_auto | quadratic_bond | 628 |
| cff91 | quartic_bond | 127 |
| cff91_auto | quadratic_angle | 330 |
| cff91 | quartic_angle | 303 |
| cff91_auto | torsion_1 | 216 |
| cff91 | torsion_3 | 492 |
| cff91 / cff91_auto | wilson_out_of_plane | 71 / 12 |
| cff91 | nonbond(9-6) | 94 |
| cff91 | bond-bond | 245 |
| cff91 | bond-bond_1_3 | 64 |
| cff91 | bond-angle | 234 |
| cff91 | angle-angle | 270 |
| cff91 | end_bond-torsion_3 | 291 |
| cff91 | middle_bond-torsion_3 | 351 |
| cff91 | angle-torsion_3 | 303 |
| cff91 | angle-angle-torsion_1 | 328 |
| cff91 | torsion-torsion_1 | 0 |

`inspect_pcff_full_source()` preserves namespaces, declarations, raw lines,
versions, reference IDs, duplicate keys, signed/zero coefficients, dimensions
and conversion factors. This additive catalog does not rewrite H1/H3 catalogs.
The empty torsion–torsion section remains uninterpreted; no coefficients or
universal zero have been invented.

`resolve_pcff_source_record()` queries one explicit family/namespace. Ordinary
lookup uses exact and family-specific equivalence. Automatic queries use the
bond, angle-end/apex, torsion-end/center and OOP-end/center columns independently.
Exact matches precede wildcard candidates; numerically conflicting wildcard
patterns remain ambiguous. Highest versions apply per source key. The query API
**does not silently select a quadratic or torsion fallback for a model**. These
forms still need explicit operational-model assembly and reference verification.

This restriction has concrete evidence: the pinned msi2lmp README excludes
bond increments and automatic equivalence. LUNAR's pinned charge fallback branch
is unreachable after an exhaustive earlier branch (H1 audit). Therefore neither
is presented as an automatic-charge oracle. Explicit automatic record queries
have numerical/role tests, but not full fallback physical-model acceptance.

## Expanded graph typing and charges

Opt in using `profile="island_pcff_source_graph_v1"`. Old default
`island_pcff_acyclic_cho_v1` is unchanged. The new profile recognizes graph-local
predicates, with explicit degree/bond-order/aromatic/formal-charge guards and
shortest-cycle membership. It never reads coordinates, repeat indices, residue
names, polymer names, or global ID ordering. Component-wide unresolved evidence
blocks successful typing of that component.

Rules include:

- Explicit-H sp3 C with specific three-/four-member-ring types before H-count
  classes; separate acetal `co/coh` predicates before generic sp3 cases.
- Aromatic `cp`, five-member `c5/cs`; nonaromatic C=C endpoint/adjacent/interior
  classes; triple-bond `ct`; carbonyl environment classes.
- Ordinary alcohol `oh/ho` retained; ether, small-ring O, aromatic O, carbonyl,
  ester and carbonate O distinguished. Water is recognized separately.
- Selected amine, aromatic N, nitrile, small-ring N, sulfur, siloxane,
  covalently bound halogen, and isolated noble-gas environments.
- Carboxylate `c-/o-` uses authoritative -1 formal charge/resonance connectivity;
  tetrahedral +1 ammonium uses `n4` with an explicit `n+` manual alias.

Each selected row includes environment evidence, source atom-row identity, and
rule/profile identity. Unsupported isotopic/radical/valence or
unresolved chemical environments are not absorbed by a generic fallback.
Assigned stereochemistry is preserved in binding metadata, not stripped. No
stereochemistry-specific charge rule is inferred. A 0.02-dalton standard-weight
guard additionally detects isotope masses when older construction adapters omit
isotope metadata; masses are never modified.

Evidence is the exact FRC plus the independently inspected
[LUNAR PCFF typer](https://github.com/CMMRLab/LUNAR/blob/67dabeda9e6bd3cc8f968aa0c6a88730886138ef/src/atom_typing/typing/PCFF.py).
No GPL implementation was copied or imported. Its assumed/fallback assignments,
formula/name-based special cases and apparent rule bugs are not treated as
chemical ground truth. Intentional differences include the historical oh/ho
ordinary-alcohol convention, authoritative bond-order/formal-charge tests, and
rejecting ambiguous ordinary amide/aromatic-amine cases rather than guessing
between `n/n_2` or `nn/nb`. The complete remaining labels appear individually in
the ledger. This phase inspected but did not execute the external LUNAR typer. Code availability is not independent typing validation.

`assign_pcff_source_types()` supplies explicit stable-ID labels with provenance,
checked against the same implemented chemical constraints. Recognized aliases
are explicit; unavailable constraints still block this route. It is **not** yet
an unrestricted manual route for all 133 labels. Automatic and explicit routes
use the same charge, Class II assignment and model representations; their
provenance identities remain distinct even when numerical content agrees.

Charges start at zero, process each authoritative bond once, retain both source
endpoint increments, and swap endpoint roles on reversal. Selection remains
exact then ordinary **bond** equivalence. Complete connected-component totals
must match their actual formal charges to `1e-12 e`, including charged components.
No projection, normalization, fitted charges or alternate backend is used.
The two tested amine fixtures lack `na–hn2` increments even after documented
ordinary equivalence; neutrality of the *partial* vector does not make them pass.

## Expanded Class II model and Wilson implementation

New records:

- `island_pcff_source_typing_v1`
- `island_pcff_source_charges_v1`
- `island_pcff_source_class2_assignment_v1`
- `island_pcff_source_model_v1`

They use the established owned/checksummed data-only containers and public
persistence methods. Validators reconstruct decisions, graphs, contributions,
parameter selection, permutations, dependencies and model definition. New
interpretation profile: `island_lammps_pcff_source_graph_v1`. Old schema payloads,
signatures and fingerprints retain the old derivation path.

Authoritative bonds generate all inventories. Wilson terms use one neighbor
triple at each degree-three center, preserving center position 2 and trying all
six arm permutations. Degree-four centers retain all angle–angle couplings.
Missing Wilson or cross-term coefficients block completion. The model currently
rejects nonzero Wilson equilibrium angles because odd-permutation sign behavior
has not been independently established for that source scope.

The executed [pinned LAMMPS improper code](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/improper_class2.cpp)
uses the mean of three signed out-of-plane arcsines, **not a single dihedral**.
For outward unit vectors u,v,w from the center, let t=u·(v×w). The mean is
`[asin(t/|v×w|)+asin(t/|w×u|)+asin(t/|u×v|)]/3`, with energy `K*(mean-chi0)^2`.
OpenMM differentiates this expression analytically. Coefficients are kJ/mol
per radian squared; angular geometry is dimensionless. No length coefficient
conversion enters this term. `SourceClass2Term` provides its independent NumPy
energy/finite-difference check; historical `Class2Term` retains its old domain.

Existing quartic/cross-term equations, torsion phase convention, kcal→kJ factor
4.184, angstrom-based expressions, sixth-power 9–6 mixing, and exact Coulomb
constant are unchanged. New evaluator settings/model identities identify the
expanded profile. Old component names, ordering and identities are unchanged;
expanded models add `wilson_out_of_plane` as a separate component.

The non-cp BB13 zero remains an explicitly named **converter compatibility
policy**, never a source row. The actual pinned GetParameters.c conditional tests
supplied cp types, not a CHO element restriction. J1 independently exercised this
branch for rings, ethene and sulfur, and cp lookup for benzene. It is not asserted
to be the universal interpretation of every PCFF source. Other missing terms
are not filled with zero. Special-pair weights remain mandatory caller choices.

## Real reference execution and a retained reference failure

Before numerical execution, `phase4j1-declared/declaration.json` fixed all 14
cases, seed 2026, unminimized coordinates, the `0.015*sin(index+0.3)` angstrom
perturbation, LJ/Coulomb `(0,0,1)`, 180-second subprocess budgets, and existing
H5 tolerances: energy/force atol `1e-5` in canonical kJ/mol and
kJ/(mol*angstrom), rtol `2e-10`. No QM, optimization or dynamics ran.

Source/type expectations were authored independently. Native charges were
checked by a separate Decimal source scanner and endpoint summation. The
reference uses msi2lmp's independently regenerated inventories/parameters and
LAMMPS, not ISLAND's compiled interaction list. Raw converter files are retained.
H5's independent angle–angle equilibrium-role correction is preserved.

| Case | Actual outcome |
|---|---|
| Cyclohexane | Typing, charge, model, two whole-system numerical comparisons pass |
| Cyclic ether (1,4-dioxane) | Pass |
| Epoxide | Pass |
| Benzene | Pass, including nonplanar perturbed Wilson forces |
| Ethene | Pass, including Wilson terms |
| Methanethiol | Pass, source sulfur types |
| Acetate (-1) | Pass, native molecular charge retained |
| Acetal | 19 unresolved model interactions; retained diagnostic |
| Acetone | 6 unresolved model interactions; retained diagnostic |
| Chloromethane | 15 unresolved model interactions; retained diagnostic |
| Ester | 39 unresolved model interactions; retained diagnostic |
| Methylamine / amino-alcohol reaction-product fixture | Missing na–hn2 charge rows |
| Tetramethylammonium (+1) | Original converter comparison failed; see below |

Across the seven original passing cases, maximum energy/component discrepancy
was `2.8422e-13 kJ/mol`; maximum force discrepancy was `1.59035e-6
kJ/(mol*angstrom)` for the near-planar initial benzene. Its perturbed geometry
agreed below `9.2e-13`. These are implementation-fidelity measurements, not
chemical accuracy thresholds.

Tetramethylammonium's original reference failed by `0.09640496198 kJ/mol` and
`1.634351616 kJ/(mol*angstrom)`. Investigation found another converter limitation:
`find_angleangle_data()` calls `find_match()`, which permits full reversal and
can move the central atom. For three identical C arms on n+, the pinned source
row **4143** gives -1.5155 kcal/mol for each coupling. The converter selected
-1.5155, -4.2781, -4.2781. The latter values correspond to another central role.
ISLAND's center-fixed selection was not changed.

A **separate declared diagnostic**, `phase4j1-aa-diagnostic`, changed only reference
AA coefficient selection. It reconstructed center-fixed rows from its own raw
source scan, converter inventory and equivalence table, retaining the H5
independent equilibrium correction. Both declared geometries then agreed below
`3.9791e-13 kJ/mol` and `1.7054e-12 kJ/(mol*angstrom)`. The original failed report,
source, coordinates and coefficients remain untouched. This is explicit reference
correction evidence, not a silently converted original pass.

## Public API and commands

```python
from island.forcefields import (
    ForceFieldRequest, PCFFOptions, prepare_forcefield, create_evaluator,
)
request = ForceFieldRequest("pcff", PCFFOptions(
    source_path="/external/pcff.frc", lj=(0, 0, 1), coulomb=(0, 0, 1),
    typing_profile="island_pcff_source_graph_v1",
))
prepared = prepare_forcefield(system, request)
with create_evaluator(system, prepared).open_session() as session:
    result = session.evaluate()
```

The facade fails on incomplete native models. Expanded records also use existing
`adopt_forcefield`, prepared-bundle save/load, fresh/session evaluation, and bound
minimization/dynamics contracts without a new numerical engine. No new simulation
acceptance or workflow schema is introduced. The validators establish internal
consistency, not authenticity of arbitrary caller-supplied records.

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY -m island.forcefields.pcff.coverage_cli --source "$FRC" --output /new/ledger.json
$PY -m island.forcefields.pcff.coverage_cli --source "$FRC" --require-full-source
$PY -m island.forcefields.pcff.coverage_cli --source "$FRC" --family quadratic_angle --namespace cff91_auto --labels hc c3 hc
$PY -m island.forcefields.pcff.coverage_cli --source "$FRC" --smiles 'c1ccccc1' --output /new/assignment
$PY examples/pcff_expanded.py --source "$FRC"
$PY examples/pcff_expanded.py --source "$FRC" --explicit
$PY scripts/validate_pcff_expanded.py --source "$FRC" \
  --converter ../island-validation/phase4h5-tools/msi2lmp \
  --lammps ../island-validation/phase4h4-upstream/build/lmp --output /new/acceptance
$PY scripts/diagnose_pcff_angle_angle.py --source "$FRC" \
  --case ../island-validation/phase4j1-declared/tetramethylammonium \
  --lammps ../island-validation/phase4h4-upstream/build/lmp --output /new/diagnostic
```

For explicit labels, add `--types file.json --provenance 'justification'` to the
assignment CLI. Format is `{"1":"c2","2":"hc",...}`, exact canonical decimal
integer IDs. Inspect the constructed system's IDs; never infer production labels
from ID order. Invalid mappings and source mismatches fail with typed errors.

Inspection exits 0 for a structurally verified inventory; `--require-full-source`
returns nonzero while full coverage is unmet. Assignment reports raw parameter
coverage, so policy-derived BB13 zeros do not falsely count as source rows.
The acceptance CLI's default gate is explicitly the six predeclared numerical
milestone cases, not full-source completion. Optional exploratory failures remain
in the report. Add `--require-full-source` to keep that full target gate nonzero.

Only graph construction requires RDKit. Offline graph typing/charge/record
validation uses core ISLAND/NumPy. OpenMM is lazy and only needed for evaluation.
No automatic downloads, Foyer, AmberTools, QM, SciPy or embedding occurs during
typing, charge assignment or record loading. Bundle loading uses explicit pinned
source resolution and performs no new preparation.

## Remaining work against the full inventory

The ledger is a continuing implementation driver, not a substitute for:

1. Independently justified predicates/fixtures for remaining protein-specific,
   ionic/charged resonance, metal, zeolite and isotopic types. Some source comments
   are ambiguous or inconsistent with the pinned upstream typer; they do not
   justify guessing a graph rule.
2. A verified automatic-charge path where source rows exist. For na–hn2 even the
   pinned automatic bond-increment labels remain na/hn2; no supported row was
   established. Changing hn2 to hn to obtain coverage is not allowed.
3. Integrating separately interpreted automatic quadratic/torsion/OOP families
   into an explicitly defined complete model, with independently justified
   cross-term applicability and equilibrium dependencies. Queries alone are not
   complete fallback support. Wildcard conflicts stay unresolved.
4. Source-missing cross terms for specific acetal/carbonyl/ester/halogen fixtures,
   listed by exact interaction IDs and types in retained parameter records.
   Lower-order base terms do not justify zeroing missing cross terms.
5. Nonzero-Wilson-equilibrium permutation semantics and additional independent
   reference capability. Torsion–torsion remains an explicit operational exclusion,
   not an invented source parameter or universal physical zero.
6. Extending numerical verification beyond the executed cases. Full typing,
   charge coverage, parameter coverage, evaluator availability and independent
   verification have separate counters. None implies scientific readiness.

No periodic simulation, packing, automated curing, or network dynamics capability
is claimed from these chemical coverage increments.

## Verification record

Executed checks:

- Ordinary suite: **1,413 passed, 10 skipped**, 490.15 seconds. Existing skips
  remain optional archive/backend checks. An earlier run caught an incomplete
  synthetic automatic-equivalence fixture; the fixture was corrected and the
  ordinary suite rerun successfully.
- Expanded regressions: **25 passed**, also included in the final ordinary run.
  These cover independently authored chemical predicates, rejection of ambiguous
  groups/isotope masses, expanded facade/bundle integrity, and the center-fixed
  converter AA regression. Synthetic tests are explicitly labelled;
  they are not substituted for real source evidence.
- Ruff passed. `pip check` passed in both the main validation and Foyer environments.
- Automatic and explicit examples passed; the historical PCFF PE DP3 single-point
  example passed. The full-source inspection CLI correctly returned **exit 1**.
- Actual declared source/LAMMPS matrix: seven original numerical passes, four
  parameter-incomplete cases, two charge-incomplete cases and one retained
  reference mismatch. The separate reference diagnostic passed for that last case.
- Fresh facade preparation, save, relocation, separate-process bundle loading,
  native identity preservation, and fresh/session forces/components passed for
  cyclohexane (`phase4j1-bundle-verified`). An initial harness used the wrong
  EvaluationResult attribute name; its output directory is retained separately
  as `phase4j1-bundle`. The corrected harness is reproducible with
  `scripts/validate_pcff_expanded_bundle.py --source FILE --output NEW_DIR`.
  That committed CLI also passed independently in `phase4j1-bundle-cli`; its
  exact command, receipt and artifact hashes are included in the evidence file.
- Read-only historical validation passed for six H3 assignments, the I2 PCFF
  bundle and the completed I3 PCFF workflow (11 frames); all **25 checked files**
  kept their original hashes. No historical preparation or trajectory was rerun.

Environment: Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1, RDKit 2026.3.6,
SciPy 1.17.1. Source and executable hashes, commands, per-case timings and
input/output hashes are in [the evidence receipt](evidence/phase_4j1.json).
The actual LAMMPS executable and msi2lmp are the retained pinned H4/H5 builds,
not an executable inferred from a documentation version. Full source and
scientific records remain external under `../island-validation/phase4j1-*`.

### Separate coverage measures

| Measure | Executed/implemented extent; full-source gaps remain |
|---|---|
| Ingestion | 133/133 type rows, all 22 data-family/namespace sections; 4,359 numerical interaction rows preserved |
| Chemical recognition | Predicates implemented for 66 labels; 23 distinct labels exercised by the 14 declared source fixtures; 38 further labels remain implementation work and 29 have unresolved specialized rule evidence |
| Native charges | 12/14 declared fixtures complete, 19 labels exercised; automatic increment fallback is not enabled |
| Parameter/model resolution | 8/14 declared fixtures have a complete named operational model; 13 labels exercised; raw BB13 policy zeros remain distinct from source coverage |
| Executable energy/force | Expanded degree-three Wilson and existing Class II terms; automatic lower-order forms are queryable but not assembled into complete executable models |
| Independent numerics | Seven original cases, plus one separate corrected-reference diagnostic; 12 distinct labels in original passes, with n4 added by the diagnostic |

[Exact unresolved interaction IDs/types/dependencies](evidence/phase_4j1_gaps.json)
are provided for the four parameter-incomplete fixtures and two charge gaps.
The full FRC is **not** implemented-and-verified. Nothing in the declared milestone
or software test count promotes that full-source gate to complete.

## Corrective follow-up: shared offline PCFF evaluation identity

This correction remains on the J1 branch based on reviewed
`2f7a89b797ebc804720acc30561a414d98741092`, with corrected I3 still an unmerged
dependency. It does not complete the full-source coverage work above.
Original J1 evidence files and retained scientific records are unchanged;
[the corrective receipt](evidence/phase_4j1_identity_correction.json) is separate.

### Reproduction and cause

Two regressions were added and executed before changing production code. The
expanded synthetic eight-site model produced a valid single point, then a
force-converged, independently verified minimum (16 iterations, 21 evaluations
for this fixture). Public workflow startup nevertheless raised
`Minimum diagnostic model mismatch`; the committed manifest remained `starting`.
The offline workflow fingerprint also differed directly from the evaluator's
fingerprint. Both tests failed on the reviewed implementation.

The evaluator selected expanded settings for `island_pcff_source_model_v1`, while
workflow inspection always hashed the historical CHO `SETTINGS`. This difference
applied even to an expanded model with no active Wilson terms. A model's schema
and validated interpretation determine these identities, not its active force list.

### Correction and compatibility

`island.evaluation.pcff_identity.pcff_evaluation_identity(specification,
system=system)` validates the native model and returns a new settings dictionary,
its parameter fingerprint and its model fingerprint. Both the evaluator and the
prepared workflow now call this single offline derivation. The historical import
`island.evaluation.pcff.SETTINGS` remains available.

The formulas and payloads are unchanged:

- Parameter fingerprint: existing identity of the complete model specification.
- Model fingerprint: existing evaluation fingerprint function applied to
  `{"specification": parameter_fingerprint, "settings": settings}`.
- Historical CHO settings remain unchanged. Expanded settings retain
  `implementation="island_pcff_source_singlepoint_v1"` and
  `compatibility_profile="island_lammps_pcff_source_graph_v1"`.

Native integrity and optional system compatibility are checked before deriving
identities. The helper does not import OpenMM, construct a System/Context, type
atoms, evaluate forces or run scientific preparation. It uses the existing offline
PCFF validation dependencies (Python/NumPy and an explicitly resolved pinned FRC).
Coordinate replacement is permitted by the native chemical binding contract.
Malformed models and rechecksummed contradictory model/settings records remain
rejected; no alternate fingerprint is accepted as a fallback.

Numerical expressions, components, charges, units, pair policies, source and
interpretation pins, evaluation identities, bundle/workflow/checkpoint schemas,
verification tolerances and I3 failed-attempt validation are unchanged.

### New bounded workflow experiment

Reproduction command (use a new output directory):

```bash
../island-validation/phase4e2-env/bin/python scripts/validate_pcff_workflow_identity.py \
  --bundle ../island-validation/phase4j1-bundle-cli/relocated \
  --source ../island-validation/phase4h1-sources/lammps/pcff.frc \
  --output ../island-validation/phase4j1-identity-correction/workflow
```

The command checks the retained bundle against its original committed J1 hashes
and publishes a declaration before propagation. It reuses the original cyclohexane
input coordinates (seed 2026), native charges, LJ/Coulomb `(0,0,1)` and source/model
identities. No molecule is rebuilt or reparameterized. This is a separate new
minimization/workflow integration experiment, not a replacement for historical
single-point acceptance or an independent force-model comparison.

Declared budgets: force tolerance 0.1 kJ/(mol*angstrom), 5,000 minimization
iterations and 10,000 evaluations; unchanged default line search and verification
criteria. BAOAB uses 300 K, 5/ps, 0.1 fs, velocity seed 78123, thermostat seed 99181,
and four steps in two two-step segments, each bounded to four evaluations and
three frames. A six-evaluation/five-frame uninterrupted comparison uses exactly
the saved minimized system and initialized velocities. Split comparisons retain
atol 1e-10 and rtol 1e-12 in canonical units, with exact RNG-state equality.

Actual result: **passed**, 47 minimization iterations and 52 evaluations;
energy 14.9015859878 → -44.7522737438 kJ/mol, final maximum force
0.0765767455 and RMS force 0.0470122279 kJ/(mol*angstrom), with
`force_converged` and independent fresh verification. Both two-step segments
completed (eight evaluations total), yielding retained steps `[0,1,2,3,4]`
and 0.0004 ps. All compared coordinates, velocities, potential/kinetic/total
energies, forces, components and times had zero observed split differences;
complete RNG state matched exactly. This observation is not a cross-platform
bitwise guarantee.

Context counts: five for start (two for minimization and three for the first
segment), three in the separate continuation process, and three for the
uninterrupted comparison (six evaluations). Completed status, frame reading and
no-op resume passed. Native/facade/model identities survived relocation and input
bundle hashes were unchanged. End-to-end wall time was 711.93 seconds, including
native reconstruction and repeated semantic inspection; this is not an optimizer
or force-kernel benchmark. Artifact paths, exact command, source hash, process IDs,
configuration, fingerprints and all artifact hashes are in the corrective receipt.

Three retained bundles (historical I2 PCFF and the J1 verified/CLI expanded bundles)
and the completed historical I3 PCFF workflow were revalidated read-only with
OpenMM, RDKit, ParmEd, Foyer and SciPy imports blocked. Historical and expanded
parameter/model fingerprints matched the originally stored evaluation records.
The historical workflow retained all 11 frames and completed no-op resume. All
24 checked files stayed unchanged. An initial receipt-check harness omitted the
existing tagged-JSON decoder; its error log is retained, and the corrected
read-only check passed. No preparation or scientific execution was repeated for
these retained-record checks.

The original two regressions now pass. The 14 focused tests additionally cover
active Wilson terms, no-Wilson expanded models, frozen historical settings,
owned settings, compatible coordinates, child-process continuation, unavailable
optional dependencies, contradictory rechecksummed records and legitimate
minimization/dynamics failures. All 50 prepared-workflow tests, including corrected
I3 rejected-attempt binding, passed. The final ordinary suite passed: **1,427 passed, 10 skipped**, in 987.90 seconds.
Ruff passed; `pip check` passed in both the main validation and separate Foyer environments.
Existing optional archive/backend skips remain unchanged. No required gate for
this bounded identity correction was unavailable.

No new independent LAMMPS matrix was run: the force model and numerical kernels
are unchanged. This correction establishes identity consistency and workflow
integration for the checked case. The broader full-source typing, charge,
fallback-model and interpretation gaps remain open, with the same conservative
readiness flags.


Verification commands:

```bash
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_pcff_workflow_identity.py
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_prepared_workflow.py
../island-validation/phase4e2-env/bin/python -m pytest -q
ruff check .
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```
