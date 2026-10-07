# Phase 4J15 — bounded urethane and benzenoid-amine source families

## Delivery status

Branch: `codex/phase-4j15-pcff-urethane-aromatic-amines`.
Base: `4fecd57d6fdd2e8eec8472627a0d9ef182a63161` (J14 squash merge, PR47).
After `git fetch origin`, the complete base tree compared equal to reviewed
`92e36427038c7174228d1537c46a444cad358929`. There is no unmerged prerequisite.
Main was neither modified nor merged.

**Partial completion.** Opt-in `island_pcff_source_graph_v6` closes the declared
neutral carbamate/urethane and benzenoid-amine hydrogen-family charge gaps.
Aniline has a complete independently verified executable model. None of the
four declared polymer cases has all required source couplings. Consequently,
**the executable polymer tranche and full-source gates remain false**. No
incomplete polymer was minimized, evaluated, bundled or propagated.

`production_validated=False`; `simulation_readiness="not_established"`.
A verified force model and local minimum do not establish scientific suitability.

## Frozen experiments and preservation

The [20-case declaration](evidence/phase_4j15_declaration.json), construction
seed 2026, source pin, reference methods, tolerances and budgets were persisted
before production edits. The external `before/` records reproduce v5 on the
complete matrix. The new [expectation contract](evidence/phase_4j15_expectations.json)
binds negatives to those original decisions and binds model diagnostics to the
independent raw-source search performed before final acceptance. An expected
missing model is a preserved failure, never an executable acceptance success.

The [receipt](evidence/phase_4j15.json) separates v5/v6 typing, charges, strict-v3
models, compatibility-v4 models and independent numerical results. Its artifact
index identifies complete external source queries, oriented charges, raw converter
files, input coordinates, declarations, failed model diagnostics and checksums.
Original J1–J14 evidence and failed attempts were not rewritten. An initial J15
inspection predating correction of three citation-only atom line numbers is also
retained; final acceptance uses lines 83/149/150. These were new-profile citations,
not changes to a historical profile, source row or numerical rule.

## Chemical evidence and limits

The sole runtime source is [LAMMPS pcff.frc](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc),
revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
The full library remains external. [Chemical evidence](evidence/phase_4j15_chemical_evidence.json)
retains the relevant atom descriptions, versions/references, ordinary and
position-specific automatic equivalence rows and external implementation hashes.

| Family | Source atom lines | Ordinary / automatic equivalence lines |
|---|---|---|
| Carbamate carbonyl `c_2` | 83, v2.1/ref8 | 225 / 351 |
| Carbonyl `o_1`, ester `o_2` | 149,150, v2.1/ref8 | 291,292 / 414,415 |
| Urethane `n_2` | 129, v1/ref1 | 271 / 394 |
| Distinct amino-H `hn2` | 106, v2.1/ref8 | 248 / 373 |
| Aromatic amine `nn` | 137, v1/ref1 | 279 / 402 |
| N-bound `hn`, broad `h*` | 105,99, v1/ref1 | `hn` maps to `h*`: 247 / 372 |
| Aromatic alias `nb` | not selected by this update | 273 / 396 |

`n_2` is explicitly the urethane N; `c_2` covers carbamate/urea carbonyl C,
with distinct `o_1`/`o_2`. Ref7 identifies the September 1993 urethane refinement;
its oriented `hn2–n_2` row belongs to that refined family. This joint chemical
and family evidence motivates the bounded N–H choice. It does **not** imply
that every `hn2` environment is a urethane, or authorize globally replacing
`hn`, `hn2` or `h*`. Aromatic `nn` belongs to the independently described ref1
amine family, whose attached `hn` reaches `h*` through existing equivalence.

Independent implementation evidence was inspected, not copied:

* [LUNAR PCFF.py](https://github.com/CMMRLab/LUNAR/blob/67dabeda9e6bd3cc8f968aa0c6a88730886138ef/src/atom_typing/typing/PCFF.py),
  revision `67dabeda9e6bd3cc8f968aa0c6a88730886138ef`, lines 861–864 recognizes
  `n_2` for the C3H7NO2 urethane example. V6 uses a bounded local motif, not that
  whole-molecule formula. Lines 905–911 and the separately inspected
  `typing_functions.is_aromatic_amine_nitrogen` corroborate non-ring N attached
  directly to aromatic C. V6 additionally restricts the ring to a six-carbon
  benzenoid component and checks state/valence.
* LUNAR's preceding sp3-amine H rule assigns `hn2`; its generic N–H rule assigns
  `hn`, including urethane H. **V6 deliberately differs for urethane H**, using
  the source's refined urethane family. For aromatic N–H the generic `hn`
  convention agrees, unlike ISLAND's inherited `hn2`. Exact LUNAR parity is
  not claimed. Its GPL implementation was not bundled or translated.
* [pysimm PCFF.py](https://github.com/polysimtools/pysimm/blob/fb33814128189a99d3c8d7c4eddac2dfe44262fe/pysimm/forcefield/pcff.py),
  revision `fb33814128189a99d3c8d7c4eddac2dfe44262fe`, lines 98,117–120,162–169,
  independently corroborates generic N-bound `h*`, carbamate `c_2`, ester `o_2`
  and carbonyl `o_1`. Its generic N rule does not establish the specialized
  urethane N/H convention. Parameter availability and charge neutrality alone
  were not used as chemical authority.

### Exact opt-in rules

`urethane_domains.py` inherits v5 and changes only the following verified H
families. Parent N, carbonyl C/O and alkoxy O labels are checked before an override:

* **Urethane:** neutral, standard-mass, nonisotopic, nonradical, acyclic N with
  three single bonds; exactly one `R–O–C(=O)–N` carbonyl neighbor. Every other
  N substituent is explicit H or saturated C with four single bonds. Carbonyl C
  and ester O are acyclic; the ester O's other neighbor is saturated C. Existing
  parent labels must be `n_2/c_2/o_1/o_2`. Only attached H becomes `hn2`.
  Tertiary carbamate retains its already charge-complete parent assignment.
* **Benzenoid amine:** neutral, acyclic N with three single bonds, one or two
  explicit H and one or two C neighbors. At least one C belongs to an entire
  six-carbon aromatic cycle; other C is saturated or in a qualifying benzenoid
  cycle. Existing N must be `nn`; attached H becomes `hn`. Fused/heteroaromatic
  rings and acylated aromatic N do not receive this override.
* Each tested atom has matching element/atomic number, formal charge, isotope,
  radical, mass, bond order and aromatic state. No coordinates, site numbering,
  repeat membership or molecule names select chemistry. Stable IDs and original
  topology/stereo/provenance are preserved; explicit supplied labels undergo the
  same predicates. Signed evidence includes parent IDs and rule precedence.

Typing v6 is **not** resolver v6. Guarded source resolver v3 and opt-in converter
compatibility v4 are unchanged, including every historical identity. V4's existing
non-cp BB13 initialization zeros retain converter provenance. They do not authorize
zero BB/BA/AA or torsional couplings, and do not make strict source coverage pass.
Nonzero Wilson equilibrium and missing torsion–torsion semantics remain unresolved.

## Charges, interactions and actual outcomes

Independently authored `pcff_j15_reference.py` supplies expected labels to the
existing separate Decimal FRC reader. The reader independently selects increment
rows, swaps endpoint contributions on reversal and sums each connected component.
Native per-site vectors are checked at **1e-12 e**, without normalization.

| Oriented pair | Source line | Endpoint contributions (e) |
|---|---:|---|
| `c–n_2` | 512, ref7 | +0.2100, −0.2100 |
| `c_2–n_2` | 690, ref7 | +0.1110, −0.1110 |
| `c_2–o_1` | 691, ref7 | +0.5850, −0.5850 |
| `c_2–o_2` | 692, ref7 | +0.1890, −0.1890 |
| `hn2–n_2` | 816, ref7 | +0.3780, −0.3780 |
| `cp–nn` | 729, ref1 | +0.0827, −0.0827 |
| `c–nn` | 515, ref2 | +0.2108, −0.2108 |
| `h*–nn` | 806, ref1 | +0.2487, −0.2487 |

Every positive case has exact raw-Decimal component total `0.0000 e`. Every
requested interaction is independently searched, including missing/ambiguous
rows, orientations and equilibrium dependencies. No successful query substitutes
for complete-model validation.

| Declared fixture | V5 charges | V6 charges | Strict v3 / compatibility v4 model | V4 blocking diagnostics |
|---|---|---|---|---:|
| Primary carbamate `COC(=O)N` | incomplete | complete | incomplete / incomplete | 45 |
| Secondary carbamate `COC(=O)NC` | incomplete | complete | incomplete / incomplete | 78 |
| Tertiary carbamate `COC(=O)N(C)C` | complete | complete | incomplete / incomplete | 113 |
| Aniline `Nc1ccccc1` | incomplete | complete | complete / complete | 0 |
| N-methylaniline `CNc1ccccc1` | incomplete | complete | incomplete / incomplete | 63 |
| Diphenylamine | incomplete | complete | incomplete / incomplete | 25 |
| Mixed aliphatic amine/urethane/ether/arylamine | incomplete | complete | incomplete / incomplete | 129 |
| Urethane DP1 / DP2 / DP3 | incomplete | complete | incomplete / incomplete | 102 / 204 / 306 |
| Arylamine polymer DP2 | incomplete | complete | incomplete / incomplete | 136 |

The urethane PSMILES is `[*:1]CCOC(=O)NCC[*:2]`. It has explicit hydrogen-terminated
alkylene ends and C–C inter-repeat bonds; urethane N remains NH. DP1/2/3 contain
19/36/53 sites. The separate aromatic series uses `[*:1]CNc1ccc([*:2])cc1`, DP2,
with hydrogen-terminated alkylene/aryl ends and aryl–alkylene inter-repeat bonds.
The mixed graph is `NCCOC(=O)NCCOc1ccc(N)cc1`; v5 `na/hn` coexists with v6
`n_2/hn2`, `nn/hn`, existing ether and aromatic labels.

**Exact blockers:** external `gaps.json` and `strict-inspection/` /
`compatibility-inspection/` retain each request's stable sites, supplied/resolved
labels, candidates, missing/ambiguity reason and independent searches. Examples:

* Urethane `c3–c2–o_2–c_2` lacks EBT, MBT, AT and AAT rows. Other ester/urethane
  orientations lack EBT/AT/AAT and required AA couplings. Compatibility BB13 zeros
  remove only their existing narrowly authorized diagnostics; the other couplings
  remain required. Tertiary carbamate additionally lacks BB/BA.
* N-methylaniline has ambiguous angle requests `c3–nn–hn` and `nn–c3–hc`, and
  missing `c3–nn–cp` BB/BA plus cp-containing BB13 and torsional/AA couplings.
  No wildcard priority or new converter zero was invented.
* Ordinary amide retains its missing BB/BA/MBT and conflicting torsion diagnostics.
  Urea, radical and isotopic carbamate remain typing-incomplete. Pyridinium,
  charged carbamate, charged aniline and heteroarylamine remain charge-incomplete.
  Guanidinium retains printed total **0.9999 e**, formal +1 and residual **−0.0001 e**.
  All nine negatives preserve the v5 rejection stage and chemistry reason.

Across all **20** fixtures: typing remains **17/20**; native charges increase
**2/20 → 12/20**; complete models increase **0/20 → 1/20** under either resolver.
Independent numerical verification is **1/20**, executable polymers **0/4**.
The current v6 source ledger retains **87/133 labels with bounded predicates**,
not globally validated labels. No new source labels are claimed: these are fixes
to local family coherence. The frozen J11 audit remains unchanged at 47 unresolved
historical obligations; it is not relabeled as a current-v6 coverage count.

## Independent executable evidence: aniline only

`aromatic_primary-declaration.json` binds the checked original v5-built geometry,
profile v6, unchanged v4 policy, explicit LJ/Coulomb `(0,0,1)` and source/tool hashes.
The independent raw-source resolver verifies inventories, rows, coefficients,
equilibria and charge contributions. `msi2lmp` receives independently authored
types/charges and optional CAR/MDF; it never receives ISLAND's compiled coefficients.
The raw converter data and commands remain retained. The existing H5 angle–angle
equilibrium-role correction is reconstructed from the reference's own angles.
No additional J14 coefficient correction was needed; no converter `-ignore`.

Pinned binaries: msi2lmp SHA256
`a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`;
LAMMPS SHA256
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`.
Both reference the pinned LAMMPS revision; execution logs retain version/build text.
Runtime native preparation and reconstruction require neither tool nor CAR/MDF.

* Original and asymmetric geometries, then retained steps 0/2/4, compare total
  potential, compatible grouped components and all forces. Maximum energy/group
  error: **3.197442310920451e-13 kJ/mol**; force error:
  **6.000477892342815e-12 kJ/(mol·angstrom)**. Declared atols: **1e-5**, rtol **2e-10**.
* Central finite-difference maximum force errors at displacements 1e-4/1e-5/1e-6 Å:
  **3.827440875e-5 / 3.879858497e-7 / 8.638292570e-8** kJ/(mol·Å).
  This supplements independent parameter selection and LAMMPS forces.
* Minimum: energy **23.560337767 → −17.707634513 kJ/mol**, **76 iterations,
  82 evaluations**, maximum force **0.0972288113**, RMS **0.0625718802** kJ/(mol·Å).
  Termination `force_converged`, independent fresh verification true. Original
  force criterion 0.1, maximum 5000 iterations / 10000 evaluations unchanged.
* Four-step BAOAB versus relocated 2+2 workflow: five shared frames, maximum
  coordinate/velocity/energy/component/force/time differences **0**, complete PCG64
  state **exactly equal**. Timestep .1 fs, 300 K, friction 5/ps, velocity seed 78123,
  thermostat seed 99181. Split tolerance 1e-10 absolute / 1e-12 relative.
* Whole dynamics uses 6 evaluations / 3 Contexts; split dynamics 8 evaluations.
  Workflow start uses 5 Contexts (minimum plus first segment); child resume uses 3.
  Completed status/frame reads/no-op resume reconstruct offline with OpenMM,
  RDKit, SciPy, Foyer and ParmEd imports blocked. Rechecksummed source/policy
  mutations are rejected; original bundle bytes remain unchanged.

Facade identity `941d709dabae5a892923ed5a47c081e9a46f084b3035f3b49e5486e2a53f35c9`;
native identity `3ef40a09c26492448de626c035e2a986512e8051d20582d640b9ad99ec4adf77`.
Observed same-host continuation equality is not a cross-platform guarantee or
independent validation of the integrator. No executable polymer is implied.

## Public API and reproduction

```python
system = build_linear_polymer("[*:1]CCOC(=O)NCC[*:2]", dp=2, random_seed=2026)
typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v6")
# Inspect typing, native charges and every source interaction separately.
# prepare_forcefield rejects this particular model's missing couplings.
options = PCFFOptions(
    source_path, (0, 0, 1), (0, 0, 1),
    typing_profile="island_pcff_source_graph_v6",
    resolution_policy="island_pcff_msi_guarded_source_v3",
)
```

Use the example for the full diagnostic report. Existing native persistence,
shared offline identities, bound evaluators, sessions and prepared workflows
support v6 via their existing validators. No new bundle/workflow schema or
scientific engine was introduced. The new rule module participates in cache
invalidation; explicit types cannot bypass its graph predicates.

Commands from the repository root (sources/tools stay external; select new output
directories for repeat experiments, never overwrite the retained ones):

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
OUT=../island-validation/phase4j15-declared

$PY scripts/validate_pcff_urethanes.py --source "$FRC" --output "$OUT/acceptance"
# Expected exit 1: all expectation checks pass, executable polymer gate fails.
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" \
  --declaration "$OUT/strict-declaration.json" --output "$OUT/strict-inspection"
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" \
  --declaration "$OUT/compatibility-declaration.json" --output "$OUT/compatibility-inspection"
$PY scripts/validate_pcff_profile.py --source "$FRC" \
  --declaration "$OUT/aromatic_primary-declaration.json" --output "$OUT/aromatic_primary-vertical"

$PY examples/pcff_urethane_typing.py --source "$FRC"
$PY examples/pcff_urethane_typing.py --source "$FRC" --converter-compatibility
# Both exit 1 with source-coupling diagnostics; historical v5 example exits 0:
$PY examples/pcff_amine_typing.py --source "$FRC"
$PY scripts/audit_pcff_completion.py --source "$FRC" \
  --declaration docs/evidence/phase_4j11_declaration.json \
  --retained-matrix ../island-validation/phase4j10-declared/parity-matrix \
  --output "$OUT/full-source-audit.json" --require-full-source
# Expected exit 1; no gate weakening.

$PY -m pytest -q tests/test_pcff_urethane_domains.py tests/test_pcff_urethane_acceptance.py \
  tests/test_pcff_amines.py tests/test_pcff_amine_acceptance.py \
  tests/test_pcff_tranche_controls.py tests/test_pcff_validation_cache.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

Baseline reproduction uses `execute_control` from
`scripts/inspect_pcff_tranche_controls.py` on every entry of the J15 declaration,
with only `typing_profile` set to `island_pcff_source_graph_v5`; it writes each
original system and result without execution of QM/dynamics. For either resolver
inspection, use the saved system path and SHA256 in `strict-declaration.json` or
`compatibility-declaration.json`. These declarations are retained and checksummed
alongside the original input files, not regenerated during numerical acceptance.

For offline historical revalidation:

```bash
$PY scripts/revalidate_pcff_tranche.py --source "$FRC" \
  --root ../island-validation/phase4j8-declared/acceptance \
  --root ../island-validation/phase4j9-declared/vertical-slice \
  --root ../island-validation/phase4j10-declared/fused-vertical \
  --root ../island-validation/phase4j12-declared/vertical-slice \
  --root ../island-validation/phase4j13-declared/peo-vertical \
  --root ../island-validation/phase4j13-declared/thioether-vertical \
  --root ../island-validation/phase4j14-declared/primary-vertical \
  --root ../island-validation/phase4j14-declared/secondary-vertical \
  --root ../island-validation/phase4j14-declared/primary_protonated-vertical \
  --root ../island-validation/phase4j14-declared/secondary_protonated-vertical \
  --root ../island-validation/phase4j14-declared/tertiary_protonated-corrected-reference \
  --root ../island-validation/phase4j14-declared/amine_oligomer-vertical \
  --root ../island-validation/phase4j14-declared/protonated_sidechain-vertical \
  --root "$OUT/aromatic_primary-vertical" --output "$OUT/historical-revalidation.json"
```

The J13 main and charged-sidechain acceptance CLIs were rerun under their original
profiles; both pass their original negative contracts. The retained J14 acceptance
was reassessed read-only under its original v5 gate. Its previous failures were
not turned into historical successes.

## Verification and remaining work

Software tests use explicitly synthetic source fixtures; actual source and
LAMMPS acceptance above are separate. New tests cover motif specificity, mixed
precedence, explicit wrong-family labels, stable-ID remapping/bond reversal,
coordinate and insertion-order independence, nonmutation, semantic rechecksummed
tampering, unexpected success, arbitrary execution errors, wrong stages/reasons
and missing/duplicate control sets. Historical PCFF/workflow tests remain intact.

Final test counts and historical hash totals are recorded below and in the receipt.
Required work for executable urethane polymers is authoritative resolution of the
listed EBT/MBT/AT/AAT/AA (and some BB/BA) gaps, not another charge adjustment.
Secondary arylamine models additionally need physically justified wildcard-angle
resolution and source couplings. Missing source information cannot be replaced by
analogy, generic aliases, dropped interactions or expanded zero conventions.

Executed verification: **91 focused tests passed** (11.33 s); **1769 ordinary tests
passed, 10 skipped** (338.49 s). Ruff passed; both primary and isolated Foyer
`pip check` commands reported no broken requirements. The ordinary suite includes
historical PCFF, bundle and corrected I3 failed-attempt workflow regressions.
All **4210** snapshotted historical files retained their hashes. Thirteen historical
bundle/workflow roots plus the new aniline root passed offline reconstruction,
status, frame and completed no-op resume checks. Neither the new diagnostic CLI
nor the full-source CLI returns zero for this partial tranche. The urethane
examples reject both policies as declared; the unchanged v5 amine example passes.
The [compact gap ledger](evidence/phase_4j15_gaps.json) records unique type patterns,
exact requests and independent searches; the external ledger retains full site roles.

Execution environment: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, RDKit 2026.3.6,
OpenMM 8.6.1 Reference. Source inspection and native reconstruction remain offline;
RDKit is needed for the example's construction, OpenMM for evaluation and SciPy
for minimization. LAMMPS/msi2lmp are reference-only dependencies. No QM was run.
