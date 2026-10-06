# Phase 4J4 — bounded source-driven chemical predicates

## Base, scope and identity

Branch: `codex/phase-4j4-pcff-source-domains`, based on reviewed J3
`b5c908babb4ce56d776f4ab9b80c1bca73e9a861`. At initial fetch, origin/main
was still `d7e30b0fb9b9e93f0c473a9b9e7098dd601558dc`; J3 was not present
there. This branch retains that dependency; main was not modified or merged.
The separate J3.1 corrective commit is documented in [phase_4j3_1.md](phase_4j3_1.md).

The opt-in profile is `island_pcff_source_graph_v4`, using the unchanged
`island_pcff_positional_fallbacks_v2` resolver. Historical v1/v2/v3 profiles,
model equations, signatures, defaults and policies are unchanged. The facade,
native persistence, bundle loader and parameter assignment explicitly recognize
v4. They still recompute semantic decisions. Explicit labels pass the same
chemical predicates; they cannot bypass an incomplete automatic result.
No new evaluator or conversion equation is introduced.

This is **partial source coverage**. The six additional label predicates do
not all produce complete models. Four hydrogen-halide molecules do; they
extend the existing `h`/halogen domain rather than adding four globally
validated atom labels. Conservation and numerical agreement are separate from
chemical accuracy. `production_validated=False` and
`simulation_readiness="not_established"` throughout.

## Source and predeclared experiment

External source: LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
No source library is bundled, edited, downloaded at runtime or substituted.
Independent raw-line inspection recounts 133 distinct atom labels, 134
ordinary equivalences, 108 automatic equivalences and 564 bond increments.

[The declaration](evidence/phase_4j4_declaration.json) was saved before J4
implementation. It fixes seed 2026, ten positive molecules, seventeen
near-miss controls, the reference, numerical tolerances, two geometries and
output locations. Scientific inputs were not tuned after failures. The
initial numerical harness used a nonexistent reporting attribute,
`calculation_fingerprint`, after a successful first-frame comparison. That
failed run is retained in `numerics/`; the corrected harness uses the existing
`evaluation_fingerprint`. `numerics-corrected/` is a separate experiment with
identical scientific settings.

## Chemical rules and limitations

Each predicate requires authoritative elements/atomic numbers, standard
masses (the existing 0.02 mass check), no undeclared isotope/radical state,
explicit hydrogens, specified formal charges, bond orders and ring/aromatic
state. A failed site blocks its component. The existing D2O exception remains
explicit and unchanged. No stereochemistry or graph is stripped or modified.

Source lines below refer to the exact external file above; the machine-readable
source audit preserves the selected atom and equivalence rows.

| Label / motif | Bounded predicate | Atom / ordinary / automatic lines |
|---|---|---|
| `h`, H–F/Cl/Br/I | Neutral isolated single-bond diatomic, standard masses, no ring/aromaticity | 98 / 240 / 366 |
| `s'` | Neutral acyclic terminal S=C with two single C/H substituents; C receives `c=` for CH2 or `c=2` otherwise | 166 / 308 / 426 |
| `nh+` | Unsubstituted isolated aromatic C5NH six-ring, exactly one N with formal +1 | 134 / 276 / 399 |
| `c+`, `nr` | Explicit guanidinium C(NH2)3 graph: one C=N+ and two C–N, formal +1 on the double-bonded N, no ring/aromaticity | 68,140 / 210,282 / 336,405 |
| `n3n`, `n4n` | Neutral three-single-bond N directly adjacent to one ordinary amide carbonyl in a three/four-membered ring | 122,125 / 264,267 / 387,390 |

Hydrogen-halide use of `h` is justified by the explicit H–halogen increment
rows (475,697,769,790) and corresponding automatic bond rows, not by a generic
hydrogen fallback. The atom-row comment lists C/Si/H and is not itself evidence
for H–halogen typing; this discrepancy is recorded. The v4 policy explicitly
limits that interpretation to these isolated neutral molecules.

For strained lactams, the source distinguishes sp2 ring N (`n3n/n4n`) from
sp3 ring N (`n3m/n4m`). We restrict the former to the adjacent-amide motif.
Pinned LUNAR's alternative interpretation of the suffix as a nitrogen-neighbor
condition is not adopted. Ring strain and electronic accuracy are not validated
by this graph predicate; these molecules still fail complete parameter coverage.

The two guanidinium source aliases `n2` and `nr` have overlapping descriptions.
This profile explicitly chooses `nr` for the stated NH2 motif and keeps `n2`
unresolved. It does not accept arbitrary interchange of those labels. No
resonance averaging or charge redistribution is applied.

Near-miss controls retain O=O, N#N, deuterated/charged HCl, sulfoxide/sulfone,
SO2, phosphonium oxide, PCl5, neutral guanidine and charged thiolate as unresolved.
Ordinary pyridine and saturated small-ring amines retain their older predicates;
they do not become `nh+` or `n3n/n4n`. The existing neutral P=O motif remains
typed but fails its original native-charge check.

## Native charges and coupling outcomes

All assignments retain the existing exact/ordinary/positional automatic
searches, oriented endpoint contributions, source rows, candidates and
component formal-charge checks (1e-12 e). Missing rows are not zero.

| Declared cases | Typing | Native charge | Complete model |
|---|---|---|---|
| HF, HCl, HBr, HI | pass | pass | pass: one source quadratic bond each |
| Thioformaldehyde, thioacetone | pass | pass | blocked: conflicting automatic angle candidates and missing required couplings |
| Pyridinium | pass | fail | missing `cp–nh+` and `nh+–h*` increments after declared searches |
| Guanidinium | pass | fail | source row 528 yields total 0.9999 e versus formal +1 |
| Three/four-membered lactams | pass | pass | conflicting angle candidates, missing bond/angle and torsional couplings, unresolved equilibrium dependencies |

Independent Decimal summation of three `c+–nr` endpoint pairs gives
`3*(0.2653+0.0680)=0.9999`; the N–H rows conserve charge pairwise. This is a
retained source-conservation failure, not a floating-point problem to repair.

Every unresolved interaction ID and its original diagnostics is in the gap
receipt. No automatic cross-term family, Wilson convention or missing
parameter was invented. Nonzero Wilson equilibrium and torsion–torsion remain
outside the established executable interpretation. The previously documented
Class II directional roles, angle–angle dependencies and reference corrections
remain unchanged.

## Independent numerical path

The local executable identifies itself as LAMMPS 30 Sep 2026, git
`patch_30Sep2026-65-ge891a3e10`; executable SHA256:
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`.
Build/version output and command/input/output hashes are retained externally.

For HF/HCl/HBr/HI the reference independently reads the **raw** automatic bond
columns and unique quadratic rows, and independently sums raw increment rows.
It checks native selected row IDs, coefficients, charges and the complete
one-bond inventory afterward. It does not export the production coefficient
list as an assignment oracle. LAMMPS uses `bond_style harmonic`:
`E = K*(r-r0)^2`, with K in kcal/(mol Å²), r in Å, no extra 1/2.
Native conversion is `K*4.184` to kJ/(mol Å²).

At this revision, [BondHarmonic::compute, lines 47–96](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/MOLECULE/bond_harmonic.cpp#L47)
implements that energy and analytical force. The entire pair inventory is one
excluded 1–2 pair under declared LJ=Coulomb=(0,0,1). The reference's zero pair
style is equivalent only for these isolated cases; source charges and LJ rows
remain present and validated in the native model. No claim about active
nonbonded interactions or larger hydrogen-halide clusters follows.

Two non-minimized coordinate frames are compared: seed 2026 input and the
predeclared 0.015 Å sinusoidal perturbation. Energy/component and every force
comparison use atol=1e-5 in canonical kJ/mol and kJ/(mol Å), rtol=2e-10.
Independent energy finite differences use 1e-4, 1e-5, 1e-6 Å displacements.
The reference does not use production conversion, interaction or pair helpers.

[msi2lmp README, lines 145–147](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README#L145)
explicitly lacks automatic-equivalence supplementation and bond increments.
Its GetParameters.c ordinary matching (bond lines 121–128, angle 177–185)
cannot establish this fallback or charge path. No converter run or `-ignore`
was used as evidence for them. Original J1/J2 converter failures/corrections
remain untouched; no converter or external typer is a runtime dependency.

## Public path and durability

```python
from island.chemistry import from_smiles
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield, create_evaluator

system = from_smiles("[H]Cl", random_seed=2026)
request = ForceFieldRequest("pcff", PCFFOptions(
    source_path="/explicit/local/pcff.frc",
    typing_profile="island_pcff_source_graph_v4",
    resolution_policy="island_pcff_positional_fallbacks_v2",
    lj=(0, 0, 1), coulomb=(0, 0, 1)))
prepared = prepare_forcefield(system, request)
result = create_evaluator(system, prepared).evaluate_fresh()
```

For each complete case, the acceptance saves the existing I2 bundle, moves its
directory, and starts a child Python process. That process reconstructs with
OpenMM/RDKit/Foyer/SciPy/ParmEd imports blocked, then restores imports and
compares fresh, reusable-session and independent-fresh evaluations. Native,
facade, parameter/model/evaluation identities and every force/component are
compared; bundle hashes must remain unchanged. No new preparation occurs in
the child. Comparison tolerances remain atol=1e-10, rtol=1e-12. This task does
not add dynamics or a new workflow schema.

## Reproduction commands

The available environment is `../island-validation/phase4e2-env/bin/python`.
Set `PY` to it, `FRC` to the explicitly verified source and `LMP` to the pinned
executable. Use new output directories: publication refuses existing paths.

```sh
$PY scripts/audit_pcff_organic_source.py --source "$FRC" --output source-audit.json
$PY scripts/validate_pcff_organic_domains.py --source "$FRC" --output matrix --require-full-source
# Expected nonzero: source coverage is partial.
$PY scripts/validate_pcff_organic_domains.py --source "$FRC" --output controls --controls-only
$PY scripts/validate_pcff_organic_numerics.py --source "$FRC" --lammps "$LMP" --output numerics
$PY examples/pcff_organic_domains.py --source "$FRC" --evaluate
$PY examples/pcff_organic_domains.py --source "$FRC" --smiles 'NC(N)=[NH2+]'
$PY -m pytest -q
ruff check .
$PY -m pip check
```

## Remaining full-source work

The full 133-label ledger retains graph, explicit-label, charge, nonbonded,
interaction-family, evaluator and independent-verification status separately.
A successful molecule never promotes a label globally. Sulfur atom rows have
connections one or two: assigning a generic S to sulfoxide/sulfone central S
would be unsupported. `sf` has a sulfonate comment but connection count one;
that conflict needs authoritative interpretation. `p=` names phosphazene but
specifies five neighbors, so PCl5 or common four-coordinate phosphazene drawings
cannot be accepted merely by element/degree. Charged phosphorus conventions,
imidazolium H equivalence, hydroxyl/metal/zeolite labels, isotope special masses
and overlapping generic aliases need further evidence.

The next useful work is an independently documented interpretation of those
specific source rows and of the ambiguous automatic angle candidates, followed
by complete interaction/force verification. Do not remove required couplings or
normalize charge residuals to expand a headline coverage count. Periodicity,
packing, curing and scientific readiness remain separate capabilities.

## Executed evidence and counts

The machine-readable [summary](evidence/phase_4j4.json),
[full label/family ledger](evidence/phase_4j4_coverage.json),
[exact failed interaction diagnostics](evidence/phase_4j4_gaps.json),
[source-row audit](evidence/phase_4j4_source_audit.json), and
[artifact hashes](evidence/phase_4j4_artifacts.json) keep the gates separate.

| Measure | Before J4 | After this bounded implementation |
|---|---:|---:|
| Source labels with bounded graph predicates | 80/133 | 86/133 |
| Labels without a justified implemented predicate | 53/133 | 47/133 |
| Universally verified source-label domains | 0/133 | 0/133 |
| Original J3 fixture typing under the selected profile | 33/37 | 34/37 (pyridinium gains typing only) |
| Historical J3 native-charge / model evidence | 27/37; 15/37 | unchanged retained records |
| J4 declared fixture typing / charges / complete models | not executed | 10/10; 8/10; 4/10 |
| J4 executable and independently verified complete models | not executed | 4/4 complete candidates (4/10 declared cases) |

The ten-case J4 matrix includes pyridinium, already an unresolved J3 fixture;
the matrix denominators must not be added as disjoint chemical sets. Counts
are bounded evidence, not full-source completion. The four complete models
are all hydrogen halides, a single narrowly supported chemical family. The
six newly recognized specialized labels still lack complete executable
interaction/charge coverage for the declared complex molecules.

Maximum discrepancies across the two declared geometries:

| Case | Energy (kJ/mol) | Force component (kJ/(mol Å)) | FD force error at 1e-6 Å | Relocated fresh/session difference |
|---|---:|---:|---:|---:|
| hydrogen_fluoride | 2.13163e-14 | 4.83169e-13 | 1.4728e-08 | 0.0 |
| hydrogen_chloride | 5.55112e-17 | 7.10543e-15 | 3.06971e-09 | 0.0 |
| hydrogen_bromide | 2.23155e-14 | 1.15108e-12 | 4.49074e-09 | 0.0 |
| hydrogen_iodide | 1.62093e-14 | 5.25802e-13 | 8.47405e-09 | 0.0 |

All four numerical/bundle candidates passed. The two finite-difference
refinements expose the expected decrease in truncation error at the perturbed
geometry, followed by floating-point subtraction error; all final errors stay
well below the declared force tolerance. No cross-platform bitwise claim follows
from observed zero reconstruction differences on this host.

The [initial harness failures](evidence/phase_4j4_initial_numerics.json) and
[near-miss outcomes](evidence/phase_4j4_controls.json) remain separate from the
successful corrected numerical evidence. The numerical/bundle summary exits 0;
`--require-full-source` exits 1, as did the complete matrix's full-source gate.

The [preservation receipt](evidence/phase_4j4_preservation.json) verifies 551
historical files without rewriting them. All 37 original typing records retain
their identities. Two retained J3 bundles and two completed workflows passed
read-only reconstruction/status/frame/no-op-resume checks with scientific
imports blocked. The original records and source hashes remain unchanged.

Remaining unimplemented or ambiguous labels:

`Ag Al Au Br Cl Cr Cu Fe K Li Mo Na Ni Pb Pd Pt Sn W az c_a ca+ cg ci cr h* h+ hb hi hoa hos n+ n1 n2 nb nho ni nz oah oas ob oe osh oss p= s- sf sz`.

Uppercase metal/ion/halogen labels are distinct source labels, not aliases for
lowercase covalent types. Their unresolved state is preserved. The ledger
provides source lines, equivalences, family-specific candidate inventories and
next actions for each label. An occurrence in a candidate list is not proof
of an interaction assignment.

Actual output root: `../island-validation/phase4j4-declared/`. Full source and
external tools remain outside Git. No new QM, optimization, dynamics, source
editing or historical acceptance rerun was needed.

### Verification executed

- Ordinary suite: **1552 passed, 10 skipped**, 1049.33 seconds.
- Historical PCFF/prepared-workflow selection: **162 passed**.
- Final isotope, new domain and synthetic aggregate-gate tests: **89 passed**.
  Additional focused tests added after ordinary collection were executed here.
- Ruff: passed. Main and isolated Foyer environment pip checks: passed.
- Real pinned-source positives and controls, independent LAMMPS/FD checks,
  child-process bundle reconstruction and read-only historical inspection: executed.
- HCl evaluation example, guanidinium diagnostic example and historical v3 D2O
  explicit-label example: executed. The guanidinium example reports failure of
  native charge conservation; its diagnostic exit is not scientific acceptance.

Environment: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, OpenMM 8.6.1,
RDKit 2026.3.6, pytest 9.1.1. The [verification receipt](evidence/phase_4j4_verification.json)
records exact commands, log/report hashes, retained failures and environment.
No large scientific matrix, QM or new trajectory was run.
