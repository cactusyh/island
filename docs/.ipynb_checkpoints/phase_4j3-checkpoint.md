# Phase 4J3 — specific PCFF domains and source-coverage accounting

## Scope and repository

Branch `codex/phase-4j3-pcff-chemical-domain-coverage` starts at
`d7e30b0fb9b9e93f0c473a9b9e7098dd601558dc`. Its tree equals corrected J2
`b963c95bd3fc8e0e1cf5bfff3e8cdd38ddb31e9c`; there is no unmerged dependency.
Main and historical records are not changed.

**Full-source coverage remains incomplete.** The bounded domain gate is separate
from the 133-label target. Neither parsing, a type predicate, native neutrality,
nor one verified molecule establishes a type's full chemical domain.
`production_validated=False`, `simulation_readiness="not_established"`.

## Versioned APIs and chemistry

Explicit opt-ins:

- Typing: `island_pcff_source_graph_v3` through `type_pcff_atoms`,
  `assign_pcff_source_types`, or `PCFFOptions.typing_profile`.
- Resolution: `island_pcff_positional_fallbacks_v2` through native assignment
  and `PCFFOptions.resolution_policy`.
- Coverage: `pcff_domain_coverage(source)` in `pcff.domain_coverage`.

Historical CHO/J1/J2 profiles and policies keep their payloads and fingerprints.
The new typing schema is `island_pcff_source_typing_v3`. Existing versioned
assignment/model containers retain an explicit new policy and interpretation
identity; no evaluator identity formula is duplicated. The shared offline helper
continues to derive fingerprints. The generic facade, bundles, sessions and
prepared workflows consume these validated native records.

`domains.py` contains the exact predicates and source evidence. They supplement
unresolved J2 sites, without weakening whole-component rejection. Explicit labels
must pass the same predicates. The 14 additional labels are `h`, `s`, `dw`, `n`,
`npc`, `nn`, `n=`, `n=1`, `n=2`, `p`, `hp`, `o=`, `o`, `ho2`.

| Domain | Bounded rule and limitation |
|---|---|
| Molecular hydrogen | Neutral explicit H–H only; no generic H/ion assignment |
| Sulfide hydride | Neutral H–S–H; sulfur oxidation states remain unresolved |
| Isotope water | Fully deuterated water: H element, isotope 2, atomic number 1, mass 2.014–2.014102; source `dw` D representation bridge |
| Ordinary amide | Neutral three-single-bond N with exactly one ordinary amide carbonyl; no urea/carbamate inference |
| Substituted aromatic N | Neutral aromatic 5/6-ring N with two aromatic bonds and a single external C, no H |
| Aromatic amine | Neutral nonring N attached to aromatic C; retains `hn2`, including missing charge coverage |
| Imines | Explicit terminal/next-to-terminal/internal C=N roles; neutral acyclic environments only |
| Phosphorus oxo | Neutral tetra-coordinate P=O with three single C/H/O neighbors; not generic phosphorus typing or an inferred oxidation state |
| Peroxide/acid OH | Peroxide O and peroxide/carboxylic-acid H conventions; ordinary alcohol `oh/ho` unchanged |

The RDKit constructor previously preserved isotope mass but dropped isotope
identity. A failing regression demonstrated this for `[2H]O[2H]`. It now retains
nonzero isotope and radical-electron metadata. Ordinary neutral metadata remains
unchanged; old isotope records are not repaired. No rule infers isotope identity
from mass alone. A remapping test initially omitted aromatic bond flags; correcting
the test preserved the authoritative aromatic graph rather than loosening typing.

Source pin remains LAMMPS `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`,
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
The source library remains external. Source atom descriptions, ordinary and
automatic equivalences and actual increment rows support each bounded rule.
The independently inspected LUNAR revision
`67dabeda9e6bd3cc8f968aa0c6a88730886138ef` supplies comparison context:
`PCFF.py` lines 513–516 (ho2 precedence), 795–807 (imine roles), 876–881
(substituted ring N), 904–907 (aromatic amine). No upstream code is copied.
Its broad element/coordination rules are not accepted as universal chemical
proof; its apparent lowercase `element == 'n'` amide branch is not an executed
reference validation. Generic `nb`, charged heterocycles and resonance aliases
remain unresolved. The committed rule evidence lists the FRC row locations.

## Native charges and parameter semantics

Every bond is still processed once. Policy v2 additionally records all attempted
lookup paths: direct, ordinary bond equivalence, automatic `bond_increment`
column, candidate orientation, and selected endpoints. Per-site zero base and
formal charge are explicit. No alias substitution, partial-vector normalization,
or component charge repair is allowed.

`na–hn2` remains missing after all three searches. `hn2` is not replaced with
`hn`. Peroxide/acid `ho2` failures remain failures; a source label existing does
not prove its increment coverage. The P=O experiment exposes a component-charge
failure from source endpoint values, rather than an adjusted neutral vector.

Policy v2 adds automatic **nonbond** column resolution into existing ordinary
9–6 rows when ordinary lookup is missing. It uses that column alone, not bond
or angle equivalence. Direct/ordinary precedence, source version selection,
wildcard ambiguity and legal center-preserving permutations remain unchanged.
There is no invented automatic Class II cross-term family. Lower-order base
substitution does not make a coupling vanish. The acetal, acetone, ester and
chloromethane gaps remain exact interaction-level failures in the new gap ledger.

Nonzero Wilson equilibrium remains unsupported: signed-angle/permutation
semantics are not independently established. The empty torsion–torsion section
remains an explicit operational scope exclusion, not a recovered zero parameter.
No CHO-specific zero rule has been extended to other chemistry.

## Bounded converter audit

The local LAMMPS checkout was verified at the exact source revision above.
[README](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README)
lines 140–148 state that automatic equivalence and bond increments are not
supported. `GetParameters.c` lines 69–74, 121–128, 177–185 apply ordinary
family equivalences; 203–245 assemble three-body cross terms; 284 onward
assembles torsional cross terms; 1009–1047 permutes out-of-plane roles;
1056–1166 implements exact/wildcard and reversal search. These routine locations
are pinned in the external audit receipt. Converter file-order precedence and
its historical angle–angle center reversal are not an ISLAND selection oracle.
The J1 raw failure and corrected reference are untouched. No `-ignore` use and
no converter runtime dependency are introduced.

The new triatomic reference uses independently declared source rows and its own
bond/angle inventory, not an export of ISLAND coefficients. H2S uses rows 1743,
2402, 3590, 3905; D2O uses 1753, 2350, 3564, 3881. The source-row selection and
multiplicities are checked separately from numerical comparison. Both systems
have two 1–2 pairs and one 1–3 pair, all excluded by the declared policy; native
charges and nonbonded records remain intact. These cases do not newly verify
active nonbonded mixing or torsional physics.

## Declared experiments and limits

`docs/evidence/phase_4j3_declaration.json` was published before implementation,
with all 133 pre-J3 labels, prior receipts, the 18 retained J2 fixtures and 19
new fixtures. Fixed seed 2026, original and asymmetric sinusoidally perturbed
coordinates, energy/force atol `1e-5`, rtol `2e-10`; independent kernel central
finite differences at `1e-4`, `1e-5`, `1e-6` angstrom.

The declaration incorrectly expected `c=1` for internal imine `CC=NC`; source
terminal-role interpretation yields `c=2`. The mismatch is retained and excluded
from independently matched-label success. It is not edited into a passing
expectation. The first matrix harness mishandled the H2S all-atom expected list;
that failed attempt remains in `matrix/`. `matrix-final/` fixes list handling
and uses the final explicit charge-search record format. The initial numerical
run completed H2S comparisons but failed JSON serialization of an owned mapping;
`numerical/` is retained. `numerical-corrected/` converts that mapping for JSON
and repeats the same bounded inputs and tolerances; it does not change physics.

Workflow settings: existing Reference session, fmax 0.1 kJ/(mol·angstrom),
5000 iterations/10000 evaluations, unchanged line-search and verification
criteria; BAOAB 4 steps, two segments of 2, 0.1 fs, 300 K, friction 5/ps, seeds
78123/99181, 4 evaluations/3 frames per segment, 6/5 uninterrupted. Pair choices
LJ and Coulomb `(0,0,1)` remain explicit model choices. The saved exact minimum
and initialized state feed both paths; child resume reconstructs checked inputs
without preparation. Status/frames/completed resume are checked with optional
scientific imports blocked. Short trajectories are interoperability evidence,
not equilibration or independent integrator validation.

## Reproduction

Set `PY` to a compatible environment, `FRC` to the verified external file,
`LMP` to the audited executable, and `OUT` to a new directory. No silent downloads.

```sh
$PY examples/pcff_domain_coverage.py --source "$FRC" --smiles S --explicit-checked
$PY scripts/validate_pcff_domains.py --source "$FRC" --output "$OUT/matrix"
$PY scripts/validate_pcff_domain_numerics.py --source "$FRC" --lammps "$LMP" --output "$OUT/numerical"
$PY scripts/validate_pcff_domain_workflow.py --source "$FRC" --bundle "$OUT/numerical/hydrogen_sulfide/bundle" --evidence "$OUT/numerical/hydrogen_sulfide/outcome.json" --output "$OUT/workflows/hydrogen_sulfide"
$PY scripts/validate_pcff_domain_workflow.py --source "$FRC" --bundle "$OUT/numerical/heavy_water/bundle" --evidence "$OUT/numerical/heavy_water/outcome.json" --output "$OUT/workflows/heavy_water"
$PY scripts/check_pcff_domain_charges.py --source "$FRC" --numerical "$OUT/numerical" --output "$OUT/charge-readback.json"
$PY scripts/inspect_pcff_domains_offline.py --source "$FRC" --bundle "$OUT/numerical/hydrogen_sulfide/bundle" --bundle "$OUT/numerical/heavy_water/bundle" --workflow "$OUT/workflows/hydrogen_sulfide/relocated" --workflow "$OUT/workflows/heavy_water/relocated" --output "$OUT/offline.json"
$PY scripts/summarize_pcff_domains.py --source "$FRC" --matrix "$OUT/matrix" --numerical "$OUT/numerical" --workflows "$OUT/workflows" --offline "$OUT/offline.json" --charge-readback "$OUT/charge-readback.json" --output "$OUT/summary"
```

The matrix-only command exits nonzero: it cannot certify numerical/workflow or
full-source gates. The summary requires the mandatory named cases, source checks,
complete numerical checks, bundles and successful matching workflow identities.
`--require-full-source` remains nonzero. Fixtures that fail are retained.

## Remaining full-source work

All 133 labels remain in the ledger, including every remaining predicate gap.
There are 80 narrowly implemented labels and 53 unresolved labels; none is
promoted to global verified coverage by a single molecule. Every family carries
source/equivalence candidates and a next action. Candidate occurrence is not a
valid interaction selection. The gap ledger records actual interaction failures.

Priority follow-ups require authoritative chemical conventions for charged
nitrogen/resonance, sulfur oxidation and ring specializations, phosphorus charge
states, coordinated metals/ions/zeolite bridging environments and specialized
mass labels. Carbonyl/ester/acetal/halogen couplings need actual applicable source
rows or a separately justified new model policy. Missing rows, ambiguous chemistry
and missing algorithms stay separate. This phase does not implement periodicity,
packing, curing, crosslinking or network simulation.

### Chemical-flag boundary regression

A final semantic check independently reproduced acceptance of impossible aromatic
flags on sulfide S, isotope H, phosphoryl O and acid OH in the new supplements.
Four tests failed before tightening the local predicates and pass afterward.
The new rules now preserve those rejection diagnostics. Valid graph payloads,
coefficient values and all retained acceptance identities remain unchanged; the
post-boundary charge/bundle readback validates them without re-signing. The
pre-fix test log is retained separately. This is validation of malformed input,
not an extension of the declared chemistry.

## Measured coverage and numerical results

| Measure | Reviewed J2 | J3 result |
|---|---:|---:|
| Labels with bounded graph predicates | 66/133 | 80/133 |
| Original 18 fixtures: typing | 18/18 | 18/18 |
| Original 18 fixtures: native charges | 16/18 | 16/18 |
| Original 18 fixtures: complete operational models | 12/18 | 12/18 |
| Additional declared fixtures: typing | not executed | 15/19 |
| Additional declared fixtures: native charges | not executed | 11/19 |
| Additional declared fixtures: complete operational models | not executed | 3/19 |
| Combined declared fixtures: typing / charges / complete models | 18 / 16 / 12 of 18 | 33 / 27 / 15 of 37 |
| Newly executed independent whole-system targets | not executed | 2/2 |
| Newly completed relocated four-step workflow targets | not executed | 2/2 |
| Globally verified full-domain source labels | not established | 0/133 claimed |

The 15 complete models are eligible for evaluation, not 15 new independently
verified force models. This phase actually executes H2S and D2O; molecular
hydrogen's new complete assignment is retained without claiming an additional
LAMMPS validation. Six new bounded domains have independently matched labels
and native-charge completion: molecular hydrogen, sulfide hydride, isotope water,
ordinary amide, substituted aromatic N and imines. Eleven new charge-complete
fixtures include an unmatched-label imine and a multifunctional ring amide without
predeclared per-site expected labels; these are not counted as independent typing
successes. Existing silicon/organic-halogen rules are not counted as new domains.

| Case | Maximum energy difference (kJ/mol) | Maximum force-component difference (kJ/(mol·Å)) | Minimum iterations/evaluations | Final fmax / RMS |
|---|---:|---:|---:|---:|
| H2S | 1.0303e-13 | 3.4817e-13 | 6 / 9 | 0.008864 / 0.006817 |
| D2O | 2.2649e-14 | 8.1002e-13 | 6 / 10 | 0.015292 / 0.011935 |

Numerical comparisons include both fixed and perturbed starting coordinates.
Independent raw-energy central differences converge to force errors below
2.7e-8 at displacement 1e-6 Å. The two minima terminate `force_converged`, with
independent fresh OpenMM verification. The independent LAMMPS comparisons are
at the declared initial/perturbed frames, not a new trajectory-force oracle.

H2S energy falls from 28.09757174 to 2.0826152e-8 kJ/mol; D2O from 2.72106286 to
3.2572920e-8. Both split/uninterrupted trajectories retain steps 0–4; every
reported coordinate, velocity, energy, force, component and time discrepancy
is zero on this host, and complete PCG64 state equality is exact. Each start uses
5 Contexts (2 minimization, 3 first-segment), child resume 3, uninterrupted
comparison 3. Whole propagation consumes 6 evaluations, split propagation 8,
including independent boundaries. Elapsed workflow harness times were 501.10 s
and 472.74 s, dominated by repeated offline integrity work; no speedup is claimed.

## Executed checks and preservation

- Ordinary suite: **1500 passed, 10 skipped**, 1039.81 s. Optional/external
  skips remain skips; no substitute success is claimed.
- Final focused J3/review/fallback/gate regressions: **79 passed**, 50.14 s,
  including the final boundary and ownership tests.
- Ruff: passed. `pip check`: passed in both retained main and Foyer environments.
- Explicit-checked domain example: executed successfully on the real source.
- Real 37-case matrix: executed; failed chemistry/charges/couplings retained.
- Independent LAMMPS comparisons, finite differences and source-row selection:
  two declared models, two frames each, passed. Executable SHA256
  `db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`;
  build/source evidence is retained separately from the pinned documentation.
- Four-step relocation/child continuation: both cases passed; offline bundle,
  status, frame and completed no-op resume checks passed with OpenMM, SciPy,
  RDKit, ParmEd and Foyer imports blocked.
- Historical read-only reconstruction: CHO and J1 bundles; I3, J1 and J2 completed
  workflows (11, 5 and 5 frames respectively), passed without mutation.
- 277 tracked external historical files and nine prior committed J1/J2 evidence
  files retain their hashes. The selected historical bundle/native identities
  are in the read-only receipt. Other unrelated directories were not accessed
  for modification, and no historical scientific record was re-signed.

Execution root: `../island-validation/phase4j3-declared/`. Its source audit,
per-case raw inputs/outputs, declarations, native records, bundles, retained
failures and separate-process receipts remain external. Committed J3 evidence
contains the manifest hashes and concise outcomes; full source libraries and
scientific bundle/checkpoint artifacts are not added to Git. The main environment
used NumPy 2.4.6, OpenMM 8.6.1 and SciPy 1.17.1. The independent charge readback
was repeated after the final chemical-flag rejection correction, preserving
both accepted bundle identities.

The initial `summary/` is a developmental aggregation. `summary-final/` also
requires the independent charge readback and dependency-blocked offline receipt;
that is the delivery gate. Original incomplete developmental outputs remain
available. These receipts establish internal consistency and bounded numerical
fidelity, not authenticity of arbitrary supplied evidence or chemical accuracy.

Evidence index: [declaration](evidence/phase_4j3_declaration.json),
[133-label ledger](evidence/phase_4j3_coverage.json),
[interaction/type gaps](evidence/phase_4j3_gaps.json),
[execution and preservation receipt](evidence/phase_4j3.json).
The delivered bounded summary exits 0; the same checks with
`--require-full-source` exit 1. No merge was performed.
