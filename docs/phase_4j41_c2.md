# Phase 4J41-C2 — pinned-source Si/P completion decisions

C2 does **not** qualify an additional mixed-silicon or phosphorus domain.
The retained candidates still lack required source rows or an independently
supported chemical type mapping. These are recorded failures, not zero-valued
coefficients. Finite neutral silane remains executable with its unchanged J41
model identity. The tranche acceptance gate remains **1**.

`production_validated=False`; `simulation_readiness="not_established"`.

The base is merged C1, `a2feb33aebfb30ceb30cf8d65edde3c6df83d0ed`.
The frozen tested implementation is
`761ae1a1f02ab710817f1b6d5b686622973e5f16`. No
historical implementation, receipt, source definition or numerical bound is
rewritten. The branch is `codex/phase-4j41-c2-source-complete-sip`.

## What C2 adds

The additive `island.forcefields.pcff.sip_completion` API produces an immutable
completion decision bound to the entire C1 assessment, final graph, atom order,
coordinates, source distribution, parameter queries and native increment
charges. It separately lists unresolved interaction contexts and phosphorus
candidate rejection reasons. Reloading reconstructs the decision from the
actual graph and pinned source; recalculating a JSON hash cannot authorize
changed source rows, missing cross terms or guessed atom types.

`assess_pcff_sip_completion(..., domain=...)` audits `mixed_silane`, `siloxane`,
`organosilicon_polymer`, `phosphorus_audit` and `mixed_sip_audit`. Domain names
are assessment scopes, not executable profiles. The only C2 execution profile
is `island_pcff_finite_silane_c2_v1`; preparation delegates to the unchanged
C1/J41 silane model. No mixed-Si or P execution profile is registered.

The separate ASCII oracle reconstructs every graph interaction, ordinary and
automatic equivalence path, charge increment and selected numerical row. For
each search it records the family/namespace, full table row IDs, row count and
SHA-256 of the newline-joined original numerical rows. Direct and permitted
equivalent searches that exhaust their source tables remain missing. A bond,
angle or automatic lower-order row cannot supply an absent Class-II coupling.
All 11 Si-partner inventories from C1 are freshly generated, including empty
families and wildcard candidates that carry no execution authority.

The source pins are unchanged:

| File | SHA-256 |
|---|---|
| `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `pcff.rlb` | `82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36` |
| `pcff_templates.dat` | `dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093` |

## Silicon source boundary

| Candidate | Source/typing result |
|---|---|
| Silane `[SiH4]` | Complete; historical finite model retained |
| Methylsilane `C[SiH3]` | 9 bond–bond-1–3 and 9 end-bond–torsion contexts unresolved |
| Disiloxane `[SiH3]O[SiH3]` | 54 required coupling contexts unresolved |
| Hexamethyldisiloxane `C[Si](C)(C)O[Si](C)(C)C` | 378 required coupling contexts unresolved |
| Sila-oligomer `C[SiH2]C[SiH2]C` | 98 missing contexts and one ambiguous angle dependency; wildcard-only Si torsions also rejected |
| Si–N/S/F/Cl/Br/I | No independently qualified complete chemical typing and coupling path |
| Mixed C/H/O/Si/P | Phosphorus mapping blocked before parameter publication |

Every coupling retains site tuples, supplied types, searched rows, equivalence
paths and unresolved dependencies. Si–C/H/O/N/S and halogen source inventories
include ordinary/automatic equivalences, increments, nonbonded rows, bonds,
angles, dihedrals, impropers, bond–bond, bond–angle, angle–angle, angle–torsion,
bond–torsion and angle–angle–torsion terms. Zeolite type `sz`, automatic Si
aliases and companion atom types are inventory evidence only; they cannot
replace the selected chemical environment.

The converter's non-`cp` bond–bond-1–3 zero initialization is an assembly
convention, not an explicit numerical FRC row. It cannot satisfy the requested
no-missing-coupling rule. Existing exact zeros in actual source rows remain
valid source values. Neither charges nor parameters are normalized,
interpolated or invented.

## Phosphorus decisions

The complete original P declarations and all three companion P template blocks
are retained. Each final-graph P site records formal charge, neighbor count,
bond-order sum, conventional oxidation-number description, direct hydrogens,
oxygen substitution and each oxygen's second shell. These graph descriptions
are not claims of supported source oxidation states.

FRC line 163 defines general `p` with four connections; line 164 defines
phosphazene `p=` with five. The template at line 1562 assigns `p=` to
`P(=O/S/N)(-*)(-*)(-*)`, which has four neighbors and bond-order sum five.
The general `(>P)` template does not independently establish chemical
precedence, oxidation-state scope or the neighboring oxygen types. Choosing
`p` merely because its neighbor count fits would silently override the more
specific conflicting template. Choosing `p=` would ignore its declared
connections. The aromatic phosphazene template does not establish PF5 typing.
No such reinterpretation is authorized.

Phosphate ester, phosphonate, phosphine oxide and the mixed Si/P fixture retain
the explicit oxo-template/connection conflict. Trivalent phosphine-like P has
no degree-three source atom definition. PF5 and charged P remain unsupported.
The C2 controls additionally distinguish trimethylphosphine, phosphoric acid
and phosphonium without selecting a P type or assigning neighboring charges.
The RLB contains `VERSION`/`elib`, not a missing chemical perception rule set.

## Executed independent comparison

The original converter is freshly compiled from the local source, with command,
source hashes, executable bytes and executable hash retained. C2 executes a
fresh C1/J41 study: the same final graphs, explicit H atoms, seed 2026,
coordinates and declared diagnostic perturbation. CAR/MDF atom labels and
charges are independently reconstructed before conversion. No ambiguous P,
Si–N/S/halogen labels or guessed charges are supplied to make conversion run.
`msi2lmp` consumes labels and charges; it does not perceive these chemical
environments. Successful assembly alone never qualifies a fixture.

The retained converter runs reproduce silane success, methylsilane/disiloxane/
hexamethyldisiloxane exit 17 and sila-oligomer exit 13, without `-ignore`.
The first missing end-bond–torsion keys are respectively `hc c3 si hsi`,
`hsi sio osi sio`, and `hc c3 sio c3`. The oligomer's `si c2 si` angle has
conflicting automatic candidates. Raw failure logs are preserved.

C2 also executes 317 distinct source searches through the pinned compiled C
leaf matcher: 79 assigned, 237 missing and one ambiguous key. Missing searches
agree on absence; the converter's ambiguous choice remains unauthorized. This supplements
the freshly repeated C1 audit of 4,359 raw numerical records and 347 direct
compiled queries. Search-table completeness does not imply model completeness.

The new adapter follows `find_angleangle_data`'s first-coefficient caller
semantics: try `a b c d`, then `d b c a`, while keeping the central and shared
sites fixed. An independently compiled regression reproduces the initial miss
for `si c h h` and finds the explicit zero row `angle-angle:cff91:4025` only
after the legal outer swap; moving the central site remains missing. This
corrects the new audit harness, not historical PCFF resolution or coefficients.
Converter file-order choices among conflicting candidates never override the
native ambiguity rejection.

Finite silane's identity remains
`a231146e09dbe71bd56baa4e3006cd97410aee8e297b6f8cd69f8dd15021486b`.
Fresh coefficient assembly, interaction counts, component energies, forces,
central finite differences, minimization, NVE/BAOAB finite software trajectories
and continuation are checked through the preserved J41 contract. The reference
bounds remain `atol=1e-5`, `rtol=2e-10`; minimization requires
`fmax <= 0.1 kJ/(mol Å)`. Whole/split final state and complete RNG equality are
distinct from checkpoint-byte equality for identical segmented lineage.
The finite silane has no active intramolecular nonbonded pairs under
`special_bonds lj/coul 0 0 1`; it does not qualify intermolecular Si forces.

## Reproduction and validation

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c2.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --converter /tmp/j41-c2-validation/converter/msi2lmp \
  --matcher /tmp/j41-c2-validation/matcher/leaf.so \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j41-c2-validation/final-study
```

A completed study exits **1** for the retained source/typing blockers. An
exception with the same exit code is not a completed scientific result: the
report, all assessments, raw audits and compiled searches must reconstruct.
The publication validator checks those artifacts as well as the receipt,
manifest, validation results and preserved file hashes:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_c2_receipt.py docs/evidence/phase_4j41_c2
```

Extract `execution.tar.xz` into another directory and pass its `study/` to
`scripts/phase_4j41_c2.py --reconstruct` for separate-process validation and
finite checkpoint/RNG reconstruction. The archives retain sources, build
commands, exact systems, converter input/output, decisions, complete source
searches and test commands/JUnit results. Test modules are partitioned into
disjoint pytest processes, with complete file lists and individual JUnit
reports retained; no test is replaced by a metadata-only check.
The supported inventory follows pytest's normal hidden-directory exclusions;
notebook checkpoint copies are excluded and left unchanged.

Verification completed with **35 focused passes**, **1,219 J22–J41-C1
regression passes**, and **3,222 supported full-suite passes with 10 existing
optional skips**. Ruff, compileall, all five applicable environment pip checks
and diff checks pass. The test commands, disjoint module lists and JUnit
reports are retained in the validation archive.

Independent raw reconstruction passes, all eight reserialized file-level
mutations are rejected, and relocation plus fresh-process checkpoint/RNG
reconstruction leaves all **226 study files byte-identical**. All **17 C1
assessment identities** and J41 silane coefficient, component-energy, force,
finite-difference and minimization results are unchanged. Both original J41
and C1 receipt validators pass. All **2,508 historical documentation/evidence
file hashes** remain unchanged, including the C1 receipts and archives.
The C2 receipt and manifest validator also passes with acceptance gate **1**.

Periodic `ewald_unvalidated`, `pressure_unavailable`, and blocked long-NPT
states are unchanged. No broader PCFF/COMPASS support, charged/ionic/metallic
coverage, AM1-BCC support, experimental agreement, scientific equilibrium or
production readiness is established.
