# Phase 4J41-C4 — full msi2lmp resolution and execution boundary

No new mixed Si/P execution profile is published. Real compiled conversions
still fail on required source couplings. C4 reconstructs the converter's actual
selection and assembly behavior, independently of the historical ISLAND
fallback resolver. The full-source acceptance gate remains **1**.

`production_validated=False`; `simulation_readiness="not_established"`.

The base was verified after `git fetch origin`: merged C3 on `origin/main`,
`7e1150a443b72d518474dc050a02a1192bc486a0` (merge PR #66). C3's final evidence
commit `8d04fa4e779cc1587ef3e71b743e2f40184880ec` is its ancestor. Work uses
`codex/phase-4j41-c4-executable-msi2lmp-resolution`, without modifying main.
The tested implementation commit is recorded in the signed C4 receipt.

## Authorized sources and input boundary

Only the authorized local `tools/msi2lmp` distribution is read. PCFF and the
separately audited published COMPASS remain isolated; no source row is copied
between them. PCFF's FRC SHA-256 remains
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
Its `pcff.rlb` and `pcff_templates.dat` retain their original hashes.
Published COMPASS retains SHA-256
`28156c9c95346a9efde17d0ae0aff8053410e1de86b05f8d284ef781c0bbadfe` and its own
embedded templates, without borrowing PCFF companion definitions.

The public PSMILES/final-graph workflow is unchanged. Audit graphs are created
with the chemistry API or reconstructed from retained final-graph records.
CAR/MDF are generated internally only to query the independent converter. No
user CAR/MDF files, proprietary downloads or network force-field searches are
required. Source inventories retain atom declarations, ordinary and automatic
equivalences, charge increments, templates, complete raw numerical rows,
companion hashes and source identities. Inventory or typing success confers
no execution authority.

## Reconstructed converter behavior

The audit compiles the original local converter and a separate diagnostic copy
with additional writes inside `find_match`. Every attempted real converter
lookup is reconstructed from raw FRC rows. Original and diagnostic builds must
have identical exit codes, stdout and stderr after removing trace lines, and
byte-identical coefficient files for successful conversions. Both builds,
original C/H source bytes and hashes, compiler/version, commands, revision,
source diff, executable hashes, LAMMPS hash/help, Python/OpenMM versions and
hardware are archived. This is a compiled converter comparison, not a claim
that its supplied type labels are chemical authority.

The following rules are directly visible in the retained source:

- `ReadCarFile` consumes type labels. `ReadMdfFile` reads connectivity and
  replaces CAR charges with MDF charges. The converter does not perceive
  atom types or compute charge increments.
- `ReadFrcFile`/`InitializeItems` load ordinary Class-II tables. They do not
  load `auto_equivalence`, automatic lower-order parameter tables or charge
  increments. ISLAND's historical automatic resolver remains separate.
- `SearchAndFill` uses the first matching section. For an identical ordered
  key it retains the original table slot, replacing values only for a strictly
  larger C-float version; an equal version retains the first entry.
- `find_match` scans all exact matches before its wildcard pass. Each pass
  uses file order and tries forward before reversed order when allowed.
  Wildcards are not ranked by specificity. Conflicting wildcard selections
  remain ambiguous for physical qualification.
- Ordinary equivalence is a single family-specific mapping, attempted after
  direct lookup fails. Wilson permutations keep the second-site center fixed.
  Angle-angle lookup searches three distinct shared-angle coefficients, each
  with two role-preserving orientations.
- The angle-angle caller accepts any of three source slots and leaves other
  slots at initialized zero. This is weaker than completeness. C4 preserves
  every missing applicable slot; converter success cannot fill it.
- `quo_cp` suppresses `bond-bond_1_3` lookup outside cp-prefixed types. The
  resulting initialized zero is recorded as `converter_zero_initialization`,
  together with a separate physical source lookup. It is never published as
  a coefficient.
- `get_r0` and `get_t0` copy equilibrium coefficient slot 0 from already
  resolved graph bond/angle types. Their exact dependencies are recorded.
  `SearchAndFill` repeats abbreviated symmetric BA/EBT/AT vectors cyclically.
  Neither operation derives a missing coupling constant.

Diagnostic states distinguish selected, absent and ambiguous rows, deterministic
converter-derived equilibrium dependencies, zero initialization, missing charge
increments, missing couplings and converter failure. Top-level blocked outcomes
use `typing_blocked`, `source_incomplete`, `parameter_incomplete`,
`converter_unavailable` or `scientifically_unqualified`. A Python exception is
not a completed acceptance result. Missing authorized sources or tools produce
an explicit classified nonzero report.

## Executed targets and chemical boundary

The priority targets are examined separately: methylsilane, Si–O–Si, then the
narrow neutral cyclic phosphazene probe. The larger siloxane and phosphorus
rejections are retained as negative controls.

| Source | Target | Actual first converter failure |
|---|---|---|
| PCFF | Methylsilane | exit 17: end-bond–torsion `hc c3 si hsi` |
| PCFF | Disiloxane | exit 17: end-bond–torsion `hsi sio osi sio` |
| PCFF | Hexamethyldisiloxane | exit 17: end-bond–torsion `hc c3 sio c3` |
| Published COMPASS | Methylsilane | exit 17: end-bond–torsion `h1 c4 si4 h1` |
| Published COMPASS | Disiloxane | exit 17: end-bond–torsion `h1 si4 o2z si4` |
| Published COMPASS | Hexamethyldisiloxane | exit 17: end-bond–torsion `h1 c4 si4c c4` |
| Published COMPASS | Neutral phosphazene probe | exit 15: bond–angle `n2= p4= n2=` |

All ordered site tuples, type tuples, raw selected rows, equivalence decisions,
missing lookup keys, charge-increment contributions, exact decimal charge
vectors and component conservation results are archived. No `-ignore` run is
used. No source values, charges, graph, force-field semantics or tolerance is
changed. Conversion failures are retained, and blocked-domain dynamics do not
run. There is consequently no C4 runtime resolver or execution profile.

PCFF's P connection/template conflict remains unchanged. Published COMPASS's
`p4=` phosphazene evidence does not authorize phosphate, phosphonate or oxo-P
mapping. Phosphate ester, phosphonate, phosphine oxide, protonated phosphorus,
PF5-like, unresolved mixed Si/P and unsupported Si-valence graphs remain
`typing_blocked`. They receive no guessed type or charge vector and are not
sent to the converter. `p`, `p=`, `p4=` and broad `>P` templates are never
transferred between these chemical states.

Finite silane is a retained control, separately checked for both sources with
source-derived energy components and central finite-difference forces against
LAMMPS. The existing PCFF execution profile and model identity
`a231146e09dbe71bd56baa4e3006cd97410aee8e297b6f8cd69f8dd15021486b` remain
unchanged. The original historical minimization, exact checkpoint and RNG
continuation evidence is preserved and reconstructed in a separate process.
The COMPASS silane diagnostic establishes no new transferable profile.

## Reproduction and validation

```bash
git fetch origin
git rev-parse origin/main
git merge-base --is-ancestor 8d04fa4e779cc1587ef3e71b743e2f40184880ec origin/main

PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c4.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --source ../lammps/lammps/tools/msi2lmp/src \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j41-c4-validation/final-study
```

The complete full-source study exits **1**, with `report.json` classifying all
results. The compiler is invoked internally; original commands and outputs are
retained. Extract `execution.tar.xz` into a fresh location, then run:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c4.py \
  --reconstruct /path/to/relocated/study

PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/pcff_j41_c4_controls.py \
  --study /path/to/relocated/study --output /tmp/c4-negative-controls

PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_c4_receipt.py docs/evidence/phase_4j41_c4
```

Fresh-process reconstruction reruns original converters at the relocated paths,
checks coefficient/failure-key equality, independently reconstructs raw source
rows and dependencies, and reconstructs retained finite-silane checkpoint/RNG
continuation. It does not write into the original or relocated study.
Negative controls reject changed source bytes, companion bytes, mixed source,
graph identity, atom order, charge vector, dependency set, missing coupling,
inserted zero, invented profile and inappropriate P label in separate processes.
Every incomplete/ambiguous target also has its own classified graph-bound
failure evidence.

Exact focused, J22–J41-C3 regression and supported full-suite commands and JUnit
results are retained in `validation.tar.xz`, with Ruff, compileall, all five
applicable pip checks and `git diff --check`. Receipt/manifest validation checks
archive membership, source/implementation hashes, numerical reconstructions,
mutation controls and historical byte hashes. Historical J41/C1/C2/C3 validators
are rerun. The signed receipt records all completed counts and commands.

Polymer `ewald_unvalidated`, `pressure_unavailable`, blocked cutoff recovery
and blocked long NPT remain unchanged. No periodic scientific closure,
scientific equilibrium, experimental agreement or production readiness is
claimed.
