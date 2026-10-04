# Phase 4H4 — explicit LAMMPS-compatible Class II model conventions

## Base and scope

Branch `codex/phase-4h4-pcff-model-conventions` starts at reviewed 4H3
`2bc40c588e22fb69ae3e2d24718f9764a1969a86`. At preparation, `origin/main`
remained `a3a480b0a5ba587ade88a36d16c3ad9120178d18`; **4H3 is an unmerged
dependency**. No prerequisite merge, main change, historical record rewrite or
checkpoint-directory change was made.

This introduces `island_lammps_pcff_acyclic_cho_v1`, an explicitly chosen
**LAMMPS-compatible operational model**, not universal/commercial PCFF semantics.
The source and automatic domain remain unchanged: neutral finite nonperiodic,
explicit-H, saturated acyclic C/H/O, with the H2 `oh/ho` ordinary-alcohol convention.

- LAMMPS implementation revision: `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- FRC: `tools/msi2lmp/frc_files/pcff.frc`, SHA256
  `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
- H3 source records and their missing-row statuses remain unchanged.
- No `ParameterizedSystem`, production energy/force evaluator, optimization,
  dynamics, periodic electrostatics or general LAMMPS exporter is exposed.
- `production_validated=False`, `simulation_readiness="not_established"`.

## Public definition and persistence

```python
from island.forcefields.pcff import (
    load_pcff_source, load_pcff_parameters, special_pair_policy,
    define_pcff_model, save_pcff_model, load_pcff_model,
)
source = load_pcff_source("/local/pinned/pcff.frc")
assignment = load_pcff_parameters("/local/H3/ethanol.json", source)
# An explicit caller choice, NOT weights recovered from the FRC:
weights = special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
model = define_pcff_model(assignment, special_pairs=weights)
print(model.payload["raw_parameter_coverage_complete"])
print(model.payload["model_definition_complete"])
save_pcff_model(model, "/new/ethanol-model.json")
restored = load_pcff_model("/new/ethanol-model.json", assignment)
restored.validate_integrity()
```

`PCFFModelSpecification` schema `island_pcff_class2_model_v1` references the
external H3 assignment identity. It stores source/typing/charge/graph identities,
full compatibility profile and equations, compact numerical terms, dependency
and source-row evidence, nonbonded records and shortest-path special-pair
inventory. Each term identifies its origin as `source_row` or
`policy_derived_zero`. The latter never has a fabricated source row.

Validation revalidates the original H3 assignment, then recomputes model policies,
terms, dependencies, pair inventory and completeness. Changes to equations,
weights, charges, values, dependencies or applicability change the identity.
Rechecksummed contradictions fail with `PCFFError`. The existing strict JSON
and exclusive atomic publication machinery is reused. Loaded records require
the compatible original assignment/source, but no scientific backend execution.
These checks establish consistency, not computational authenticity.

Payloads are owned copies. `model.numerical_terms()` validates at the boundary
and returns immutable **bonded** `Class2Term` objects containing only small
coefficient/role/equilibrium tuples. Their repeated energy/finite-difference
calls neither reparse FRC nor reconstruct H3 catalogs. Nonbonded records remain
separate. This is a validation interface, not a whole-system simulation API.
Only NumPy/core imports are required; RDKit, OpenMM, Foyer, ParmEd, SciPy and
AmberTools are not needed for loading, definition or term kernels.

## Equations, roles and units

The numerical convention is **kJ/mol, angstrom, radians and elementary charge**.
H3's conversions are preserved: energy-bearing coefficients multiply by 4.184;
lengths do not change; only equilibrium angles/phases multiply by pi/180.
Angular coefficients act on radian displacements. No one-half factor is added.

Let `d_r = r-r0`, `d_theta = theta-theta0`, and `C_n = cos(n phi)`.

| Family | Operational equation |
|---|---|
| Quartic bond | `K2*d_r² + K3*d_r³ + K4*d_r⁴` |
| Quartic angle | `K2*d_theta² + K3*d_theta³ + K4*d_theta⁴` |
| Adjacent bond–bond | `K*d_r12*d_r23` |
| Bond–angle | `(Kleft*d_r12 + Kright*d_r23)*d_theta123` |
| Proper torsion | `sum_n Vn*(1-cos(n*phi-phase_n))`, n=1,2,3 |
| Middle-bond–torsion | `d_r23*sum_n Fn*C_n` |
| End-bond–torsion | `d_r12*sum_n Fleft_n*C_n + d_r34*sum_n Fright_n*C_n` |
| Angle–torsion | `d_theta123*sum_n Fleft_n*C_n + d_theta234*sum_n Fright_n*C_n` |
| Angle–angle–torsion | `K*d_theta123*d_theta234*cos(phi)` |
| Bond–bond 1–3 | `K*d_r12*d_r34` |
| Angle–angle | `K*d_theta123*d_theta324`, site 2 central, site 3 shared arm |

Coefficient ordering is H3 order: quartic equilibrium first, then powers 2–4;
torsion `(V1,phase1,V2,phase2,V3,phase3)`; directional cross terms are the three
left coefficients followed by three right coefficients. Equilibrium arrays
follow the corresponding physical roles. AAT now explicitly has energy/rad²
units in this model; H3's unresolved-unit annotation is **not rewritten**.

### Signed torsion geometry and the FRC discrepancies

The kernel uses projected vectors, an independent expression of the geometry
in pinned `src/CLASS2/dihedral_class2.cpp`:

```
v1 = x1-x2; t = (x3-x2)/|x3-x2|; v3 = x4-x3
a = v1-(v1.t)*t; b = v3-(v3.t)*t
phi = atan2((t cross a).b, a.b)
```

Cis is zero; trans is pi. Full proper reversal preserves signed phi, and exchanges
left/right coefficient and equilibrium blocks. No coefficient negation or phase
shift is applied. Nonzero phases and asymmetric blocks were exercised against
the executable, so a sign error cannot be hidden by zero-phase fixtures.

The pinned executable computes `1-cos(n*phi-phase)` and AAT's `cos(phi)`.
These are **explicit interpretation choices** over the FRC's conflicting
plus-cos and linear-Phi comments. Both discrepancies remain in model provenance.
The named operational compatibility target resolves the choice; it does not
establish what every historical PCFF implementation intended.

Primary evidence: [pinned executable torsion code](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/dihedral_class2.cpp),
[pinned angle code](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/angle_class2.cpp),
and [class2 dihedral documentation](https://docs.lammps.org/dihedral_class2.html).
Actual comparison uses the pinned compiled code, not a moving documentation page.

### BB1–3, Wilson, angle–angle and torsion–torsion scope

Pinned `tools/msi2lmp/src/GetParameters.c` initializes BB1–3 coefficients to zero,
sets the two equilibrium distances, and only looks up a source coefficient when
one of the **supplied** torsion types is exactly `cp`. The supported automatic
CHO profile never assigns `cp`. `msi2lmp_non_cp_zero_v1` therefore produces a
**policy-derived zero** for these propers, only with resolved equilibrium
bond dependencies. This reproduces the audited converter behavior, without
claiming a missing FRC row was a zero-valued row. No unrelated missing,
ambiguous or unresolved assignment is filled. Incomplete definitions retain
explicit diagnostics and cannot provide numerical terms.

Wilson terms stay structurally outside this profile: C is four-connected, O
two-connected and H one-connected. Angle–angle remains active: each unordered
neighbor triple yields three distinct angle pairs. A tetrahedral center has
four triples and twelve couplings. Site 2 is always central; swapping all four
sites as an unordered improper is invalid. The standalone five-site executable
fixture verifies the twelve-coupling energy and force sum.

Torsion–torsion is **excluded by the operational model's declared scope**. An
empty FRC section is not offered as proof of physical inapplicability.

Evidence: [pinned converter](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c),
[pinned inventory construction](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/MakeLists.c),
and [pinned improper implementation](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/improper_class2.cpp).

### Nonbonded definition and explicit special-pair choice

```
rmin_ij = ((ri^6+rj^6)/2)^(1/6)
epsilon_ij = 2*sqrt(ei*ej)*ri^3*rj^3/(ri^6+rj^6)
U_LJ = w_LJ * epsilon_ij * [2*(rmin_ij/r)^9 - 3*(rmin_ij/r)^6]
U_C  = w_C * (332.06371*4.184) * qi*qj/r
```

`rmin` is the minimum-energy distance, not a 12–6 sigma. The electrostatic
constant is the pinned `units real` value from `src/update.cpp`, converted to
kJ·angstrom/(mol·e²). Dielectric is one. Zero epsilon does not suppress Coulomb.

LJ and Coulomb weights are independent explicit triples for shortest-path
1–2, 1–3 and 1–4 pairs. No default is provided by the model API. Other pairs
have weight one. The source does not establish a trustworthy special-pair
choice here; weights are **caller model choices**, not source facts. The seven
chemical definitions use `(0,0,1)` for each; term verification deliberately
also uses LJ `(.25,.5,.75)` versus Coulomb `(.2,.4,.6)` to expose accidental
coupling or inherited defaults.

The intended finite model has no cutoff, switching, shift, long-range correction
or reciprocal electrostatics. The isolated executable oracle uses a 30 Å cutoff
in a fixed 100 Å box, with every tested pair far inside the cutoff, no shifting,
no tails and no periodic boundaries. This is equivalent for these fixtures;
it is not a proposed finite-cutoff production model.

Evidence: [pinned pair implementation](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/CLASS2/pair_lj_class2_coul_cut.cpp),
[pinned unit constants](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/src/update.cpp),
and [pair class2 documentation](https://docs.lammps.org/pair_class2.html).

## Executable verification and reproduction

No preinstalled `lmp` was found. A separate upstream checkout was built at the
exact pinned commit, with a clean source working tree:

```bash
# External directory; no vendored upstream implementation/data in ISLAND.
git clone --filter=blob:none --no-checkout https://github.com/lammps/lammps.git \
  ../island-validation/phase4h4-upstream/lammps
git -C ../island-validation/phase4h4-upstream/lammps sparse-checkout init --cone
git -C ../island-validation/phase4h4-upstream/lammps sparse-checkout set \
  src cmake tools/msi2lmp doc/src
git -C ../island-validation/phase4h4-upstream/lammps checkout \
  e891a3e10973c1a729e391a0aefaa02fd70f8c0f
cmake -S ../island-validation/phase4h4-upstream/lammps/cmake \
  -B ../island-validation/phase4h4-upstream/build -D BUILD_MPI=off \
  -D BUILD_OMP=off -D PKG_CLASS2=on -D PKG_MOLECULE=on \
  -D CMAKE_BUILD_TYPE=Release
cmake --build ../island-validation/phase4h4-upstream/build -j 4
```

The executable reports `30 Sep 2026 - Development`, Git
`patch_30Sep2026-65-ge891a3e10`. Compiler, configuration, executable/source hashes
and retained paths are in `docs/evidence/phase_4h4.json`. This is an actual pinned
executable comparison, not merely a pinned-source inspection.

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY scripts/validate_pcff_terms.py \
  --lammps ../island-validation/phase4h4-upstream/build/lmp \
  --output ../island-validation/phase4h4-term-acceptance
$PY scripts/validate_pcff_models.py --source "$FRC" \
  --assignments ../island-validation/phase4h3-acceptance-final \
  --output ../island-validation/phase4h4-model-acceptance \
  --lj 0 0 1 --coulomb 0 0 1
$PY examples/pcff_model_specification.py --source "$FRC" \
  --assignment ../island-validation/phase4h3-acceptance-final/ethanol.json \
  --output /new/ethanol-model.json --lj 0 0 1 --coulomb 0 0 1
```

Each acceptance command publishes a declaration before calculations, refuses an
existing destination and records durable outcomes. Term CLI exit 0 requires all
24 declared comparisons; missing executables and failed comparisons return
nonzero. Inputs, stdout/stderr, logs, energies and force dumps are retained and
checksummed. No retry, seed search, tolerance increase or coefficient change was
needed. LAMMPS inputs are handwritten isolated-term constructions from synthetic
fixtures; they do not call H3/model conversion or pair-inventory helpers. This
harness is not a general-purpose exporter.

The energy kernels use kJ/mol; LAMMPS real-unit energies and forces are multiplied
by 4.184 independently. Central differences use the **negative** energy gradient
and displacements `1e-3, 3e-4, 1e-4` Å. Predeclared tolerances:

- Energy: absolute `1e-10` kJ/mol, relative `1e-12`.
- Force: absolute `3e-5` kJ/(mol·Å), relative `1e-7` at the finest displacement.
- Finite-difference convergence: finest error < 0.2 of coarsest error, or below
  an absolute `1e-8` roundoff floor. Exact-zero fixtures are allowed.

All 24 comparisons passed. Maximum energy discrepancy was
**1.7763568394002505e-15 kJ/mol**. Maximum force errors at the three displacements
were **5.6780245618615055e-6**, **5.11023642291164e-7**, and
**5.678181280721617e-8 kJ/(mol·Å)**, demonstrating second-order convergence.
The test suite also contains small retained numerical outputs from the actual
executable (`tests/fixtures/pcff_class2_terms.json`). Those are synthetic software
fixtures, not independently validated chemical PCFF parameters.

## Seven-case definition and integrity results

| Case | Source-derived bonded terms | Policy-derived BB1–3 zeros | Special pairs | Model bytes | Audit time (s) |
|---|---:|---:|---:|---:|---:|
| butane | 268 | 27 | 64 | 268774 | 15.227 |
| neopentane | 346 | 36 | 82 | 332631 | 16.356 |
| ethanol | 131 | 12 | 33 | 129013 | 13.450 |
| dimethyl_ether | 101 | 6 | 27 | 97183 | 12.809 |
| PE_DP3 | 424 | 45 | 100 | 423714 | 17.674 |
| PEO_DP3 | 421 | 42 | 103 | 403603 | 18.142 |
| PE_DP50 | 7756 | 891 | 1792 | 7814135 | 178.850 |

All seven definitions are complete under the explicitly named policy. All seven
H3 records still have `raw_parameter_coverage_complete=False`. These outcomes
are definition/integrity checks; no molecular energy or force was calculated.
PE DP50 is a 302-site construction/integrity smoke test, not dynamics evidence.
Times include loading/revalidating H3, model definition, repeated boundary
validation, identity calculation and persistence; they are not kernel timings.

The retained H1/H2/H3 file hashes were all rechecked: **5, 42 and 9 files**,
respectively, unchanged. Original source bytes, charges and assignment identities
were preserved. Full source/data remain external; only small generated synthetic
numerical fixtures and evidence manifests are committed.

## Verification and remaining gates

- Full ordinary suite: **1,193 passed, 10 expected optional skips**.
- Focused model/kernel regressions: **13 passed**, including the final strict
  stable-ID input check, source versus policy zeros, missing-dependency blocking,
  rechecksummed contradictions, policy identities, exclusive persistence,
  optional-dependency isolation, reversal/remapping and retained executable data.
- Ruff and `python -m pip check`: passed.
- New model specification example and existing H3 parameter example: executed,
  exit 0.
- Actual pinned LAMMPS: **24/24 isolated comparisons passed**; no mocked result
  stands in for this gate. All inputs/outputs and checksums are retained.
- Seven historical-source model definitions: **7/7 complete** with the declared
  caller pair policy; no historical raw coverage was relabeled as complete.

**Ready for the next bounded single-point implementation phase:** the selected
operational equations, applicability and caller-selected pair policy now have
executable term-level evidence. The next phase must implement complete system
aggregation and analytic/OpenMM forces, independently verify inventories and
special-pair application, and compare molecular energies/forces at varied
geometries against a separately constructed pinned reference. Term agreement
does not establish whole-system implementation fidelity or scientific accuracy.

No universal PCFF special-pair rule was recovered. Callers must continue to
select it explicitly. Wilson and torsion–torsion extensions, unsupported chemical
environments, periodic models and scientific validation remain outside this
profile. No simulation-ready snapshot or production workflow is exposed.
