# Phase 4H5 — finite PCFF single points

## Repository and scope

Branch: `codex/phase-4h5-pcff-single-point`. Base:
`159774a2141dd4b6481e089470c92e810e499e41` (`origin/main`, merged H3/H4).
Its tree is identical to reviewed H4 `436d940ed06fed662c034bf40d3aa4d1aa3d1f32`:
`d28bdfea1db7231d9cbfc630890a0fc9cf81d5be`. No prerequisite merge was needed.

The new evaluator implements **island_lammps_pcff_acyclic_cho_v1**, preserving
H1–H4 records and interpretation. It supports the existing neutral, saturated,
acyclic, explicit-H C/H/O domain, including the documented ordinary-alcohol
`oh`/`ho` choice. It does not establish universal/commercial PCFF equivalence.
No optimizer, dynamics, reusable session, periodic model, or simulation snapshot
is added. `production_validated=False`, `simulation_readiness="not_established"`.

## API, ownership, and binding

```python
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import define_pcff_model, special_pair_policy

specification = define_pcff_model(
    validated_assignment,
    special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
)
evaluator = PCFFSinglePointEvaluator(system, specification)
result = evaluator.evaluate()  # owned default coordinates
other = evaluator.evaluate_fresh(coordinates, coordinate_unit="angstrom")
```

The constructor validates the complete specification against the authoritative
system, owns the chemical binding/default coordinates, and compiles owned
numerical records into private OpenMM XML. Evaluation neither parses FRC nor
rebuilds the H3 catalog. Later caller mutations cannot change the bound model.

`OpenMMBoundPotential` supplies the existing binding boundary:
`validate_system()` rejects changed IDs, elements, masses, bonds, formal charge,
periodicity, representation, or chemical stereo/isotope/radical metadata.
Coordinate changes and unrelated coordinate-provenance metadata are permitted.
Contract compatibility is not a claim that PCFF minimization/dynamics have been
accepted; those remain later work.

Every `evaluate()` and `evaluate_fresh()` deserializes its own System and creates
its own Integrator/Reference Context. The integrator is never stepped. Complete
positions are installed on each call. Resources are disposed deterministically;
cleanup failure does not obscure the original backend exception. No mutable
backend object is exposed through the evaluator API.

Outputs use the existing owned `EvaluationResult`: kJ/mol energies,
kJ/(mol·angstrom) stable-ID forces, all 13 components, coordinate/model/parameter/
calculation fingerprints, and backend/version/platform identity. The parameter
fingerprint binds the complete H4 model, including caller-selected pair weights.

Raw coordinate mappings require exact integer site coverage and finite real
nonboolean triples; NumPy booleans are rejected before coercion. Only explicit
angstrom input is accepted. Coincident sites and singular active angular
geometries, including overflow during geometry checks, raise
`EvaluationInputError` before Context construction. The coincidence threshold is
1e-10 Å; the cross-product angular threshold is 1e-12 times the product of bond
lengths. These are numerical boundaries, not changed force-field parameters.
Malformed/incomplete specifications also raise `EvaluationInputError`.
Backend failure/nonfinite energy, any component, or any force raises
`EvaluationError`, never a partial result. OpenMM remains a lazy optional import.

## Equations and units

The equations and physical roles remain those in [H4](phase_4h4.md). OpenMM
custom forces analytically differentiate all terms; production forces never use
finite differences. Backend distances are multiplied by 10 **inside** the
expressions, so all stored angstrom-based coefficients remain unchanged.
Backend angular geometry is already in radians. Returned kJ/(mol·nm) forces are
multiplied by 0.1.

| Component | Compilation convention |
|---|---|
| Quartic bond | K2·dr² + K3·dr³ + K4·dr⁴; dr in Å |
| Quartic angle | K2·dθ² + K3·dθ³ + K4·dθ⁴; dθ in radians |
| Bond–bond / bond–angle | Preserve signed coefficients, left/right roles, and Å/radian dependencies |
| Proper torsion | Σ Vn·[1−cos(nφ−phase_n)], n=1,2,3 |
| Middle/end bond–torsion | H4 cosine blocks times the corresponding Å displacement |
| Angle–torsion | Directional cosine blocks times the corresponding radian displacement |
| Angle–angle–torsion | K·dθ_left·dθ_right·cos(φ) |
| Bond–bond 1–3 | K·dr_left·dr_right, including explicitly marked H4 policy zeros |
| Angle–angle | Each H4 coupling once; three couplings per neighbor triple, twelve per tetrahedral center |
| LJ 9–6 | εij·[2(rmin,ij/r)^9−3(rmin,ij/r)^6] |
| Coulomb | 138.935456264·qi·qj/r_nm kJ/mol |

OpenMM's signed `dihedral(p1,p2,p3,p4)` agrees with H4's projected-vector
atan2 convention (cis=0, trans=π), including a nonzero-phase test. Coefficients
are not negated or phase-shifted. Angular polynomial/cross coefficients are not
converted merely because source equilibrium angles were written in degrees;
H3/H4 already store those equilibrium angles/phases in radians.

LJ is a separate 9–6 custom force, not OpenMM's built-in 12–6 nonbonded force.
For each unordered pair, once:

- rmin,ij = ((ri⁶+rj⁶)/2)^(1/6);
- εij = 2√(εiεj)·ri³rj³/(ri⁶+rj⁶);
- apply the independent LJ and Coulomb weights from H4's shortest-path pair list;
- absent special-pair entries, including disconnected pairs, have unit weights.

Zero epsilon never deletes Coulomb. The Coulomb constant is precisely the H4
332.06371×4.184 kJ·Å/(mol·e²) convention, expressed in nm above. There is no
cutoff, switching, shift, tail, periodic image, or reciprocal term. Pair
compilation and geometry checks are O(N²); the bounded large case was 302 sites.

## Independent reference and a converter discrepancy

Source pins remain:

- LAMMPS `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`;
- `pcff.frc` SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

The validation-only CLI writes CAR/MDF from authoritative bonds and independently
authored H2 SMARTS labels/native-charge row sums. The pinned `msi2lmp` executable
then regenerates bonds, angles, proper torsions, tetrahedral neighbor triples,
source selection, equilibrium dependencies, and cross coefficients from raw
FRC. It does not receive ISLAND's compiled interaction lists. There is no
`-ignore`/missing-parameter zero fallback. Verbosity 3 is deliberate: the pinned
converter's AA equivalence lookup is nested in that verbosity condition.
Independent graph combinatorics check every converted interaction and charge;
the existing separate source-line/Decimal oracle checks H3 selected rows and
normalized coefficients. Inventories and force comparisons therefore do not
share the production compiler.

**Raw converter output is not used blindly.** Inspection and execution exposed
an upstream equilibrium-column discrepancy:

- `tools/msi2lmp/src/GetParameters.c` writes AA angles `ABC, CBD, ABD`;
- `src/CLASS2/improper_class2.cpp` consumes `ABC, ABD, CBD`;
- the K1/K2/K3 coefficient roles themselves agree.

The reference independently reconstructs these three dependencies from its
converted angle inventory, verifies the raw converter ordering, and supplies
explicit `improper_coeff ... aa` overrides with `ABC, ABD, CBD`. Raw data files
are retained unchanged. These overrides use neither H3 dependencies nor the
production model helper. A regression rejects unexpected converter behavior.
For the retained butane baseline, raw-versus-corrected reference energy differs
by **0.03624289956 kJ/mol**, entirely in AA, and maximum force difference is
**0.7263134762 kJ/(mol·Å)**. This is a documented reference adapter correction,
not a change to H4 records or a tolerance relaxation. Source file hashes and
both inputs/outputs are in the evidence manifest.

The LAMMPS executable is the actual H4 build of the pinned revision (30 Sep 2026
Development; GCC 13.3, Release, CLASS2/MOLECULE, MPI/OMP disabled). SHA256:
`ee945ead72f4997905ea42636bc987f3cc4a2a9cf89ab2af30d8e02c4c1cfebc`.
`msi2lmp` was compiled separately from the same unchanged checkout. Its command,
binary hash, source hashes, and version 3.9.11 are retained.

LAMMPS uses finite `boundary f f f`, `lj/class2/coul/cut`, an all-pair-containing
cutoff, no shift/tail, and explicit `special_bonds`. Full-precision declared
coordinates are restored after converter formatting/recentering. LAMMPS real
kcal/mol and kcal/(mol·Å) are multiplied by 4.184. Its grouped components are:
`ebond`=quartic bond; `eangle`=angle+BB+BA; `edihed`=proper+all proper cross terms;
`eimp`=AA; `evdwl`=9–6; `ecoul`=Coulomb. Wilson terms remain structurally excluded.

## Declared and executed acceptance

A separate declaration was published before each matrix execution. Seed 2026;
ETKDG small molecules, local templates with template/assembly seed 2026 for
polymers; no minimization or dynamics. Six small cases also use the fixed
perturbation 0.015·sin(arange(3N)+0.3) Å. Main pair weights are LJ=Coulomb=(0,0,1).
Ethanol additionally uses LJ=(0,0.2,0.5), Coulomb=(0,0.3,0.7), for both frames.
Tolerance: energy 1e-5 kJ/mol, force 1e-5 kJ/(mol·Å), rtol=2e-10. No seed,
parameter, or tolerance changed following a numerical comparison.

| Case | Frames/policies | Maximum energy discrepancy, kJ/mol | Maximum force discrepancy, kJ/(mol·Å) |
|---|---:|---:|---:|
| Butane | 2 | 1.35e-13 | 1.52e-12 |
| Neopentane | 2 | 8.53e-14 | 1.49e-12 |
| Ethanol | 4 | 5.68e-14 | 8.53e-13 |
| Dimethyl ether | 2 | 2.84e-14 | 1.59e-12 |
| PE DP3 | 2 | 4.62e-14 | 1.95e-12 |
| PEO DP3 | 2 | 2.98e-13 | 2.16e-12 |
| PE DP50 | 1 | 1.05e-9 | 4.86e-10 |

All 15 comparisons passed totals, grouped components, and every force component.
Coefficient checking maximum was below 5e-13 in normalized units. PE DP50 has
301 bonds, 600 angles, 891 propers, 400 neighbor triples / 1,200 AA couplings.
Its construction/assignment/source checking took 63.4 s, evaluator boundary
validation/compilation 32.7 s, and fresh evaluation including cleanup 0.734 s.
These are individual measurements, not speedup benchmarks. Large signed record
validation/compilation dominates this single-point workflow, not force execution.

Artifacts: `../island-validation/phase4h5-whole-system-final/`; hashes, full
per-frame measurements, commands, and provenance in
[evidence/phase_4h5.json](evidence/phase_4h5.json). The initial CAR/MDF harness
attempt failed seven cases because of a missing required blank separator; its
exit-4 logs/report remain in `phase4h5-whole-system/`. The converter-discrepancy
baseline remains in `phase4h5-harness-check/`. Neither is relabeled as passed.
Original H3 (9 files), H4 model (9 files), and H4 term (171 files including nested
artifacts) hashes were rechecked and unchanged.

## Reproduction

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
UPSTREAM=../island-validation/phase4h4-upstream/lammps
mkdir -p ../island-validation/phase4h5-tools
gcc -O2 "$UPSTREAM"/tools/msi2lmp/src/*.c -lm \
  -o ../island-validation/phase4h5-tools/msi2lmp
$PY scripts/validate_pcff_singlepoint.py --source "$FRC" \
  --converter ../island-validation/phase4h5-tools/msi2lmp \
  --lammps ../island-validation/phase4h4-upstream/build/lmp \
  --output /new/exclusive/acceptance-directory
$PY examples/pcff_singlepoint.py --source "$FRC" --lj 0 0 1 --coulomb 0 0 1
```

The acceptance command refuses an existing output directory and exits nonzero
for any failed required case. FRC and executables must already be explicitly
installed/supplied. No runtime download. Production evaluation needs NumPy and
OpenMM; construction/reference labels additionally need RDKit. The reference
executable/compiler are acceptance-only dependencies. Actual environment:
Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1, RDKit 2026.03.6.

## Verification and remaining limits

The focused suite includes retained real LAMMPS isolated-term outputs for all
11 bonded families; synthetic nonzero-phase and unequal-sigma/zero-LJ fixtures;
finite-difference convergence at 1e-3, 1e-4, 1e-5 Å; remapped/reordered IDs;
translation/rotation covariance; ownership and graph/mass/stereo rejection;
malformed inputs and rechecksummed model tampering; fresh Context counts; no
integrator stepping; injected nonfinite backend output and cleanup failures;
and optional-import isolation. Synthetic fixtures test software, not chemistry.

Final ordinary suite: **1,221 passed, 10 opt-in skipped** (156.32 s); focused
regressions: **28 passed**. Ruff and `python -m pip check` passed. The skips
are the existing Amber archive/QM/workflow and pinned Foyer opt-ins; PCFF
real-source acceptance was executed separately. Additional post-run checks
verified exact native labels/charges, independently enumerated special-pair
walks, and BB13 policy provenance for all eight model variants (including
1,792 special pairs and 891 policy zeros in PE50). Their results are retained
in `supplemental-integrity.json`; original records were not rewritten. The new
single-point, H4 model-specification, and H3 parameter examples were executed.
No QM, minimization, dynamics, source edits, or historical schema changes.

Within the named operational model, the whole-system numerical gate is met and
a **bounded independently verified minimization phase can begin**. It must add
its own convergence, budget, failure, application, and final-reference evidence.
Reusable sessions remain deferred. No scientific transferability, commercial
PCFF equivalence, production readiness, or equilibration follows from these
single-point comparisons. Bare pinned `msi2lmp` output retains the documented AA
ordering discrepancy and must not be mistaken for the role-corrected oracle.
