# Phase 4J60 — executable Class-II residual and executed identifiable pilot fit

J60 binds an executable research Hamiltonian, executes **87 real Psi4 energy/
full-gradient calculations**, and fits all twenty joint coefficients. Corrected
basis and numerical gates pass. The actual fit fails untouched-conformer energy,
force and uncertainty criteria, so **zero coefficients are accepted**. The
remaining numerical obstacle is a method-matched Si–O/Si–H baseline constraint,
not an unbound engine setting or a missing executable.

J60 starts from fetched `origin/main`
`b4aecb4c601ea4b7d24cec9eba17d92747280297`, which contains J59
`9681e4926d92f5a38ea9dbc329c15c1207d819cd`. The branch is
`codex/phase-4j60-pcff-hamiltonian-identifiable-pilot-fit`. No merge or push to
main is performed. The complete base tree is the preservation baseline.

The scope is the exact **9-atom Si–H disiloxane control**, graph
`7aa04b55c75e48d32abd159a51e6f6e39fe51244521118c94df49484ba1a738f`,
atom order Si/O/Si/H/H/H/H/H/H, neutral singlet. It does not establish
methyl-PDMS transfer. The ordinary `forcefield="pcff"` API remains unchanged.
No final PCFF-PDMS profile is registered and no PCFF MD is run.

## Exact selected Hamiltonian

The separately identified research model is
`island_research_class2_sih_disiloxane_j60_v1`. It is not claimed native PCFF
or a faithful implementation of every statement in the FRC. The recursive
source snapshots preserve the pinned FRC, converter and official LAMMPS
statements, their original line locations, SHA-256 hashes and citations.

* Main torsion: the FRC prints `Vn*(1+cos(n*Phi-delta_n))`; LAMMPS uses
  `Kn*(1-cos(n*phi-phi_n))`. The pinned converter copies all six cells without
  a phase transformation. J60 explicitly adopts LAMMPS's cis-zero atom-order
  geometry for this research model and transforms `Kn=Vn`,
  `phi_n=(delta_n+pi) mod 2pi`. The identity
  `1-cos(x-pi)=1+cos(x)` proves exact energy, constant and force agreement for
  this declared geometry. It does not prove an absent original author definition.
  More generally, if a separately evidenced source geometry were
  `Phi=phi+alpha`, the mapping would be `phase=delta-n*alpha+pi`.
  A trans-zero source convention could therefore explain unchanged odd-harmonic
  converter phases (the pilot's nonzero source amplitude is n=3). J60 does not
  call the converter wrong or claim an independently established author angle
  definition. Its explicit research choice is alpha=0; no alternative is chosen
  after inspecting held-out errors. The conditional phase identity is tested for
  arbitrary alpha as well as the selected value.
* Angle–angle–torsion: the FRC literally prints
  `K*(Theta-Theta0)*(Theta'-Theta0')*(Phi-Phi1(0))`; LAMMPS uses
  `M*(theta1-theta10)*(theta2-theta20)*cos(phi)`. A nonzero linear-Phi term
  gains `2*pi*K*dt1*dt2` per turn; a cosine term is periodic. There is no global
  phase shift or constant coefficient conversion between them. Both original
  statements are preserved. The research variant selects the engine cosine
  term with a new explicit unknown coefficient. No linear-FRC coefficient is
  copied, relabeled or invented.
* Torsion–torsion is outside this selected finite LAMMPS class2 Hamiltonian.
  The nine five-site paths are marked model-inapplicable, not assigned zero
  coefficients. The empty FRC section is retained without a native absence claim.
* The pinned converter creates Wilson terms only at degree-three centers.
  This graph has Si degree four, O degree two and H degree one. All degree-four
  angle–angle products remain explicit; no Wilson coefficient is synthesized.
* Research-owned pair settings are `special_bonds lj/coul 0 0 1 angle no
  dihedral no`, sixth-power LJ9–6 mixing, dielectric 1 and the pinned LAMMPS
  real-unit Coulomb constant 332.06371 kcal Å/(mol e²). The 21 pairs at graph
  distances one/two are excluded; all 15 pairs at distances three/four have
  unit weight. Pair weights are fixed implementation choices, never fit variables.
  Source charges remain unchanged and sum to zero without normalization.
* The molecule is isolated: direct included pairs, 100 Å cutoff, no switching,
  shifting, tail correction, periodic images or KSPACE. All sampled separations
  lie inside the cutoff. A containing fixed-boundary box has no energetic volume
  role. This is not a periodic electrostatics or pressure convention claim.

`interaction-accounting.json` assigns every historical record to verified fixed,
explicit fitted unknown, structural/model inapplicability, or a precise dependency.
The selected functional accounting is complete. Numerical parameter completeness
is **incomplete**: all 20 independent unknown coefficients are required explicitly
for a full energy evaluation. No unknown is silently subtracted as zero.
The physical relationship remains

```text
QM = verified fixed MM + explicit joint unknown contributions + model discrepancy
```

The 28 confirmations retain their actual user-attested provenance and
`original_author_artifact_status="not_yet_supplied"`. No author email, document,
coefficient or missing numerical interpretation is manufactured.

For an ordered path `i-j-k-l`, let `dr1,dr2,dr3` be bond
length deviations from the verified source equilibrium rows and `dt1,dt2`
be the two angle deviations in radians. The explicit unknown part is the sum
over the inventory sites of the following functions, tied only by the recorded
atom-type/path symmetries:

| Family | Energy feature | Independent coefficients |
| --- | --- | ---: |
| Angle–angle | `dt_a * dt_b`, with the declared shared center/neighbor | 3 |
| First/third bond–bond | `dr1 * dr3` | 1 |
| End bond–torsion | `dr1 * sum(Bn*cos(n*phi)) + dr3 * sum(Cn*cos(n*phi))` | 6 |
| Middle bond–torsion | `dr2 * sum(Dn*cos(n*phi))` | 3 |
| Angle–torsion | `dt1 * sum(En*cos(n*phi)) + dt2 * sum(Fn*cos(n*phi))` | 6 |
| Cosine angle–angle–torsion | `M * dt1 * dt2 * cos(phi)` | 1 |

The harmonic index is 1–3. The fixed part includes the complete source quartic
bond/angle polynomials, adjacent bond–bond and bond–angle terms, transformed
main torsions, and every included LJ9–6/Coulomb pair. The source equilibrium
values are retained; the QM center does not replace them. Forces are minus the
Cartesian energy derivative in kcal/mol/Å. The 197 accounting records comprise
125 fixed records (including typing/mass/charge records), 54 unknown interaction
sites, nine structural inapplicabilities and nine model inapplicabilities.
Graph-distance three **and greater** pairs have unit LJ/Coulomb weight.

## Numerical kernels and joint design

The evaluator reuses J59's verified quartic bond/angle, adjacent bond–bond and
bond–angle contributions. Forward Cartesian differentiation adds source-phase
transformed main torsions, the complete selected unknown-column Jacobian and
fixed LJ9–6/Coulomb pairs. Independent scalar equations use separate atan2
geometry and finite differences at 1e-4/5e-5 Å. Controls include rigid translation
and rotation, net force/torque, atom-order reversal, torsion derivatives and
continuity/periodicity across ±pi. Synthetic nonzero coefficient controls remain
separate from physical parameter values.

Four actual LAMMPS `run 0` synthetic controls exercise all dihedral class2 terms
at distinct torsions and across the branch boundary. They verify the selected
energy/force kernels; they do not qualify the incomplete physical profile or run MD.

Valid source/type symmetries leave **20 independent coefficients**: three
angle–angle, one first/third bond–bond, six end-bond–torsion, three middle-bond–
torsion, six angle–torsion and one cosine angle–angle–torsion. No arbitrary
nonbonded weights or extra phase nuisance parameters remain.

Complete energy/force Jacobian screening uses 100 deterministic candidate
geometries with independent Si–H/Si–O bond, local/central angle and left/right
silyl torsion deformations. Selecting 18 training geometries reduces the scaled
condition from **478.8312 to 5.7171** and the correlated end/middle n=3 column
cosine from **0.9999526 to 0.0068458**. The full rank is 20. This addresses the
previous correlation with real geometric information, rather than hiding
unbound nuisance columns. Two separately reserved torsional conformers lie
outside the training torsion range; neither participates in design selection,
weighting, bounds or parameter fitting. Rank is a design result, not evidence
that QM values or transferable coefficients already exist.

The predeclared objective uses relative energy/0.05 kcal/mol and Cartesian
force/0.5 kcal/mol/Å, one training-reference energy offset, natural-unit bounds
±1000, no regularization and no zero prior. The conditional solver rejects
missing coefficients, rank loss, a failed reference gate and conformer/geometry
leakage. Manufactured controls recover all 20 nonzero input coefficients; they
are not QM data or physical fits.

## Bounded reference study

`plan/reference-plan.json` was saved before execution. It permits 35 new full
energy/gradient jobs, two concurrent jobs, four threads and 4 GiB per job,
3600 seconds per job, 18000 seconds campaign wall time and a 20 GiB scratch cap.
Native requests, outputs, lossless gradients, failure records and hashes remain
separate from normalized ISLAND results. Exact completed requests are verified
and skipped on resume; failed/interrupted attempts are preserved and never
silently overwritten. Shared scientific environments are unchanged.

Psi4 1.9.1 runs in the existing `phase4j36-native-env`. Its installed basis
composition gives 377 functions for both def2-QZVP and def2-QZVPP, with identical
Si/O/H molecular shells. That would not be a useful added-cardinal comparison.
J60 therefore compares actual cc-pVQZ and cc-pV5Z compositions at common
geometries. It does not call QZVPP a higher-cardinal basis or claim a CBS limit.
The saved inventory retains every shell, exponent, contraction and basis-file hash.

The new J59 QZ/dense optimized center is reused without changing its identity.
Both 0.08-radian angle corners and 0.04-Å bond corners are computed at the same
center, with PBE0, no dispersion, DF and explicit def2-universal-jkfit auxiliary
basis. Separate comparisons address orbital-basis cardinality/family, unpruned
Treutler 99×590 versus 150×974 grids, SCF 1e-10/1e-8 versus 1e-12/1e-10 and
PBE0 versus B3LYP. No optimization or protocol blending occurs in these contrasts.
The old center effects and failed J59 gates remain immutable historical results.

The original numerical criteria are not widened: cardinal-basis relative
energy ≤0.05 kcal/mol, full-gradient RMS difference ≤0.5 kcal/mol/Å, mixed
energy/gradient derivative change ≤10%; grid energy ≤0.005 and gradient ≤0.02;
SCF energy ≤0.0001 and gradient ≤0.005. Relative deformation energies, gradients
and total-QM mixed scan derivatives are retained separately. SCF convergence
alone is not a QM adequacy claim. A conditional supplementary plan predeclares
22 fitting/held-out basis-control jobs, activated only if the complete reference
gates permit that work.

## Chemical scope and independent statuses

`scope-map.json` lists all 19 directed missing queries. Only four Si–H control
queries can be constrained here; the other 15 require separate methyl environments.
HMDS needs `sio` with `osi,c3,c3,c3`, methyl `c3/hc`, Si–C/C–H charge/equilibrium
records and distinct terminal/methyl torsions. DP2 adds interior `sio` with
`osi,osi,c3,c3`, `osi-sio-osi-sio`/`sio-osi-sio-osi` chain couplings, exact
end/interior typing and charges, and independent chain conformers. Si–H values
are never transferred automatically to Si–C. No cyclic, crosslinked, Si–OH,
Si–O–C or phosphorus coverage is inferred.

Parameter completeness, kernel correctness, physical backend preparation,
QM execution/integrity, fit completion/validation, optional cross-engine
equivalence and scientific validation are separate results under J55.
`production_validated=False` and `simulation_readiness="not_established"`.
No universal PCFF/COMPASS, complete FRC/PDMS coverage, pressure/Ewald/PPPM
closure, production equilibrium or long-NPT readiness is claimed.

The 72 graph automorphisms (equivalent H permutations and whole-end exchange)
plus proper rotation/translation alignment give nearest training RMS separations
of 0.1338 and 0.1954 Å for the two held-outs. They are not relabeled copies of
training points. The exact per-coefficient units, finite bounds, optimizer and
uncertainty limitations are in `parameter-specifications.json`.

`plan/grid-derivative-amendment.json` separately caps six missing grid corners
for the original reference, conditionally required before any fit. That path is
inactive when its original cardinal gate fails. The additive
`plan/tight-d-refinement-plan.json` instead tests silicon tight-d polarization
with 18 matched cc-pV(Q+d)Z/cc-pV(5+d)Z jobs. Their actual molecular compositions
are 363/621 functions, compared with 353/611 without the tight-d shell.
This separates fixed-cardinality polarization effects from corrected QZ/5Z
cardinality. If the unchanged gate passes, twelve corrected-basis grid/SCF
controls and then 22 reserved fitting/held-out basis-control jobs are allowed.
The maximum for this primary-failure/refinement path is **87 gradients**.
Neither the old failure nor the primary 35-job plan is rewritten. No 6Z or
augmented calculation is counted as executed merely because its basis was
discovered.

The pinned converter has an additional BB13 `cp`-type guard in
`GetParameters.c:506–529`, with its `quo_cp()` definition at lines 1222–1240.
The selected graph has no `cp`. J60 explicitly retains all six requested Si–H
BB13 paths under the general source equation and LAMMPS class2 form; it does
not reproduce or treat the converter's omission as a physical zero. This is
another recorded research-variant difference, with the coefficient still unknown.


## Reproduction and retention

Run from the repository root with the existing interpreters (their exact resolved
paths, versions and dependency inventories are in `engine/inventory.json`):

```bash
export PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
J60_CHECK_PY=../island-validation/phase4g1-env/bin/python
J60_QM_PY=../island-validation/phase4j36-native-env/bin/python
"$J60_CHECK_PY" -m island.pcff_qm_pilot.hamiltonian_campaign \
  --work ../island-validation/phase4j60-qm --python "$J60_QM_PY"
"$J60_CHECK_PY" scripts/finalize_phase_4j60.py \
  --work ../island-validation/phase4j60-qm --python "$J60_QM_PY"
"$J60_CHECK_PY" scripts/refine_phase_4j60.py \
  --work ../island-validation/phase4j60-tight-d --python "$J60_QM_PY"
"$J60_CHECK_PY" scripts/phase_4j60.py replay
"$J60_CHECK_PY" scripts/validate_phase_4j60_receipt.py docs/evidence/phase_4j60
```

These are bounded, conditional commands; they do not launch the 889-point J57
inventory. A completed output is skipped only after its request, native output
hashes and completion record verify. A failed or interrupted attempt is retained
and rejects an in-place retry: use a separately identified additive attempt,
retaining its parent requests and failure. Never remove an output directory to
circumvent those checks. Re-execution on a fresh work directory creates real new
QM calculations; replay reconstructs retained data without running a QM engine.

The native archives preserve exact requests, atom order, Psi4 text, lossless
Cartesian gradients, convergence/energy records, file hashes, execution logs and
attempt summaries. Normalized arrays, independent native/unit reconstruction and
QM-minus-fixed-MM results remain separate. Large scratch integrals stay outside
Git. Only Psi4 is executed; other QM adapters remain unexecuted.

`fit-design-metrics.json` eliminates the single training-reference energy offset
exactly as the solver does: the resulting 20-column rank is 20 and condition is
**5.9609**. The earlier selection metric 5.7171 uses absolute feature rows; both
are retained with their distinct definitions. The solver reports the actual
objective conditioning rather than substituting the screening metric.

J60 pilot execution contains no MD. Historical regression tests retain their own
original scope, evidence and failure/retry records. Software development tooling
iterations are described separately from final checks; a passing final Ruff
result is not described as an unqualified clean first-pass tooling run.


The additional complete nine-atom synthetic `run 0` control exports all eight
bonds, thirteen angles, six proper paths and eight AA-only improper carriers.
All twenty synthetic unknown coefficients and all 24 AA products are nonzero.
It matches the assembled evaluator to 7.11e-15 kcal/mol in energy and 2.17e-13
kcal/mol/Å in maximum Cartesian force. LAMMPS's AA-only carrier explicitly
encodes structurally absent Wilson energy with `K=0`; this is documented backend
encoding of model inapplicability, not a supplied or silently filled source row.
The initial export failed on a Python atom-ID mapping before native launch;
its source and traceback are retained separately from the successful correction.
Raw LAMMPS padding is retained byte-for-byte, with J60-local whitespace attributes;
the original diff-check failure and its scoped resolution remain in validation.


## Actual reference-accuracy results

| Contrast | Maximum relative energy difference, kcal/mol | Maximum full-gradient RMS difference, kcal/mol/Å | Maximum mixed-derivative change | Declared gate |
| --- | ---: | ---: | ---: | --- |
| Original bare QZ/5Z | 0.140644896 | 0.863700566 | 0.817482% | failed |
| Corrected QZ/5Z | 0.0318209447 | 0.2303391 | 1.0025% | passed |
| Corrected grid | 0.000156689086 | 0.00456348127 | 0.0200581% | passed |
| Corrected SCF | 2.61934474e-09 | 8.65161089e-06 | not sampled | passed |

The matched-family def2-QZVP/corrected-5Z ancillary contrast gives 0.0358048
kcal/mol and 0.204335 kcal/mol/Å, with mixed-derivative changes below 0.149%.
Its exact raw/normalized inputs are retained; it adds no QM jobs and introduces
no new acceptance threshold. The original bare-5Z PBE0/B3LYP contrast gives
0.0953608 kcal/mol and 1.87484 kcal/mol/Å. It is method sensitivity at that
exact basis, not a conflated corrected-basis method test. Any fitted result
therefore remains PBE0-specific research evidence. General/final QM adequacy
and scientific validation remain not established.

`plan/reference-gate-before-fit.json` binds the successful complete basis/grid/
SCF analysis and its input hashes before any fitting-reference native result
existed. The old bare-basis and J59 failures remain unchanged. For a later
campaign, use an additive evidence-directory copy with `--evidence`; the
execution commands above document the original campaign.


## Executed joint fitting decision

The actual calculations total 35 primary contrasts, 18 corrected-basis contrasts,
12 corrected numerical controls and 22 fitting/held-out basis-control jobs. All
87 SCF energy/gradient operations completed; no synthetic QM values, new QM
Hessians, optimizations or MD runs are counted. The old J59 optimized-center
result is reused separately. There are 29 distinct input geometries. Native
outputs, lossless gradients, exact requests, hashes and completed attempts are
retained in the four native archives. Six additional native LAMMPS static
synthetic controls are separate from the QM and physical-profile execution counts.

The held-out QZ/5Z reference controls pass: relative-energy differences are
+0.0152110 and -0.0207386 kcal/mol, with full-gradient RMS differences 0.195462
and 0.227906 kcal/mol/Å. The joint fit therefore executes against QM minus the
verified fixed contribution, using all twenty unknown columns and the declared
energy offset, weights and bounds.

| Data | Geometries | Energy RMSE, kcal/mol | Force RMS, kcal/mol/Å |
| --- | ---: | ---: | ---: |
| Training | 18 | 0.440747 | 2.592741 |
| Untouched distinct conformers | 2 | 0.941158 | 5.454228 |
| Declared held-out acceptance | — | ≤0.1 | ≤0.5 |

Maximum held-out force error is 17.5094 kcal/mol/Å. Rank is 20, scaled condition
5.9609, and no bound is active. All 100 training-geometry bootstrap replicates
are valid, but **nine of twenty** coefficients fail the predeclared 50% relative
sigma gate. The full covariance/correlation and bootstrap samples are retained.
Maximum absolute off-diagonal residual-covariance correlation is 0.766239 (two AA
coefficients); bootstrap correlation reaches -0.882620 between end-angle n=1
AT and middle n=1 BT. Residual-only Cartesian covariance is a linearized proxy,
not a formal confidence interval based on independent force components.
Bootstrap excludes source, method, basis and transfer-systematic uncertainty.

`experimental-joint-fit.json` retains all twenty trial values as
`fitted_external` research-extension candidates, together with per-geometry
energy and full Cartesian-force residuals. `fitted-candidate-provenance.json`
binds each value, unit and sigma to the source model, target, actual native
training/validation jobs, dataset, optimizer/version and coefficient checksum.
None is a source-derived/native PCFF row, an accepted parameter, or a validated
parameter even for this control. The original nineteen directed queries remain
incomplete and no methyl-PDMS transfer is established.

A post-fit **training-only** diagnostic projects source-fixed polynomial
stretch/bend baseline directions after removing all twenty coupling columns.
The remaining weighted residual sum of squares is 14466.8, entirely outside the
coupling Jacobian to numerical precision. Linear Si–O and Si–H stretching probes
individually project 47.8311% and 41.7646%. These are overlapping, nonadditive
projections; no probe coefficient is inferred or inserted. They identify a
baseline incompatibility with the chosen PBE0 target at this accuracy, not a
proof that the original FRC rows are wrong. More sampling of the same twenty-
column model cannot eliminate its already demonstrated training residual.

The smallest next numerical constraint is a **PBE0/cc-pV(5+d)Z stationary-center
geometry and projected stretching Hessian for the two Si–O and six Si–H bonds**,
or independently justified method-matched numerical baseline rows. Its purpose
is to constrain equilibrium forces and local curvature jointly with the cross
terms, not equate a molecular Hessian element to a bare bond coefficient.
Project translational/rotational modes and never demand an unprojected molecular
Hessian be strictly positive. Any baseline refit needs a new research-model/
target identity, propagation of equilibrium changes through every Class-II
cross term, and fresh validation conformers; the original hold-outs have now
been examined. The current weights, tolerances and source rows are unchanged.
A methyl-PDMS program separately requires the exact fifteen HMDS/DP2-directed
query environments listed in `scope-map.json`.

## Verification and preservation

The supported ordinary suite covers 4,900 historical nodes: 4,839 passed and 61
skipped, with all eight original shard logs retained and no shard failure.
The final J60 suite contains 52 focused nodes. The combined unique collection
is 4,952 nodes (4,891 passed, 61 skipped), including 430 J54–J59 regressions.
Ruff, compileall, five applicable pip checks, source-object reconstruction,
diff checks, native/unit reconstruction, manifest/receipt validation, relocation
and fresh-process replay pass. The historical validator preserves **5,849**
files, including its additional legacy J40 archive, byte-for-byte. The initial
tracked-tree snapshot has 5,848 entries; that distinction does not drop the archive.

Original scientific failures and tooling/export attempts are retained separately.
Final passing checks are not described as an unqualified first-pass phase result.
The additive final interpretation preserves the initial compute fitting decision;
no fitted value or raw output is revised. Parameter completeness stays incomplete,
physical LAMMPS/OpenMM preparation stays not run, optional cross-engine equivalence
stays not claimed, and scientific/production readiness stays not established.
`production_validated=False` and `simulation_readiness="not_established"`.
