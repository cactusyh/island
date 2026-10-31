# Phase 4J63 — HMDS stationary reference and coupled Class-II research model

This experiment addresses the exact explicit-H hexamethyldisiloxane graph,
Me3Si–O–SiMe3: 27 atoms, 26 single bonds, neutral singlet, stable atom IDs 1–27.
Its graph identity is
`d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd`.
The Si–H investigation ended in J62; no new Si–H calculation or coefficient
transfer is part of J63. The ordinary API remains `forcefield="pcff"`.

## Base and frozen experiment

Fetched `origin/main` and the exact branch base are
`251e17db8a19b652196e559409ddd601fc20ceef`. This contains J62 prerequisite
`197a4b8921cf8180124485d9f0ba840c2853f7af`. The branch is
`codex/phase-4j63-hmds-stationary-coupled-mm`. Main is neither merged nor pushed.
Historical evidence, implementations, source hashes and input identities remain
unchanged. Unrelated local work remains outside this change.

All 21 retained HMDS gradients were independently reconstructed before new QM;
see [the reconstruction](evidence/phase_4j63/J62-HMDS-reconstruction.json).
The last unfinished J62 SVP geometry is only an initial guess. It is not
relabeled as a method-matched stationary point.

The [frozen scientific plan](evidence/phase_4j63/plan/scientific-plan.json)
SHA-256 is `ad92b223c2acb4fd3265a5c7800107f133fe8f648d7b3a2503cd6a40c92218a8`.
It caps all attempted new energy/gradient evaluations at 80, including failures
and optimizer callbacks; optimization receives at most 40. The staged allocation
is 40 optimizer callbacks, 16 primary path points, 19 precision controls and three
conditional validation points, leaving two unallocated. There is no implicit
full or numerical QM Hessian and no full J57 scan campaign.

The existing pinned Psi4 1.9.1 environment runs PBE0/def2-TZVP without dispersion,
DF/`def2-universal-jkfit`, unpruned 150×974 quadrature, SCF energy/density
thresholds 1e−10/1e−8 and integral threshold 1e−12. Actual contracted orbital
composition is 399 functions / 175 shells; the auxiliary has 1107 functions / 359 shells.
The 99×590 grid and tighter 1e−12/1e−10 SCF contrasts are separate controls.
TZVPP has 553 functions / 213 shells: this tests additional polarization, not a
higher-cardinal or complete-basis limit. QZVP's 1079 functions / 343 shells are
inventoried but not activated automatically. Successful SCF or optimization
does not establish reference adequacy.

Limits are four threads and 4 GiB per process, at most two simultaneous workers,
24 GiB scratch, 1,800 seconds per standalone gradient, 24,000 seconds for the capped
optimizer and 36,000 seconds for the campaign. Prior native timings informed these
limits. The ledger counts every native gradient attempt, rather than counting
a numerical derivative or optimizer as one cheap logical job.

Native GAU_TIGHT force/displacement convergence is retained alongside independently
checked final Cartesian-gradient maximum ≤1.5e−5 and RMS ≤1e−5 hartree/bohr.
Optking's displacement vector mixes internal bond lengths in bohr and angles in
radians; the frozen plan's displacement field names are clarified additively in
[the coordinate note](evidence/phase_4j63/plan/native-convergence-coordinate-clarification.json).
No Cartesian displacement proof or full-minimum claim follows.

The separate TZVP optimization converged on callback 32, with process exit 0.
Its final Cartesian gradient maximum is 2.4398081e−6 and RMS is 6.7518276e−7
hartree/bohr; both independently checked stationary gates pass. The final
native internal maximum/RMS displacements are 3.94e−5/1.13e−5 in Optking's
mixed internal coordinates. The new geometry identity is
`d17b1ac5465d9f30c5904ffe4cbba26c9b5ad0c969b7c2be3e705bfbe5841cfd`,
and its target identity is
`63def2061d91e8eff8f0a8b05dec61af2f37b3e0b976e846972a299323493361`.
The retained callback geometry and its matching gradient define this stationary
candidate; no full molecular Hessian or verified-minimum claim is made.

## Graph-controlled paths and reference tests

Si–C(2–1) and C–H(1–10) stretches translate the complete connected fragment
on the moving end of the cut acyclic bond. Only that bond length changes among
graph bond/angle/proper coordinates. Steps are ±0.006 and ±0.003Å.
Methyl rotation cuts Si–C(2–1); siloxane rotation cuts O–Si(5–6). Each rotates
the complete connected moving fragment rigidly, with steps ±0.04 and ±0.02rad.
All graph bond lengths and angles remain unchanged; only proper dihedrals
spanning the selected axis change. Nonbonded distances across a cut can change
and remain fully included in the research Hamiltonian.

Every point retains exact coordinates, stable IDs, tangent, acceleration,
energy, all 81 Cartesian gradient components and derived forces. Rotational
curvature uses `d²E/ds²=tᵀHt+g·a`; it is not identified with `tᵀHt` alone.
The full Cartesian response, central energy curvature, path-gradient curvature,
acceleration contribution and both displacement sizes are retained. These
are total molecular responses, not isolated bond or torsion coefficients.

Precision gates use the frozen absolute-plus-relative allowances separately
for each response. Common-geometry TZVPP, grid and SCF controls compare relative
deformation energies, every Cartesian force component, maximum local vector
force errors, path forces/torques and curvature. A large C–H response cannot
conceal a soft siloxane failure. The two original J62 SVP/TZVP results and
failed or unfinished historical work remain unchanged.

All 16 primary gradients succeeded. The gradient-derived path curvatures are:

| Path | Coarse step | Fine step | Curvature units |
| --- | ---: | ---: | --- |
| Si–C | 425.911226 | 425.890158 | kcal/mol/Å² |
| C–H | 762.596070 | 762.506606 | kcal/mol/Å² |
| Methyl rotation | 6.624301 | 6.635788 | kcal/mol/rad² |
| Siloxane rotation | 1.420111 | 1.420258 | kcal/mol/rad² |

The frozen two-step gradient-curvature tests pass, but the primary siloxane
energy/gradient consistency gate **fails**. Energy second differences are
1.337506 and 1.309160 kcal/mol/rad² at 0.04 and 0.02 rad; their differences
from gradient-derived curvature exceed the unchanged allowance
`0.02 + 0.02*abs(gradient curvature)`. The curved-path Cartesian chain rule
includes the center-gradient/acceleration contribution; it is not replaced
by a straight Cartesian Hessian contraction. These successful calculations
therefore do not establish all-objective reference precision.

The campaign completed **67/67 gradients on their original attempts**:
32 optimization, 16 primary paths and 19 controls, below the 80-call cap.
Individual calls took 80.40–314.29 seconds; their summed wall times are
12,489.12 seconds (concurrent calls overlap). The raw archive contains 553
hash-bound files. Other J57 engine adapters remain unexecuted.

| Common-geometry contrast | Si–C | C–H | Methyl | Siloxane |
| --- | --- | --- | --- | --- |
| TZVP/TZVPP full precision gate | fail | fail | fail | fail |
| 150×974 / 99×590 grid gate | untested | untested | pass | fail |
| SCF 1e−10/1e−8 / 1e−12/1e−10 gate | untested | untested | pass | pass |

Every basis contrast passes the declared relative-energy and path-curvature
tests, but local Cartesian-force errors fail for every path. Maximum local
vector-force differences are 0.571012, 0.580656, 0.570880 and 0.570747
kcal/mol/Å, respectively. The C–H path-force differences also fail. At the
common center the force RMS shift is 0.289559 kcal/mol/Å. TZVPP curvatures
are 425.513372, 759.152781 kcal/mol/Å² and 6.701630, 1.319840 kcal/mol/rad².
This is an additional-polarization sensitivity study, not a complete-basis claim.

The grid/siloxane failure is specifically path torque: differences are
0.008069 and 0.006735 kcal/mol/rad, exceeding their frozen allowances.
Its Cartesian-force and curvature comparisons pass. Tighter SCF changes
the center force RMS by only 5.3754e−5 kcal/mol/Å. Its comparison gates pass,
yet its own siloxane energy/gradient curvature discrepancy remains
0.105309 kcal/mol/rad², 2.176 times the original allowance. Tight SCF therefore
does not repair the failed consistency objective. No tolerance was widened.

## Complete research Hamiltonian and missing coefficients

The internal model is `island_research_class2_hmds_j63_v1`, with binding identity
`d9ae6279735fb219c7cdaf447971dde5995980292f2c019b71710213d5b3851c`.
It is a source-bound research implementation, not a native-PCFF equivalence
claim. [The definition](evidence/phase_4j63/Hamiltonian-definition.json) pins
J54 source assignments and J60's documented execution-equation authority.
[The interaction table](evidence/phase_4j63/interaction-accounting.json)
classifies every graph interaction, numerical source row, equilibrium
dependency, role, pair convention and missing coefficient.

Independent graph enumeration verifies 26 quartic bonds, 49 quartic angles,
49 adjacent bond-bond and 49 bond-angle interactions, 96 angle-angle roles,
and 60 proper paths for each main torsion and applicable dihedral coupling
family. All 351 atom pairs are classified. Source-complete terms are executable;
missing terms require explicit caller-supplied coefficients for any complete
evaluation. No missing physical coefficient defaults to zero.

The 26 directed missing queries reduce to **56 independent coefficients**:

| Family | Coefficients |
| --- | ---: |
| Angle-angle |5|
| First/third bond-bond |3|
| End bond-torsion, three harmonics at each end |18|
| Middle bond-torsion, three harmonics |9|
| Angle-torsion, three harmonics at each end |18|
| Angle-angle-torsion |3|

For angle-angle, J is the center and K is the shared neighbor: only I/L
exchange is allowed. Proper reversal exchanges the end blocks. Identical source
types share coefficients under these role rules; there is no invented equality
between distinct environments. Eighteen angle-angle sites already have source
coefficients. Missing cross-term rows remain missing numerical information,
regardless of confirmation of their interpretation.

Source plus-cosine main torsions, zero source phases and signed amplitudes are
preserved; the equivalent engine minus-cosine convention has phase π.
Angle-angle-torsion uses the separately identified cosine research form,
explicitly distinct from the FRC's printed linear-Phi expression. BB1_3 terms
on all declared paths remain unknowns; converter `cp` lookup behavior is not
authority for physical zero. Wilson out-of-plane terms use the selected
degree-three converter applicability; HMDS has no degree-three centers.
All degree-four angle-angle roles remain included. Torsion-torsion lies outside
the selected LAMMPS class2 research sum. That model boundary is documented,
not presented as a numerical zero or original-PCFF omission.

LJ9-6 uses sixth-power mixing of source r-minima and well depths. Graph-distance
one/two pairs are excluded and distance≥three pairs have full LJ/Coulomb
weight. Isolated electrostatics has dielectric 1, prefactor 332.06371
kcalÅ/mol/e², no periodic summation, shifts or tail correction, and a 100 Å finite
molecular bound. No arbitrary pair fitting weights are introduced.

Unchanged PCFF source bond-increment charges are Si +0.52, O −0.44, C −0.259,
H +0.053 e. They conserve molecular charge without normalization. No OPLS/GAFF
rows, Si–H trial values or fabricated author documents are used. The 28
convention confirmations remain user-attested; original author artifacts are
still `not_yet_supplied`.

Independent NumPy scalar equations verify the production energy/Cartesian-force
and 56-column parameter Jacobian at two finite-difference steps, together with
rigid invariance, role permutations, path reversal, periodicity and all
cross-term equilibrium dependencies. Explicitly manufactured nonzero values
test complete recovery; they are software controls, never physical fitted rows.

## Joint design and scope

Geometry-only screening uses the complete 56-column relative-energy/Cartesian-force
Jacobian. A frozen pool of 48 joint bond/angle/torsion geometries is screened
by rank and column-scaled condition, selecting at most six greedy blocks and
exhausting subsets of that selected set. Any minimum-cardinality claim is
restricted to that bounded pool; it is not global design optimality.
Diagnostic energy/force scales 0.05 kcal/mol and 0.5 kcal/mol/Å are not an adopted
physical fitting objective or validation-derived weights.

The source-baseline compatibility calculation projects QM minus all verified
fixed contributions against the complete unknown column space. Its residual
is an optimistic unconstrained joint lower bound. It does not assign omitted
terms zero, select physical coefficients, prove which source row is wrong,
or automatically release source baselines.

An additional development diagnostic uses every successfully retained new TZVP
optimizer snapshot under one exact protocol, with the same frozen diagnostic
energy/force scales. Its orthogonal residual is independently reconstructed by
pivoted QR when the full Jacobian has rank 56. Per-configuration residuals and
maximum atom-force errors accompany RMS errors. These optimizer points are
correlated development data; they do not supply independent conformer validation,
coefficient uncertainty or a physically constrained fit. This analysis requires
no new QM and does not change stationary or response-stage activation rules.

For all 32 terminal optimizer snapshots, rank is 56 and scaled condition is
13,649.89. The optimistic joint projection leaves relative-energy RMS
0.051570 kcal/mol, Cartesian-force RMS 0.175899 kcal/mol/Å and maximum
atom-force residual 1.168810 kcal/mol/Å. The SVD residual is independently
checked with pivoted QR. These errors describe an unlimited-coefficient
representation bound on correlated development data, not a validated fit or
a statistical uncertainty estimate.

At the actual stationary geometry, the center plus 16 primary path geometries
have rank 56 but scaled condition 73,330.57. The frozen screen selects
`joint-040`, `joint-028`, `joint-038`, `joint-044`, `joint-020` from its 48-point
pool; together they reduce condition to 967.46. This requires **five additional
energy/gradient calculations**, whose exact coordinates and actions are retained.
Those QM values have not been supplied. It is the smallest passing subset
within the bounded greedy-selected set, not a global minimum-design claim.

Pure Si–C and C–H fragment stretches provide a distinct baseline diagnostic.
Only one bond changes, linearly, while all angles and proper torsions remain
fixed. Each of the 56 missing coupling kernels is at most linear in that
displacement, so its second path derivative is identically zero. Independent
scalar equations and the full force/parameter Jacobian check this null relation.
All verified fixed contributions, including nonbonded changes, remain included.

| Path | Step (Å) | Total QM curvature | Total fixed-model curvature |
| --- | ---: | ---: | ---: |
| Si–C | 0.006 | 425.911226 | 374.754336 |
| Si–C | 0.003 | 425.890158 | 374.715691 |
| C–H | 0.006 | 762.596070 | 720.733700 |
| C–H | 0.003 | 762.506606 | 720.642484 |

Curvature units are kcal/mol/Å². The fine-step differences are 51.174467
(12.016%) and 41.864122 (5.490%). Missing couplings cannot repair these
controlled total-model curvatures. This does not isolate a uniquely incorrect
bond row: the selected fixed source/nonbonded model and reference protocol
are assessed together. No source baseline is automatically released.

Three genuinely new deformed configurations are reserved before their QM,
fitting weights or bounds. HMDS has 3,359,232 graph automorphisms. Distance-colored
graph backtracking plus proper Kabsch checks methyl, H and whole-end mappings;
search exhaustion fails closed. The nine-atom 72-map shortcut is not reused.
These are distinct deformations, not independently sampled conformational
basins. Historical J60–J62 validation is development data.

HMDS covers terminal `sio(osi,c3,c3,c3)` and `c3(sio,hc,hc,hc)` environments.
Capped-DP2's interior `sio(osi,osi,c3,c3)`, O–Si–O and Si–O–Si–O remain a
separate graph/source/reference dependency. HMDS alone cannot establish full
PDMS parameterization, cyclic or crosslinked domains, Si–H/Si–OH/Si–O–C or
phosphorus support. No complete profile, automatic assignment promotion or
PCFF MD is created. Parameter completeness, backend execution, optional
cross-engine equivalence and scientific validation remain separate:
`production_validated=False`, `simulation_readiness="not_established"`.

## Fitting decision and smallest next requirements

**No physical coefficient was fitted or accepted.** There are no rejected
physical candidates hidden in a successful optimizer status: the prerequisites
for a preregistered physical fit were not met. The complete primary-data
column-space projection leaves energy RMS 0.004428 kcal/mol, Cartesian-force
RMS 0.055431 kcal/mol/Å and maximum atom-force residual 0.571056 kcal/mol/Å.
Its good aggregate errors cannot establish adequate local accuracy, repair
the independently measured unknown-null curvature gaps or remove conditioning.
All 56 physical coefficients remain numerically missing.

[The specific next-data decision](evidence/phase_4j63/next-data-requirements.json)
first proposes **three** common-geometry gradients: the retained center and
siloxane ±0.02rad, PBE0/def2-TZVP, unpruned 250×1454 grid, SCF 1e−12/1e−10.
Comparing these with the retained tighter-SCF 150×974 triple isolates the grid
effect without changing centers or thresholds. A separate **one-gradient**
QZVP center-force comparison would discriminate the observed basis-force shift;
it needs a new resource plan and would test only the center, not close all
deformation or complete-basis gates. Neither proposal was executed here.

Before further coupling-fitting QM, the total-model Si–C/C–H baseline gaps
require independently justified source/baseline accounting or a separately
identified, preregistered HMDS research extension. Source rows remain fixed
and unchanged in J63. The five exact joint geometries address conditioning
only after reference and baseline prerequisites are resolved; the three reserved
validation deformations must remain untouched by parameter/weight selection.
Capped-DP2's interior graph, typing/charges, source dependencies and O–Si–O /
Si–O–Si–O response data remain a separate methyl-PDMS entry requirement.

## Reproduction and exact-job resume

Offline native reconstruction unpacks into a fresh temporary directory, checks
the complete native inventory, independently reconstructs ASCII gradients and
text energies, and replays center, path, design and baseline diagnostics:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  scripts/replay_phase_4j63.py replay
```

Numerical reruns use a separate source checkout, an isolated evidence copy and
a new external work directory. For exact v1 producer reconstruction, restore
the two archived execution sources **in that separate checkout only**:

```python
import lzma
from pathlib import Path

root = Path("docs/evidence/phase_4j63/plan")
for archive, target in (
    ("executed-hmds-engine.py.xz", "src/island/pcff_qm_pilot/hmds_engine.py"),
    ("executed-phase-coordinator.py.xz", "scripts/phase_4j63.py"),
):
    Path(target).write_bytes(lzma.decompress((root / archive).read_bytes()))
```

Copy the frozen scientific plan into the new work directory first. Run the
staged commands from that checkout with the pinned interpreter. Rerun outputs
have new execution identities; they must not overwrite the retained J63 archive:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/phase_4j63.py optimize \
  --evidence /absolute/path/to/copied-j63-evidence \
  --work /absolute/path/to/new-j63-work \
  --python ../island-validation/phase4j36-native-env/bin/python

# The deformation command independently checks the frozen stationary gate.
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/phase_4j63.py deform \
  --evidence /absolute/path/to/copied-j63-evidence \
  --work /absolute/path/to/new-j63-work \
  --python ../island-validation/phase4j36-native-env/bin/python
```

Repeating a command on the same work directory reuses only exact hash-verified
completed jobs. An incomplete, failed or interrupted attempt is retained and
requires a separately identified retry; it is not silently overwritten. Native
optimization returns0 only on native convergence,3 for exhausted callback budget
and1 for unexpected failure. Subprocess timeouts remain124. Archive collection
and scientific acceptance are separate from these numerical process statuses.
Commands, selected environment settings, return codes, stdout/stderr and native
failure artifacts are preserved before scratch cleanup. Large scratch and
wavefunctions stay outside Git.

After all 67 native processes exited, the explicit
[maintenance amendment](evidence/phase_4j63/plan/engine-maintenance-amendment.json)
removed an unused import and tightened prospective budget-exception handling.
The original v1 source SHA-256
`a82f5c12a57b4cd01268b028d24f27ce2ec7513e7789b8d9da7faedadde0bb8e`
and original coordinator remain archived. The maintained adapter is v1.1,
SHA-256 `069df0819effa5f6d1510ed1d4be8ea486b0e1357d9b5b413e172fd524014e6a`.
Strict AST reconstruction permits only the declared import/status changes;
numerical settings, callbacks, output capture and kernels are unchanged.
Unrelated metadata failures at callback40 can no longer be mislabeled as
budget exhaustion. No such failure occurred in the actual 32-callback run.
The current coordinator accepts the explicit amendment, and exact optimize /
deform resume verified every completed job with an unchanged ledger and
**zero new gradients**. Every retained native file remains bound to v1.

## Verification and preservation

The supported historical suite passed on its first eight-shard run:
4,997 passed, 61 skipped. The terminal J63 suite passed all 115 tests, including
all energy/force/parameter derivatives, graph role closure and symmetry,
controlled paths, equilibrium dependencies, optimizer status/budget controls,
native mapping/unit mutations and relocated fresh-process reconstruction.
The combined accounting contains 5,173 unique cases: **5,112 passed, 61 skipped**,
with no duplicated cases or failed/error cases in that terminal coverage.

Earlier development attempts are retained separately in
[the attempt ledger](evidence/phase_4j63/validation/early-attempts.json).
One manufactured-analysis test initially failed because its fixture aliased
the baseline and mutated force objects; the deepcopy correction and retry
are both retained. Initial lint findings and command errors also remain
visible. The first receipt-sealing attempt also exposed a validator tuple/list
comparison bug; its unchanged accounting hash and exact JSON content were
checked before correcting the comparison. That failure and retry are separate
from the successful numerical reconstruction. This is not an unqualified clean
first-pass aggregate.

J54–J62 receipt/regression checks, five applicable `pip check` environments,
Ruff, compileall and diff checks passed. The receipt validates its manifest,
complete test-node accounting, source pins, native/source reconstruction,
scope, compute limits and every historical file. All **6,785 base files** remain
byte-for-byte unchanged. Native producer hashes and failed scientific gates
are preserved alongside successful source/status maintenance.

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python \
  scripts/validate_phase_4j63_receipt.py docs/evidence/phase_4j63 --repo .
```
