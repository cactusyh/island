# Phase 4J40-C12 — explicit strained-box special-pair comparison completion

C12 reproduces and isolates the special-pair cutoff discrepancy without changing
charges, force-field parameters, topology, ordinary pair cutoffs, or numerical
acceptance bounds. It adds an **explicit analytical single-frame comparison
adapter**. Native LAMMPS output is retained verbatim; a separate composite row
contains the missing special-pair energy, force and virial contributions.
This is not a patch to LAMMPS, a modified production force field, or a claim that
unmodified native LAMMPS now implements OpenMM's untruncated exceptions.

The tested implementation is `67d237a2045db2259fd95c0a2d7b39ed16dd4ef8`, based on
C11 `ac0c38185e7b316b68109cd8be24e00c897f3fa8`. The production exporter retains
its existing rejection of scaled special pairs outside cutoff. Scientific
acceptance, backend pressure qualification, cutoff recovery and long NPT remain
blocked. Historical evidence through C11 is preserved in place.

## Exact reproduction and identity

The first retained strained case contains two neutral, hydrogen-complete C8H18
chains, 52 atoms, exact-graph RESP/Psi4 charges, and a 1.5× isotropic strain of the
C11 coordinates and orthorhombic box. Its box lengths are
1.6220832216509522 nm each. The real-space cutoff is 0.5 nm, switching is disabled,
and the native LAMMPS skin is 0.2 nm. Eighteen special pairs lie beyond cutoff.
The unstrained box is 1.081388814433968 nm per axis; no special pair lies beyond
cutoff there or in the 1.25× box.

Frame records bind graph, numerical parameters, charge identity, stable atom IDs,
coordinates and box. Pair membership is independently reconstructed from shortest
covalent paths and checked against the authoritative explicit exception list.
The adapter rejects altered scaling, cutoff/native-setting mismatches, changed
boxes, changed charges, stale native frames and reordered atom IDs. The native
LAMMPS dump supplies independent coordinates modulo periodic images, charges,
atom ordering and force-unit checks. Coordinates are exported with 17 significant
digits in angstroms; the native execution checks their nanometer round trip
modulo the box with absolute tolerance 3e-15 nm and zero relative tolerance.
The standalone parser uses 1e-13 nm for its independent mapping check. Neither
round-trip audit changes the energy or force acceptance bounds. The original charge vectors are unchanged. A standalone raw-source audit reads
the tracked J40 tar archive directly and verifies all eight retained backend/data
members against the bytes used here and the historical inventory. All eight
match byte-for-byte; `validation/raw-source-audit.json` records their SHA-256
values and the historical archive hash. This check does not import the production
backend loader.

| Convention | Retained setting |
| --- | --- |
| GAFF2 1–4 LJ / Coulomb | 0.5 / 5/6 |
| OPLS-compatible 1–4 LJ / Coulomb | 0.5 / 0.5 |
| 1–2 and 1–3 | excluded |
| GAFF2 LJ mixing | Lorentz–Berthelot |
| OPLS-compatible LJ mixing | geometric |
| Boundary | three-dimensional periodic; conducting Ewald/PME boundary |
| OpenMM platform | Reference, double precision |
| Native LAMMPS | KSPACE Ewald and PPPM, double FFT build |
| Energy / force units | kJ/mol / kJ/(mol nm) |
| Native real-unit force conversion | kcal/(mol Å) × 4.184 × 10 |

The pinned native executable SHA-256 is
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
Its raw help/build output and the Python/OpenMM versions, library hash,
platform inventory and precision limitations are retained in
`study/c11-rerun/environment.json` and `lammps-help.txt` inside the artifact archive.
OpenMM is `8.6.1.dev-b399af4`; LAMMPS is the retained 10 Sep 2025 Development
`g36f7c81ced` KSPACE build. No new quantum charge generation is claimed.

## Demonstrated cause and correction equation

OpenMM's Coulomb exceptions and LJ14 bond force are not truncated at the ordinary
nonbonded cutoff. The retained LAMMPS `lj/cut/coul/long` pair evaluation applies
its real-space neighbor cutoff to special pairs. Outside that cutoff, it omits
both the scaled LJ term and the direct-space special Coulomb expression, while
its reciprocal contribution remains. Native real/reciprocal pieces use each
engine's own Ewald splitting parameter; their individual values are retained
without treating them as equal physical observables. The unchanged component
acceptance test uses total Coulomb energy. `special_bonds` factors themselves agree;
atom ordering, units, charge origin and minimum-image conventions also agree.

For a special pair outside cutoff, let `r` be its orthorhombic minimum-image
distance, `fL` and `fC` its source scaling factors, and `alpha` the native LAMMPS
Ewald splitting parameter. The explicit missing terms are

```text
Delta U_LJ = 4 fL sqrt(epsilon_i epsilon_j) [(sigma_ij/r)^12 - (sigma_ij/r)^6]
Delta U_C  = C q_i q_j [erfc(alpha r) - (1 - fC)] / r
C          = 138.935456264 kJ nm mol^-1 e^-2 (pinned LAMMPS prefactor)
F_i        = (d Delta U/dr) (r_j - r_i)_minimum_image / r
F_j        = -F_i
W_pair     = -(r_j - r_i)_minimum_image outer F_i
Delta P    = trace(W_pair) / (3 V)
```

The small erfc tail is evaluated, not discarded. Reciprocal output is untouched.
The formulas also retain the missing Coulomb exclusion correction for 1–2 and
1–3 pairs (`fL=fC=0`); they do not invent LJ terms for excluded pairs. A separate
control exercises both exclusion classes beyond cutoff. LJ and Coulomb energies,
forces and virial increments are persisted separately for each affected pair.
The native log rounds alpha to about eight digits; its value and this limitation
are recorded. No fitted residual, atom-charge repair or force rescaling is used.

The completion applies only to the stated single-frame orthorhombic, unconstrained
GAFF2/OPLS-compatible models with no switching. Ambiguous half-box images are
rejected. It does not make a cutoff-crossing trajectory differentiable, supply a
production integrator, or qualify constrained/triclinic pressure.

## Executed results and remaining failures

The full unchanged C11 numerical plan was rerun: 38 paired single points plus 24
strain evaluations (100 native evaluations). An additional 48 paired points
cover point charges, explicit counterions, pairwise controls, a complete C8H18
component, retained GAFF2 PE, OPLS PEO, aromatic OPLS-compatible RESP/Psi4 and the
packed OPLS network. Each is evaluated under Ewald and PME/PPPM in the original
box, the 1.5× box and its wrapped rigid translation. Twelve new native evaluations
supply a three-step finite-volume pressure sweep at the failing strained frame.
The final frozen matrix therefore executes **208 native engine evaluations**.
A separate coordinate-only excluded-pair control executes two more native
evaluations, for 210 native evaluations including that additional control.

Requested diagnostic tolerances are 1e-9 for the additional direct Ewald points
and 1e-7 for PME/PPPM, matching C11's retained PME diagnostic accuracy. The full
C11 tolerance/mesh/cutoff/order/alpha sweeps remain unchanged. A preliminary
1e-9 expanded-box PME run was stopped for excessive memory use and is recorded
as interrupted, not successful. These are reciprocal solver requests; they do
not replace the unchanged acceptance bounds: total/component energy 1e-5 kJ/mol,
force absolute 1e-5 kJ/(mol Å), force relative 2e-10.

| Exact GAFF2 PE comparison | Native energy error (kJ/mol) | Composite energy error (kJ/mol) | Composite max force-component error (kJ/mol/nm) |
| --- | ---: | ---: | ---: |
| Ewald, unstrained, 1e-9 | 1.90109e-5 | 1.90109e-5 | 4.36304e-5 |
| Ewald, 1.5× box, 1e-9 | 1.754276325868 | 1.75198e-7 | 1.71331e-5 |
| Ewald, 1.5× box, translated | 1.754276325864 | 1.75194e-7 | 1.71330e-5 |
| PME/PPPM, unstrained, 1e-7 | 2.90897e-4 | 2.90897e-4 | 4.58447e-3 |
| PME/PPPM, 1.5× box, 1e-7 | 1.755187128478 | 9.17901e-4 | 3.38084e-3 |
| PME/PPPM, 1.5× box, translated | 1.754393471878 | 1.24244e-4 | 5.61263e-3 |

The Ewald strained-frame special-pair mismatch is resolved by the labelled
completion under the existing energy/force bounds. The unstrained residual is
unchanged, and PME/PPPM mesh-sensitive residuals remain above the bounds. Aromatic,
PEO and network rows also retain explicit residual failures. No setting is selected
as a universal repair. The native long-range decomposition is preserved separately
from the completed totals, and both native and composite pass/fail states appear
in the matrix.

The strained direct-Ewald finite-volume derivative is approximately
-184800.103 bar in OpenMM, versus -184795.357 bar natively in LAMMPS and
-184800.103 bar after the explicit special-pair completion. These are zero-kinetic,
high-energy stretched-frame configurational diagnostics, not measured NPT
pressures or evidence of equilibration. The pair virial increment passes independent
energy-volume derivatives, but OpenMM's full force/reciprocal virial and general
backend pressure remain unqualified. No new pressure acceptance bound is invented.

## Independent audit and persistence

The standalone archived `independent_audit.py` imports no ISLAND/C11/C12 importer,
metric or completion helper. It reads native OpenMM XML particle charges,
exceptions and LJ14 parameters and evaluates an independently constructed
`CustomBondForce` using `fC-erf(alpha*r)`. It checks all 48 additional frames plus
the exact C11 strained case: energy, per-atom force, pair virial versus finite-volume
derivative, raw LAMMPS atom/coordinate/charge mapping and real-unit force conversion.
All 49 main controls passed. A separately executed coordinate-only experiment
moves atom 1 by +0.6 nm along x, retaining the exact graph, charges and parameters.
It places three 1–2 exclusions, one 1–3 exclusion and seventeen 1–4 pairs beyond
cutoff. Its native/composite energy errors are 1.619411 kJ/mol and
1.72586e-8 kJ/mol, respectively; the composite maximum force-component error is
1.99800e-5 kJ/(mol nm), within the unchanged bound. The same standalone XML audit
passes this geometry and the repeated original strained control. These two audit
rows include one duplicate frame and are not two new independent physical samples.
The extra native outputs and independent report are under `validation/excluded-control`
and `validation/exclusion-independent-audit.json` in the validation archive.
Its numerical bounds apply to the isolated correction
identity, not to full-engine acceptance. Independent multiplication of force-unit
factors is allowed a documented four-epsilon float64 rounding bound; charges and
atom ordering are compared directly. The standalone audit output records its
own source hash, numerical errors and bounds.

The receipt validator reconstructs every completion and discrepancy from raw
records, checks native dumps and XML/source hashes, validates both receipt and
manifest identities, verifies the tested Git source commit and preserves all
historical evidence byte hashes. Mutation tests include renewed outer digests;
changed vectors, cutoff, scaling, box, ordering, native outputs and readiness
claims are rejected. Truncated force arrays cannot broadcast across atoms.

The archived original failures and native outputs remain reviewable. The
analytical adapter does not authorize exporting or simulating these strained
frames with unmodified native LAMMPS. `ewald_unvalidated`, `pressure_unavailable`,
cutoff recovery `not_executed` and long NPT `not_executed` remain explicit.
RESP/Psi4, native OPLS-AA, analytic controls, unsupported PCFF and unavailable
AM1-BCC remain distinct; no unsupported path is promoted.

## Reproduction commands

The exact final command and source/executable hashes are in `study/invocation.json`;
all native LAMMPS commands, raw logs and dumps are retained. In the declared
five-environment workspace, the executed command is:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python \
  -m scripts.phase_4j40_c12_study \
  --source /tmp/j40-validation/final/run \
  --reference /tmp/j40-c2-validation/run/reference \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j40-c12-validation/study
```

Its exit code is intentionally 1 because scientific acceptance is blocked. For a
portable rerun, extract `artifacts.tar.xz` to an empty `EXTRACTED` directory and
`validation.tar.xz` to `CHECKS`, then reconstruct the retained input layout:

```bash
python CHECKS/validation/reconstruct_inputs.py EXTRACTED NEW-INPUTS
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python \
  -m scripts.phase_4j40_c12_study \
  --source NEW-INPUTS/source --reference NEW-INPUTS/reference \
  --lammps ../island-validation/phase4j22/lmp-kspace --output NEW-STUDY
```

Offline evidence validation and independent audit require no original temporary
input paths:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j40_c12_receipt \
  docs/evidence/phase_4j40_c12/receipt.json --diagnostic-only
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python \
  CHECKS/validation/independent_audit.py EXTRACTED/study independent-audit.json
```

Integrity returns 0; scientific acceptance without `--diagnostic-only` returns 1.
The validation archive records the exact focused, regression, full-suite, Ruff,
compileall and five environment-specific pip commands, environment variables,
source commit, elapsed times, logs and JUnit reports.

## Final verification

The frozen implementation passed 20 focused C12 tests and all 1,099 J22–J40
regression tests. The full supported suite reported 3,067 passed and 10 skipped,
with no failures or errors. The skips require optional AmberTools integration
configuration or archived reference manifests, retained analysis bundles, or an
OPLS source in the pinned Foyer environment; they are not counted as executed
integration passes. Exact reasons are retained in `validation/full.xml`.
Ruff, compileall, all five execution-environment pip checks and the diff check
passed. All 2,310 historical evidence entries through C11 remain byte-identical.

Two fresh-process receipt/manifest checks returned 0 for diagnostic integrity
and 1 for scientific acceptance, as intended. Their commands, receipt SHA-256,
elapsed times and exit codes are in `receipt-validation.json`. The nonzero
scientific gate preserves the remaining PME/Ewald and pressure limitations;
it does not indicate corrupt evidence. Offline reconstruction after extraction
to a different temporary directory also passed, without claiming a second MD
execution or checkpoint/RNG continuation.

`production_validated=False`; `simulation_readiness="not_established"`.
No experimental agreement, production equilibrium or universal force-field
coverage is claimed.
