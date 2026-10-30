# Phase 4J57 — user-attested conventions and PDMS QM input preparation

J57 binds the user's explicit attestation for 28 previously reported author
confirmations and prepares reference calculations for the remaining 19 directed
PCFF/COMPASS Class-II source queries. **No QM/reference datasets have been
supplied. No final PCFF profile has been created. No PCFF MD has been run.**

The 28 records now have `status="author_confirmed_audited"`, with:

```text
evidence_kind="user_attested_author_confirmation"
authority="user_reported_author_confirmation"
source_artifact_status="not_yet_supplied"
audit_scope="the 28 conventions and parameter interpretations explicitly confirmed by the user"
```

“Audited” describes the explicit user attestation, query membership, scope and
immutable links. It does not describe inspection of an original author document.
Original author artifacts have not yet been supplied. No author email, author
confirmation date, filename, document text, author-file hash or missing numerical
coefficient is invented. Hashes of the user-attestation records are hashes of
those records, never hashes of absent author files. Separately pinned historical
FRC files remain source-file evidence with their original, narrower authority.

## Base and preservation

The fetched exact `origin/main` is `5fd06ec0f25532ed164ba534b04a52a7d7a0aedb`
(the J54 merge following J55). J56 was available on its completed branch, not
merged into main. This branch starts at that exact main and carries the two
existing J56 commits as prerequisites, preserving their file contents. The
prerequisite snapshot is `24ad7698f64d62aeb08ed893eb8372bdbbf33819`; its complete
Git tree is the J17–J56 preservation baseline. The J57 receipt records both
identities. This phase does not merge main or rewrite historical J54/J56 claims.

All implementation is additive in `island.pcff_qm`. In particular, the J56
package and its source-derived code identities remain unchanged. Historical
receipts still describe their original reported-confirmation policy.

## Effective convention contract and numerical targets

The 28 confirmations comprise seven directed environments in each of:
`end_bond-torsion_3`, `middle_bond-torsion_3`, `angle-torsion_3`, and
`angle-angle-torsion_1`. They are bound to the effective convention contract and
are excluded from the missing numerical blocker list. Original-artifact absence
is not a convention blocker under the user's explicit provenance decision.
This does not create any missing numerical row or redistribute charges.

The remaining 19 directed queries are **12 angle-angle** and **7 first/third
bond-bond (`bond-bond_1_3`)** queries. These are the existing J54/J56 queries,
not 19 new transferable physical coefficients. They account for 1,068 numerical
site requests; the confirmed queries account for 1,848 of the 2,916 original
requests. The complete inventory retains the directed query, exact sites,
source equation/line/hash, dependencies, applicable graphs and previous target
identity. No unrelated targets are combined into one fit. Symmetry or a future
joint residual must be justified explicitly before fitting.

Every immutable target includes the final graph and bond orders, atom ordering
and labels, source typing and unmodified charge vector, neutral singlet state,
starting-geometry checksum, scan definition, functional form and fitting fields, bounds policy, constraints,
units, QM protocol, required observables, held-out split and uncertainty plan.
There is no reproducible equation that numerically determines an independent
missing mixed coefficient from the verified equilibrium parameters alone.

The initial proposed protocol is PBE0/def2-SVP, no dispersion, isolated neutral
singlet single-point gradients without geometry optimization. These are explicit
**J57 preparation choices**, not author-confirmed values or a claim of adequate
QM accuracy. Method/basis convergence, numerical bounds and held-out error
thresholds must be established before a fit can be validated. Unknown coefficient
bounds and uncertainty values are not invented.

## Reproducible geometries and engine templates

Each numerical target has a 5×5 independent paired scan, with bond offsets
`[-0.04, -0.02, 0, 0.02, 0.04]` Å or angle offsets
`[-0.08, -0.04, 0, 0.04, 0.08]` radians relative to its exact J54 input geometry.
The input geometry is not asserted to be an optimized equilibrium. Tree-fragment
translations and rotations preserve the exact graph and order. Achieved internal
coordinates are measured and recorded; failed or degenerate reconstruction
rejects the bundle. No live source system is changed.

There are 475 numerical-target geometries. There are also 52 bounded supporting
scans containing 414 geometries. These scans
cover available Si–O/Si–C bonds, Si–O–Si angles, Si–O–Si–O torsions, terminal
C–Si–O–Si torsions and joint bond-bond, bond-angle, angle-angle, bond-torsion,
angle-torsion and angle-angle-torsion coordinates. Supporting scans do not add
missing targets, authorize combined fitting or generate coefficients. No
Si–C–O or Si–O–C environment is invented.

`audit/geometry-input-bundle.tar.xz` is a deterministic regular-file-only archive.
Each target/supporting directory contains the initial XYZ, all scan XYZs,
coordinates and measured residuals, atom labels/order, full graph, target metadata,
expected native output filenames and a hash manifest. There are five explicit
input adapters: Psi4, PySCF, CP2K, ORCA and Gaussian. Central points request a
Hessian; every point requests energy/gradient information. Gradients must be
negated to become forces with an explicitly documented unit conversion.

These are input templates, not results or assertions of engine compatibility.
CP2K visibly requires a pinned all-electron def2-SVP basis library and element
KIND declarations; its cell/grid/vacuum convergence must also be checked before
use. No engine executable is assumed, called or installed. Executable inspection
returns `asset_unavailable` for absent engines, without manufacturing outputs.
Adapter syntax and actual QM accuracy remain unverified until real runs exist.

## Commands

Run with the repository and `src` on `PYTHONPATH`. `--distribution` points to the
pinned FRC distribution used by J54. It is an advanced source-asset path, not a
profile or charge-method selection.

```bash
python -m island.pcff_qm confirmations --json confirmations.json
python -m island.pcff_qm artifacts --json original-artifact-status.json
python -m island.pcff_qm assets --json engines.json
python -m island.pcff_qm targets --distribution /path/to/pinned/frc_files --json targets.json
python -m island.pcff_qm generate --distribution /path/to/pinned/frc_files --output qm-inputs --json inputs.json
python -m island.pcff_qm closure --distribution /path/to/pinned/frc_files --json closure.json
python -m island.pcff_qm readiness --distribution /path/to/pinned/frc_files --json readiness.json
```

Use `--fixture hmds --query 'angle-angle:c3,sio,c3,c3'` to select an exact target;
`--engine psi4` restricts generated templates. The default target inventory
selects the smallest verified exact fixture for each directed query and retains
all other applicable graph bindings. Selection never creates transfer authority.

The additive family-only Python preparation API is:

```python
from island.pcff_qm import prepare_polymer_forcefield
report = prepare_polymer_forcefield(system, forcefield="pcff")
```

Source assets can be supplied explicitly when needed. No ordinary profile or
charge-method argument exists. With no completed parameter domain this returns
an auditable incomplete preparation report, not an executable partial profile.
The historical API/registry is preserved.

For returned data, first provide a complete request conforming to
`audit/reference-data.schema.json`, normalized samples and hashed native outputs
with a pinned importer/verification record. Then run:

```bash
python -m island.pcff_qm validate-reference --distribution /path/to/pinned/frc_files \
  --fixture hmds --query 'angle-angle:c3,sio,c3,c3' \
  --request returned-request.json --request-sha256 REQUEST_SHA256 \
  --data-root returned-files --verification-sha256 VERIFICATION_SHA256 --json intake.json
```

The recursive inventory must exactly match every file. Every sample must have
its exact graph/bond orders/order, charge/spin, geometry, planned displacement,
train/held-out split, method/basis/dispersion/software/version/settings, finite
energy and full force vector. A full Cartesian Hessian is required at the central
point. All planned points must be covered; optimization or reordered outputs do
not silently substitute for frozen scan geometries. Missing units, duplicates,
mixed charge/spin states, undeclared mixed methods, incompatible functional forms,
missing permission, unverified outputs or incomplete coverage reject intake.
The importer record must also preserve native energies, gradients and central
Hessians. Intake independently reconstructs their normalized values, including
the gradient sign and CODATA/SI unit factors; altered factors or signs reject.
Importer/output hashes establish integrity of a supplied attestation; they do
not independently prove that a calculation occurred or its scientific accuracy.

A later original author submission uses `attach-author --query ... --request ...
--request-sha256 ... --data-root ... --json ...`. It must supply the artifact
filename/hash/date/source identity, exact text/location, permission and independent
reconstruction. It creates a separate parent-linked original-artifact record;
the current user-attestation record remains immutable. The future fitted-parameter
schema likewise binds dataset/target/convention identities and requires
`origin="fitted_external"`; J57 creates no such parameter.

## Readiness and independent claims

Target fitting readiness has exactly six possible states: `qm_data_missing`,
`qm_data_available`, `fit_ready`, `fit_completed`, `fit_validated`, `fit_rejected`.
All 19 current targets are `qm_data_missing` and retain numerical classification
`qm_fit_required`. A successful integrity intake alone does not fit or validate
a coefficient. Residual isolation, bounds, objective/weights, identifiability,
uncertainty and held-out criteria remain explicit prerequisites.

Completeness requires the 28 bound attestations, all 19 supplied/fitted/validated
terms, complete Class-II dependencies, typing, charges, units, functional forms,
fit provenance/uncertainty and independent reconstruction. No final profile is
registered here and no zero couplings, borrowed OPLS/GAFF/GAFF2 rows, charge
normalization or harmonic simplification is permitted.

Under J55, parameter completeness, LAMMPS export/runtime, OpenMM build,
cross-engine equivalence, scientific validation and production readiness are
separate fields. The claim flags `--require-parameter-complete`,
`--require-lammps-runtime`, `--require-openmm-runtime`,
`--require-cross-engine-equivalence`, and `--require-scientific-validation` fail
independently when their claims have not passed. Cross-engine diagnostics do not
block future parameter completeness unless requested explicitly.

The bounded domain is exact disiloxane diagnostic control, HMDS, and capped
DP2/DP4/DP6 graphs. Disiloxane's Si–H is an exact supplied control, not transferable
Si–H support. No extension to cyclic/crosslinked, Si–OH, Si–O–C, phosphorus or
other terminal groups is authorized. `production_validated=False` and
`simulation_readiness="not_established"` remain fixed. No universal PCFF/COMPASS,
complete `.frc`/PDMS coverage, production MD, pressure/Ewald/PPPM validation or
long-NPT readiness is claimed.

## Evidence and verification

The receipt and reproducibility manifest pin every new implementation/evidence
file, the original main base and complete J17–J56 prerequisite tree, including
literal symlink contents. The validator reconstructs every target and generated
archive from the pinned J54 bundles/source files. Relocated sources and all five
fixture graphs and geometry checksums are reconstructed in separate processes
under two hash seeds.
The standalone `validation/independent_geometry.py` checker uses separate
angle/determinant formulas without importing the production geometry functions.
It verifies 1,823 coordinates across all 889 generated points, with maximum
residual about 4.4e-16 Å/radians. This is geometry construction evidence only,
not QM or force-field scientific validation.

Focused tests exercise graph/charge/source/convention/atom-order mutations,
geometry reconstruction, adapters, strict returned-data rejection and additive
artifact upgrades. Software sentinel outputs occur only in tests and never enter
the reference-dataset inventory.

Validation records contain focused J57 and J54–J56 regressions, the ordinary full
suite with exact collected-node accounting, Ruff, compileall, all five environment
pip checks, diff checks, historical receipt validation and complete historical
hash comparison. The final receipt reports actual counts and statuses.
