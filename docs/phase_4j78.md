# Phase 4J78 — bounded PCFF-compatible PDMS extension qualification

The decision is **derived_profile_blocked**. This phase creates an offline,
non-executable proposal named `island_pcff_pdms_derived_v1`, labeled
**PCFF-compatible PDMS extension**, with `profile_origin="qm_fitted_extension"`.
That origin describes the intended fitting route; no physical coefficient has
been fitted. It does not recover, complete or amend native PCFF.

## Base and preservation

The branch is `codex/phase-4j78-pcff-compatible-pdms-extension`. Both the exact
base and fetched `origin/main` are
`b4354693d6f888fda442f6944d6af7700860e639`. Main is untouched. The receipt records
the tested implementation commit separately from the final evidence commit.
The frozen base inventory preserves **10,154 historical tracked files**,
including all J1–J76 evidence and native PCFF source/registry implementation.
Unrelated untracked files and the external literature directory are untouched.
No historical identity glob is broadened. New library code is a subpackage of
`island.pcff_pdms`, so historical `pcff_pdms/*.py` hashes remain unchanged.

J76's conclusions remain: zero exact native rows resolved, 378 native Class-II
requests unresolved, zero accepted native coefficients, no universal executable
PCFF profile. J72's four unintended SQM invocations and 46 verified / 48
possible MD steps remain violations with their receipts and C1 quarantine
unchanged. The 28 author confirmations remain user-attested metadata; original
author artifacts remain `not_yet_supplied` and audited confirmations remain zero.

## Exact chemical scope

The candidate scope is the four retained J54 finite, linear, trimethylsilyl-capped
graphs: HMDS and capped DP2, DP4 and DP6 (27, 37, 57 and 77 explicit atoms).
Their stable atom IDs, original serialized input hashes, graph identities,
coordinates and neighborhood mappings are frozen as data. The DP labels retain
their historical graph meanings; no homologous-length extrapolation is claimed.
Only neutral closed-shell singlets with explicit H, single Si–O–Si bridges,
Si–C methyl attachments and the exact declared graph environments are eligible
for this proposal. The systems are nonperiodic finite molecules.

Si–H, Si–OH, Si–O–C, rings, crosslinks, arbitrary silicon polymers, P, Br, Se,
charged/radical systems and periodic models fail scope checks. Si–C–C is absent
from these graphs: every carbon is methyl bonded to Si. It is reported as
structurally absent rather than supplied a guessed angle. Wilson impropers
have no three-connected center. Torsion–torsion applicability remains a required
unresolved policy; an empty historical section does not establish a zero term.

## Existing code reused

`implementation-reuse.json` inventories the inspected original modules and their
hashes, including `pcff_pdms`, `pcff_qm`, `pcff_qm_pilot`, and the native PCFF
operational profile, registry, resolution and coverage modules.

The extension reuses `pcff_pdms.domain_diagnostics`, core topology enumeration,
`pcff_qm.records.Record`/canonical identities, the proposed J57 protocol metadata,
`pcff_qm_pilot.hmds_design.metrics`, and the established pilot's
`scipy.optimize.lsq_linear` backend. The relocation test reuses the pilot's
data-only archive pack/unpack functions. No molecular evaluator is added.

Historical fitting entry points bind a 58-column/27-atom HMDS problem or a
9-atom source-bound Hamiltonian. They cannot accept this inventory without
changing historical contracts. The new staged adapter accepts supplied numeric
matrices and reuses the existing solver and diagnostics. It does not duplicate
the molecular fitting/evaluation framework or create a second parameter model.

The explicit new implementation inventory is five library modules under
`src/island/pcff_pdms/extension/`, two scripts, the focused test module and this
document. The receipt lists every path and SHA-256.

## Terms, parameters and source separation

`unique-term-inventory.json` separately enumerates parameter identities and
stable site-instance identities over all four graphs. It includes quartic bonds
and angles, proper torsions, bond–bond, bond–angle, angle–angle, first/third-bond
couplings, end/middle-bond–torsion, angle–torsion, AAT, partial charges and 9–6 LJ
targets. Atom typing, charge origin, mixing, exclusion/scaling, 1–2/1–3/1–4,
long-range, cutoff, tail and torsion–torsion applicability policies are explicit.

New `ext_*` types identify terminal/interior Si and their methyl environments.
They are extension types, not author-authorized FRC types. Reversal-equivalent
terms share one identity. Symmetric left/right coefficient aliases map to one
variable; they are not fitted independently. Angle–angle retains its center and
shared-arm roles, allowing only the documented outer-site exchange.
Every coupling instance retains its required reference bond/angle sites.
Reference values and coefficients remain unresolved, not zero-filled.

The scalar equations, units, phases, polynomial prefactors and ordering in this
inventory are **proposed derived-model contracts**. For example, a cosine AAT
proposal is explicitly separate from unresolved native Phi/AAT authority. It
does not resolve a native source convention. All candidate values are `null`.
Native PCFF, COMPASS, CFF91, OPLS-AA, GAFF/GAFF2, external literature PDMS and
analog values remain separate. No literature or analog numerical value is
imported, averaged or substituted. Source hashes are provenance context only.

## QM target and split contracts

`qm-target-manifest.json` contains deterministic *requests*, not QM observations
or runnable engine decks. Each binds its graph, stable atom/coordinate order,
atom mapping, seed checksum, method, basis, dispersion, charge/multiplicity,
convergence policy, constraint design, role and input-specification hash.
Output hashes are explicitly absent because nothing was executed. Software
version and actual independent constrained geometries remain required before
output intake. A missing version or original native input/output hash blocks
intake. No engine is probed or launched.

The plan supports equilibrium geometries, Si–O/Si–C/C–H bond scans, O–Si–O and
Si–O–Si angles, Si–O–Si–O/C–Si–O–Si and methyl rotations, paired signed
deformations, stationary gradients and Hessian/projected curvature targets.
Charge/LJ reference design is separately required; monomer energies alone do
not identify a transferable charge/nonbonded model.

PBE0/def2-SVP without dispersion is retained as a **proposed** historical input
policy, with method adequacy explicitly unestablished. It is not accepted as a
scientifically converged reference protocol. A software/version/execution
profile and adequacy study are still dependencies, not silent defaults.

Training, validation and test request IDs are assigned before fitting. Their
seed copies are the same historical unoptimized inputs, not independent data.
The split manifest therefore reports no completed samples and
`independence_verified=false`. Completed data must prove disjoint geometry and
conformer blocks. J68's previously used 26 observations are not relabeled as
independent J78 held-out data. Reserved validation is not activated.

## Fitting and acceptance

The staged contract freezes primary terms before torsions, those terms before
orthogonal Class-II cross-term fitting, and treats charges/LJ separately.
Optional refinement requires all earlier gates. There is no joint all-term fit.
Physical parameter ordering is immutable. Finite numeric physical bounds must
be supplied and frozen before real fitting; this phase invents no such bounds.
An unfilled bound blocks fitting, and a bound/order/hash mutation fails intake.
Nonlinear reference geometry, phase and charge/LJ fits are not silently treated
as affine columns. They require an identified design and declared future stage.

Synthetic matrix controls use explicit finite bounds and the existing bounded
solver, with deterministic column scaling, SVD rank/singular values, condition,
correlation and training-only target perturbation diagnostics. They preserve
solver states and all perturbation attempts in their returned records. These
controls test software only, never supply PDMS coefficients or claim physical
qualification. A marker that says “synthetic” cannot create a registered profile.

The conservative numerical limits are declared research gates, not widened to
force acceptance: energy RMSE 0.1 kcal/mol, Cartesian force RMS 0.5
kcal/mol/angstrom, maximum energy error 0.3 kcal/mol and maximum local vector
force error 2.0 kcal/mol/angstrom. Proposed geometry gates are 0.02 angstrom and
0.04 radian; stationarity RMS 0.05 kcal/mol/angstrom; torsion RMSE 0.1 kcal/mol;
relative curvature error 0.1; scaled condition at most 1,000 with full rank;
95% conditional sensitivity at most 0.1 of each bound width. Threshold changes
are semantic mutations. Model/QM/source uncertainty requires separate evidence.

The gate report additionally requires complete unique terms/dependencies,
independent training/held-out observations, stable geometries and curvature,
charge conservation, full Class-II/exclusion conventions, immutable native
files, provenance and independent replay. The real candidate has no training
error, held-out error, design matrix, rank, condition or uncertainty estimate:
those fields are `null` with deterministic blockers. They are not numerical
passes. No fitted or executable coefficients are accepted.

## API, verification and replay

The extension's `prepare_polymer_forcefield(final_system_payload,
forcefield="pcff", ...)` accepts the existing family spelling without requiring
an internal profile name. It validates scope and provenance and then fails
closed because no derived bundle is qualified or registered. Historical
`prepare_polymer`, ordinary selection and native registry behavior are unchanged.
No PCFF runtime adapter or automatic assignment path is added. A checksum or
re-signed passing report cannot grant registry authority.

All verification runs through the inherited J72-C1 process/kernel fence with
J73–J75 evidence-only restrictions. The fence is verified before evidence
access. Synthetic child processes use fake marker executables only. Real QM,
MD, LAMMPS, OpenMM, msi2lmp, converters, physical fitting, parameter generation,
scientific evaluation and reserved validation counts are zero. Explicit live
controls are scope skips, never successful scientific checks.

The terminal receipt separates focused coverage from regressions and records
unique node IDs without adding overlapping reruns. Ruff, compileall, supported
pip check, historical receipt validation, source/evidence byte preservation,
relocation and strict replay are recorded separately. Broader historical live
scientific suites are not run. Failed attempts and corrections are retained:
an initial metadata inspection mistook a packed system envelope for plain JSON,
and the initial Ruff pass required binding loop variables in a local collector.
The first focused run rejected periodic input through the loader instead of
the structured scope diagnostic; this was corrected before terminal testing.
None of these attempts launched scientific execution. Historical PDF-parser
environment pip conflicts remain disclosed in J75 and are not changed by this phase.

Use a fresh fence-output directory for every command:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j78-replay -- -m scripts.phase_4j78_extension replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j78-strict -- -m scripts.validate_phase_4j78_receipt
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j78-tests -- -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j73.py tests/test_phase_4j74.py tests/test_phase_4j75.py tests/test_phase_4j76.py tests/test_phase_4j78.py
```

Offline replay uses frozen graph/provenance data and explicit source-code
dependencies, without external source assets, network or scientific backends.
The relocated bundle contains metadata/code only and replays in a fresh fenced
process. Strict receipt validation additionally checks the historical checkout
and tested Git blobs. Synthetic controls require SciPy as the established
pilot's optional numerical dependency; real manifest replay does not.

The next dependency is an explicitly authorized, protocol-pinned independent
QM/reference campaign and complete derived-model bounds/conventions. This phase
establishes software contracts and reports missing evidence. It establishes no
chemical transfer or scientific/experimental/condensed-phase validation.
`native_pcff_complete=false`, `executable=false`,
`automatic_assignment_enabled=false`, `accepted_coefficients=0`,
`production_validated=false`, and `simulation_readiness="not_established"`.
