# Phase 4J56: author-confirmed Class-II PDMS planning and QM intake

This corrected J56 uses exact fetched `origin/main`
`5fd06ec0f25532ed164ba534b04a52a7d7a0aedb` (J54 after J55), on
`codex/phase-4j56-pcff-pdms-author-confirmed-fitting`. J17–J55 tracked files,
receipts, bundles and scientific claims remain byte-for-byte unchanged.

The active 47-query inventory now reports **28 `author_confirmed` and 19
`qm_fit_required`**. The confirmed items appear in the effective execution
contract and are excluded from the missing-parameter blocker list. No item is
marked permanently unresolved. Historical J54 reports are retained as historical
observations; they are not rewritten.

| Interaction family | Directed queries | Corrected planning status |
| --- | ---: | --- |
| End-bond–torsion | 7 | `author_confirmed` |
| Middle-bond–torsion | 7 | `author_confirmed` |
| Angle–torsion | 7 | `author_confirmed` |
| Angle–angle–torsion | 7 | `author_confirmed` |
| Angle–angle | 12 | `qm_fit_required` |
| First/third bond–bond | 7 | `qm_fit_required` |

These are directed source queries spanning 2,916 historical graph-site requests,
not 47 independent fitted numerical coefficients. The original 19 numerical
queries consist of angle–angle and first/third bond–bond couplings. Existing
bond–angle, Si–O–Si proper-torsion and other present rows are not refitted merely
to populate a requested family list. No coupling is silently set to zero.

## Confirmation provenance and its limits

The supplied correction states that the author confirmed the 28 interpretations.
J56 preserves the exact relevant user message excerpt, its SHA-256, recording
date (2026-10-30), scope, affected families and evidence identity. Each confirmed
query binds that evidence. The status is explicitly
`user_reported_author_confirmation`; this is not represented as an original
email or signed author document.

**No separate primary author artifact, author confirmation date, confirmed
coefficient values or precise interpreted equations were present in the supplied
files.** Those fields remain null. The recording date is not substituted for the
author's date. The confirmation validator can bind a subsequently supplied primary
artifact using a recursive inventory, independently pinned record, exact source
hashes, query coverage, source location, verbatim text, interpreted value, units
and evidence checksum. It does not authenticate an email sender or confer a
redistribution license.

The 28 confirmations authorize planning under the user's correction. They do not
create missing numerical rows or silently select between conflicting equations.
The effective contract exposes the absence of the primary numerical
interpretations. A new executable Hamiltonian cannot be reconstructed from the
user's status assertion alone. The 19 actual QM targets independently keep
parameter completeness incomplete. No final PCFF profile is registered and no
PCFF MD or QM calculation is run.

## Immutable QM targets

`qm-targets.json` contains exactly 19 query targets, with separate immutable
bindings for every applicable bounded fixture: disiloxane, HMDS and TMS-capped
DP2, DP4 and DP6. Each target includes graph, source, type and charge identities,
atom ordering, exact sites, equilibrium dependencies, raw source equation with
file hash and lines, coefficient order/dimensions, required paired deformation,
observables, identifiability and held-out requirements, and uncertainty fields.
`QMTarget` stores canonical JSON in a frozen record; reconstruction regenerates
the expected source-bound target and rejects even re-signed mutations.

The provided sources and user message do not prescribe a QM method, basis,
dispersion treatment, convergence threshold, deformation range, held-out error
threshold or expected numerical uncertainty. Each is an explicit null planning
gap. A configured plan is a new immutable revision requiring all these settings;
no common method or arbitrary range is silently selected. Uncertainty must
separate statistical, QM, model and identifiability contributions. No uncertainty
number or reference output is invented.

The six reporting states are `author_confirmed`, `qm_fit_required`,
`derived_external`, `fit_completed`, `fit_validated`, and `parameter_incomplete`.
Only the first two item states and the incomplete domain state apply to this
release. A future candidate fit cannot claim `fit_validated` merely because an
optimizer converged. There is no defensible coefficient derivation equation in
the supplied rows: equilibrium dependencies do not determine mixed derivatives.
The derivation route continues to return exact blocked evidence without a value.

## Reference intake and fitting

The new reference schema and intake require a recursive SHA-256 inventory,
license/permission, graph/bond-order identity, exact atom order, source-bound
charges, charge state, coordinates, units, method/basis/dispersion, software and
version, and raw output verification. Optional forces must be complete Nx3;
Hessians must be 3N by 3N; frequencies have declared units. Nonfinite values,
duplicate geometries, undeclared mixed QM methods, mixed charge states, missing
graphs/units and unverified or changed output files reject.

An independently pinned output-verification file binds every normalized sample
to its raw engine output hash and the importer code hash/version. This verifies
the integrity of user-supplied evidence, not that a physical calculation occurred
or that its results are scientifically accurate. The exact normalized data and
all raw files remain retained. A declared mixed-method protocol must be explicit;
it does not implicitly make disparate calculations interchangeable.

The existing bounded fitting adapter accepts only identifiable, independently
isolated component energies for supported non-torsional coupling equations. It
retains objective, weights, bounds, constraints, initials, SVD rank/convergence,
statistical covariance, held-out residuals and a separate reconstruction.
QM requests must additionally pass the new raw-output intake and match its data.
Total QM energies still require a verified joint residual/decomposition adapter;
forces and Hessians are retained but no unsupported inverse solver is claimed.
Fitted candidates use `origin="fitted_external"`; permissible future derivations
use `origin="derived_external"`. Neither creates native PCFF/COMPASS values.
There are no physical datasets or fitted parameter outputs in this phase.

## Input templates and commands

The explicit Psi4, PySCF, ORCA, CP2K and Gaussian adapters emit XYZ geometries and
**visibly unfilled, non-executable templates**, together with the exact supplied
plan and deformation manifest. Bond, angle and torsion scans move graph-separated
rigid fragments using explicit displacements/units. Cyclic cuts and degenerate
axes reject. Missing executables report `asset_unavailable`; templates do not
claim outputs or execute programs. Engine-specific method, dispersion and
convergence syntax still needs expert completion.

The corrected family-only planning API is additive, so historical J54 code hashes
remain intact:

```python
from island.pcff_completion.readiness import prepare_polymer_forcefield

# The pinned distribution is installed/configured once by the administrator.
result = prepare_polymer_forcefield(system, forcefield="pcff")
report = result.validate_integrity(system)
```

The result retains its original J54 assessment and reports the corrected effective
contract and active numerical blockers. It does not register a parameter profile.
The readiness CLI also supports J55's independent `--require-parameter-complete`,
`--require-lammps-runtime`, `--require-openmm-runtime`,
`--require-cross-engine-equivalence` and `--require-scientific-validation` flags.
An explicit claim failure is listed independently in `claim_failures`.

Inspect readiness with only the public family (no profile or charge method):

```bash
python -m island.pcff_completion.readiness readiness \
  --forcefield pcff --fixture hmds \
  --bundle docs/evidence/phase_4j54/audit/hmds \
  --manifest-sha256 <exact-J54-manifest-sha256> \
  --distribution /path/to/pinned/frc_files --json readiness.json
```

An incomplete readiness report is written and returns exit code 2. The same CLI
supports `confirmations`, `reference-intake` and `qm-templates`. A primary
confirmation uses `--request`, `--request-sha256` and `--data-root`; its query set
must match the selected fixture. Reference intake additionally needs `--query`
and an independently pinned `--verification-sha256`. Templates need `--query`,
`--engine`, `--settings`, `--scan` and a new `--templates-dir`. Inputs and graph
bundles cannot be overwritten. The lower-level `island.pcff_completion` CLI
retains `classify`, `template`, `derive`, `intake`, `fit` and `validate`.

Reproduce all fixture plans and ten relocated fresh-process classification checks:

```bash
PYTHONPATH=.:src python -m scripts.phase_4j56 \
  --distribution /path/to/pinned/frc_files --output /tmp/new-j56-audit
```

Validate the sealed receipt, generated plans, implementation membership, exact
suite membership and all historical bytes:

```bash
PYTHONPATH=.:src python scripts/validate_phase_4j56_receipt.py \
  docs/evidence/phase_4j56 --repo . --distribution /path/to/pinned/frc_files
```

The ordinary source/typing/charge inspection remains `forcefield="pcff"`; source
installation is administrator configuration. J54/J55 code and registry identities
are untouched. The isolated disiloxane Si–H control grants no transferable Si–H
domain. No cyclic, crosslinked, Si–OH, Si–O–C or phosphorus coverage is added.

## Completeness and verification

The J55 independent status fields remain separate: typing and candidate charges
are complete; parameter assignment, backend export/build/runtime remain blocked
by the 19 numerical targets and lack of an executable complete term set.
Cross-engine equivalence is `not_claimed`. An optional comparison failure cannot
invalidate an otherwise complete source-bound profile, and a comparison cannot
supply missing coefficients or conventions.

Registration requires the full author-confirmed numerical contract, all 19 terms
with validated fit/derivation provenance, full Class-II closure, exact typing and
charges, units/forms and independent reconstruction. None is bypassed here.
Focused tests use clearly labelled software sentinel records only. They do not
produce QM data or scientific parameters. Validation counts and all requested
check results are recorded in the receipt.

`production_validated=False` and `simulation_readiness="not_established"`.
No universal PCFF/COMPASS or complete `.frc` coverage, production PDMS MD, pressure
validation, Ewald/PPPM closure or long-NPT readiness is claimed.

Validation: 116 focused tests and 129 J54/J55 regression tests passed. The full
suite union contains 4,715 distinct nodes: 4,654 passed and 61 documented skips.
Ruff, compileall, all five pip checks and the historical receipt checks passed.
All 5,483 baseline tracked files remain byte-for-byte unchanged.
