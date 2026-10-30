# Phase 4J56: controlled external Class-II PDMS completion

J56 starts from exact fetched main `5fd06ec0f25532ed164ba534b04a52a7d7a0aedb`,
the J54 merge after J55. Work is on
`codex/phase-4j56-pcff-pdms-missing-term-fitting`. All baseline tracked files,
including J17–J55 receipts, bundles, source hashes and code identities, are
preserved byte-for-byte. The implementation is additive in
`island.pcff_completion`; it does not edit native PCFF or any OPLS/Amber records.

**No reference fitting dataset was supplied. No coefficient, executable PCFF
profile or new MD run is produced.** Parameter completeness remains incomplete.
The existing family-only PCFF inspection and canonical `forcefield="opls"` API
retain their existing behavior. Completion requests are an expert workflow;
ordinary users do not need profile or charge-method arguments.

## Classification and dependencies

The exact J54 missing-row report contains 47 directed source queries covering
2,916 interaction requests in disiloxane, HMDS and TMS-capped DP2/DP4/DP6. These
are source queries and graph sites, not 2,916 independent fitted coefficients.
The new audit replays J54 source, final-graph, type and charge identities; it does
not repeat the PDF inventory or overlay different parameter versions.

| Source family | Directed queries | Classification |
| --- | ---: | --- |
| angle–angle | 12 | `requires_qm_fit` |
| bond–bond 1–3 | 7 | `requires_qm_fit` |
| end-bond–torsion | 7 | `requires_author_confirmation` |
| middle-bond–torsion | 7 | `requires_author_confirmation` |
| angle–torsion | 7 | `requires_author_confirmation` |
| angle–angle–torsion | 7 | `requires_author_confirmation` |

Every query also retains `source_status="source_missing"`. Each individual site
retains its directed labels, equilibrium dependencies, classification, reason and
source equation identity. All six requested classification values are supported;
none is forced into a category merely to populate the list. `requires_qm_fit`
is a proposed route requiring identifiable reference data, not a claim that data
exist. Alternative experimental fits would need an explicitly supported inverse
model. A rank-deficient fit rejects with `not_identifiable`.

The 28 torsion-related queries first require a source-authorized torsion geometry.
Seven additionally have the documented linear-Phi versus cos(Phi) functional-form
conflict and unresolved angular dimensions. A fit cannot resolve which physical
Hamiltonian the author intended. No backend convention replaces source authority.

There is no source equation determining any of these missing constants from the
available equilibrium lengths, angles or single-coordinate force constants.
Those values specify geometry dependencies, not mixed derivatives. Consequently
**zero current queries are `analytically_derivable`**. The derivation command
returns a checksummed blocked record with exact equation evidence, missing
coefficient relation, input dependencies and unavailable uncertainty. There is no
arbitrary expression evaluator. A future authorized derivation adapter must bind
verified input rows/units, an explicit coefficient equation, conversion/version,
uncertainty and independent reconstruction, and use `origin="derived_external"`.

The dependency closure report includes all five fixtures, exact missing sites,
source conventions, charge component totals and graph/type/charge identities.
Even a successful individual candidate fit cannot close the remaining Class-II
inventory, improper applicability, torsion–torsion applicability, authorization,
exclusion, 1–4, cutoff or long-range gaps. No silent zero, interpolation, borrowing,
charge redistribution or harmonic downgrade is permitted.

## Fitting contract and bounded adapter

`missing-term-fitting.schema.json` documents the request, reference samples and
completion records. Requests bind the candidate identity, single FRC version,
exact graph, atom ordering, typing, charge vector/origin, chemical domain, target
site/query and source equation. A recursive payload inventory pins every file,
including nested files, license text and reference producer assets. Symlinks,
path traversal, unlisted files and changed hashes reject. The contract itself is
outside the payload directory and is independently pinned by the CLI.

Reference intake accepts Cartesian geometries, energies, forces, Hessians,
vibrational frequencies with units, experimental observables with uncertainties,
QM method/basis/software and version, or reference-engine provenance. All samples
use the exact source-bound graph and atom order. Train and held-out assignments,
positive weights, explicit coefficient bounds/initial values, constraints,
objective and optimizer are recorded. Unfilled templates use null placeholders;
no initial coefficients or charges are invented.

The v1 numerical adapter is deliberately bounded: documented **isolated component
energies** for bond–bond, bond–angle, angle–angle or bond–bond 1–3 forms. Columns
are reconstructed from Cartesian geometry and verified equilibrium rows; a user
cannot supply an arbitrary design matrix. Fixed equilibria are not refitted.
Bond–bond and bond–angle rows already present in the current source are protected
from refitting; adapter support does not permit overwriting those rows.
The objective is `sum(w_i*(E_model_i-E_reference_i)^2)` on training samples.
Full-rank weighted SVD solves the linear coefficients. There must be positive
residual degrees of freedom, held-out data and explicit bounds; an optimum outside
bounds rejects rather than being clipped. Separate normal-equation reconstruction
checks the result. Initial values are recorded as supplied provenance; this direct
linear solver does not use an iterative initial guess.

A reference producer's pinned isolation description, code file, raw component
records, command, graph, exact target-site and equation identities are required. Reference values
and geometries must match those raw records exactly. ISLAND imports those as
**user-supplied evidence**; it does not certify that the referenced QM/engine
calculation occurred or that its Hamiltonian is physically accurate. Total QM
energies cannot be relabelled isolated terms by setting a Boolean flag. They need
a separately verified joint residual/decomposition adapter. Force, Hessian,
vibrational and experimental channels are accepted and retained, but have no
numerical inversion adapter in v1. Their absence does not become fabricated data;
force/Hessian-only fits remain blocked. Torsional forms have explicit record
support and source evidence but stay blocked by the J54 convention conflicts.

A fitted candidate, when sufficient external data are actually supplied, is marked
`origin="fitted_external"`, never native PCFF/COMPASS. Its immutable record contains
all dataset hashes/provenance, target/form/units, objective and weights, constraints,
bounds/initials, optimizer/version/convergence, residual-based coefficient
covariance, degrees of freedom, held-out residuals, independent reconstruction,
code identity and final parameter checksum. Statistical covariance is conditional
on the supplied data; it excludes source, QM and model error. Held-out residuals
are recorded without inventing a scientific accuracy threshold. Such a candidate
has no profile or execution authority and cannot mark the whole domain complete.

The supported ordinary domain remains finite, linear, explicit-H methyl PDMS with
TMS ends. The exact disiloxane fixture is an explicit Si–H diagnostic control,
not a transferable Si–H domain. There is no extrapolation to cyclic siloxanes,
networks, Si–OH, Si–O–C, phosphorus or unsupported end groups. Si–C–O/Si–O–C
requests need new source typing, charge evidence and reference data; none is
supplied here.

## Reproduction

Use the project validation Python or an installation with the chemistry extras:

```bash
export PYTHONPATH=.:src
python -m scripts.phase_4j56 \
  --distribution /path/to/pinned/frc_files \
  --output /tmp/j56-audit
```

The output directory must be new. This command reconstructs all five historical
bundles, classifies all 47 queries, emits dependency/provenance/uncertainty reports,
and checks relocated reconstruction in ten fresh processes with two hash seeds.
Original source paths are not part of the serialized completion identity.

For one fixture, without supplying a profile or charge method:

```bash
python -m island.pcff_completion classify \
  --bundle docs/evidence/phase_4j54/audit/hmds \
  --manifest-sha256 <pinned-J54-bundle-manifest-hash> \
  --distribution /path/to/pinned/frc_files --json /tmp/classification.json
```

The same CLI supports `template`, `derive --target-id ...`, `intake`, `fit` and
`validate`. `intake`/`fit` require `--request`, `--request-sha256` and `--data-root`.
`validate` requires `--record`, `--record-sha256` and, for data-bound records,
`--data-root`. A blocked derivation/intake or rejected fit returns exit code 2.
`template` emits an unfilled contract, not a valid dataset or fitted term. Output
cannot overwrite graph bundles, the source distribution, payloads or requests.
Mutations are checked against regenerated evidence even if a record is re-signed.

Validate the sealed receipt and all unchanged historical bytes with:

```bash
python scripts/validate_phase_4j56_receipt.py docs/evidence/phase_4j56 \
  --repo . --distribution /path/to/pinned/frc_files
```

## Validation and claims

Focused checks cover the 47-query inventory, all current missing families,
source equations, refusal to derive constants, absence/insufficiency of reference
data, rank deficiency, provenance/units/weights/bounds/license/source/graph/type/
charge mutation, all supported reference channels, schema, relocation and fresh
processes. Test sentinel arrays are software contract inputs only; no physical
fit is performed and no parameter output is retained. J54/J55 regressions and the
ordinary suite run separately; their historical integration tests retain their
original scope. Ruff, compileall, five pip checks, diff checks and historical
receipt/hash checks are recorded in the J56 receipt.

J55 remains unchanged: parameter completeness is independent of optional
cross-engine equivalence. A completed external domain would additionally need
full Class-II closure, provenance/charges and deterministic serialization before
registration. Backend static, finite-difference, minimization, export and
relocation checks apply only after such a domain genuinely exists. J56 does not
skip missing source requirements by substituting a numerical comparison.

`production_validated=False` and `simulation_readiness="not_established"` remain
unchanged. No universal PCFF/COMPASS, complete `.frc` coverage, production PDMS MD,
pressure validation, Ewald/PPPM closure, equilibrium or long-NPT readiness is claimed.

Validation: 57 focused tests and 129 J54/J55 regression tests passed. The ordinary
union contains 4,656 distinct nodes: 4,595 passed and 61 documented skips. All
5,483 baseline tracked files remain byte-for-byte unchanged.
