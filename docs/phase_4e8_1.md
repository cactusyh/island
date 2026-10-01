# Phase 4E8.1 — validated reports and consistent input errors

## Repository and reproduced defects

The fetched remote still had Phase 4E8 unmerged. This correction remains on
`codex/phase-4e8-trajectory-analysis`, based on reviewed commit
`60b5abec4860bba94e6891ede1a562ae41e59183`. Remote main was still
`a563d7f27451c17659a68bac8686bdae95f2c475`. Main and unrelated notebook checkpoint
directories were not modified.

Before the fix, the new focused regression module produced **nine failures**:

- Independently changing sample count to 999, Rg to -10, or readiness to True
  still published the report and completion marker (three cases).
- Invalid JSON and missing fields leaked JSONDecodeError/KeyError (two cases).
- Mixed Python and NumPy boolean coordinates became numeric values (two cases).
- A 10**400 step bound overflowed in `math.isfinite()` (one case). Time bounds
  used the same unsafe conversion path.
- The valid reconstruction test exposed the absent `validate_integrity()` API.

The cause was an unchecked JSON-string wrapper at the publication boundary,
combined with checking coordinate dtype only after coercion and applying a
floating-point finite check indiscriminately to integers.

## Report validation contract

`AnalysisReport.validate_integrity()` checks the existing
`island.trajectory-analysis.v1` / implementation version `1` payload.
`AnalysisReport.payload` uses the same validator and returns a fresh owned copy.
`analyze_workflow()` validates its result before returning. `export_analysis()`
obtains a validated payload **before path creation or publication**.
Malformed records raise `AnalysisError`; implementation exceptions retain their
useful chained causes. Constructing the JSON wrapper itself remains inexpensive;
validation is explicit or occurs on payload access/export.

The validator checks:

- Strict JSON, including duplicate keys, NaN/Infinity and exponent overflow such
  as 1e999; exact required fields and supported schema/conventions.
- Nonempty frames and exact sample count; integer increasing steps, finite
  increasing times, source progress, inclusive selection windows and necessary
  stride spacing. No assumption of uniform retained sampling is added.
- Source artifact descriptors, relative names and identities; unique listed
  segment names; frame references in source order; checkpoint membership;
  trajectory-origin frame content, units and identities via the existing
  dynamics validators; agreement of source model/parameter/backend identities.
- Ascending unique selected integer IDs belonging to the recorded origin
  inventory; complete all-site coverage or matching explicit selection;
  positive finite masses and weights, equal array lengths, normalized mass or
  uniform weights, and unchanged selection/masses/endpoints across frames.
- Finite center/tensor/eigenvalues and dimensions/units, tensor symmetry,
  nonnegative ascending eigenvalues, nonnegative Rg/Re, anisotropy in [0,1],
  and consistent tensor/eigenvalue/Rg/shape/ratio arithmetic.
- Zero-spread nulls and diagnostics; single-site zero spread; unavailable-endpoint
  nulls and explanations; explicit or provenance endpoint conventions; equal-ID
  and coincident-endpoint Re=0 diagnostics.
- U+K and 3N kinetic-temperature arithmetic, origin energy agreement, and origin
  kinetic energy versus saved velocities when the selection supplies all masses.
- Both report and source readiness remain `production_validated=False` and
  `simulation_readiness="not_established"`.

### Numerical policy

Report validation does not change kernel definitions or workflow verification
tolerances. Internal geometry comparisons use **256 float64 eps** relative to
`max(abs(G))`, with tensor comparisons performed after scaling. This allows
roundoff from constructing, diagonalizing, serializing and independently
checking the tensor without a fixed angstrom-scale absolute floor. The existing
64-eps bound for tiny negative computed eigenvalues remains; substantial
indefiniteness fails. Rg and Re/Rg comparisons use 256-eps relative tolerance;
anisotropy and normalized-weight summation allow 256-eps absolute roundoff.
Stored eigenvalues, Rg/Re and anisotropy still must satisfy their sign/range rules.
Zero tensor/spread and required nulls are checked explicitly.

Energy and kinetic-temperature arithmetic retain the existing dynamics policy:
relative tolerance 1e-12, absolute tolerance 1e-10 in their respective units.
Fresh-evaluation energy/force tolerances and boundary deduplication are unchanged.

### Limits of validation

This is **internal consistency, not computational authenticity**. No source
workflow files are opened by report validation, so an already generated report
can be exported after the source is moved or deleted. Recorded source hashes
are checked for structure and cross-reference consistency; unavailable source
bytes cannot be authenticated. File-byte and signed-payload hashes are distinct
identities and are not compared as if they used the same normalization.

The v1 report includes the original trajectory state, but not each analyzed
frame's coordinates, all unselected atom masses/elements, omitted retained
frames or segment payloads. Consequently standalone validation cannot recompute
all geometry, establish heavy-atom membership, prove exact stride selection from
the omitted samples, reconstruct segment lineage, or prove the original
trajectory fingerprint from absent physical settings. The existing deep workflow
validation performs its checks during generation. No new data was invented to
extend those guarantees during standalone export. A mutually consistent caller
fabrication cannot be distinguished from computed data by these checks.

## Input and publication boundaries

Raw coordinate sequence elements are inspected for both Python `bool` and
`numpy.bool_` before NumPy conversion. Ordinary numeric sequences and numeric
arrays continue to work. An already floating-point array is interpreted as
numeric data; its earlier construction history cannot be inferred.

Step bounds use integer-only sign/type validation. Arbitrarily large nonnegative
Python integers, including 10**400, are supported without conversion to float.
A window beyond available steps simply selects no frames and produces the
existing domain error. Time bounds require representable finite numerical
values; huge time integers raise `AnalysisError` with the overflow chained.

Report and path validation are inside the public export error boundary.
Exclusive output-directory ownership, atomic per-file publication and publishing
`COMPLETE.json` last remain unchanged. Existing destinations are refused and
left untouched. If publication fails, cleanup is attempted only for the newly
owned output; a cleanup failure is attached as a note to the original failure,
which remains the chained cause of `AnalysisError`. No successful marker is
fabricated. As before, an interrupted/incomplete export requires manual inspection
or a new destination; source workflow records remain unchanged.

## Compatibility and verification

There is **no schema migration, signature change, readiness promotion, geometry
change, tolerance relaxation or scientific-engine change**. Previously generated
valid v1 reports remain exportable. Previously accepted malformed or contradictory
reports now fail. `payload` access now performs validation rather than exposing
unchecked parsed content. Optional scientific dependencies remain lazy; standalone
report validation uses core data validators and NumPy without scientific engines.

The two retained Phase 4E7.1 capped-PE bundles were analyzed read-only. Both their
old Phase 4E8 report reconstructions and newly generated reports passed validation;
the entire decoded new payload was **exactly equal** to its old Phase 4E8 payload.
The existing independent acceptance script retained atol=rtol=1e-12:

| Saved source | Samples | Maximum Rg² error (Å²) | Tensor error (Å²) | Re error (Å) |
| --- | ---: | ---: | ---: | ---: |
| Live-origin GAFF2/AM1-BCC PE | 7 | 2.990e-15 | 1.193e-15 | 9.558e-16 |
| Archived GAFF2/provided-charge PE | 7 | 3.175e-15 | 1.068e-15 | 1.192e-15 |

These are existing steps 0,10,20,30,40,50,60. No PE dynamics, minimization,
initialization or parameterization was rerun. Outputs and comparison records are
outside the repository in `island-validation/phase4e8-1/`. Short saved trajectories
still do not establish equilibration or equilibrium statistics.

Actual delivery checks:

- Complete ordinary suite: **828 passed, 9 skipped** (107.94 s), including 14
  new regression tests. The nine skips are the existing opt-in live/archive
  checks, including the retained-analysis module; the two PE analysis acceptance
  paths were executed separately as described above.
- Targeted geometry/workflow/initial-regression run: 49 passed. Expanded focused
  report/input regression module: 14 passed, including 22 independent malformed
  field mutations, strict JSON, degenerate reports, removed sources, existing
  destinations, and injected publication/cleanup failures.
- Ruff passed; `python -m pip check` reported no broken requirements.
- Actual analysis CLI on a validated saved synthetic workflow: three retained
  samples exported successfully. The regression suite also exercises the CLI
  in a separate process.
- Cross-process checkpoint example passed: NVE and Langevin each had zero
  split-state difference, matching RNG state, and 21 split / 17 uninterrupted
  evaluator calls. These are analytical example trajectories, not regenerated
  PE acceptance runs. Existing conformation example passed unchanged.
- Standalone old-v1 report validation succeeded in a separate process with
  OpenMM, SciPy, RDKit and ParmEd imports explicitly blocked. Kernel optional
  isolation, coherent snapshots and tolerated-boundary regressions also passed.

Environment: Python 3.11.16, NumPy 2.4.6, ParmEd 4.3.1, RDKit 2026.03.6,
OpenMM 8.6.1 and SciPy 1.17.1. No scientific dependencies were added. Injected
malformed records and failures are software tests; the read-only PE comparisons
use actual previously retained scientific records. No new equilibration or
production-readiness claim is made.
