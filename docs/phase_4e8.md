# Phase 4E8 — traceable retained-trajectory analysis

## Scope and APIs

`island.analysis` provides `GeometryOptions`, `GeometryResult`,
`geometry_metrics()`, `AnalysisOptions`, `AnalysisReport`, `analyze_workflow()` and
`export_analysis()`. Invalid input, unsafe numerical arithmetic, invalid saved
records and export failures raise `island.exceptions.AnalysisError`, with the
underlying exception chained where applicable.

The numerical kernel requires only NumPy. It accepts explicit integer stable-ID
coordinate and mass mappings with **identical complete inventories**, and optional
atomic numbers. Coordinates are angstrom and masses dalton; alternative units
are rejected rather than guessed. Every mass must be positive and finite, even
under uniform weighting. Boolean IDs, duplicates, unknown selections, empty
selections, nonfinite coordinates and malformed vectors are rejected.

```python
from island.analysis import GeometryOptions, geometry_metrics

geometry = geometry_metrics(
    {7: (0, 0, 0), 90: (4, 0, 0)}, {7: 1.0, 90: 3.0},
    options=GeometryOptions(weighting="mass", endpoints=(7, 90)),
)
assert abs(geometry.radius_of_gyration**2 - 3.0) < 1e-12
```

Selections are `all`, `heavy` (explicit atomic number > 1), or `explicit` with
`site_ids`. Deterministic ascending IDs determine the calculation order; neither
insertion order nor coordinate array position establishes atom identity. Results
record selected IDs, their masses and normalized weights. The default is **mass
weighting**; `weighting="uniform"` assigns equal weight to each selected site.
Existing conformation diagnostics remain **uniform-site** Rg and are unchanged.
Results contain immutable tuples; `GeometryResult.to_dict()` and
`AnalysisReport.payload` return owned copies.

## Definitions and numerical policy

For selected weights w and positions r:

- center = sum(w r) / sum(w)
- G = sum(w (r-center)(r-center)^T) / sum(w)
- Rg = sqrt(trace(G))
- lambda_1 <= lambda_2 <= lambda_3 are the symmetric tensor eigenvalues
- kappa_squared = 1.5 sum(lambda_i^2) / sum(lambda_i)^2 - 0.5

Center and Rg are angstrom. G and eigenvalues are angstrom squared. Relative
anisotropy and Re/Rg are dimensionless. Collinear nonzero-spread configurations
have anisotropy 1; an isotropic tensor has anisotropy 0. Definitions follow the
official [LAMMPS gyration](https://docs.lammps.org/compute_gyration.html) and
[gyration/shape](https://docs.lammps.org/compute_gyration_shape.html) documentation.
This implementation has no periodic unwrapping.

A reference shift reduces cancellation from translation. We scale masses before
normalization and the centered coordinates before forming a symmetric tensor;
`numpy.linalg.eigvalsh` supplies real ordered eigenvalues. Only negative
eigenvalues within `64 * float64_eps * max(abs(G))` are treated as roundoff and
set to zero, with a diagnostic. The same equivalent criterion is applied before
rescaling. Substantial negative values fail. Anisotropy roundoff outside [0,1]
is corrected only within 64 eps and recorded. Positive eigenvalues are not
arbitrarily clipped. Unrepresentable mass ranges, overflow and tensor underflow
raise a domain error; scaling does not promise support for every finite input.

A single site or zero-spread selection has Rg=0. Normalized shape and Re/Rg are
**undefined**, represented as JSON null/CSV blank with a zero-spread diagnostic.
They are never reported as NaN, infinity or silently assigned zero.

## Endpoints

Explicit `endpoints=(head_id, tail_id)` are validated against the complete
inventory, independently of the Rg selection. Equal IDs and coincident distinct
endpoints are legitimate Re=0 cases, with distinct diagnostics.

For workflow analysis without an explicit pair, the adapter reuses the existing
linear-polymer provenance validator: connected finite atomistic topology,
contiguous repeat provenance, matching repeat count, terminal-repeat head/tail
membership, and the supported linear inter-repeat connections. Endpoint atoms
must additionally carry original heavy-atom repeat provenance, not generated
hydrogen provenance. It uses the recorded `head_site_id` and `tail_site_id`.
Missing, invalid or ambiguous provenance yields an unavailable Re and an explicit
reason while preserving valid Rg/shape results. No graph-path, minimum-ID or
maximum-ID fallback is used. The resolved pair and convention are recorded
separately from the selected Rg IDs. Re/Rg uses that selection's Rg.

## Validated workflow snapshot

```python
from island.analysis import AnalysisOptions, analyze_workflow, export_analysis

report = analyze_workflow("my-run", AnalysisOptions(start_step=10, stride=2))
export_analysis(report, "my-analysis")
```

The adapter reads one manifest byte sequence and verifies its content checksum.
It copies exactly that manifest's referenced immutable bytes into a private
system-temporary directory, checking relative paths, lengths, SHA-256 checksums
and the declared artifact budget. It then applies the existing deep workflow
validation and combined-frame reader to that fixed snapshot. A concurrent newer
manifest is irrelevant to this analysis. Changed/deleted old artifacts fail
validation; files that merely exist but are not referenced are ignored. No
source lock or source-directory file is written. The private copy is cleaned up.

Deep validation reconstructs the saved parameter records and minimum without
performing scientific calculations. It requires ParmEd and, where assigned
stereochemistry needs checking, RDKit. It does not need an OpenMM Context, SciPy
optimization, AmberTools executables, or velocity initialization. Missing required
validation dependencies cause an explicit error; checks are not skipped.

Phase 4E7.1 consistency rules remain authoritative. Shared synchronized boundary
coordinates, velocities and identities must agree exactly; independently
recomputed energy/component/force values use the existing tolerances (no changes).
The combined reader **retains the earlier accepted frame**. Original distinct
signed evaluation/segment records are neither rewritten nor re-signed.

Paused, budget-limited or failed workflows with valid accepted history may be
analyzed and keep their original status. Rejected diagnostic trials are excluded.
An incomplete setup with no accepted frames gives `No accepted retained frames
available`. Analysis does not invent an empty successful trajectory.

Step and time windows are inclusive and intersect. Filtering occurs first;
`stride` then selects every stride-th **retained sample**, beginning with the
first qualifying frame. Unstored frames are never interpolated. Exact selected
steps, times, coordinate fingerprints and source segments are recorded. Unequal
retained intervals are supported. Per-frame U, K, U+K and instantaneous kinetic
temperature are read from saved records (3N degrees of freedom). Langevin energy
changes are not classified as NVE conservation failures.

These are statistics/measurements of retained samples, not time-weighted or
equilibrium ensemble estimates. No automatic burn-in, standard errors,
independent-sample assumption or convergence certification is provided.

## Reports, CLI and publication

```sh
python -m island.analysis my-run --output my-analysis --selection heavy
python examples/analyze_trajectory.py my-run --output selected-analysis \
  --selection explicit --site-ids 7 90 --endpoints 7 90 \
  --weighting uniform --start-step 10 --end-step 50 --stride 2
```

Output must be separate from the source run (neither inside it nor its ancestor).
An existing output directory is refused. Exclusive directory creation reserves
the destination. Each file is fsynced and atomically linked using existing safe
publication utilities; the completion marker is published last. An exception
cleans up only this newly owned directory. A process interruption may leave an
incomplete directory: absence of a valid `COMPLETE.json` means it is not a
published export. Inspect/remove that output or choose a new destination; source
workflow recovery is unaffected. There is no overwrite or merge operation.

Files:

- `report.json`: schema `island.trajectory-analysis.v1`, implementation version 1;
  original manifest content/file hashes, run/evidence/status/readiness, original
  trajectory origin, model/parameter/backend identities, all referenced artifact
  hashes, ordered segment content identities and selected checkpoint reference;
  options, exact frame selection, geometry, selected masses/weights, endpoints,
  formulas, units and diagnostics.
- `frames.csv`: one row per included retained frame. `step`, `time_ps`,
  `source_segment`, coordinate/evaluation fingerprints, `weighting`, selected and
  endpoint IDs (space-separated), endpoint convention; U/K/total columns in
  kJ/mol and temperature in K; center xyz in angstrom; six G components and three
  ascending eigenvalues in angstrom squared; Rg/Re in angstrom, Re/Rg and
  kappa_squared dimensionless, and semicolon-separated diagnostics. Undefined
  values are **empty fields**. Full masses/provenance are in the paired JSON.
- `COMPLETE.json`: schema `island.analysis-export.v1`, SHA-256 of both report files.
  These are content-integrity checks, not computational authenticity signatures.

Exports are analysis products, **not restart files**. Original chemical,
preparation, parameter and coordinate-provenance records remain untouched.

## Acceptance

Analytical tests independently check the unequal-mass identity
Rg^2=m1*m2*d^2/(m1+m2)^2, equal weighting, isotropic and collinear limits,
translation/rotation scalar invariance and tensor covariance. Tests also cover
selection, degenerate endpoints, ownership, overflow, negative-eigenvalue bounds,
optional dependency isolation and validation failures.

Real validated synthetic saved workflows exercise partial histories, unequal
spacing, CLI exports, source non-mutation, injected positive/negative 5e-9 kJ/mol
startup offsets, and concurrent publication while a manifest is pinned. Engine
entry points are forbidden during a regression analysis to detect accidental
scientific execution. Controlled offsets/publication failures are injected
software tests, not claims of naturally observed numerical discrepancies.

`scripts/validate_trajectory_analysis.py RUN --output NEW_DIRECTORY` independently
checks selected frames using the weighted pair-distance Rg identity and
long-double direct tensor/Re calculations. Predetermined atol=rtol=1e-12 are used.
The opt-in regression uses `ISLAND_ANALYSIS_REFERENCE_ROOT` containing the two
retained Phase 4E7.1 relocated PE bundles.

Both existing capped PE DP=3 bundles were actually analyzed: previously live
GAFF2/AM1-BCC and archived GAFF2/provided-charge workflows. Each supplied seven
retained frames at steps 0,10,20,30,40,50,60 (0.006 ps). No new dynamics,
minimization, initialization or AmberTools execution was performed.

| Existing source | Maximum Rg^2 error (A^2) | Tensor error (A^2) | Re error (A) |
| --- | ---: | ---: | ---: |
| Live-origin PE | 2.990e-15 | 1.193e-15 | 9.558e-16 |
| Archived PE | 3.175e-15 | 1.068e-15 | 1.192e-15 |

Source manifest content checksums:

- Live: `5548991ca896ce1ccc22d698e74c5d0b47e575d9ded7ff2cd3b3826ce50a081f`
- Archive: `d64804ad9402d8367624fcf2b2cbb18dfd2dc1e21968dc8ffe17d58a085f17e3`

Reports remain outside the repository in `island-validation/phase4e8/`.
Historical charges, charge tolerances, signatures and force-field provenance
were preserved. Short trajectories do not establish equilibrium statistics.

## Delivery and limitations

The branch starts from merged `origin/main` at
`a563d7f27451c17659a68bac8686bdae95f2c475` (squash merge PR #15). Its
repository tree is identical to reviewed correction
`58f76762f4049a56cf9cd82df3d5f06ca21c235f`; the original feature commit is
not an ancestor because of the squash. No unmerged dependency remains.

Delivery validation:

- Complete ordinary suite: **814 passed, 9 skipped** (100.52 s). Forty new
  ordinary analysis regressions passed. Skips are the pre-existing live/archive
  opt-ins plus the new retained-bundle analysis opt-in.
- Separately enabled retained-PE analysis acceptance: **1 passed** (5.37 s),
  exercising both real saved bundles. The dedicated direct-comparison script
  also passed for each bundle and produced the errors listed above.
- Ruff: passed. `python -m pip check`: no broken requirements.
- Analysis example with heavy selection: seven saved live-origin PE samples,
  successfully exported. Cross-process checkpoint example: zero split-state
  discrepancy for NVE and Langevin, matching RNG state, 21 split versus 17
  uninterrupted calls. Existing random-walk conformation example passed with
  its unchanged uniform-site diagnostics.
- A separate process analyzed the saved live-origin bundle with OpenMM and
  SciPy imports blocked; seven samples passed deep validation. The kernel
  regression blocks RDKit, ParmEd, OpenMM and SciPy together.
- No new live AmberTools or archived dynamics acceptance was run in this
  read-only phase. Historical numerical records were analyzed, not regenerated.

Environment: Python 3.11.16, NumPy 2.4.6, ParmEd 4.3.1, RDKit 2026.03.6;
existing source backend OpenMM 8.6.1 Reference. The full test environment also
contains SciPy 1.17.1 and OpenMM 8.6.1. The kernel adds no optional dependency. Snapshotting uses temporary disk space proportional to
the referenced workflow artifacts and repeats existing deep integrity checks;
this deliberately prioritizes consistency within the current bounded workflow.
Finite nonperiodic systems only; no unwrapping, new dynamics, force fields,
packing or AmberTools scope expansion. `production_validated=False` and
`simulation_readiness="not_established"` remain mandatory.
