# Phase 4E2: bounded nonperiodic local minimization

## Delivery and scope

Branch: `codex/phase-4e2-local-minimization`. Base:
`a2dece53fbd7e4bb73682deae979ba8cb472592c` (`origin/main`, Phase 4E1 merged
through PR #9). Its tree matches reviewed Phase 4E1 `b3a4547` exactly; the
squash merge does not retain that commit as an ancestor. Main is not modified
or merged by this phase.

`island.minimization` exports `MinimizationOptions`, `MinimizationStep`,
`MinimizationResult`, and `minimize_geometry(system, evaluator, options=None)`.
The optimizer consumes `PotentialEvaluator`; parameter assignment, Amber parsing,
and OpenMM force construction remain in their existing modules. The production
entry point is a previously bound `OpenMMSinglePointEvaluator`. Its new
`validate_system()` checks actual imported parameter integrity, graph, and stable
IDs before optimization. A mutable parameter snapshot is not an authority.

```python
from island.evaluation import OpenMMSinglePointEvaluator
from island.minimization import MinimizationOptions, minimize_geometry

# Bind preparation against its original system, before changing coordinates.
evaluator = OpenMMSinglePointEvaluator.from_preparation(original, preparation)
result = minimize_geometry(starting_system, evaluator,
                           MinimizationOptions(force_tolerance=0.1))
print(result.initial_energy, result.final_energy, result.final_fmax)
print(result.termination_reason, result.iterations, result.evaluations)
if result.converged:
    optimized_system = result.to_system(starting_system)
```

The starting system must have the same authoritative graph and stable IDs as the
bound evaluator. Coordinates may differ. This phase supports finite, atomistic,
nonperiodic, unconstrained systems and the existing restricted Amber physical
model only. No AmberTools operation occurs inside minimization. Its existing
100-site preparation scope is unchanged.

## Units and convergence

Cartesian variables use deterministic ascending stable-ID order and **angstroms**.
Energy is **kJ/mol** and force is **kJ/(mol*angstrom)**. The analytical gradient
passed to SciPy is `-force`; no additional nm conversion occurs here.

For N sites:

- `fmax = max_i sqrt(F_ix² + F_iy² + F_iz²)`.
- `rms_force = sqrt(sum_i (F_ix² + F_iy² + F_iz²) / N)`.
  This is the RMS atomic vector norm, not the RMS of 3N scalar components.

The public force threshold defaults to 0.1 kJ/(mol*angstrom). SciPy's Cartesian
component `gtol` is set to this threshold divided by sqrt(3), while accepted
iterations are checked directly using the public norm. See the official
[SciPy L-BFGS-B options](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html).
SciPy energy stopping (`energy_change_tolerance`, default `1e-12`, relative `ftol`)
is not force convergence. A result is converged only if independently evaluated
returned coordinates have finite energy/forces, acceptable assigned stereo,
fmax at or below tolerance, and energy no higher than the initial energy plus
`energy_increase_tolerance` (default absolute `1e-8` kJ/mol).

`optimizer_message` retains the native termination message when SciPy returns.
When ISLAND stops via its callback or budget guard, it explicitly records that
message and `optimizer_success=None`; it does not invent a SciPy success flag.
The final `converged` flag and `termination_reason` are independently determined.

## Budgets, failures, and ownership

Options require finite positive integer `max_iterations` (500),
`max_evaluations` (2000), and `max_line_search_steps` (20). Every actual evaluator
call, including failed calls, initial checks, and the final independent call,
is counted. A cache avoids repeating the initial objective call; it adds no
imaginary evaluations. One slot is reserved for final verification. With a budget
of one, only initial diagnostic evaluation is possible and convergence is false.
No native SciPy call may overrun the wrapper's evaluator budget.

Termination reasons include `force_converged`, `maximum_iterations`,
`maximum_evaluations`, `energy_stagnation`, `line_search_failed`,
`invalid_geometry`, `evaluation_failed`, and `stereochemistry_changed`.
Malformed inputs raise `MinimizationInputError`; missing required dependencies
raise `MinimizationUnavailableError`. Failed evaluations never become zero
gradients or fabricated energies. If the initial evaluation fails, energy/force
metrics are absent and coordinates cannot be applied as an evaluated result.

Successful callback termination returns the accepted iterate. Other termination
paths return the **lowest-energy admissible evaluated trial**, which need not
have been accepted by SciPy. `returned_state` explicitly identifies this choice.
The independent final call is at exactly these returned coordinates. If it fails,
the cached best admissible evaluation remains diagnostic data and
`final_evaluation_verified=False`. The compact history contains the initial state
and accepted iterations only, each with energy, force metrics, call count, and
coordinate fingerprint; it excludes unaccepted line-search trials.

Inputs and parameter/preparation records are never mutated. Results own immutable
coordinate mappings and immutable evaluation records. Model and parameter
fingerprints, units, exact force-site coverage, coordinate fingerprints, and
finite output are checked on every call. The evaluator retains its existing
per-call private Context lifecycle; no shared mutable Context optimization was
introduced.

`result.to_system(system)` checks chemistry, site IDs, repeat metadata, and
provenance against the input identity, then deep-copies the system. Applying an
unconverged diagnostic requires `allow_unconverged=True`. This choice cannot
convert its status to convergence. Current coordinate provenance is synchronized across root and polymer source/
generation fields, with superseded records retained in `coordinate_history`.
Application validates result content and independently checks coordinate stereo;
see [Phase 4E2.1](phase_4e2_1.md). Historical AmberTools preparation records remain
unchanged. Rebinding preparation from optimized coordinates is
incorrect: validate it against the original preparation system first.

## Stereochemistry and geometry

The AmberTools coordinate CIP validator has been extracted into shared
`island.chemistry.coordinate_stereo` logic. Initial, trial, and returned coordinates
are checked for explicitly assigned tetrahedral centers: stored tags are removed
before RDKit assigns CIP from 3D geometry. Assigned planar/degenerate tetrahedra
are rejected using the existing absolute determinant threshold of `1e-3 Å³`.
A rejected trial terminates honestly and returns the best prior admissible state;
`stereochemistry` describes those returned coordinates. No graph is inferred
from distances. OpenMM's existing finite-coordinate and singular-geometry guards
remain authoritative for the evaluated model.

RDKit is required if assigned centers exist; missing RDKit never silently skips
these checks. Core, evaluation, and minimization imports do not require SciPy,
OpenMM, RDKit, or ParmEd. SciPy is loaded only when minimization is requested:
`pip install '.[minimization,evaluation,amber]'`, adding `chemistry` for assigned
stereo. Graph preservation and tetrahedral validation do **not** establish absence
of bond-through-ring intersections.

## Numerical acceptance

The retained Phase 4D2.1 archive was loaded through an explicit manifest argument.
Every retained artifact checksum and historical preparation/import signature was
verified, with schema-defined integer keys restored. No historical record was
re-signed. All five cases used deterministic seed 20260929 Gaussian perturbations
of 0.025 Å per component, force tolerance 0.1, 500 iterations, 2000 evaluations,
and 20 line-search steps. All converged at the originally declared tolerance.

| Actual case | Initial → final E (kJ/mol) | Final fmax | Iterations / calls |
| --- | ---: | ---: | ---: |
| Phenol GAFF AM1-BCC | 32.292667 → -51.660061 | 0.082590 | 68 / 75 |
| Phenol GAFF2 AM1-BCC | 27.376523 → -54.405157 | 0.081282 | 67 / 74 |
| Phenol GAFF2 provided | 93.415463 → 10.763386 | 0.069584 | 65 / 70 |
| Capped PE DP=3 provided | 84.754409 → 4.088709 | 0.079889 | 135 / 143 |
| Assigned Cl/Br provided | 13.202373 → 0.173058 | 0.058959 | 11 / 15 |

All force metrics use kJ/(mol*angstrom). Full initial/final fmax and RMS values,
coordinates, fingerprints, counts, options, versions, source checksums, and
charge residuals are in [the acceptance report](references/phase_4e2/acceptance.json).
The assigned Cl/Br case passed independent coordinate stereo validation.
Each returned frame was evaluated again through both the bound evaluator and
OpenMM's independent original-prmtop reader. Energy tolerances are
`rtol=2e-10, atol=2e-7 kJ/mol`; forces use
`rtol=2e-9, atol=2e-6 kJ/(mol*angstrom)`. The largest observed absolute
energy discrepancy was `2.11e-8 kJ/mol`; the largest force-component discrepancy
was `6.80e-9 kJ/(mol*angstrom)`.

The historical phenol AM1-BCC **0.002 e** tolerance and approximately **-0.001 e**
residual remain unchanged. No charge neutralization or unsupported attribution
of that residual was made. Provided-charge examples are software tests, not
scientific validation. The older external phenol fixture retains its unknown
exact force-field/charge provenance.

Reproduce with an externally supplied path (no server path in public APIs):

```sh
python scripts/validate_minimization_references.py \
  --preparation-manifest "$ISLAND_AMBERTOOLS_REFERENCE_MANIFEST" \
  --output /new/path/minimization-acceptance.json
# Opt-in archive assertions:
pytest -q tests/test_minimization_archive.py
```

The new standalone `examples/minimize_geometry.py` evaluates synthetic Amber
bond records through real OpenMM: E 67.827126 → 0 kJ/mol, fmax 336.920582 → 0,
one accepted iteration and four evaluator calls. This is a software demonstration.
Tests additionally cover an anisotropic analytical harmonic minimum, gradient
sign/angstrom units, energy stagnation with residual force, initial convergence,
exact budgets, bounded line-search failure, ownership, stable ordering, model and
unit changes, failed evaluations, inverted/planar centers, missing dependencies,
and diagnostic application rules.

## Verification environment and limits

Verification uses a separate environment with Python 3.11.16, SciPy 1.17.1,
OpenMM 8.6.1 Reference, NumPy 2.4.6, ParmEd 4.3.1, and RDKit 2026.03.6. The base scientific
environments were not altered. The Reference platform uses its normal double
precision; GPU behavior is not an acceptance gate.

Final verification:

- Complete suite with both archive checks and live AmberTools regeneration
  enabled: **486 passed, zero skipped**. The live test regenerated all five
  AmberTools cases; the archive minimization test evaluated all five retained
  cases independently of regeneration.
- The earlier ordinary network-free run passed **481 tests, 3 opt-in skips**;
  two further real OpenMM regressions (singular start and ordered improper
  minimization) were added and passed in the final complete run. Default skips
  are live AmberTools generation, archived single-point acceptance, and archived
  minimization acceptance, each with an explicit enabling instruction.
- Base environment without SciPy: optional-import check **1 passed** and
  optimizer/OpenMM minimization modules **2 dependency skips**. A subprocess
  additionally blocked SciPy, OpenMM, RDKit, and ParmEd together and successfully
  imported core, evaluation, and minimization.
- Ruff and `pip check` passed. All **14 examples** passed, including actual
  short-polymer preparation and the new minimization example.
- The separate acceptance script reproduced all five report rows; it refuses to
  overwrite an existing report. Raw historical archives remain external and
  checksum-bound; the committed report contains numerical evidence and provenance.

Numerical agreement establishes implementation consistency only. A converged
local minimization is not a global minimum, equilibrium conformation, or thermal
equilibration. `production_validated=False` and
`simulation_readiness="not_established"` remain fixed. MD, thermostatting,
periodicity, packing, constraints, restraints, export, other force-field families,
MLIP, long-chain charge transfer, and expansion of AmberTools scope are excluded.
