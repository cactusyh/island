# Phase 4F4 — explicit whole-oligomer conservation experiment

## Base and scope

Branch `codex/phase-4f4-charge-conservation-audit` starts at reviewed
`da8258ace238059e310b308ee16aa3fd50250ab8`. At fetch, origin/main remained
`796acf7e791b9fb6fdd08e5672ddc40bc4b5b968`; Phase 4F3/4F3.1 is an unmerged
prerequisite. Main was not modified or merged.

This is an **experimental transformation of post-BCC charges** on the eight
existing neutral, achiral, connected, explicit-H, hydrogen-terminated PE/PEO
oligomers. The existing oriented repeat definitions, correspondence restrictions,
standard masses, odd DP >=3 and 100-site reference ceiling apply. No default
AmberTools/charge-assignment engine changes, new QM, parameterization, dynamics,
seed search, or target-chain atom charges are introduced.

The original six-pass/two-fail experiment and its 0.002 e tolerance remain
unchanged. Failed imports provide raw observations, not successful historical
preparations or ChargeReferences. All products retain
`production_validated=False` and `simulation_readiness="not_established"`.

## Explicit policy and numerical convention

The caller must select **`uniform_molecular_l2_v1`**. For every explicit atom,
including H, with `N = number_of_sites`:

```
Q = sum of authoritative formal charges (zero in this scope)
delta = (Q - fsum(raw charges in ascending stable-ID order)) / N
projected[i] = raw[i] + delta
```

This uniquely minimizes `sum_i (projected[i]-raw[i])**2` subject to the molecular
charge constraint. It is an unweighted molecular L2 projection, not a fitted
charge model. The fixed **1e-12 e conservation tolerance** was declared before
processing the matrix. There is no final-atom residue repair, second normalization,
retry with different settings or relaxed tolerance. Metrics report the actual
floating-point corrections; tiny rounding differences from the nominal uniform
offset are possible. Equal input charges receive equal output charges. `math.fsum`
is used for molecular totals and RMS accumulation. Undefined/malformed, boolean,
nonfinite, missing or extra charge entries are rejected with `ChargeReferenceError`.
The numerical kernel accepts native integer/float charge scalars and integer IDs.

Whole-oligomer neutrality neither implies zero central-repeat charge nor neutral
head+tail charge. For DP >3, other interior repeats also contribute to neutrality;
using the central repeat for every interior repeat is an additional hypothesis.

## APIs and data-only contracts

From `island.charge_references`:

* `observe_retained_case(root, manifest_path, case, historical_reference=None)`
  creates a `RawChargeObservation` from exact, unique manifest paths.
* `project_molecular_charges(system, raw_charges, policy=...)` is the small
  numerical kernel. Its standalone mapping has no claim of calculation provenance.
* `project_charge_observation(observation, policy=...)` creates a validated,
  owned `ChargeProjection`, binding the transformation to artifact evidence.
* `audit_charge_projections(projections)` creates a `ConservationAudit` using
  the existing verified heavy-atom and parent-H correspondence/comparison kernels.
* `load_raw_observation`, `load_charge_projection`, `load_conservation_audit` and
  the existing exclusive `save_record` load/publish validated records.

New schemas are `island_post_bcc_observation_v1`,
`island_experimental_charge_projection_v1`, and
`island_conditional_conservation_audit_v1`. Existing schemas/signatures remain
unchanged. Each new record has strict, tagged JSON plus a content checksum,
`validate_integrity()`, and fresh owned `payload` access. Observation/projection
identities identify the recorded serialized evidence; they are not new physical
model fingerprints. Numerical projection is independent of insertion order and
uses noncontiguous stable IDs without renumbering.

Observation evidence embeds the original manifest, all its listed case artifacts,
and (for the six historical successes only) the original validated reference
content. Creation verifies every listed hash before scientific parsing. All
operational paths come from exact manifest entries; ambiguous basenames/backend
directories, unlisted directory substitutions, absent files and checksum changes
are rejected. No wildcard selects the first available prmtop or MOL2.

The loader validates source identity, repeat provenance, generated-name/index
lineage, elements, exact atom coverage, full single-bond connectivity, AM1/neutral
SQM inputs, GAFF2/BCC command evidence, SQM completion/convergence and LEaP zero
errors. Post-BCC AC and typed MOL2 charges must agree. Prmtop charges (Amber scale
18.2223), their totals and per-site differences from typed MOL2 remain separate.
Source coordinates are checked against the six-decimal input serialization.
The pre-QM AC file has three-decimal coordinates; across all eight retained cases,
`sqm.in` contains exactly those values, padded to four decimals. Input-to-AC
changes reach 0.0005 angstrom. Validation checks this actual intermediate
(the half-unit bound of three-decimal serialization), rather than assuming
direct input-MOL2-to-SQM four-decimal rounding. A synthetic 0.0003 angstrom
regression covers this distinction. Final SQM coordinates use its explicit
indices resolved through named `sqm.in`;
input, final SQM and typed-MOL2 coordinate identities are recorded separately.

For historical successes, creation can reconstruct through the existing validated
Amber importer using the **listed** prmtop, then verify the original reference,
preparation/import identities and charges. This may require ParmEd. It does not
run AmberTools or parameter assignment. A supplied already validated historical
reference avoids that reconstruction. Failed cases never take this route and
never acquire a fabricated preparation signature. Their completion evidence is
explicitly log-based; missing signed stage return codes are not invented.

After creation, offline loading/comparison is NumPy/core-only: no RDKit, ParmEd,
OpenMM, SciPy, executable, or source directory is needed. Records embed the
required evidence. Validation recomputes the projection, uniform offset, totals,
maximum/RMS corrections, repeat changes, comparison metrics and geometry evidence.
Canonical, type-sensitive comparisons reject rechecksummed contradictory derived
fields. Checksums and these consistency checks are **not computational
authentication**; a wholly fabricated coherent evidence set cannot be proved real
by its checksum. These charges are not unmodified AM1-BCC or recovered unrounded
SQM charges. No conversion to a production parameterized system is provided.

## Reproducible CLI and publication

The retained source is `../island-validation/phase4f3-references`. The final new
experiment is `../island-validation/phase4f4-final`. Use new output directories:

```sh
python -m scripts.audit_charge_conservation \
  --source ../island-validation/phase4f3-references \
  --policy uniform_molecular_l2_v1 --output /user/new-experiment --declare-only
# Inspect declaration.json; then execute precisely that declaration:
python -m scripts.audit_charge_conservation \
  --source ../island-validation/phase4f3-references \
  --policy uniform_molecular_l2_v1 --output /user/new-experiment --execute-declared
python examples/inspect_charge_projection.py \
  /user/new-experiment/projections/pe-dp5-seed80317.json
```

Without either execution flag the command exclusively creates the directory,
publishes the declaration, then executes. `--declare-only` exit 0 means only that
the declaration was published. Execution exit 0 requires all eight validated
projections and the eight declared comparisons, plus unchanged source hashes.
A failed case does not count. Incomplete execution writes truthful outcomes or a
global `failure.json`; no replacement charges are synthesized. Output is separate
from the source. A writer lock rejects concurrent execution. Existing records are
not overwritten; abandoned unpublished/partial runs require a new experiment
directory. Per-record publication uses the existing fsync/exclusive utility.

Outputs: `declaration.json`, `observations/<case>.json`,
`projections/<case>.json`, `audit.json`, `outcomes.json`. Records contain local
artifact text; keep those in the user directory and consider external-data
redistribution separately. The committed
[acceptance manifest](phase_4f4_acceptance.json) contains derived metrics,
source/declaration identities and output checksums, not the Amber binaries or
parameter artifact bundle. It includes every heavy/H-group and repeat comparison,
raw/projected values, correction contributions, and target DP20/50/100 scalar
extrapolations. These reports are not restart files or target charge assignments.

## Executed eight-case diagnostic

All **212 original artifact hashes** matched before and after execution. All six
historical reference identities and the old audit checksum were rechecked and
remain unchanged. Eight raw observations and eight projected cases validated;
the original raw acceptance remains six passes/two failures. No original file was
rewritten. The complete projected comparison is a **new conditional diagnostic**,
not completion of the original raw conformation gate.

| Case | Sites | Raw typed total (e) | Uniform offset per site (e) | Original status |
| --- | ---: | ---: | ---: | --- |
| PE DP3, 2026 | 20 | -0.000002 | +1.000000e-7 | pass |
| PE DP5, 2026 | 32 | +0.001998 | -6.243750e-5 | pass |
| PE DP7, 2026 | 44 | +0.000998 | -2.268182e-5 | pass |
| PE DP5, 80317 | 32 | +0.004000 | -1.250000e-4 | **fail** |
| PEO DP3, 2026 | 23 | -0.000001 | +4.347826e-8 | pass |
| PEO DP5, 2026 | 37 | +0.003001 | -8.110811e-5 | **fail** |
| PEO DP7, 2026 | 51 | +0.001001 | -1.962745e-5 | pass |
| PEO DP5, 80317 | 37 | -0.000999 | +2.700000e-5 | pass |

Largest projected molecular residual: **1.5959456e-16 e**, below predeclared
1e-12 e. Largest maximum/RMS per-site correction is approximately **0.000125 e**.
This is conservation by construction, not independent scientific evidence.
Each repeat changes by its explicit-atom count times the common offset, within
floating-point roundoff. H-group sums change by group count times the offset;
means by the offset; equal charges and group spreads are preserved to roundoff.

### Length and seed variation

Maximum absolute corresponding-heavy-site differences, across head, central and
tail groups (e; full RMS and H-group metrics are in the manifest):

| Chemistry / comparison | Raw | Projected |
| --- | ---: | ---: |
| PE DP3–5, seed2026 | 0.000500 | 0.0004374625 |
| PE DP3–7, seed2026 | 0.001000 | 0.0009772182 |
| PE DP5–7, seed2026 | 0.000500 | 0.0005397557 |
| PE DP5, two seeds | 0.000500 | 0.0005625625 |
| PEO DP3–5, seed2026 | 0.019000 | 0.0190811516 |
| PEO DP3–7, seed2026 | 0.035000 | 0.0349803291 |
| PEO DP5–7, seed2026 | 0.043000 | 0.0430614807 |
| PEO DP5, two seeds | 0.043000 | 0.0431081081 |

DP5 heavy-site RMS differences change from 0.0002886751 to 0.0003287880 e
(PE), and 0.0195306028 to 0.0195622431 e (PEO). DP5 maximum H-group mean
differences change from 0.000500 to 0.0004374375 e (PE) and 0.027000 to
0.0268918919 e (PEO). Uniform projection is small relative to PEO's observed
variation; it neither isolates nor eliminates conformation effects. Raw and
projected comparisons use exactly the same verified groups. The report explicitly
subtracts raw differences from projected differences to show the projection's
contribution, including head/central/tail net-charge differences.

Distinct seeds are not taken as evidence of distinct optimized conformations.
For the verified heavy-atom correspondence, RMS differences between all pairwise
internal distances (angstrom; invariant to rigid motion) are:

| DP5 seed pair | Input geometry | Final SQM geometry | Typed MOL2 geometry |
| --- | ---: | ---: | ---: |
| PE | 0.3464520 | 0.5257345 | 0.3465885 |
| PEO | 1.9481474 | 2.0049100 | 1.9480475 |

Thus the retained final SQM geometries do differ, well beyond their printed
precision. This is not a local-basin classification or an ensemble sample.
Typed MOL2 retains near-input geometry here; it must not be mislabeled as the
final SQM geometry. Hydrogen permutations are excluded from this geometry metric.
The raw charges still originate from the declared post-BCC typed MOL2 stage.

### Neutrality does not establish naive transfer

`Q(N) = Q_head + (N-2)*Q_central + Q_tail` is computed separately for each
raw/projected source, without any further correction. Examples at DP100:

| Source | Raw Q(100), e | Projected Q(100), e |
| --- | ---: | ---: |
| PE DP3, 2026 | -0.194002 | -0.1939418 |
| PE DP5, 2026 | -0.095002 | -0.1325894 |
| PE DP7, 2026 | +0.000998 | -0.0126565 |
| PE DP5, 80317 | +0.005000 | -0.0702500 |
| PEO DP3, 2026 | -1.843001 | -1.8429705 |
| PEO DP5, 2026 | +0.394001 | +0.3370631 |
| PEO DP7, 2026 | -0.284999 | -0.2987775 |
| PEO DP5, 80317 | +0.691001 | +0.7099550 |

Projected central-repeat charges remain nonzero: e.g. PE DP5/2026
-0.001374625 e, and PEO DP5/2026 +0.0034322432 e. A neutral diagnostic oligomer
therefore still gives substantial source-dependent extrapolation drift. No
universal scientific pass threshold was selected after observing these results.

## Next narrowly scoped experiment and limitations

PE's small variation makes a restricted head/interior/tail hypothesis worth a
**held-out oligomer** test, not production transfer. PEO's much larger length/seed
sensitivity argues for examining an environment extending beyond one repeat and
local torsional geometry first. A next bounded experiment should predeclare
held-out oligomers/conformations within 100 sites, compare single-repeat versus
neighbor-aware descriptions against separate calculations, and evaluate an
explicit conservation constraint at the model level. It must retain original and
projected targets separately and avoid fitting to both training and test cases.
Neither an atom template nor that conservation rule is implemented here.

Unrounded SQM charges remain unavailable. This experiment cannot establish the
correct physical distribution of the missing charge, rescue the original failed
raw reference gate, or establish scientifically transferable polymer charges.
Eight calculated geometries are not an equilibrium ensemble.

## Verification

Complete ordinary suite: **942 passed, 9 opt-in skips**. Ruff and
`python -m pip check` passed. Projection/reference inspection examples ran.
Tests include analytical projection, neutrality/idempotence, equal-charge
preservation, stable-ID/insertion-order independence, scope/coverage/value errors,
rechecksummed tampering, failed-status preservation, H/repeat changes, nonneutral
extrapolation from a neutral oligomer, exact artifact path failures, exclusive
publication, CLI partial gates and a subprocess blocking RDKit/ParmEd/OpenMM/SciPy.

The eight-case retained-artifact execution is separate from synthetic software
fixtures/injected failures. The nine skips are the existing opt-in live QM and
archived analysis/checkpoint/NVE/session/Langevin/minimization/single-point/workflow
gates; their environment switches were unset. No new QM/MD was requested or run.
The relevant eight retained oligomer calculations were all available and used.
Environment: Python 3.11.16, NumPy 2.4.6, ParmEd 4.3.1 for historical-success
reconstruction; the original tool/package identities remain those in Phase 4F3.1.
