# Phase 4F5 — frozen experimental PE charge template

## Repository and scope

Branch `codex/phase-4f5-pe-charge-transfer` starts at
`a7b57da05dce90f485c887f42f27f9fe47a48e33`. Fetch confirmed its tree equals
reviewed `220dcf4def05a98275dc849255507361dd89b7c2` exactly. The squash includes
Phase 4F3/4F3.1 and 4F4/4F4.1; old branches were not merged again.

Only literal `[*:1]CC[*:2]`, neutral, achiral, standard-mass, explicit-H,
hydrogen-terminated linear PE is supported. Provenance and the full chemical
path must agree. No PEO, isotopes, other end groups, copolymers, periodic systems,
implicit hydrogens or stereochemical assignments. Historical odd-DP reference
validation and the 100-site AM1-BCC limit are unchanged. The separate target
contract permits odd **and even DP >=3**, up to **1,000 explicit sites**.

This is an explicit experimental model of conservation-projected post-BCC
observations. It is not unmodified AM1-BCC, recovered unrounded QM, or a
scientifically validated production charge assignment. No production engine
selects it automatically. No new parameterization, minimization or MD backend
is introduced; larger-chain demonstrations are prediction only.

## Frozen fit

Exactly PE DP3/2026, DP5/2026, DP7/2026 and DP5/80317 are used. Embedded original
observations retain raw charges and historical statuses, including the failed
DP5/80317 import. All fit targets are `uniform_molecular_l2_v1` projected
charges. PEO and held-out observations cannot enter fitting.

Verified `(repeat index, source_repeat_atom_index)` and parent-H membership
identify each atom. The oriented CC path reverses as `(r,j) -> (DP-1-r,3-j)`.
Thus head1 equals tail2, head2 equals tail1, and both interior heavy sites are
equivalent; their attached-H classes follow the same verified correspondence.
This gives six coefficients in order:

1. terminal carbon, terminal-parent hydrogen (three H per parent);
2. near-end carbon, near-end-parent hydrogen (two H per parent);
3. interior carbon, interior-parent hydrogen (two H per parent).

Every atom in every interior repeat contributes. Individual reference H charges
are preserved; only predictions share a parent-class parameter. Each design row
has weight `sqrt(1/(4*N_case))`, giving exactly
`L = sum_cases(sum_atoms(error**2)/N_case)/4`. H multiplicities therefore contribute
to the objective, not just to reported net charges.

Model A is this symmetry-reduced unconstrained least-squares fit. Model B adds:

```
Q_head = terminal_C + 3*terminal_H + near_end_C + 2*near_end_H = 0
Q_interior = 2*interior_C + 4*interior_H = 0
```

Reversal already implies `Q_head = Q_tail`, so the first condition is equivalent
to the requested `Q_head + Q_tail = 0`. Both end repeats individually become
neutral. This is an explicit hypothesis, not a consequence of the references.
No target-dependent normalization or last-atom correction is used.

The deterministic solver uses an SVD constraint null space followed by
`numpy.linalg.lstsq(rcond=None)`. Dependent constraint rows or a deficient reduced
design are errors; no regularization/arbitrary coefficient choice. Diagnostics
include rank, singular values, objective (e²), and constraint residuals (e).
Numerical conservation is checked at **1e-12 e**, declared before execution.

## APIs and records

From `island.charge_references`:

- `fit_pe_template(projections, mode="baseline" | "conserving")`
- `pe_target_correspondence(system)`
- `predict_pe_charges(model, system)`
- `validate_pe_templates(models, targets, frozen_model_identities=...)`
- `load_pe_record(path)` and `save_record(record, path)`

Owned models are `PETemplateModel`, `PEChargePrediction` and
`PETemplateValidation`. Schemas are respectively
`island_experimental_pe_template_v1`, `island_experimental_pe_prediction_v1`,
and `island_experimental_pe_validation_v1`. Strict tagged-JSON checksum envelopes
preserve integer IDs and reject malformed or rechecksummed contradictory content.
Validation rechecks embedded original evidence, fit definitions/coefficients,
prediction coverage, and comparison arithmetic before publication. Source records
are never re-signed or changed. Data-only loading/comparison needs core ISLAND and
NumPy; actual generation uses the established optional AmberTools stack.

Predicted stable-ID mappings and chemical identities depend on verified topology,
masses and repeat provenance, not coordinates, seeds, insertion order or ID
magnitude. The complete prediction record also retains the caller's source system
as provenance, so its record identity can differ between coordinate frames even
when its chemical identity and assignments are identical. Source objects and
previously returned payloads are not mutated. Returned payloads are owned copies.

Model fit recomputation uses the recorded deterministic algorithm; numerical
records are generated and checked with the same NumPy environment. This does not
promise cross-build bitwise model identities. Content checksums and semantic
validation establish internal consistency, not computational authentication.

## Reproducible experiment

```
python -m scripts.validate_pe_charge_transfer --freeze \
  --training /path/to/phase4f4-final \
  --original-source /path/to/phase4f3-references \
  --output /path/to/NEW-phase4f5
python -m scripts.validate_pe_charge_transfer --execute \
  --output /path/to/NEW-phase4f5 --amberhome /path/to/ambertools
python examples/predict_pe_template.py /path/to/NEW-phase4f5/conserving.json --dp 20
```

Freeze verifies the committed acceptance hashes and exact original source files,
then exclusively publishes both models, their coefficients/identities, training
identities and the held-out plan. Execute accepts only that untouched frozen
directory, locks its writer and checks model hashes before and after evaluation.
A run cannot be silently retried in place. Failed cases have durable diagnostics;
no seed substitution, method fallback, tolerance change or refit is performed.

The declared held-out matrix is DP5/314159, DP9/314159, DP9/271828 and DP11/314159
(template and assembly seeds both equal the listed seed), with 32, 56, 56 and
68 expected explicit sites. Actual inventories are checked. GAFF2/AM1-BCC uses
`max_atoms=100`, `charge_tolerance=0.002 e`, 600 seconds **per stage**, zero outer
retries. These new artifacts are separate from all eight original cases.

Raw import acceptance, verified completed-QM observation eligibility, and projected
target availability are distinct gates. A raw charge-total failure can provide
a diagnostic observation only through the unchanged strict observation validator;
it never becomes a successful historical preparation. Exit 0 means all four
verified projected held-out comparisons and larger-target coverage checks ran,
not successful raw import for every case or a scientific accuracy threshold.

Reports compare baseline/conserving predictions to raw and projected targets
separately: heavy maximum error/RMSE, parent-H mean/sum errors, all group errors
with role and distance from nearest end, and repeat/molecular totals. There is
no individual-H correspondence across conformations. DP9 geometry uses mapped
heavy-atom pair distances for input and final SQM frames, as in the previous
audit; neither seed difference nor nonzero geometry distance proves distinct
optimized basins. No scientific pass threshold is selected.

## Executed evidence

Local output: `island-validation/phase4f5-frozen/`. The declaration was published
before executing any held-out calculation, with SHA256
`e39cbad6b06ad6b7e75be6962d346465b8cc04aafaf3ab46cde9aa6a712ead87`.
[phase_4f5_acceptance.json](phase_4f5_acceptance.json) records the exact frozen
models/training identities, all new source paths/checksums, stage outcomes,
comparison profiles and output hashes. External Amber files remain in the user
output directory; only derived numerical evidence/manifests are committed.

Frozen model identities:

- Baseline: `698867c028e8587d2ba7ca7cabac62b0e67de557e6d7b73556e6a08edf5b1642`
- Conserving: `87bf6d2c13c10de53b1473e9ba0f5c05dce62601c27111e9db8e323ff5db833c`

| Class | Baseline charge (e) | Conserving charge (e) |
| --- | ---: | ---: |
| terminal C | -0.0922311240975936 | -0.09238080252100883 |
| terminal H | 0.03252450195282658 | 0.03237482352941173 |
| near-end C | -0.07966347703877003 | -0.07981315546218466 |
| near-end H | 0.03768442212089381 | 0.037534743697479125 |
| interior C | -0.07882792030928214 | -0.07869341963322546 |
| interior H | 0.039212209140556044 | 0.03934670981661272 |

Baseline rank is 6; conserving constraint rank is 2 and reduced rank is 4.
Weighted objectives are respectively `7.569249401200903e-8` and
`9.582434315240151e-8` e². Conserving head/interior constraint residuals are
`-5.551115123125783e-17` and `-2.7755575615628914e-17` e. All interior repeats
were used. Neither model was changed after held-out results became available.

### Actual GAFF2/AM1-BCC calculations

AmberTools **24.8**, conda build `cuda_None_nompi_py310h834fefc_101`, was used
through the existing backend. GAFF2 data header is version **2.2.20 (March 2021)**,
SHA256 `14ad62c8e532c47e2e400e2ca6ad8052b33bda4c4b9bc6f49b6a528e527512be`.
The evidence manifest records executable, package, leaprc and data hashes from
successful preparation records. Executable help banners do not expose individual
version strings; those remain explicitly unavailable. No protocol modification,
extra QM execution or charge tolerance change was made.

| PE case | Sites | Raw import | Typed MOL2 total (e) | Build + preparation wall time (s) |
| --- | ---: | --- | ---: | ---: |
| DP5 / 314159 | 32 | failed charge check | +0.002002 | 8.554 |
| DP9 / 314159 | 56 | passed | +0.001002 | 49.655 |
| DP9 / 271828 | 56 | failed charge check | +0.002000 | 32.784 |
| DP11 / 314159 | 68 | passed | -0.001002 | 130.503 |

All four completed QM/observation validation and have projected targets. Only two
are successful raw imports. The DP9/271828 case sits on the declared boundary:
a read-only ParmEd 4.3.1 parse of the checksummed prmtop sums to
`0.0020000000000000157 e`, strictly above 0.002. The data-only decimal-conversion
observation sum is `0.001999999999999988 e`. We preserve the actual failed import
and its numerical convention; we do not relabel it as passing, alter the tolerance,
or claim that six-decimal display equality overrides the importer. DP5/314159's
ParmEd sum is `0.0020019999670733146 e`. Both failed observations have null
successful historical identities and retain the failure reason.

### Held-out comparisons

These are conditional errors against projected post-BCC observations, **not errors
against independent high-precision QM truth**. All errors below are in e. Complete
raw-target and projected-target results, group means/sums, end-distance profiles,
repeat totals and molecular totals are in the evidence manifest.

| Case | Baseline heavy max / RMSE | Conserving heavy max / RMSE | Conserving H-group mean max / RMSE |
| --- | --- | --- | --- |
| DP5 / 314159 | 0.00036536 / 0.00029851 | 0.00035059 / 0.00025707 | 0.00020927 / 0.00017140 |
| DP9 / 314159 | 0.00091003 / 0.00039921 | 0.00077553 / 0.00034326 | 0.00041460 / 0.00024663 |
| DP9 / 271828 | 0.00189221 / 0.00112003 | 0.00175771 / 0.00109337 | 0.00143242 / 0.00084594 |
| DP11 / 314159 | 0.00094266 / 0.00044410 | 0.00080815 / 0.00037742 | 0.00063197 / 0.00042421 |

Constrained heavy errors decrease for these four cases, but H-group errors do
not uniformly improve. For example DP9/271828 projected H-sum RMSE increases
from 0.00167991 to 0.00169915 e; constraints impose neutrality, not best accuracy
for every observable. The input heavy pair-distance RMS difference for the DP9
seed pair is **2.437416711275074 angstrom** and the final-SQM difference is
**2.925495278357244 angstrom**. Thus retained final geometries do differ under
this metric. Conformation effects and printed-charge precision are not isolated
by these few observations.

### Graph-only larger-chain prediction

| DP | Explicit sites / covered sites | Baseline total (e) | Conserving total (e) |
| --- | ---: | ---: | ---: |
| 20 | 122 / 122 | -0.012430575086314316 | -6.245004513516506e-16 |
| 50 | 302 / 302 | -0.036640696776517295 | -1.457167719820518e-15 |
| 100 | 602 / 602 | -0.07699089959352226 | -2.8449465006019636e-15 |

Repeated assignments were identical and complete. No whole-chain QM,
parameterization, relaxation or dynamics was performed for these targets. The
baseline drift is reported without normalization. The constrained totals meet
the predeclared 1e-12 e software tolerance; this is coverage/conservation evidence,
not independently validated larger-chain accuracy.

## Verification and recommendation

Ordinary suite: **978 passed, 9 existing opt-in skips** (98.84 s). The focused
PE tests passed **21 tests**, covering independent constrained algebra and
weighted means, rank deficiency, even/odd DP, reversal, noncontiguous IDs,
insertion/coordinate independence, H multiplicities, scope rejection,
rechecksummed tampering, ownership, frozen split, exclusive publication, injected
CLI gate failures and offline dependency isolation. Ruff and pip check passed.
The prediction and existing projection-inspection examples executed. Software
failure injection is separate from the four actual AmberTools calculations above.
Read-only post-run verification confirmed all 212 original source-file hashes and
19 Phase 4F4 output hashes unchanged, and checked all 104 new source-file hashes.
The nine skipped gates concern prior live/archived dynamics and workflow suites;
no new dynamics acceptance is claimed here.

The evidence supports a further **predeclared PE-only conformation validation**,
with a larger independent set of conformations and explicit end-distance error
stratification, while holding this candidate fixed. The larger DP9 seed-pair
sensitivity warrants testing whether end/neighbor environments need refinement
on a separate training set, followed by a new untouched test set. We have not
implemented such refinement. Precision-preserving charge references or an
independently justified charge method are still needed to evaluate accuracy
beyond the rounded/projected targets.

Extending to PEO requires new oriented oxygen/environment classes and independent
training/validation; PE reversal/terminal neutrality assumptions cannot simply
be transferred. Automatic force-field integration additionally needs an explicit
scientifically reviewed charge policy and broader independent physical validation.
No scientific threshold, equilibration claim or production suitability is inferred
from neutrality or these short oligomers. All products keep
`production_validated=False` and `simulation_readiness="not_established"`.
