# Phase 4F3 — whole-oligomer charge references and transferability audit

## Scope and repository base

Branch `codex/phase-4f3-oligomer-charge-references` starts from `796acf7`
(`origin/main`). Main contains the squash-merged Phase 4F1/F2 changes; its tree
matches reviewed `c2d953b`. No merge, main modification, or historical charge
replacement was performed. Phase 4F2's synthetic charges and workflow remain
synthetic and unchanged.

This phase provides **whole-oligomer AM1-BCC reference charges**, using GAFF2
through the existing AmberTools engine. It does not provide RESP, a long-chain
charge assignment model, transferable charge templates, or new dynamics.
The initial chemistry contract deliberately recognizes only the literal mapped
repeat definitions `[*:1]CC[*:2]` and `[*:1]CCO[*:2]`, with builder-generated
hydrogen terminations, neutral atomic formal charges, standard masses, achiral
nonperiodic atomistic chains, and odd DP >=3. All AM1-BCC references remain
within 100 explicit sites. Copolymers, other end groups, isotope/stereo cases,
reversed or alternative repeat definitions, and ambiguous provenance are rejected.

## APIs and data contracts

`island.charge_references` exposes:

- `repeat_correspondence(system)`: verified repeat-local groups, or structured
  `compatible=False` diagnostics. No guess or coordinate matching is attempted.
- `create_charge_reference(system, preparation)`: first validates the actual
  `AmberToolsPreparationResult`, then creates an owned `ChargeReference`.
- `load_charge_reference(path)` and `save_record(record, path)`: strict loading
  and exclusive, durable publication of validated references/audits.
- `audit_charge_references(references)`: returns an owned `ChargeAudit`.
- Both record types expose `validate_integrity()` and a fresh owned `payload`.
  References additionally expose their content `identity`.

Schemas are `island_oligomer_charge_reference_v1` and
`island_oligomer_charge_audit_v1`. Their strict JSON envelopes contain a payload
and SHA256; existing tagged workflow serialization preserves integer keys and
tuples. Duplicate keys, NaN/Infinity, unsupported schema, stale signatures,
contradictory source/charge content, failed preparation stages and readiness
promotion are rejected with `ChargeReferenceError`. Publication validates before
writing, refuses overwrite, and uses the existing durable publication utility.

Each reference stores the complete authoritative source graph/metadata and
coordinates, seeds, verified correspondence, unmodified imported per-site charge
record, original preparation record/signature, and the original imported-result
signature payload. The signed imported payload links graph, mapping, charge
signature, preparation provenance and artifact checksums. Source coordinates,
charge coverage/residuals, calculation outcome, construction seeds and local
correspondence are cross-checked. Imported charges use the existing static-site
charge-assignment representation; the signed preparation identifies the actual
AM1-BCC calculation, rather than inferring its method from the assignment
container's name.

Loading/comparing saved records requires no RDKit, ParmEd, OpenMM or SciPy.
Creation from a preparation uses its established integrity checks. Reconstructing
from retained prmtop artifacts in the acceptance command additionally uses
ParmEd and existing import validation. Generation needs the existing construction
and AmberTools dependencies. No comparison invokes QM or creates an OpenMM Context.

These are internal consistency checks, not cryptographic proof that a calculation
was performed. The creation path validates original engine results; retained
artifact reconstruction checks their original hashes and signatures. Data-only
loading verifies the stored evidence without rerunning the physical calculation.
The full imported signature payload is retained as evidence, not as an additional
parameter parser or force-construction implementation. Historical signatures
are never regenerated to accept changed content.

## Correspondence and audit conventions

Each of the two literal repeat definitions has explicit source atom positions.
Heavy-atom correspondence uses `(repeat_unit_index, source_repeat_atom_index)`,
verified against the definition's elements, zero formal charges, oriented
inter-repeat connectivity, within-repeat bonds, and terminal valences. Builder
head/tail IDs must agree. The entire graph must match the expected hydrogen-capped
linear chain; site coverage is exact. Stable-ID values and insertion order do not
select correspondence.

Generated hydrogens have no unique cross-oligomer atom correspondence. They are
validated as singly bonded members of their verified parent's group, with the
correct repeat provenance and valence count. Stored summaries preserve every raw
site charge and report group count, sum, mean and **spread = maximum minus minimum**.
No raw value is overwritten by a group mean. Interior ether oxygens have zero
hydrogens: sum=0, mean/spread=null. Terminal hydroxyl hydrogens remain distinct
from these empty groups.

Every repeat has a head/interior/tail role; the central repeat is identified
separately. The audit reports every repeat's net charge and molecular residual.
Pairwise length comparisons require identical construction seeds and chemistry;
DP=5 conformation comparisons require distinct seeds. Head, central and tail
heavy atoms and parent hydrogen groups are compared only when their verified
local environments and group counts agree. Differences are right minus left;
maximum absolute and RMS differences are over explicitly listed matched groups,
not over arbitrarily matched hydrogens. Missing reference lengths or seed pairs
produce diagnostics. Length and conformational effects are not perfectly isolated.

For **each identified reference**, the diagnostic extrapolation is:

`Q_predicted(N) = Q_head + (N-2)*Q_central + Q_tail`, for N=20,50,100.

This is an audit of naive repeat transfer. It is not a target atom charge set,
neutralization policy or scientific acceptance threshold. No extrapolated value
can be sent directly to a parameterization entry point by this API.

## Commands and declared experiment

```bash
# A new, fixed matrix; no retries or method fallback.
python scripts/validate_oligomer_charges.py \
  --amberhome /path/to/ambertools --output /path/to/new-references

# Validate retained signatures/artifacts and publish a new audit; no recalculation.
python scripts/validate_oligomer_charges.py \
  --audit-existing /path/to/references --output /path/to/new-audit

# Offline inspection of one validated record.
python examples/inspect_charge_reference.py /path/to/new-audit/pe-dp3-seed2026.json
```

Generation publishes `declared-settings.json` before any calculation and updates
`matrix.json` after each attempted case. The fixed matrix is PE and PEO, each
DP=3,5,7 with template/assembly seed 2026 and a second DP=5 with both seeds 80317.
Per-stage timeout is 600 s; retries=0; GAFF2/AM1-BCC; maximum sites=100.
The explicitly predeclared molecular charge tolerance is **0.002 e** (not the
backend default). Raw residuals are recorded; the tolerance was not increased
following failures. Local-template internal construction attempts remain the
builder's recorded behavior; no outer seed search is performed.

The read-only command creates a separate audit directory, reconstructs only
successful preparations, and checks source artifact hashes before import. It
writes reference files, a self-contained `audit.json` and `outcomes.json`.
The latter records per-case file hashes, failures and validated coverage. CLI exit
0 requires all eight declared references validated; incomplete coverage exits 1
while preserving a valid partial audit. Command-level failures produce
`failure.json`. Existing output directories are refused.

## Actual executed matrix

All eight calculations ran using the installed AmberTools **24.8** conda-forge
build `cuda_None_nompi_py310h834fefc_101`. GAFF2 data header is version **2.2.20,
March 2021**, SHA256
`14ad62c8e532c47e2e400e2ca6ad8052b33bda4c4b9bc6f49b6a528e527512be`.
Executable hashes, package metadata hash, source file hashes and original
preparation/import signatures are in `phase_4f3_acceptance.json` and the full
local references. Help output does not give independent executable version
banners; that limitation remains explicitly recorded.

| Chemistry | DP | Seeds | Explicit sites | Validated reference | Original molecular residual, e |
| --- | ---: | ---: | ---: | --- | ---: |
| PE | 3 | 2026 | 20 | yes | -0.000001999967 |
| PE | 5 | 2026 | 32 | yes | +0.001998000033 |
| PE | 7 | 2026 | 44 | yes | +0.000998000033 |
| PE | 5 | 80317 | 32 | **failed import** | +0.004000000000 (diagnostic only) |
| PEO | 3 | 2026 | 23 | yes | -0.000000997788 |
| PEO | 5 | 2026 | 37 | **failed import** | +0.003001002179 (diagnostic only) |
| PEO | 7 | 2026 | 51 | yes | +0.001001002179 |
| PEO | 5 | 80317 | 37 | yes | -0.000998997821 |

Both rejected cases have SQM completion in retained artifacts, but the existing
importer rejects their charge residuals above the declared 0.002 e tolerance.
Read-only diagnostic inspection found sums of +0.004000 e and +0.003001 e already
in their `ANTECHAMBER_AM1BCC.AC` and `typed.mol2` files. The prmtop retains those
residuals within serialization precision. This evidence places the discrepancy
before tleap; it does not establish its earlier computational cause. Failed
charges are not promoted to validated references or substituted into comparisons.

Six references passed. No validated DP=5 seed pair exists for either chemistry,
so **conformation-sensitivity acceptance remains unmet**. PEO's DP=5 baseline
is also unavailable. All eight bounded attempts are preserved; none were retried.
The aggregate acceptance command correctly exited **1** for this incomplete matrix.

Raw artifacts are retained under `island-validation/phase4f3-references/`, with
the original generation driver and declaration. Validated outputs are under
`island-validation/phase4f3-audit-final/`. The initial driver used exactly the
same existing builder/parameterization engine and declaration; the committed
command subsequently reconstructed and validated these artifacts without reruns.
Only computed summaries and checksum manifests are committed, not Amber force-field
files or tool artifacts whose redistribution permissions have not been established.

## Measured discrepancies and naive-transfer drift

Matched head/central/tail heavy-atom differences, in e:

| Chemistry / lengths | Maximum absolute | RMS |
| --- | ---: | ---: |
| PE 3 → 5 | 0.000500 | 0.000289 |
| PE 3 → 7 | 0.001000 | 0.000645 |
| PE 5 → 7 | 0.000500 | 0.000408 |
| PEO 3 → 7 | 0.035000 | 0.016180 |

For PEO 3 → 7, hydrogen-group **sum** differences have maximum 0.048 e and
RMS 0.023357 e. Hydrogen-group mean differences have maximum 0.024 e and
RMS 0.014237 e over the seven nonempty groups. The complete per-group values,
counts, raw hydrogen values and repeat sums are retained in the audit.

Hypothetical charge totals for naive transfer, in e:

| Reference providing head/central/tail | N=20 | N=50 | N=100 |
| --- | ---: | ---: | ---: |
| PE DP3, seed 2026 | -0.034002 | -0.094002 | -0.194002 |
| PE DP5, seed 2026 | -0.015002 | -0.045002 | -0.095002 |
| PE DP7, seed 2026 | +0.000998 | +0.000998 | +0.000998 |
| PEO DP3, seed 2026 | -0.323001 | -0.893001 | -1.843001 |
| PEO DP7, seed 2026 | -0.044999 | -0.134999 | -0.284999 |
| PEO DP5, seed 80317 | +0.131001 | +0.341001 | +0.691001 |

These totals expose sensitivity to the chosen reference and accumulation of a
small nonzero central-repeat charge. They are not predicted measured polymer
charges. Even the nearly neutral PE DP7 extrapolation is not evidence of universal
transferability; terminal/environment and finite-precision effects remain.

## Recommendation and scientific boundary

First investigate the upstream AM1-BCC charge residual using the retained
intermediate artifacts and a separately declared precision/charge-conservation
experiment. Do not accept the missing seed comparisons by raising this run's
tolerance or correcting its charges. A fresh declared validation should restore
complete DP5 conformation coverage before fitting a transfer rule.

PE's small observed differences make a restricted head/interior/tail hypothesis
worth a future test, not deployment. PEO's larger variation and reference-dependent
drift motivate testing a larger local environment and more than one conformation.
A future neutral-repeat model would need to scientifically assess explicit
constraints such as `Q_interior=0` and `Q_head+Q_tail=0`, with a justified rule for
how any adjustment is fitted and validated. This phase implements none of those
constraints, residual redistributions, or assignments.

RadonPy's published approach uses RESP/HF charges for an optimized repeat and
adds removed capping-H charges to their bonded atoms when constructing chains.
That is relevant background for an explicit future conservation treatment, but
ISLAND's whole-oligomer AM1-BCC audit does **not** reproduce that RESP protocol or
its capping-charge redistribution. See [Hayashi et al., Methods, 2022](https://www.nature.com/articles/s41524-022-00906-4).

`production_validated=False` and `simulation_readiness="not_established"` remain
mandatory. No long-chain MD, equilibration claim, new force field or scope increase
was introduced.

## Verification

Synthetic software regressions cover ordering/remapped IDs, hydrogen groups and
terminal counts, invalid provenance/orientation/isotopes/stereo, signed-content
tampering, failed outcomes, analytical sums/drift, ownership, strict/exclusive
publication and missing-dependency isolation. They explicitly label fabricated
charge evidence and do not count it as QM acceptance. The live matrix and six
reconstructed records are separate evidence.

The retained SQM inputs explicitly request semi-empirical minimization (AM1,
`grms_tol=0.0005`, `scfconv=1e-10`, `ndiis_attempts=700`). Thus seed comparisons
would describe the backend's full charge workflow, including its geometry
relaxation, not an isolated charge evaluation on unchanged input geometries.
The missing valid seed pairs prevent measuring even that sensitivity here.

The reference-inspection, local-template builder and cross-process checkpoint
examples executed successfully. Ruff and `python -m pip check` passed.
The ordinary suite retains nine opt-in skips for separate historical/live
acceptance commands. The eight new AM1-BCC runs described above actually executed
outside that ordinary suite; their two validation failures remain unmet acceptance,
not skips or successful references.

Final ordinary result: **885 passed, 9 expected opt-in skips** in 71.89 s.
This includes 28 new synthetic/contract tests. The six actual saved references
and their audit were revalidated after implementation; optional-dependency
isolation also passed in a fresh process with RDKit, ParmEd, OpenMM and SciPy
blocked. The historical Phase 4F1 source inventory still matches every checksum
recorded by Phase 4F2. No historical workflow or preparation was rewritten.
