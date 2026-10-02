# Phase 4F3.1 — consistent charge evidence and residual diagnosis

## Scope and repository

Correction on `codex/phase-4f3-oligomer-charge-references`, based on reviewed
`8f8de8ffdb1db6d63e3a856876ba686a0941580e`. At inspection, origin/main was
`796acf7e791b9fb6fdd08e5672ddc40bc4b5b968`; Phase 4F3 remained unmerged.
No historical reference, parameter, charge, tolerance or signature is changed.
`production_validated=False`, `simulation_readiness="not_established"`.

## Reproductions and correction

Eight focused regressions failed on the reviewed code before correction:

* Five independently rechecksummed preparation contradictions: `-c rc`,
  `-cf charges.txt` with AM1-BCC, duplicate `-c`, conflicting duplicate `-at`,
  and missing `-fi`. Saving accepted these despite declared AM1-BCC success.
* Matrix acceptance admitted actual tolerance 0.01 e/residual 0.004 e under
  the declared 0.002 e experiment, and incorrect reported site count or residual.

`validate_preparation_record()` now supplies the same data-only consistency
checks to `AmberToolsPreparationResult` and `ChargeReference`. The former still
validates its actual ImportedAmberResult first. The latter validates its stored
source, charge and imported-content signatures, then uses the shared validator.
Checks include schema-specific policies, canonical inner/outer payloads,
coordinates, generated-name lineage, source/artifact hashes, executable/data
provenance, stage results, required option/value pairs, duplicate/conflicting
options, and charge outcome. Malformed references raise `ChargeReferenceError`
before publication. Synthetic evidence fixtures now include full commands and
provenance and remain explicitly synthetic; they do not claim an executed QM run.

No schema or signature algorithm changed. Legitimate preparation v2/v3 and
reference v1 contracts remain readable. Previously accepted contradictory or
incomplete evidence is rejected. Offline loading still needs only core ISLAND
and NumPy, without scientific backends or original artifact files. Checksums and
consistency checks do not authenticate an arbitrary caller's computational claim.

The matrix reader validates the summary and row types/statuses, case membership,
signatures, actual method/size/tolerance and authoritative count/residual. It
rejects contradictory reported metrics before counting or saving that reference.
After comparison it derives outcome metrics from the reconstructed reference.
Exit 0 requires all eight cases; failed and absent cases remain visible.
Timeout (600 s/stage) and retry policy (zero) are **declared only**, not newly
claimed as verified by the historical signed preparation records.

## Read-only evidence and reproduction

Original inputs are in `island-validation/phase4f3-references/`; original accepted
records are in `island-validation/phase4f3-audit-final/`. New diagnostics and
re-audit outputs are in `phase4f3-1-diagnosis/` and `phase4f3-1-reaudit/` alongside
those directories. All 212 original case-file checksums in
[phase_4f3_acceptance.json](phase_4f3_acceptance.json) were verified before charge
parsing and again afterward. This includes input systems, SQM inputs/outputs,
lineage, AC files, MOL2, prmtop and logs, including both failed cases.

[phase_4f3_1_evidence.json](phase_4f3_1_evidence.json) records those exact checksums,
per-stable-site decimal charges and adjacent-stage changes, sums, observed decimal
precision, installed binary identities, and re-audit identity checks. It contains
our derived numerical evidence and manifests, not redistributed Amber source,
binaries or a replacement preparation. No QM or parameterization was rerun.

From the repository, with the validation Python environment:

```sh
python -m scripts.diagnose_oligomer_residuals \
  --source ../island-validation/phase4f3-references \
  --manifest docs/phase_4f3_acceptance.json \
  --output ../island-validation/phase4f3-1-diagnosis/stages-verified.json
python -m scripts.validate_oligomer_charges \
  --audit-existing ../island-validation/phase4f3-references \
  --output ../island-validation/phase4f3-1-reaudit
```

Outputs are exclusive: use a new output path to repeat. Re-audit intentionally
exits 1: six validated references, two original failures. All six reference
identities and raw charges are unchanged. Neither chemistry has a complete DP5
seed pair, so the conformation comparisons remain missing.

### Mapping and formats

The diagnostic resolves names through saved `input_name_to_site_id` against the
saved authoritative topology. It checks exact coverage and uniqueness. SQM's
unnamed output indices resolve through named `sqm.in` entries with checked atomic
numbers, then output element/index coverage. AC generated names and textual atom
types are checked; numeric BCC types are **not** treated as atomic numbers. AC
has no independent element column. MOL2 type elements and prmtop atomic numbers
are checked. This diagnostic intentionally supports the declared H/C/O cases.
It does not infer correspondence from charge values or spatial proximity.

`Decimal` sums preserve printed decimal tokens (prmtop division uses 28-digit
Decimal arithmetic). SQM atoms have three fractional decimals; AC and typed MOL2
have six. Prmtop uses E16.8 scaled charges, divided by **18.2223** to obtain e.
Tiny differences from the historical binary-float residuals are arithmetic
representation differences, not altered charges. An independent ParmEd 4.3.1
read of all eight prmtops agreed within 1.1102230246251565e-16 e per site
(predeclared 1e-15 e tolerance). The evidence includes input
MOL2 and pre-QM bond/atom-typing AC files: their charges are all zero and are not
AM1-BCC reference values.

## Stage-by-stage findings

All values below are molecular charge sums in e; formal charge is zero.
`PRE`, `BCC` and `MOL2` sums coincide for each row.

| Case | SQM atomic sum | PRE/BCC/MOL2 sum | prmtop sum (rounded here) | Historical outcome |
| --- | ---: | ---: | ---: | --- |
| PE DP3, 2026 | 0 | -0.000002 | -0.0000019999670733 | pass |
| PE DP5, 2026 | 0.002 | 0.001998 | 0.0019980000329267 | pass |
| PE DP7, 2026 | 0.001 | 0.000998 | 0.0009980000329267 | pass |
| PE DP5, 80317 | **0.004** | **0.004000** | **0.004** | fail |
| PEO DP3, 2026 | 0 | -0.000001 | -0.0000009977884241 | pass |
| PEO DP5, 2026 | **0.003** | **0.003001** | **0.0030010021786492** | fail |
| PEO DP7, 2026 | 0.001 | 0.001001 | 0.0010010021786492 | pass |
| PEO DP5, 80317 | -0.001 | -0.000999 | -0.0009989978213508 | pass |

SQM's separately printed **total** is `0.000` or `-0.000` in every case. It is
not the sum of the printed atomic tokens. For the two failed cases, the first
retained stage demonstrably containing the residual is the three-decimal SQM
atomic table, before BCC. The underlying unrounded per-atom values are absent;
we cannot reconstruct their exact sum from the rounded total either.

The observed sequence, supported by logs and installed implementation inspection,
is SQM -> parse printed charges -> equivalent-atom averaging -> write PRE.AC ->
AM1BCC -> read BCC.AC -> typed MOL2 -> tleap prmtop. The initial AC/AC0 files
precede SQM. The nested `atomtype -p bcc` log line belongs to am1bcc processing;
filenames alone were not used to establish ordering.

Averaging is directly visible. In failed PE, source parent 1 has hydrogen sites
11/12/13 with SQM values 0.071/0.072/0.073 and PRE values 0.072000 each. The
mirror-equivalent carbon names C001/C008 have -0.158/-0.159 in SQM and -0.158500
in PRE. In failed PEO, parent 1's H sites 16/17/18 have 0.090/0.079/0.076:
mean 0.081666666..., written as 0.081667 each. That group's six-decimal
serialization contributes +0.000001 e, accounting for the observed sum change.
These examples use verified bonded groups; raw per-site values are never replaced
in a reference by our calculations.

Across all eight, SQM->PRE maximum per-site changes range 0.0015–0.0275 e,
with total changes -0.000002, -0.000001, zero or +0.000001 e. BCC changes individual
sites by up to 0.1179 e (PE) or 0.2728 e (PEO), but preserves the printed molecular
sum in every retained case. BCC->MOL2 changes are zero at every site. The maximum
MOL2->prmtop per-site change is 2.195112582e-9 e; total changes are at most
2.211575926e-9 e. Those final serialization effects cannot explain the failed
0.003–0.004 e residuals. No post-BCC neutralizing adjustment is observed.

### Installed implementation, not an upstream assumption

The installed conda-forge package is AmberTools **24.8**,
`cuda_None_nompi_py310h834fefc_101`, package SHA256
`737d262831945c85c591d46d9c9a83acda0d049de2f985886f536aefdcf46710`.
The `antechamber` launcher is a shell wrapper; its historical wrapper checksum
alone would not identify the implementation. We additionally hashed the actual
ELFs and matched both to the cached package `info/paths.json`:

* `bin/wrapped_progs/antechamber`:
  `c679f0b7435b80e4d5f874e3fa8f6dba2800f7bb8fe782754871cd65cec55b4d`
* `bin/sqm`:
  `c82c78d22801d2899c0a8d622d5fe7f421eee4e541c67f6b6ca354008c56eb57`

With `AMBERHOME` pointing to the retained `md_rg` environment, exact inspection
commands were:

```sh
objdump -d --disassemble=rsqmcharge "$AMBERHOME/bin/wrapped_progs/antechamber"
objdump -d --disassemble=bccharge "$AMBERHOME/bin/wrapped_progs/antechamber"
objdump -d --disassemble=bcc "$AMBERHOME/bin/wrapped_progs/antechamber"
objdump -d --disassemble=qm2_print_charges_ "$AMBERHOME/bin/sqm"
objdump -s --start-address=0x31ba1 --stop-address=0x31bc0 "$AMBERHOME/bin/wrapped_progs/antechamber"
objdump -s --start-address=0x3242f --stop-address=0x32490 "$AMBERHOME/bin/wrapped_progs/antechamber"
objdump -s --start-address=0xd16d8 --stop-address=0xd17c0 "$AMBERHOME/bin/sqm"
```

`qm2_print_charges_` selects `F14.3` for atoms (format at 0xd16d8), accumulates
unformatted doubles (`addsd` at 0x5d447), and prints the total with `F12.3`
(format at 0xd1790). `rsqmcharge` parses `%d%s%lf` (0x31ba1) and stores the parsed
double, explaining why lost printed precision becomes input to BCC. `bccharge`
calls `identify_equatom`, sums/divides equivalent charges (0x24abd/0x24ae2),
then writes PRE.AC (0x2475c), runs am1bcc and reads BCC.AC (0x24964).
Disassembly hashes and format dumps are included in the evidence manifest.
This explains the observed data path without pretending the rounded zero total
proves an exactly neutral unrounded charge vector.

The cached recipe names AmberTools24_rc5.tar.bz2, SHA256
`52fb4fb3370a89b7ce738a2dc3e513c2fc1943fde4b4381846d9e75cc48d840f`, followed by
updates through AmberTools.8. An attempted source download was stopped because
of slow transfer; its incomplete bytes were not used as evidence. The exact
patched source tree is not available. No current upstream checkout was used as
proof of the installed implementation. These newly measured binary hashes are
not retroactively inserted into historical signed records; their agreement with
the installed/cached package and historical output supports this diagnosis, not
a new historical authenticity claim.

## Next experiment and limits

No supported precision-preserving charge-export option was established for this
installed path. Re-reading the same `sqm.out` cannot restore lost digits. Do not
raise the historical 0.002 e tolerance or convert failed references into passes.

A concrete next experiment, for a future phase with separate versioned records,
is an explicit **whole-oligomer charge-conservation projection hypothesis**:
`q'_i = q_i + (Q_formal - sum(q))/N`, on a copied charge vector with original
values retained. This is a proposed minimum-unweighted-squared-change correction,
not recovered QM precision. Apply only in a separately named diagnostic matrix,
never silently in preparation or long-chain transfer. Predeclare all eight cases;
measure maximum/RMS site changes, equivalent-group equality, length/conformation
audit changes and naive-transfer drift. Evaluate against independently obtained
higher-precision charge evidence when available before considering scientific
use. Neutrality alone would not validate that hypothesis, and none is implemented
here. The missing DP5 conformation pairs and scientific transferability gate
remain unresolved. No long-chain MD or charge transfer was performed.

## Verification

The complete ordinary suite: **908 passed, 9 expected opt-in skips**. Focused
new matrix/command/parser tests include rechecksummed contradictions, exact decimal
sums, element/coverage failures, partial/complete CLI gates and failed publication.
Existing offline subprocess isolation tests still block RDKit, ParmEd, OpenMM and
SciPy while loading/auditing records. Ruff and `python -m pip check` passed.
Reference inspection, local-template builder and cross-process checkpoint examples
ran successfully. Python 3.11.16, NumPy 2.4.6, RDKit 2026.03.6, ParmEd 4.3.1,
OpenMM 8.6.1, SciPy 1.17.1 were used. Opt-in live simulations were not rerun:
this correction executed retained-artifact reconstruction/diagnosis separately,
not a new eight-case QM acceptance. Historical charge tolerances and outcomes
remain unchanged.

The nine ordinary-suite skips were explicit opt-in gates: live AmberTools QM;
archived analysis; and archived checkpoint, NVE, session, Langevin, minimization,
single-point and workflow acceptance. Their environment switches were unset;
this does not imply those older archives are unavailable. The relevant Phase
4F3 retained charge artifacts were available and all eight were inspected.
