# Phase 4F4.1 — final geometry and historical identity validation

## Base and scope

Correction on `codex/phase-4f4-charge-conservation-audit`, based on reviewed
`18bbcc7c9d7dedd1482639e0ddb9efbfbcd8bb82`. Fetch confirmed origin/main at
`796acf7e791b9fb6fdd08e5672ddc40bc4b5b968`; Phase 4F4 and its Phase 4F3/4F3.1
prerequisite remain unmerged. Main and historical evidence were not changed.

## Reproductions and correction

Seven focused regressions failed before the correction: initial-only geometry,
missing final marker, concatenated complete runs, and each of the three forbidden
successful historical identities on a failed observation, individually and together.
The tests recompute artifact/manifest checksums and affected derived fields.
These failures concern contradictory evidence, not stale checksums.

Final geometry previously used `rsplit("Final Structure", 1)[-1]`, which scanned
all output when the marker was absent. Validation now requires exactly one
standalone `Final Structure` marker and exactly one completion marker, in that
order. Only the bounded Cartesian section between them is parsed. The supported
AmberTools Cartesian headers are recognized explicitly; headerless v1 synthetic
fixtures remain supported. Every row must have the expected index pair, element,
unique site and finite coordinates, with exact atom coverage. Unknown intervening
content, incomplete tables, repeated markers, concatenated runs and content after
terminal completion are rejected. Completion-line dash decoration is allowed.
Initial tables before the final section are not used. Existing fatal/convergence
checks remain mandatory. Valid coordinate values and fingerprints are unchanged.

Historical identity rules are status dependent:

- Failed import: `reference_identity`, `record_signature`, and `import_signature`
  must be absent or JSON null. A successful reference object or preparation file
  is also forbidden. The actual historical failure reason remains required.
- Passed reference: failure must be absent or null; the original reference,
  preparation and imported-result signatures continue to be validated.
- Observation and projection identities describe derived evidence and remain
  legitimate for failed historical imports.

Both creation and offline loading use the same validator. Nested projection/audit
validation and exclusive publication consequently reject invalid observations
with `ChargeReferenceError`. There is no schema change, signature migration,
charge-policy change or numerical-tolerance change. Checksums and semantic checks
establish internal consistency, not computational authenticity.

## Verification

Focused tests additionally cover truncated final tables after complete initial
ones, duplicate/inconsistent indices, wrong elements, nonfinite values, a distinct
valid final table, null failed identities, and rechecksummed nested records rejected
before publication. Existing dependency-isolation tests continue to load and
compare evidence without RDKit, ParmEd, OpenMM or SciPy.

Retained-artifact validation was executed read-only using:

- `island-validation/phase4f3-references`, resolved by exact paths and hashes in
  `docs/phase_4f3_acceptance.json` (212 original source files).
- `island-validation/phase4f4-final`, checked against all 19 output hashes in
  `docs/phase_4f4_acceptance.json` before and after validation.

All eight observations, eight projections and the audit validate. Observation and
projection identities match the committed acceptance manifest exactly. Raw and
projected charges, final-SQM fingerprints and all audit geometry comparisons are
unchanged. Six passed historical references retain their original identities;
PE DP5 seed80317 and PEO DP5 seed2026 retain failed status and null successful
historical identities. No retained case lacks final-geometry evidence. No files
were rewritten, no QM was rerun, and no numerical evidence was replaced.

The projection and reference inspection examples ran on the retained PE DP3
seed2026 records. Tests use the existing phase4e2 validation environment.

Final verification: **957 passed, 9 opt-in skips** (75.74 seconds); the focused
observation/conservation tests passed **49 tests**. Ruff and `python -m pip check`
passed. Python 3.11.16 and NumPy 2.4.6 were used. The ordinary skips are the
pre-existing live AmberTools and archived single-point, minimization, NVE,
session, Langevin, checkpoint, workflow and analysis gates; those unrelated
execution gates were not requested or rerun. The retained charge-evidence checks
above were executed separately, not counted as synthetic test evidence.

## Unchanged limitations

The original six-pass/two-fail raw experiment and incomplete raw conformation gate
remain unchanged. Conservation projection is experimental and does not recover
unrounded QM charges or validate a transferable charge model. No new calculation,
parameterization or dynamics was performed. All records retain
`production_validated=False` and `simulation_readiness="not_established"`.
