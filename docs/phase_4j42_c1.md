# Phase 4J42-C1 — strict external PDMS blocked outcomes

C1 repairs the raw exception leak when an invalid external-PDMS graph is used
at the evaluator boundary. Preparation still returns its immutable diagnostic
record. Evaluator creation now raises `ExternalPDMSBlockedError` for all three
blocked statuses, preserving the full preparation result and diagnostics.
No PDMS parameters, numerical execution profile or MD capability is introduced.

The fetched base is the published J42 branch at
`2c48a3eaa69e111279869ff6da12640225f78421`. The C1 branch is
`codex/phase-4j42-c1-blocked-diagnostics`; the tested implementation is
`c2dc8e3c6440a1a1998bf9bece98f3d507fe0946`. Main is not modified or merged.

`production_validated=False`; `simulation_readiness="not_established"`.

## Reproduced defect

Before editing, four real final graphs reproduced eight raw `ValueError`
exceptions: one evaluator attempt and one evaluator-before-minimization attempt
per graph. `prepare_forcefield()` had already returned `typing_blocked`,
`requirements=None`, and a `final_graph_invalid` diagnostic in every case.

| Original input | Recorded preparation diagnostic | Old evaluator result |
|---|---|---|
| `[SiH3][SiH3].O` | `silicon_oxygen_bond_required` | raw `ValueError` |
| `[SiH2]O` | `wrong_valence_or_missing_H:1` | raw `ValueError` |
| `[SiH3][O-]` | `wrong_valence_or_missing_H:2` | raw `ValueError` |
| `P` | `outside_pdms_domain` | raw `ValueError` |

The old `create_evaluator()` called `validate_integrity(system)`, which always
recomputed `requirements(system)`. Repeating the rejected graph validation
raised the original raw error before `require_executable()` could classify it.
The archive retains the original reproducer, exact output and original commit.
The defect was independently reexecuted in a detached checkout of that commit.

## Correction and validation policy

For `typing_blocked`, integrity validation now requires `requirements=None`, the
single recorded v1 typing stage `final_graph_invalid`, and its nonempty context.
It requires empty interaction/dependency inventories and no bundle/import
evidence. It then returns without validating the rejected graph again. A
typing-blocked record is a diagnostic, not an assignment to a valid graph.

For `asset_unavailable` and `parameter_incomplete`, integrity validation still
requires a graph inventory and recomputes it when a system is supplied.
Mismatched graph, molecule or component identity remains rejected. These statuses
cannot contain the typing-failure stage. Re-signed mutations that promote a
typing result to `prepared`, `asset_unavailable` or `parameter_incomplete`, insert
a fake inventory, or alter the typing stage are rejected.

`create_evaluator()` delegates to `require_executable(system)` before any native
evaluator dispatch. This guard raises `ExternalPDMSBlockedError`, including when
a stale or tampered record fails integrity validation. The exception retains:

- `result`: the exact preparation object;
- `status`: the classified blocked status;
- `diagnostics`: a separate copy of the original diagnostics;
- `validation_diagnostic`: the integrity error, when present, recorded separately
  as `blocked_record_invalid` without replacing the original diagnostics.

Invalid status data introduced by bypassing the frozen constructor is classified
as `parameter_incomplete`. Invalid system attributes encountered during initial
graph preparation are also captured as `typing_blocked`. No numerical engine is
created to diagnose these failures.

Only the external adapter, its evaluator dispatch, and the relevant stale-graph
test expectation change. Native PCFF, the C5 bundle schema, historical receipts,
minimization, workflow, dynamics and export implementations are unchanged.
The prepared-bundle serializer still refuses a blocked result because it is not
a `PreparedForceField`. Tests install failure sentinels at numerical execution,
minimization, initialization, checkpoint and dynamics boundaries. The public
evaluator-before-minimization path raises the classified exception first; no
checkpoint, trajectory or prepared-bundle directory is published.

## Executed controls

| Case | Result | Evaluator and minimization/evaluator boundaries |
|---|---|---|
| Invalid Si valence | `typing_blocked` | classified |
| Missing Si–O bond | `typing_blocked` | classified |
| Formal charge on an otherwise four-coordinate Si graph | `typing_blocked` | classified |
| Radical metadata on an otherwise four-coordinate Si graph | `typing_blocked` | classified |
| Unsupported phosphorus graph | `typing_blocked` | classified |
| No supplied external bundle | `asset_unavailable` | classified |
| Missing independent bundle pin | `parameter_incomplete` | classified |
| Malformed bundle | `parameter_incomplete` | classified |
| Empty C5 template | `parameter_incomplete` | classified |

The empty template is a negative schema control with no coefficients or charge
vector. No positive external bundle was supplied or discovered. Graph-bound
preparation records, selections and all 18 boundary outcomes reconstruct
identically after relocation in a fresh Python process. Focused tests also cover
status/stage/inventory mutations, stale valid graphs, constructor bypass,
malformed envelopes and preservation of exception diagnostics.

The C1 software gate is zero for correct blocking. The scientific gate remains
one. Native PCFF is not completed. No source fallback, missing-to-zero assignment,
charge normalization, external execution authority, LAMMPS/OpenMM execution or
equilibration result is introduced. Polymer Ewald/PME remains unvalidated,
pressure remains unavailable, and cutoff recovery and long NPT remain blocked.

## Verification and historical evidence

The focused C1 suite passes **32 tests**. The pre-publication combined J42/C1
check passes **80 tests**. The supported full suite passes **3,425 tests**, with
**10 existing optional skips** and no failures or errors. Its exact J22–J42
regression subset passes **1,425 tests**; this is a subset of the full execution,
not a second run. Counts and command records are retained in
`docs/evidence/phase_4j42_c1/receipt.json`. Ruff, compileall, five pip checks,
diff checks and fresh-process relocation pass.

An independent Git-blob comparison preserves all **2,533 historical documentation
and evidence paths**, including symlink text and target bytes, plus **645 source
and test files** outside the three intentionally modified files. This includes
native PCFF, the C5 schema and unrelated workflow/minimization/dynamics code.

The unchanged J42 receipt binds its original implementation hashes, including
the defective adapter. Its validator is executed with `repo` pointing to an
unchanged detached J42 checkout for those source hashes. Study reconstruction
runs with the current repaired runtime and retains J42's exact records. This
preserves historical evidence without claiming the corrected source equals its
old hash. No historical receipt or validator is rewritten.

## Reproduction

```bash
git fetch origin
git rev-parse origin/codex/phase-4j42-external-pdms-execution
# Verified base: 2c48a3eaa69e111279869ff6da12640225f78421

export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PYTHON=../island-validation/phase4e2-env/bin/python

# Fresh output directory; exit 0 means correct software blocking, not MD readiness.
"$PYTHON" scripts/phase_4j42_c1.py --output /tmp/island-j42-c1-study
cp -a /tmp/island-j42-c1-study /tmp/island-j42-c1-relocated
"$PYTHON" scripts/phase_4j42_c1.py --reconstruct /tmp/island-j42-c1-relocated
"$PYTHON" -m pytest -q tests/test_phase_4j42_c1.py tests/test_phase_4j42.py

"$PYTHON" -m ruff check .
"$PYTHON" -m compileall -q src scripts tests
git diff --check
"$PYTHON" scripts/validate_phase_4j42_c1_receipt.py docs/evidence/phase_4j42_c1

# Validate historical source pins against their original checkout while using
# the current runtime for J42 study reconstruction.
git worktree add --detach /tmp/island-j42-original 2c48a3eaa69e111279869ff6da12640225f78421
"$PYTHON" -c "from scripts.validate_phase_4j42_receipt import validate; print(validate('docs/evidence/phase_4j42', repo='/tmp/island-j42-original'))"
```

The execution archive contains the nine-case study. The validation archive
contains `reproduction/reproduce.py`, `reproduction/before.json`, exact command
records and environments, JUnit outputs, relocation hashes, preservation
inventories and orchestration scripts. Reexecute the original reproducer from
the original checkout, using its interpreter and `PYTHONPATH=.:src`, to observe
the old raw errors. Full-suite commands retain all supported optional asset
paths and disjoint shard node IDs; regression results are the exact verified
subset of that full execution, not a second claimed run. The five pip checks
use the supported, Foyer, analysis, AmberTools and quantum environments recorded
there. No simulation parameters are needed for these blocked-boundary controls.
