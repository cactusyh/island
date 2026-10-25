# Phase 4J40-C14: backend-scoped readiness and cross-engine boundary

C14 adds a fail-closed backend readiness API and a reproducible comparison
report. It preserves C13's **300 native evaluations, 108 comparisons and 13
passing comparisons**, including every numerical residual, frame identity,
setting and inherited bound. It does not execute new molecular dynamics,
change charges, parameters or electrostatic semantics, or qualify a new
PME/Ewald numerical-equivalence policy.

`production_validated=False`; `simulation_readiness="not_established"`.
Scientific equilibrium and scientific NPT acceptance remain unestablished.

## Separate backend states

The API is `island.analysis.backend_readiness`. Its versioned record binds
checks to the exact final graph, parameter assignment, charge identity,
backend and protocol. It reports these states independently for each retained
polymer cell:

| field | retained C14 state | meaning |
|---|---|---|
| `openmm_periodic_backend_status` | `software_tests_only` | required OpenMM evidence is incomplete |
| `lammps_periodic_export_status` | `exact_frame_exported_only` | retained data decks executed for particular frames |
| `lammps_kspace_runtime_status` | `native_kspace_executed_only` | real KSPACE execution occurred; not physical equivalence |
| `cross_engine_equivalence_status` | `ewald_unvalidated` | actual-polymer component, force, virial, pressure and policy gates unresolved |
| `pressure_validation_status` | `pressure_unavailable` | independent full periodic reciprocal virial pressure validation missing |
| `scientific_equilibrium_status` | `not_established` | no new sufficient sampling or physical validation |

Export success never satisfies a runtime or equivalence gate. KSPACE runtime
success never implies physical agreement. OpenMM-only results cannot satisfy
cross-engine checks. Native LAMMPS output and the explicitly labelled C12
single-frame analytical special-pair composite remain distinct; the production
LAMMPS exporter is not repaired or promoted by C14.

## OpenMM permission boundary

`assess_periodic_backend_readiness(bindings=..., checks=...)` requires seven
OpenMM-specific checks before setting `openmm_periodic_workflow_permitted=True`:
exact graph/parameter identity, validated charge origin, safe cutoff/box
preflight, internal energy/force checks, independently validated backend
pressure controls, exact checkpoint/continuation and sufficient sampling.
Every passed check requires SHA-256 artifact references, a reason, canonical
check name, correct backend scope and identical system/protocol bindings.
Missing checks fail closed. Unknown checks, invalid statuses, empty successful
evidence, or mismatched bindings raise `AnalysisError`.

`require_openmm_periodic_workflow(record, expected_bindings=...)` validates
persistence and blocks an opt-in qualified workflow unless all seven OpenMM
requirements pass. LAMMPS failures do not substitute for OpenMM's own
requirements or automatically fail a fully evidenced OpenMM-only decision.
Conversely, no LAMMPS success can compensate for a missing OpenMM requirement.
Even a complete backend decision is scoped to the bound system/protocol; it
never establishes production validation, scientific NPT acceptance or
thermodynamic equilibrium.

This is an additive policy/guard API. Existing low-level software testing
entry points and historical force-field behavior are unchanged. The policy
validates scope, bindings, completeness and persisted decisions; callers must
independently verify artifacts before supplying successful witnesses. A hash
is an integrity reference, not proof of scientific validity. Synthetic passing
unit-test witnesses exercise the contract and are not publication evidence.
The C14 report builder accepts only the pinned C13 source bundle and cannot
promote these inputs through caller-supplied successful assertions.

## Actual retained systems and blocked stages

The report retains all tested point-charge, counterion, pairwise, C8H18
fragment, GAFF2 PE, OPLS PEO, aromatic RESP/Psi4 and network comparisons,
including strained/unstrained boxes and convergence sweeps. Four actual
polymer decisions are reported separately. Exact backend bytes and graph,
assignment and charge identities are bound to the report. Original charge
metadata and provenance strings are preserved verbatim, including historical
adapter-level `provided` metadata; C14 does not relabel them as native OPLS,
RESP or AM1-BCC or independently requalify quantum generation evidence.
The validated-charge-origin gate therefore remains `not_measured` in this
new exact-protocol qualification. This does not rewrite prior charge receipts.

The retained initial boxes are checked with the existing periodic preflight:
cutoff plus a declared 0.05 nm safety margin and 0.2 nm neighbor skin must be
strictly below half the shortest box length. The 0.2 nm skin corresponds to
the retained 2-Angstrom neighbor policy. This is an explicit conservative
assessment, not a packing transformation, cutoff reduction or density change.
No minimum safe NPT box bound is established; an initially safe frame alone
would still not qualify the long protocol. Per-cell box lengths, required
margin and exact diagnostic are in the report.

OpenMM internal force/energy qualification, independent reciprocal virial
pressure, exact continuation for the proposed qualified protocol and sufficient
sampling are not established by C13's single-point outputs. Historical software
continuation tests are not silently promoted into evidence for a new long NPT
protocol. Backend pressure remains unavailable even where finite-volume
cross-engine derivatives are close. Full independent tensor/force-virial
validation remains required. Cutoff recovery and scientific long NPT are
blocked. Software-level workflow tests remain permitted as software tests.

## Report, reproduction and integrity

`docs/evidence/phase_4j40_c14/comparison-report.json` includes:

- all 108 C13 comparison rows unchanged, with system sizes, boxes, energy and
  force residuals, percentiles, settings, bounds and unavailable pressure fields;
- exact engine versions/build output, precision, executable/source hashes and
  source references, preserving native/composite comparison semantics;
- each polymer's graph/parameter/charge/backend/protocol bindings, original
  charge metadata, explicit gates, blocked stages and preflight diagnostics;
- C13's unchanged candidate policy and convergence limitations, with no newly
  qualified equivalence policy.

From the repository root:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python -m scripts.phase_4j40_c14 --source docs/evidence/phase_4j40_c13 --output /tmp/j40-c14-comparison-report.json
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python -m scripts.validate_phase_4j40_c14_receipt docs/evidence/phase_4j40_c14
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_phase_4j40_c14.py
```

The report command writes the report then exits **1** for blocked scientific
acceptance. The receipt validator exits **0** for intact blocked evidence and
reports scientific acceptance exit code **1** separately. It reconstructs the
entire report from pinned C13 artifacts and verifies the manifest, current
implementation, backend decisions and historical bytes. Re-signing a mutated
report does not make it valid. The C13 receipt SHA-256 anchor is
`4d394a73291a535f9c03f6770376780f060edb49be6e132ec0198638efae06e2`.

Focused tests exercise each OpenMM requirement, backend separation, missing
and stale evidence, persistence mutations, and exact preservation of C13
residuals and counts. Required regression/full-suite results, Ruff, compileall,
five environment pip checks, diff checks and historical preservation are
recorded with exact commands, source commit, environment manifests and logs
in the receipt's validation archive. All historical evidence remains unchanged.
No experimental, production-equilibrium or universal force-field claim follows
from these software-level results.
