# Phase 4F2 — bounded retained-input relaxation acceptance

This phase reuses the Phase 4F1 oxygen-containing `[*]CCO[*]`, DP=20,
142-explicit-site GAFF2 preparation. Its alternating nonzero charges are
**synthetic software-test inputs**, not a scientifically validated charge model.
No construction, AmberTools execution, charge generation, optimizer replacement,
or parameter reassignment occurs in the new acceptance command.

## Repository prerequisites and CLI correction

Phase 4E8/4E8.1 was squash merged into `origin/main` at `5be5960`.
Its tree exactly matches reviewed `1d216234`. Normal merge `e9c3347` brought
that main into Phase 4F1 without content changes. The 17-file difference
against main contains the Phase 4F1 additions. Correction `495ff7a` on Phase
4F1 fixes the old acceptance CLI's successful exit despite failed cases.
Phase 4F2 branches from that corrected tip and still depends on unmerged
Phase 4F1. Main was not modified.

The Phase 4F1 command now records two distinct aggregate gates:

- `parameterization_numerical`: every declared case passes parameterization,
  charge-preservation and independent numerical comparisons.
- `workflow_dynamics_analysis`: additionally the four-step workflow completes
  and its five retained samples have validated analysis/export.

Default exit status tests the first gate; `--require-workflow` tests the second.
Missing dependencies, missing cases and failed required comparisons fail the
requested gate. Later workflow failure does not erase parameterization success.
Injected outcome tests cover these exit semantics; they are not live scientific
integration tests.

## Reusable retained-input command

```bash
python scripts/validate_long_chain_relaxation.py \
  --source /path/to/island-validation/phase4f1-final \
  --output /path/to/new-phase4f2-experiment
```

The output must be new and separate from the source directory. Inputs are the
original `peo20/input-system.json`, `charges.json`, `preparation.json`, original
prmtop and all signed preparation artifacts. Reconstruction validates the
original preparation/import signatures, chemistry, coordinates, charges, mapping
and artifact checksums. Recorded historical paths remain historical; operational
reads resolve the checked artifact directory under the retained case directory.
The command compares the saved charge mapping and source label to signed content
and anchors preparation/import identities to the original Phase 4F1 report.

`report.json` records checksums of every retained source file and verifies they
remain unchanged, including the original failed 500-iteration workflow. New
experiments never overwrite that diagnostic. Content checksums are integrity
checks, not computational authenticity guarantees.

## Settings declared before execution

`declared-settings.json` is published before the first minimization call:

| Setting | Value |
| --- | --- |
| Platform | OpenMM Reference |
| Min iterations / evaluations | 5,000 / 10,000 |
| Force criterion | maximum atomic force-vector norm <= 0.1 kJ/(mol*angstrom) |
| Line search budget | 20 |
| Energy-change / increase tolerances | 1e-12 / 1e-8 (existing conventions) |
| Temperature / friction | 300 K / 5 per ps |
| Timestep | 0.1 fs |
| Steps / segments | 20 / two segments of 10 |
| Recording interval | 5 |
| Velocity / thermostat seeds | 78123 / 99181 |
| Per-segment calls / frames | 12 / 3 |
| Uninterrupted calls / frames | 22 / 5 |
| Continuation comparisons | absolute 1e-10, relative 1e-12 |

There is no automatic budget escalation, force clipping, potential modification,
or relaxation of charge, force, stereo or verification tolerances. The existing
workflow requires both force convergence and independent final verification
before initializing velocities or beginning propagation.

After the first ten-step segment, the directory is relocated and the existing
workflow CLI resumes it in a separate Python process. The uninterrupted
comparison uses exactly the saved minimized system and initialized velocities;
neither minimization nor initialization is repeated. Comparisons cover retained
coordinates, full-step velocities, energies, components and forces, and verify
complete final PCG64 state, trajectory origin, absolute step/time and draw counts.
Extra startup/final checks explain 24 split calls versus 22 uninterrupted calls.

Expected retained steps are 0, 5, 10, 15 and 20, with the earlier accepted shared
boundary retained. Analysis validates/export these samples without modifying the
workflow. Restart data and CSV/JSON analysis output remain separate.

## Failure and evidence policy

The command writes machine-readable outcomes and exits nonzero on unmet
acceptance. Minimization histories and independently verified failure diagnostics
remain in the workflow directory. If relaxation does not converge, no dynamics,
continuation or analysis acceptance is claimed. Script-level failures retain an
explanatory diagnostic; existing output directories are refused.

The reported setup wall time includes artifact reconstruction within workflow
setup, minimization, initialization and the first segment when reached. It is
not an isolated optimizer benchmark. Historical Phase 4F1 timing and its
500-iteration failure are preserved independently.

Production code, physical-model fingerprints, integration/checkpoint contracts
and public APIs are unchanged. This is one bounded isolated-chain software test;
even successful 20-step propagation does not establish equilibration, meaningful
polymer properties or production suitability. `production_validated=False` and
`simulation_readiness="not_established"` remain unchanged.

## Actual verification

The retained-input experiment actually executed on OpenMM 8.6.1 Reference,
Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, RDKit 2026.03.6 and ParmEd 4.3.1.
No AmberTools program was rerun.

| Quantity | Initial | Returned minimum |
| --- | ---: | ---: |
| Potential energy, kJ/mol | 297652.80804778053 | 418.08727432871683 |
| Maximum atomic force, kJ/(mol*angstrom) | 1910787.250569282 | 0.09673819782172746 |
| RMS atomic force, kJ/(mol*angstrom) | 222024.20268919342 | 0.0452036430578729 |

Termination was `force_converged`: **3,113 iterations / 3,193 evaluations**
within the declared 5,000 / 10,000 budgets. The independently evaluated final
state passed verification. The complete accepted-iteration history is retained
in the immutable minimization record. Setup plus first segment took 194.2555 s;
separate-process resume took 3.95384 s. These measurements were taken with other
verification work running and are not controlled performance benchmarks.

Twenty steps completed (0.002 ps), with 24 split evaluator calls versus 22
uninterrupted calls. Both runs consumed 8,520 Gaussian scalar draws across
20 steps. The full PCG64 state and original trajectory identity matched.
At each retained step **0, 5, 10, 15, 20**, maximum absolute discrepancies in
coordinates, synchronized velocities, forces, potential, kinetic and total
energy were **0.0**. Components also passed the predeclared 1e-10 absolute /
1e-12 relative comparisons. This exact agreement is observed on this environment,
not a cross-platform bitwise guarantee.

Five retained samples passed report integrity and JSON/CSV export. The analysis
left workflow files unchanged; every retained Phase 4F1 source file was unchanged,
including its original failed workflow. The increased declared work budget
resolved this case without an optimizer or scientific-model modification.

Results are retained outside the repository under
`island-validation/phase4f2-bounded/`: `declared-settings.json`, `report.json`,
`workflow-relocated/`, `uninterrupted.json`, `resume.stdout`, `resume.stderr`, and
`analysis/`. Operational artifact paths are configurable through `--source`.

### Original provenance anchors

- Preparation signature:
  `230177cc467b088eb61544806caa94affb6194d12eb978e4348035e8632445f5`.
- Imported result signature:
  `6d5952387b0820511336744cef4e004ab7b17d2e25a2e0d4df35302075a4eeb3`.
- Authoritative input JSON SHA256:
  `130adff5e634cc199926df04e127b9655c2dc463a86c82efc89d321d3b31dbde`.
- Supplied charge JSON SHA256:
  `73e47ebe33c7d9614b6719bd3365b5d3a00dae3af1ade05273669aab444169b8`.
- Preparation JSON SHA256:
  `9dd33435650a4776a32724b54656d1dcc2607d9e4f4ae2e85cf92ad47fc2d215`.
- Original prmtop SHA256:
  `39d00fcd07be4246f24a2362934f21388f337af81a7bd6d22b5dc9f7925beca5`.

All other artifact hashes, including the old failure records, are in the new
report's source inventory. Charges remain labelled “Synthetic alternating
sorted-site +/-0.01 e; software acceptance only; no scientific charge model.”

### Regression and compatibility checks

The merge-related subset passed 75 tests. Six injected CLI-gate tests and four
acceptance-tool guard tests passed. The latter cover missing retained inputs,
existing-output preservation, out-of-tolerance comparisons and rejection of
corrupt artifacts before reconstruction. These tests supplement the actual
142-site numerical execution; they do not replace it.

Four opt-in historical compatibility tests executed and passed: archived
single-point references, archived PE workflow, checkpoint continuation and
read-only analysis of retained live-origin/archived PE bundles. Historical
parameters and tolerances were unchanged. Checkpoint, Langevin, NVE, evaluation
session, minimization and conformation examples executed, as did the workflow
status CLI and analysis example on the new relocated 142-site workflow.


Final ordinary suite: **857 passed, 9 expected opt-in skips**, in 69.24 s.
Ruff and `python -m pip check` passed. Four of those opt-in checks were run
separately as described above (4 passed in 40.76 s). The standalone archived
NVE, session, Langevin and minimization acceptance commands were not additionally
rerun; their ordinary regressions ran, and the new retained workflow exercised
actual session-backed minimization and Langevin continuation. Live AmberTools/QM
acceptance was not run because this phase explicitly reuses retained parameters.
There were no missing retained inputs for the requested 142-site experiment.

Remaining limitations: only this declared above-100-site synthetic-charge
workflow was established here. Larger systems or other starting geometries may
still exhaust their budgets. The old 500-iteration failure remains valid evidence
of insufficient work for that experiment; it is not reclassified as success.
