# Phase 4J31-C1 — strict network acceptance and failure classification

Actual base: `6606895d643806c731b52cf0fafa00e338ce3adb`, the reviewed J31
commit, confirmed after fetching origin and with `git ls-remote`.
Branch: `codex/phase-4j31-c1-strict-network-acceptance`. No main merge or
alternative base is used. This correction targets the established remote Linux
environments; the existing Linux atomic-publication implementation is unchanged.

## Reproduced defects and correction

Before changing production code, a fresh subprocess monkeypatched only
`prepare_network_backend` to raise an explicitly injected GAFF `OSError`. J31
wrote `status="blocked"` and exited zero. A separate temporary copy of the real
receipt changed all eight cases to blocked, supplied unrelated failure reasons,
recomputed its identity and matching manifest section; `validate_structure`
accepted it. Neither reproduction modified the retained receipt. The two new
regression tests were also run against exact original scripts extracted with
`git show` into a temporary import overlay: both failed on J31. They pass on C1.
Injected failures, mutated receipts and dispatcher fixtures are **software
negative controls, never scientific evidence**.

The shared `scripts/phase_4j31_contract.py` freezes eight cases. OPLS-AA, GAFF and
GAFF2 must each complete both branched and crosslinked fixtures. Each successful
case must contain the existing complete interaction inventory, source and
assignment bindings, independent verification, five retained frames, exact
split/continuous frame and checkpoint/RNG comparison, relocation, separate-process
reconstruction, and all required network observables. Blocked cases cannot
satisfy these six successes. Case keys must agree with row family/fixture labels,
and final-graph, J20 plan and J21 packing identities are frozen to reviewed J31.
Unexpected PCFF success also fails this contract.

## Exact PCFF negative controls

The pinned source is SHA-256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
The original default profile is `island_pcff_acyclic_cho_v1`. Its native read-only
typing result is incomplete for each exact final chemical graph:

- Branched: site 3 reports `multiple_oxygen_carbon`.
- Crosslinked: site 3 reports `multiple_oxygen_carbon`; the component containing
  sites 1–28 additionally reports `cyclic_component`.

The documented rejection is exactly `PCFFError: Incomplete automatic typing
cannot assign charges`, raised by the typing-completeness guard in
`island.forcefields.pcff.automatic.bridge_data`, before charge assignment and
parameterization can continue. The driver independently verifies the pinned
source, exact final graph and native diagnostic payload, then requires the exact
exception class/message and native traceback origin at this guard. It does not
recognize an unrelated exception containing the word “typing”. Missing/corrupt
sources, missing dependencies/executables, permissions, other stages and other
diagnostics fail acceptance. No PCFF dynamics run. Expanded PCFF support requires
a separately reviewed change to the frozen negative-control expectation.

## Diagnostics and exit codes

New diagnostic classifications distinguish `expected_typing_rejection`,
`asset_unavailable`, `unexpected_preparation_failure`, `execution_failure`, and
`completed`. The existing `complete` and `blocked` status meanings remain;
`failed` records unexpected preparation or numerical/workflow failures. Partial
success evidence is separated from success claims. An execution failure can
record `dynamics_executed=null` when the failure does not establish how much
execution occurred; this can never meet acceptance.

A single-case CLI exits zero only when its own declared outcome is met. Every
other outcome preserves `result.json`, available subprocess output and a
traceback in `failure.log`, then exits nonzero. Matrix runs continue through
independent failed cases where safe, append available logs, persist results after
each case, and preserve a partial matrix even on interruption. Completion of the
loop is insufficient: all eight outcomes and all child exit codes must pass.
There is no inventory-success mode and no production injection/bypass option.
Tests patch callables only inside controlled fresh interpreters.

## Historical compatibility and additive evidence

The original J31 receipt, manifest, all J17–J31 evidence and serialized identities
remain byte-for-byte unchanged. The full receipt entry point enforces the shared
contract, verifies each retained row against its archived `result.json`, and
rechecks workflow/backend/source/assignment and analysis bindings. Original PCFF
rows lack additive stage fields: compatibility is restricted to their exact
frozen row digests and a required independent replay of source verification and
native typing. No generic legacy blocked row is accepted.

J31 implementation-file hashes continue to be checked against its actual
validated implementation commit, `ba0e2abd951278958a852964da790126f49bbf65`.
Four acceptance-only files necessarily differ in C1: the two corrected validators
and their two existing test files. Their current hashes, along with the new
contract, tests, documentation and validation scripts, are bound by the additive
C1 receipt to the corrective implementation commit actually tested. All other
J31 implementation hashes are still checked against current bytes. A later
receipt/evidence commit is distinct from the validated implementation commit.

The additive record is `docs/evidence/phase_4j31_c1/receipt.json`; it does not alter
the historical reproducibility manifest. It archives exact commands, environment
versions, actual JUnit counts, check logs, before/after controls and the read-only
revalidation report. The C1 validator checks source against the recorded commit,
artifact hashes, historical bytes and the corrected full J31 receipt validator.

Six retained analysis bundles are copied only into temporary relocation paths,
loaded with their expected identities, compared with archived observables and
workflows, and reconstructed in fresh subprocesses. Originals are only read.
The revalidation script runs no dynamics. Existing integration tests still run:
each focused/full J31 run includes six executable four-step uninterrupted versus
two-plus-two-step fixtures; the two PCFF controls reject before dynamics.
J22/J23 and historical analysis suites retain their existing short integration
runs. No new scientific campaign, fitting, equilibration study or property scan
is introduced.

## Reproduction and verification

Use `../island-validation/phase4e2-env/bin/python` with
`PYTHONPATH=src:.`, `PYTHONDONTWRITEBYTECODE=1`, `OPENMM_CPU_THREADS=1` and
`OMP_NUM_THREADS=1`. Real OPLS tests select the pinned
`../island-validation/phase4g1-env/bin/python`; AmberTools is
`../miniforge3/envs/md_rg`. No required test is disabled for this correction.
Exact expanded commands and versions are recorded in the additive receipt.

```sh
python -m pytest -p no:cacheprovider -q --tb=short \
  tests/test_phase_4j31_network.py tests/test_phase_4j31_execution.py \
  tests/test_phase_4j31_receipt.py tests/test_phase_4j31_c1_outcomes.py
python scripts/validate_phase_4j31_c1_evidence.py --output /tmp/j31-c1-retained.json
python scripts/validate_phase_4j31_receipt.py \
  docs/evidence/phase_4j31.json docs/paper/reproducibility_manifest.json
python scripts/phase_4j31_c1_receipt.py docs/evidence/phase_4j31_c1/receipt.json
```

Required gates also include J26–J30 analysis regressions, J22/J23 compatibility,
the full supported suite, Ruff, compileall, both applicable environment pip
checks and `git diff --check`. Actual final counts and any remaining unmet gates
belong to the additive receipt; negative-control failures from the before run
are not failed correction gates.

No force-field parameters, charges, typing coverage, numerical tolerances,
network geometry definitions or Linux publication logic are changed. No AI
prediction models or experimental-property or production claims are added.

```text
production_validated=False
simulation_readiness="not_established"
```
