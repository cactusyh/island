# Phase 4J80 — authorized PDMS source and QM-pilot gates

The decision is **qm_pilot_blocked**. Gate A failed because the original
Fritz–Hofmann parameter source is unavailable. Gates B and C are **not evaluated**.
Scientific work stopped at Gate A. No QM request, input deck, geometry,
scientific design matrix, fitted parameter or accepted coefficient was created.

## Base and preservation

Branch: `codex/phase-4j80-authorized-pdms-source-and-qm-pilot`.
The exact base is verified J79 commit
`1676a780d1f5aa1ce46e707a13aadeff14ccf39d`; its origin branch matched before J80
was created. Fetched `origin/main` was
`8b483278c17bbd6da8123fe5e346857e12a46394`. Main is untouched; it is not merged.
The terminal receipt records the tested implementation commit independently
of the final evidence commit.

All **10,476 tracked files at the base** are hash-inventoried and preserved,
including J79's result, source hashes, parameter values, receipts, manifests,
native PCFF registry and all historical evidence. J72's four unintended SQM
invocations and 46 verified / 48 possible MD steps remain historical violations;
the C1 quarantine is unchanged. Unrelated untracked files and source assets
remain untouched. Explicit implementation inventories prevent new modules
from broadening historical receipt globs.

## Gate A: source recovery

The inherited J72-C1 kernel/process fence and J73 evidence-only restrictions
were installed before source inspection. A fresh hash check of the exact local
source root
`/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/island_ffref/pcff`
found the same **27 files**, byte-identical to J75, with no additions, changes or
missing files. Their unchanged citations do not identify the requested primary
paper. Existing FRC/COMPASS/external files are catalogued as context and are
not substituted for the reported rows' source.

The exact publisher PDF endpoint for PII `S0032386196006003` was retried and
returned HTTP 403. No full text was recovered or read. No other directory,
paper or model is silently substituted. Publisher metadata remains citation
context, never parameter evidence. No original PDF, FRC or archive is committed.

`source-provenance-manifest.json` preserves the bibliographic identity, all
27 context-file hashes/families/versions/license statuses and all 20 candidate
identities. Every candidate has explicit missing primary-artifact hash,
source version and exact printed location. The mode is
`unverified_user_reported`, distinct from `directly_reported`, `reconstructed`
and `inferred`. J79's exact reported strings remain unchanged, including the
unverified C–Si–O quartic coefficient.

The reusable source gate requires original-author/publisher/native-release
identity, immutable filename, SHA-256, citation, year, version, access/license
status and an independently reviewed, hash-bound authority audit. Each candidate
must match the independently audited raw values, units and exact source
location. Duplicate candidate IDs or source identities and conflicting
convention assertions block the gate. A hash proves bytes, not authorship.
User attestations and fabricated/re-signed records cannot supply authority.

The production `trusted-authority-audits.json` is empty. Synthetic audit fixtures
live only in tests. They cannot enter the frozen production evidence package.

## Gate B: inactive convention contract

The report defines twelve required groups, each with explicit fields and
source-evidence references. Their current values are null, with
`not_evaluated_Gate_A_blocked` status:

1. Source and target van der Waals equations and powers.
2. Energy, length and angular-coordinate units.
3. Sigma, epsilon and minimum-distance meanings.
4. Distance and energy combination equations.
5. 1–2/1–3 exclusions and separate LJ/Coulomb 1–4 scaling.
6. Charge origin, source-type mapping, repeat and capped-chain neutrality.
7. Bond/angle polynomials, prefactors, ordering and equilibrium dependencies.
8. Torsion equation, sign, phase, periodicity, normalization, origin and order.
9. Improper/out-of-plane equations or documented nonapplicability, with order.
10. Cross-term equations, ordering, equilibrium dependencies and coverage.
11. Source atom types, chemical environments, exact graph mapping and caps.
12. Finite/periodic scope, electrostatics, cutoff and tail policies.

A filled field must exactly match independently audited source assertions.
Source 9–6 parameters cannot be coerced to 12–6. The 9–6 distance must retain
its minimum-energy meaning. Element symbols cannot establish exact source
atom types. Charge neutrality requires explicit source-type charges and
integer compositions for the repeat and each pilot topology, with independent
charge-sum arithmetic. Repeat neutrality alone does not qualify caps.

J79's read-only comparisons remain relevant but unqualified: generic ISLAND LJ
is 12–6; isolated PCFF code and the J78 proposal describe 9–6. The paper's mixing,
exclusions and Fourier conventions are still unknown. Transferring the reported
repeat charges to the capped molecules would yield +0.240 e. J80 neither fixes
that charge nor selects a missing convention from another force field.

## Candidate and target classifications

The six classes are `exact_native_row`, `exact_native_row_convention_blocked`,
`literature_candidate_verified`, `literature_candidate_convention_blocked`,
`inferred_or_analog_value_blocked` and `end_group_scope_blocked`.

All **20 reported candidate rows** remain
`literature_candidate_convention_blocked`; **zero** literature rows are verified
and **zero** native rows are resolved. Literature verification would still not
be native authority. Only an independently exact native source/type/row,
environment, equation/unit/phase and policy identity could receive
`exact_native_row`. This is a prerequisite for consideration in a future phase,
not runtime authorization in J80.

The separate 129-target audit preserves original J78/J79 identities. End-group
scope takes precedence for 69 targets; 20 non-cap targets have convention-blocked
literature candidates; 40 other targets block any inferred missing-value fill.
The latter class describes prohibited inference, not newly supplied analog
values. Orthogonal flags preserve **87 unresolved cross terms** and
**69 unresolved end-group terms**, including 48 overlapping targets. All
129 extension terms, 365 scalars and 11 policies remain unresolved. Native
PCFF's 378 outstanding instances and J79's decision remain unchanged.

## Gate C: inactive minimal-pilot request contract

`qm-pilot-request-manifest.json` contains **empty requests and completions**,
`execution_enabled=false`, no execution profile and zero QM decks. A prospective
contract describes the prerequisites without generating configurations or
launching the historical 423-request campaign.

The bounded prospective scope is explicit-H, neutral, closed-shell singlet,
finite, linear trimethylsilyl-capped HMDS and DP2, with DP4 as an already-supported
J78 topology held out from training. DP6 and periodic calculations are excluded
from this pilot. Prospective requests are capped at 36; reaching that cap
does not establish sufficient independent data or identifiability.

The contract requires independent equilibrium references, Si–O/Si–C/C–H
distortions, relevant angle distortions, representative torsions, coupled
cross-term deformations and separate intermolecular nonbonded/electrostatic
cases. One independent coupled geometry family is validation-only, and DP4
is the held-out topology. Sample IDs, geometry hashes and conformer groups
must be disjoint across roles. The historical reserved-validation dataset is
not activated or relabeled as pilot validation.

Software/version, method, basis, dispersion, convergence, constraints, actual
independent coordinates, mapping and input/output hashes remain unspecified
and required before executable requests. No proposed J78 protocol silently
becomes authorized or scientifically adequate. There is no backend launcher
in this intake implementation. A later source-complete intake must establish
the actual execution profile and independent data before any scientific work.

## Identifiability and fitting policy

Six separate blocks are required: bonded, angles, torsions, nonbonded, charges
and cross terms. No joint all-term fitting is permitted. Every physical rank,
condition number, null-space dimension, parameter-correlation, design identity,
uncertainty and bound identity is currently **null / not computed**. No zero
or synthetic diagnostic is substituted for missing scientific evidence.

The software reuses J78's `check_bounds` and `split_check`, the established
pilot's `hmds_design.metrics` and J78's unchanged diagnostic limits. Manufactured
matrices test rank/nullity, singular/infinite conditioning, column correlation,
weak columns and finite bounds. The diagnostic adapter calls no fitting solver
and returns no coefficients. Rank deficiency, poor conditioning, correlation,
missing physical bounds or insufficient independent data keeps a block
unresolved. No GAFF, OPLS, CVFF, CFF91, COMPASS or other analog values are imported.

## Evidence, verification and deterministic replay

Evidence under `docs/evidence/phase_4j80/` includes:

- `source-provenance-manifest.json`, `local-source-check.json`, `source-access.json`.
- `gate-report.json`, `convention-contract.json`, `candidate-row-classification.json`.
- `qm-pilot-request-manifest.json`, `identifiability-and-conditioning.json`.
- `result.json`, sealed `receipt.json` and `manifest.json`.
- Explicit parent/input/historical inventories and archived validation records.

New code is confined to the `pcff_pdms.pilot_intake` subpackage, the J80 builder
and validator, focused tests and this document. No historical module or runtime
registry is changed. Validation pins source inputs and independently reconstructs
every report; merely re-signing scientific metadata cannot grant a gate pass.

All tests and checks run under the inherited C1 fence and evidence-only policy.
Focused tests include provenance/hash mutations, duplicate identities,
metadata/attestation rejection, all convention fields, 9–6/12–6 separation,
exact typing, cap charge arithmetic, family separation, classification,
training leakage, rank/condition/bound failures and no runtime activation.
Fake native executables and a fake backend module prove child-process denial;
read-only module discovery remains permitted. No real backend is invoked.
Scientific controls are scope skips. The archived terminal collection and JUnit
identify unique focused/regression counts; development failures are retained
separately. An initial synthetic test exposed shared mutable audit references,
which were isolated before terminal verification. The import-denial fixture was
also corrected to use an available fake module rather than an absent backend.

J73–J79 evidence-only regressions, Ruff, compileall, supported pip checks,
J79/J78/J76/C1 strict receipt checks, historical preservation, offline replay
and relocated fresh-process replay are recorded. A full scientific suite is
not run because Gate A prohibits scientific execution.

Use fresh output directories:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j80-replay -- -m scripts.phase_4j80_pilot replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j80-validate -- -m scripts.validate_phase_4j80_receipt
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j80-tests -- -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j73.py tests/test_phase_4j74.py tests/test_phase_4j75.py tests/test_phase_4j76.py tests/test_phase_4j78.py tests/test_phase_4j79.py tests/test_phase_4j80.py
```

Replay uses frozen evidence only: no network, local source directory, scientific
backend, PDF parser or live QM data is required. A scientific `qm_pilot_blocked`
result is a successful integrity replay, not a passed scientific gate.

Final counts are zero completed QM requests, fitted parameters and accepted
parameters. Automatic typing, assignment, registration, production export/MD
and analog fallback remain disabled. Production validation is false and
simulation readiness remains `not_established`. Source completeness, convention
closure, chemical transfer and scientific/production validation are not claimed.
