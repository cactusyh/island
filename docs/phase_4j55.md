# Phase 4J55: parameter completeness and independent claims

Base: exact fetched `origin/main`,
`6715a9d742b9bf85de37acabb512fbacfa74b9fb` (merged J53).
Branch: `codex/phase-4j55-parameter-completeness-policy`.
There is no J54 on this base. Regression/migration accounting covers actual
J48–J53 records and explicitly records J54 as absent.

This additive policy separates parameter completeness from backend preparation,
runtime, optional cross-engine equivalence and scientific validation. Source
parameters, charge literals, acceptance tolerances, canonical profile selection
and every tracked base file remain unchanged. No new physical profile is created.
The policy and readiness matrix are in [parameter_completeness_policy.md](parameter_completeness_policy.md)
and [forcefield_readiness_matrix_j55.md](forcefield_readiness_matrix_j55.md).

## Public use

```python
from island.forcefield_policy import (
    prepare_polymer_forcefield, ClaimRequirements, build_openmm, export_lammps,
)

# Use installed source assets; no profile or charge-method choice is required.
p = prepare_polymer_forcefield(system, forcefield="opls")
report = p.payload
# Optional independent claim, never implicit cross-engine verification:
p.require(ClaimRequirements(parameter_complete=True))
# Independently requested backend preparation:
p = export_lammps(system, p, "lammps-output")
p, evaluator = build_openmm(system, p)
```

`prepare_polymer` supports the existing PSMILES builder with explicit chemistry
(DP and terminal groups), supplied coordinates/box, and family-only OPLS selection.
`prepare_forcefield(system, ForceFieldRequest(...))` supports existing final graph,
type and charge constraints. Native, literature-derived and user-provided identities
remain distinct. The registered PDMS model stays
`island_external_opls_pdms_xiang2025_dft_local_v1`, external literature-derived
OPLS-compatible for the declared TMS-capped domain, not native universal OPLS-AA
or PCFF. The Xiang paper's Qeq material results are not transferred to DFT/NPA.

The new CLI keeps claims independent:

```bash
python -m island.inspect_parameters --smiles CC --forcefield opls \
  --json report.json --require-parameter-complete
# Administrator installation: ISLAND_OPLS_ASSETS plus its independent SHA256,
# or --assets assets.json --assets-sha256 <pin>.
# Independently requested flags:
# --require-lammps-runtime (also pinned --lammps/--lammps-sha256/--runtime-output)
# --require-openmm-runtime
# --require-cross-engine-equivalence
# --require-scientific-validation
# Optional --build-openmm, --export-lammps, --diagnostic/--diagnostic-sha256.
# PDMS source-audit installation: --paper ... --supplement ...; explicit periodic --box.
```

Inspection returns success when a valid source-bound report is generated unless
a requested claim fails. `--require-parameter-complete` makes incomplete coverage
a nonzero exit. Runtime flags never imply equivalence; equivalence/science flags
cannot promote a complete assignment into an unverified scientific result.

## Current bounded PDMS result

HMDS, capped DP2/DP4/DP6 retain all source interaction rows, exact source charges,
typing and independent interaction enumeration. Necessary source dielectric and
periodic-image/surface conventions are still absent, so their overall parameter
completeness remains `incomplete`, with exactly those two model-convention gaps.
J55 does not cure missing physical conventions by selecting defaults. Seven other
historical gaps remain numerical implementation or pressure/virial qualifications,
separate from the interaction inventory. Their old MD restrictions remain intact.

Diagnostic LAMMPS export/static execution and OpenMM build/static checks are
independently requested in the audit. They use the existing explicit backend
policy with dielectric/image choices labelled backend-owned, without author
promotion. The current source replay and tighter OpenMM solver are not claimed
equivalent. No MD/NPT, pressure qualification, new cross-engine comparison or new
trajectory/RNG equality is claimed. Exact historical component comparisons remain
hashed references in `audit/prior-evidence.json`.

A native finite OPLS ethane control demonstrates parameter completeness without
cross-engine comparison. An existing explicitly analytical external XML fixture
proves that a failed comparison observation does not invalidate complete external
assignment. These software controls do not add PDMS source parameters.

## Historical migration

J47–J49 static/FD/trajectory checks retain their original scope, bytes and claims.
Consumers must not reinterpret `executable`, `numerically_verified` or an MD gate
as a synonym for parameter completeness. J50–J52 four authority layers, nine
permanent author gaps and execution modes remain immutable. J53 intake source
completeness is a source-package contract, distinct from J55 graph-specific
parameter completeness. J54 is absent. Old APIs/CLIs remain historical facades;
new applications opt into `island.forcefield_policy` and `island.inspect_parameters`.
The new report includes its original assignment and source gaps for auditing.
No serialized historical record is upgraded or rewritten.

## Reproduction and validation

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J55_PY=../island-validation/phase4g1-env/bin/python
"$J55_PY" -m scripts.phase_4j55 audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --native ../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --supplement island_ffref/pcff/jp4c08471_si_001.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --output /tmp/j55-audit
"$J55_PY" -m pytest -q tests/test_phase_4j55.py
"$J55_PY" -m pytest -q tests/test_phase_4j{48,49,50,51,52,53}.py
"$J55_PY" -m pytest -q
ruff check src tests scripts
"$J55_PY" -m compileall -q src tests scripts
# pip check in the five environments listed in validation/checks.json
"$J55_PY" scripts/seal_phase_4j55_receipt.py docs/evidence/phase_4j55 --repo .
"$J55_PY" scripts/validate_phase_4j55_receipt.py docs/evidence/phase_4j55 --repo .
git diff --check
```

Receipt/manifest, complete test membership, checks, source hashes, relocated
separate-process reconstruction and unchanged historical file hashes are retained
under `docs/evidence/phase_4j55`. Readiness remains
`production_validated=False`, `simulation_readiness="not_established"`.
There is no cross-engine equivalence, pressure validation, production equilibrium,
Ewald/PPPM closure, universal coverage or long-NPT readiness claim.

Completed validation: **74 focused tests passed**, **486 J48–J53 regressions
passed**, and **4,483 tests passed with 61 existing skips** in the ordinary-suite
union (4,544 distinct nodes, no failures/errors). The union includes all 4,470
base nodes and all 74 J55 nodes exactly once; the separate regression run is not
counted again. Ruff, compileall, five pip checks, diff checks, historical receipt
validators and all **5,237 historical file hashes** passed. Five workflow bundles
were reconstructed in ten relocated fresh processes. Re-signed receipt controls
and final relocated receipt validation are retained alongside the manifest.
