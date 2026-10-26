# Phase 4J42 — explicit external PDMS preparation boundary

No authorized external PDMS bundle path or independent out-of-band SHA-256 was
supplied. J42 therefore stops after executing the integration boundary. All five
fixtures return `asset_unavailable`; there is no executable PDMS profile,
parameter assignment, LAMMPS export, OpenMM Context, minimization or MD result.
No numerical parameter, charge vector or source hash is fabricated.

`production_validated=False`; `simulation_readiness="not_established"`.

After `git fetch origin`, the verified base was
`3825c838246ab2de75bb0c0a494c457c3153fee9`, merging C5 in PR #68.
Published C5 commit `810145f630fd90f63fa69845fc8d3c2751345fa4` is an ancestor.
The branch is `codex/phase-4j42-external-pdms-execution`. The tested implementation
is `2302d5275c924aefd840756b1ac4256fce95f038`. Main is not modified or merged.

## Public preparation contract

`ExternalPDMSAdapter`, `ExternalPDMSOptions` and `ExternalPDMSPreparation` are
exported from `island.forcefields`. The explicit selector is
`external_pdms_class2`, separate from `pcff`, `oplsaa`, `gaff` and `gaff2`.
The reserved requested profile remains
`island_pcff_external_pdms_class2_v1`. This string is not a registered execution
profile or a claim of native PCFF, published COMPASS, OPLS-AA, GAFF/GAFF2,
RESP/Psi4 or universal msi2lmp coverage.

The public input remains a chemically validated final graph constructed from
SMILES or PSMILES. No user CAR/MDF is requested. For example, the absent-asset
boundary can be reproduced directly:

```python
from island.builders import build_linear_polymer
from island.forcefields import (
    ExternalPDMSOptions, ForceFieldRequest, prepare_forcefield,
)

system = build_linear_polymer(
    "[*:1]O[Si](C)(C)[*:2]", dp=2, random_seed=2026,
    add_hydrogens=True, head_end_group="[*]C", tail_end_group="[*]O",
)
result = prepare_forcefield(
    system, ForceFieldRequest("external_pdms_class2", ExternalPDMSOptions()),
)
assert result.payload["status"] == "asset_unavailable"
assert result.payload["execution_profile"] is None
```

An import request instead supplies
`ExternalPDMSOptions(bundle_path=explicit_path, expected_sha256=independent_pin)`.
The pin must be 64 lowercase hexadecimal characters and is never inferred from
the candidate. A missing or malformed pin blocks the request before the bundle
is opened. Source/evidence files resolve only relative to the explicitly
supplied bundle. No environment variable, local parameter search, network
download, native-source fallback or CAR/MDF input is used.

The C5 bundle schema, validators and historical identities are unchanged.
The adapter delegates graph/type/charge/source/row/unit validation to C5. That
includes explicit Si/C/O/H valence and hydrogen completion, stable atom IDs,
atom ordering, molecule/component membership, provenance, exact decimal charges,
all Class-II couplings and equilibrium dependencies, raw-line evidence, mixing
and 1–4 scaling. Native typed/charge records cannot substitute for external
records. Missing coefficients are never assigned zero; charge normalization,
redistribution and cross-family borrowing remain forbidden.

`ExternalPDMSPreparation` is a frozen, JSON-owned, checksummed diagnostic record;
payload access returns a copy. It preserves the complete graph and interaction
inventory, all absent interaction keys, unresolved equilibrium dependency edges,
and C5 import diagnostics when a candidate is supplied. It never impersonates a
`PreparedForceField`. Re-signed execution-authority mutations and stale graph or
component identities are rejected. `create_evaluator` raises a classified
`ExternalPDMSBlockedError` carrying the result before creating an engine.
Existing prepared workflows consequently receive no executable external model.

| Condition | Classified result |
|---|---|
| No explicit path or missing external asset | `asset_unavailable` |
| Missing independent pin, malformed/tampered/incomplete bundle | `parameter_incomplete` |
| Invalid valence, chemical state or unsupported final graph | `typing_blocked` |
| Hypothetical schema-complete data without executed adapter qualification | `asset_unavailable`, stage `backend_qualification_unavailable` |

The last condition distinguishes schema completeness from execution authority.
It has no positive PDMS acceptance evidence in this tranche. A future supplied,
authorized complete bundle still requires independent row reconstruction,
immutable assignment construction, deterministic LAMMPS export, supported
OpenMM construction, energy/component/force and finite-difference checks,
minimization, permitted short dynamics and exact continuation qualification.
No such numerical adapter or qualification is published from absent assets.

## Executed fixture matrix

| Final graph | Atoms | Required/missing interactions | Unresolved dependency edges | Result |
|---|---:|---:|---:|---|
| `[SiH3]O[SiH3]` | 9 | 116 | 167 | `asset_unavailable` |
| Hexamethyldisiloxane | 27 | 656 | 977 | `asset_unavailable` |
| Explicit-H capped PDMS repeat, DP1 | 16 | 319 | 469 | `asset_unavailable` |
| PSMILES PDMS DP2 | 26 | 576 | 852 | `asset_unavailable` |
| PSMILES PDMS DP4 | 46 | 1,090 | 1,618 | `asset_unavailable` |

The repeat is `[*:1]O[Si](C)(C)[*:2]`. Polymer fixtures retain the explicit
methoxy/hydroxyl terminal groups `[*]C` and `[*]O`, hydrogen completion,
coordinates, component identity and polymer provenance. They are finite capped
graphs. Each requirement includes exact sites and equilibrium keys. Inapplicable
impropers retain explicit structural reasons.

Every fixture also imports an intentionally empty C5 template as a negative
control. These controls contain no types, charges or coefficients and return
`parameter_incomplete`, retaining the detailed missing-source and interaction
diagnostics. Their locally computed hashes test integrity only; they are never
represented as user-authorized positive assets. Original and relocated copies
reconstruct the same preparation records in a fresh Python process.

The execution archive retains `system.json`, `preparation.json`,
`incomplete-template.json` and `incomplete-preparation.json` for each fixture,
plus the study report. No external source/parameter-row reconstruction is
claimed: its status is `asset_unavailable`. Unchanged C5 receipt validation
independently reconstructs its historical native raw rows and compiled converter
failures. Those native failures remain failures, and finite-silane identities
remain unchanged.

## Scientific and native-source boundary

The native pinned-PCFF acceptance and external execution acceptance both remain
**1**. C5's source incompleteness is not repaired by this adapter. PCFF and
COMPASS remain isolated; source-state-specific Si and P boundaries are unchanged.
The original native PCFF source identities remain:

| Preserved asset | SHA-256 |
|---|---|
| `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `pcff.rlb` | `82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36` |
| `pcff_templates.dat` | `dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093` |

These are historical native assets, not external PDMS parameter provenance.
There is no supplied external bundle SHA-256 to report. No converter run is
promoted to chemical authority or native msi2lmp parity.

Polymer Ewald/PME remains `ewald_unvalidated`, backend pressure remains
`pressure_unavailable`, cutoff recovery and long NPT remain blocked. No
equilibrium, experimental or production-readiness claim is made. Energy/force,
finite differences, minimization, NVT/NPT, checkpoint and RNG acceptance for an
external PDMS model are unexecuted because the required asset is absent.

## Verification and reproduction

The focused J42 integration suite (the request's focused C6 scope) passed **48**
tests. The supported full suite passed **3,393** tests with **10 existing optional
skips**, no errors and no failures. Its exact J22–J41-C5 subset contains **1,377
passing tests**. Fresh-process relocation preserves all **21 study files**
byte-for-byte. Ruff, compileall, all five applicable pip checks and
`git diff --check` pass. The receipt/manifest validators cover J41 and C1–C5.
The independent historical byte audit compares **2,528** paths directly
with the fetched Git base, including archive-symlink text and target bytes.
J41-C1 through C5 evidence and the C5 schema code are unchanged.

The execution and validation archives contain exact commands, environment
settings, disjoint full-suite shard node IDs, JUnit output, the exact regression
subset of that full execution, reconstruction results and historical hashes.
Regression results are not claimed as a second execution. Test temporary
directories are excluded from the archives.

```bash
git fetch origin
git rev-parse origin/main
git merge-base --is-ancestor 810145f630fd90f63fa69845fc8d3c2751345fa4 origin/main

export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PYTHON=../island-validation/phase4e2-env/bin/python

# Use a fresh output directory. Expected classified acceptance exit: 1.
"$PYTHON" scripts/phase_4j42.py --output /tmp/island-j42-study

# Relocate the data; reconstruct in another Python process. Expected exit: 0.
cp -a /tmp/island-j42-study /tmp/island-j42-relocated
"$PYTHON" scripts/phase_4j42.py --reconstruct /tmp/island-j42-relocated

"$PYTHON" -m pytest -q tests/test_phase_4j42.py
"$PYTHON" -m ruff check .
"$PYTHON" -m compileall -q src scripts tests
git diff --check

"$PYTHON" scripts/validate_phase_4j41_c5_receipt.py docs/evidence/phase_4j41_c5
"$PYTHON" scripts/validate_phase_4j42_receipt.py docs/evidence/phase_4j42

# Only an explicitly supplied bundle and independent pin may be requested:
# "$PYTHON" scripts/phase_4j42.py --case linear_pdms_dp2 \
#   --bundle /explicit/supplier/bundle.json --bundle-sha256 "$INDEPENDENT_PIN" \
#   --output /tmp/island-j42-supplied
```

The validation archive's `run_gates.py`, `run_extra.py` and
`derive_regressions.py` reproduce full, focused, checks, receipt preservation,
relocation and regression derivation. Their command records retain every
interpreter, optional retained-asset environment variable and test node ID.
The five pip checks use the supported, Foyer, analysis, AmberTools and quantum
interpreters listed there. Reproducing optional historical runtime tests requires
the same retained local assets; no external PDMS parameters are substituted.
