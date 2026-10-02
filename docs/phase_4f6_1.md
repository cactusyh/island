# Phase 4F6.1 — aromatic oxygen and serialization acceptance

## Repository and scope

Correction on `codex/phase-4f6-generic-fragment-charges`, based on reviewed
`1339d580fbac558f679fb92405f940e24f30bb6a`. Fetched main remained `a7b57da`;
Phase 4F6 was unmerged. Phase 4F5, historical records and checkpoint directories
were preserved. No historical schema, engine default or charge policy changed.

## Reproduction and chemical correction

`prepare_capped_fragment("[*:1]Cc1ccoc1C[*:2]")` failed with
`FragmentChargeError("Invalid aromatic valence")` before the fix. The new focused
suite reproduced two failures at that construction boundary while C/N/S controls
passed. Neutral furan-like O now requires exactly two bonds, both aromatic with
order 1.5, to heavy neighbours. Thus there is no O–H or third substituent; the
formal-charge/radical restrictions continue to apply. The aromatic 1.5 convention
must not be interpreted as ordinary oxygen valence three. Other aromatic element
rules are unchanged.

Provided-charge construction and mapping do not require Amber support. A backend
can retain narrower capabilities. Tests use explicitly synthetic charges for
furan DP1/DP3, coverage/conservation/ownership, rechecksummed malformed oxygen,
C/N/S controls and offline mapping with RDKit/ParmEd/OpenMM/SciPy imports blocked.
No furan QM was run.

## Two distinct numerical contracts

Projection and cap-transfer algebra remain checked at **1e-12 e**. The original
Phase 4F6 PE DP3 interoperability experiment remains failed under its declared
1e-12 e molecular import tolerance. Its declaration, records, charges, artifacts
and outcome were not overwritten or re-signed.

Read-only inspection of its actual formats established:

- `charges.txt`: fixed ten decimal places (`.10f`).
- Typed MOL2: six decimal places in the charge column.
- Prmtop CHARGE: `%FORMAT(5E16.8)`, Amber-scaled charge, divided by **18.2223**
  on ParmEd readback. Direct parsing independently agreed with ParmEd.

The reused vector has 20 sites, each with magnitude below 1 e. Its values already
lie on the six-decimal grid to floating-point accuracy: summed quantization
error is `5.5511e-17 e`. **An arbitrary 20-site vector could lose up to 1e-5 e
through six-decimal rounding alone; the following bound is case-specific.**
For scaled magnitudes below 18.2223, half an E16.8 last-place unit is at most
`5e-8` scaled units. Adding input rounding gives the conservative bound
`20*(5e-11 + 5e-8/18.2223) + MOL2_quantization = 5.5878e-8 e`.
This assumes ordinary rounding to the inspected formats; actual per-site values
are also checked, so the bound is not a substitute for execution validation.
The declaration sets **1e-6 e molecular** and **1e-8 e per-site** acceptance.
No final normalization or global engine tolerance change is made.

Tests accept small serialization changes and independently reject both neutral
compensating per-site corruption and accumulated molecular error.

## Separately declared real execution

A declaration was published before running in
`island-validation/phase4f6-1-interoperability/`. That launch failed dependency
preflight because executables were absent from PATH; no parameterization stage
ran. Its diagnostic remains intact. The environment was corrected explicitly,
and an identical scientific declaration was published in a new directory
`island-validation/phase4f6-1-interoperability-live/` before its single execution.
There were no tool-stage retries or changes to scientific settings.

Reproduce with the existing audited environment:

```sh
python -m scripts.validate_fragment_interoperability \
  --source /path/to/phase4f6-fragments --output /path/to/NEW
PATH=/path/to/ambertools/bin:$PATH python -m scripts.validate_fragment_interoperability \
  --source /path/to/phase4f6-fragments --output /path/to/NEW \
  --execute-declared --amberhome /path/to/ambertools
```

The command verifies every original PE artifact checksum and both retained
records, including the template-to-assignment relationship. It reuses the saved
PE DP3 system/coordinates and supplied vector without rebuilding or embedding.
It refuses execution in a directory with previous execution artifacts.
A failed required gate returns nonzero and preserves diagnostics.

**Live acceptance passed**, using AmberTools 24.8, GAFF2 data SHA256
`14ad62c8e532c47e2e400e2ca6ad8052b33bda4c4b9bc6f49b6a528e527512be`.
Antechamber `-c rc -cf charges.txt`, parmchk2 and tleap ran once. No new fragment
QM or SQM invocation occurred. The existing backend checked atom lineage, exact
coverage, connectivity and per-site charges. Preparation integrity passed and
an owned 20-site `ParameterizedSystem` snapshot was created. Structured outcome
provenance binds the retained assignment/template and transformations; the
signed preparation charge source links its assignment identity.

| Measurement | Original failed run | New run |
| --- | ---: | ---: |
| Supplied total, e | -6.9389e-17 | -6.9389e-17 |
| Typed MOL2 total, e | -1.3878e-17 | -1.3878e-17 |
| Prmtop/ParmEd total, e | -5.4877814792e-10 | -5.4877814792e-10 |
| Maximum supplied/readback difference, e | 2.7438906702e-10 | 2.7438906702e-10 |

[Machine-readable evidence](phase_4f6_1_acceptance.json) contains the predeclared
policy, commands, tool/data hashes, artifact checksums, original identities,
per-site comparisons, preparation/import signatures and preflight diagnostic.
These are actual AmberTools results, separate from injected software tests.

## Verification and limits

The full ordinary suite passed **991 tests**, with **9 existing opt-in skips**,
in 80.02 seconds. All **34** focused fragment tests passed. Ruff, pip check and
the cached fragment-charge example passed. The real interoperability command
exited 0, separately from these software tests. Original source checksums were
rechecked after execution. No source workflow or failed experiment was repaired.

This completes a bounded GAFF2 interoperability check, not independent charge
accuracy validation. Isolated fragment charges still approximate a polymer's
chemical environment. `production_validated=False` and
`simulation_readiness="not_established"` remain mandatory.

The planned next product milestone is force-field-specific selection and
compatibility followed by an explicitly versioned OPLS-AA backend. Further GAFF
charge research, RESP or context-containing fragments are not prerequisites for
starting it; none is implemented in this correction.
