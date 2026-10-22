# Phase 4J35: explicit provided-charge GAFF2 workflows

J35 adds `island.gaff2-provided-charge.v1` and an additive
`island.gaff2-provided-preparation.v1` native preparation wrapper. The reviewed
base is `81e60033d33e979983aee418359986b7951c9b49`. Historical native AM1-BCC
preparation, tolerances, parameters and J17–J34 evidence are preserved.

The charge record binds the exact periodic final graph, canonical stable atom
IDs and elements, original numeric charges, authoritative decimal strings,
component membership and formal totals, elementary-charge units, source and
method identities, provenance, and an independently established reference
identity. A self-signed vector cannot establish source verification: the caller
must supply the independently verified reference identity. Mutations including
ID reordering, booleans, nonfinite values, missing/extra IDs, rounding,
normalization, changed sources/methods or stale graphs are rejected. No charge
is redistributed or normalized. Decimal component totals must equal formal
charges exactly; existing native charge tolerances are unchanged.

The acceptance fixture uses **external OPLS-AA source electrostatics with GAFF2
bonded, Lennard-Jones and 1–4 conventions**. This is an explicit software
acceptance model, not standard AM1-BCC GAFF2, a validated physical hybrid model,
or universal GAFF2 coverage. The published XML charge decimals are used at their
full source precision; short decimal source constants are not padded and called
new high-precision quantum charges.

`scripts/phase_4j35_charge_reference.py` independently reconstructs each vector
from the pinned XML using upstream Foyer SMARTS and a separate bounded chemical
neighbourhood oracle. It reads neither frozen Amber charges nor production
provided-charge/typing/assignment helpers. Each source receipt binds graph,
atom/type mapping, source revision/hash, generator hash, implementation hashes,
environment versions and exact decimal component sums. Both graphs are the
retained complete J34 linear and J33 connected-network fixtures. Evidence uses
a fresh J35 namespace and directory.

AmberTools `antechamber -c rc` is used for native typing and parameterization,
without SQM. Its limited-precision bridge is not the authoritative charge
vector. The final Amber native `CHARGE` section uses `3E24.16` and the unchanged
original provided values. Independent import must agree within **1e-14 e**
(binary Amber charge scaling roundtrip); this tightens precision checking and
does not relax any charge-conservation tolerance. Native component records and
all intermediate hashes are retained separately. A synthetic serialization test
uses more than ten decimal places; it is software test data, not scientific
charge-source evidence.

Diagnostics distinguish `provided_charge_validated`,
`provided_charge_source_unverified`, and `parameterization_incomplete`.
`native_charge_generation_unavailable` remains a separate blocked outcome.
Provided preparation always declares `native_am1bcc_success=False`.
The C2 native acceptance CLI remains nonzero, with its original −0.003 e
rounded-SQM failure and lack of a supported native high-precision export.

Each real case runs independent topology and AmberPrmtopFile energy/force
checks, minimization, 8 NVT steps and 32 NPT steps at 0.5 fs, with four density
blocks. A separate workflow stops at NPT step 16, relocates, and continues in a
fresh Python process. Existing J33-C1 exact comparisons cover retained
coordinates, velocities, boxes, volume/density, energies, temperature, pressure
availability, stages, checkpoint bytes, RNG and barostat state. Pressure is
explicitly unavailable when the existing engine does not measure it; target
pressure is not substituted. These very short runs test software continuation;
equilibrium and density convergence are not established.

Reproduce from repository root using the declared Foyer/OpenMM environment:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
PY=../island-validation/phase4g1-env/bin/python
"$PY" -m scripts.validate_phase_4j35 \
  --runtime docs/evidence/phase_4j34/runtime.json \
  --runtime-identity cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab \
  --output /tmp/j35-new-empty-run
"$PY" -m scripts.validate_phase_4j34_c2 --output /tmp/j35-new-native-control
# Expected exit 1: native fresh acceptance remains blocked.
"$PY" -m scripts.validate_phase_4j35_receipt docs/evidence/phase_4j35/receipt.json
```

Activation of pinned AmberTools 24.8 is scoped to a child process. The runtime
manifest checks executable resolution, package/library and source fingerprints.
Offline receipt validation (`--offline`) verifies archived bindings without
claiming current asset availability. Original J23/J32 serialized schemas and
Linux atomic publication are unchanged; only additive J35 schema dispatch is
introduced.

The final receipt and archived command logs provide actual commits, counts,
environments, results and remaining limitations. Readiness remains exactly
`production_validated=False`, `simulation_readiness="not_established"`.
No experimental validity, equilibrium, Tg, mechanical properties, production
readiness, universal force-field or msi2lmp coverage is established.

## Recorded validation

Tested implementation: `c38cc40055c661f28238d0aa49f895140540dccc`.
The subsequent evidence commit adds validation results without changing that
implementation. Its exact final SHA is verified on the pushed branch.

| Gate | Actual outcome |
| --- | --- |
| J35 focused | 89 passed |
| J22–J35 regressions | 548 passed |
| Full supported suite | 2,516 passed, 10 skipped, no failures/errors |
| Ruff / compileall / diff checks | passed |
| Four execution-environment pip checks | passed |
| Historical protected files | 855 preserved byte-for-byte |
| Re-signed receipt negative controls | 11 rejected |
| Linear / network provided-charge workflows | both complete |
| PCFF same-network control | verified expected typing rejection, no NPT |
| Native AM1-BCC fresh acceptance | blocked, CLI exit 1 |

The 10 skips are existing optional QM/reference-asset tests: one opt-in native
QM integration, one analysis-reference root, seven archived Amber reference
manifest cases, and one OPLS XML opt-in. They are recorded individually in
`artifacts/skipped-tests.json` and the full JUnit file. The two mandatory J35
native-tool provided-charge workflows, relocation and fresh-process checks
actually executed; they were not skipped or replaced by synthetic trajectories.

The additive [receipt](evidence/phase_4j35/receipt.json) binds 228 archived run
files, both independent source references, native inputs/outputs and hashes,
complete continuation comparisons, four Python environments, source/executable
fingerprints, exact commands and JUnit counts. Historical receipt validators
passed their intended modes; C1/C2 native acceptance remains nonzero while
diagnostic validation passes. The original J34 frozen-input result remains a
separate, unchanged historical control.

## J35-C1 provenance correction

The retained J35 successful linear and network trajectories used
`charge_origin="external_oplsaa"`, `charge_method="oplsaa_provided"` under the
additive C1 classification. Their historical serialized method strings and
reference identities are preserved. **They are not RESP/Psi4 results and are
not evidence of fresh RESP or native AM1-BCC generation.** No charge origin is
inferred from use of GAFF2 bonded/LJ parameters.

[Phase 4J35-C1](phase_4j35_c1.md) adds a separate, strict RESP/Psi4 evidence
loader. RESP/Psi4 acceptance is blocked in the declared environments: neither
Psi4 nor a RESP fitter/result is available. No RESP-based trajectory has been
executed. Historical J35 external-charge workflow success remains a separate
control. Only a genuinely executed, independently reviewed Psi4 ESP/RESP fit
can enter the new provider; changing the labels on an OPLS-AA vector is rejected.
