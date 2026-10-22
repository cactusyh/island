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
