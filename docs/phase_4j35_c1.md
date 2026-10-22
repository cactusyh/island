# Phase 4J35-C1: RESP/Psi4 provenance separation

Base: `7f1a4201c3e8bb82628acf5f01860005874827dd`, verified by fetch and
`git ls-remote origin refs/heads/codex/phase-4j35-provided-charge-gaff2-workflow`.
Corrective branch: `codex/phase-4j35-c1-resp-psi4-provenance`.

**RESP/Psi4 acceptance is blocked (`asset_unavailable`).** Each of the four
declared execution environments lacks the `psi4` and `resp` Python modules;
`psi4` is absent from PATH. No pinned Psi4 runtime, RESP fitter, or independently
validated Psi4 ESP/RESP-fit result for either retained graph is declared.
Installing an arbitrary unpinned package would not establish that provenance.
The receipt archives executable paths, Python versions and exact audit commands.

The original J35 receipt was revalidated before production edits. It retains
actual linear/network external OPLS-AA electrostatics with GAFF2 parameterization
and NPT continuation. The original references and artifacts are unchanged.
C1 classifies that historical path as `external_oplsaa` / `oplsaa_provided`.
It establishes neither RESP/Psi4 nor native AM1-BCC reproducibility.

## Validated-result provider

`island.forcefields.ambertools.resp.RESPCharges.load(directory, graph,
expected_identity=independently_reviewed_manifest_identity)` accepts the additive
`island.resp-psi4-result.v1` schema. It is a **loader**, not a Psi4 runner or RESP
fitter. A trusted identity must come from independent review of the actual
quantum/fit calculation; self-signed declarations alone are insufficient.

The manifest binds `charge_origin=resp_psi4`, `charge_method=RESP`,
`quantum_software=Psi4`, Psi4 and RESP versions/implementation, quantum method
and basis, explicit geometry-optimization execution/settings/input geometry,
formal charge, multiplicity, final graph and ordered stable atom mapping,
normalization/redistribution prohibitions, provenance and evidence references.
Every required evidence role has a distinct local file and SHA-256:
Psi4 input/output, ESP grid, RESP input/output, full-precision charge-vector JSON,
independent verification, Psi4 executable, fitter source and environment manifest.
Paths must remain inside the bundle. The manifest identity binds all fields.

An independently reviewed `island.resp-independent-verification.v1` certificate
must confirm actual quantum and fit execution and bind the exact settings,
source/intermediate hashes and reference-vector identity. This review is an
external scientific trust boundary: ISLAND does not prove an ESP fit by looking
for a word in a log. The high-precision vector exporter must be part of that
review, with original native fit output retained; a rounded MOL2/SQM vector
cannot become authoritative by being written to JSON.

The nested independent reference uses the existing provided-charge validator:
canonical stable atom IDs/elements, finite nonboolean charges, elementary-charge
units, authoritative decimal strings, exact component and total formal charge,
no repair. File/vector equality rejects reordering, rounding or truncation.
`provided(directory, graph, expected_identity=...)` reloads and revalidates all
assets before bridging to `GAFF2ProvidedCharges`; native AM1-BCC success remains
false. The existing generic external-charge path and its serialized identities
are unchanged. No quantum settings, parameter/tolerance or publication changes
are made to historical workflows.

## Acceptance and reproduction

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
PY=../island-validation/phase4e2-env/bin/python
"$PY" -m scripts.validate_phase_4j35_receipt docs/evidence/phase_4j35/receipt.json
"$PY" -m scripts.validate_phase_4j35_c1 --output /tmp/j35-c1-new-empty
# Exit 1, preserving result.json; this is not a successful acceptance run.
"$PY" -m pytest tests/test_phase_4j35_c1.py -q
"$PY" -m scripts.validate_phase_4j35_c1_receipt docs/evidence/phase_4j35_c1/receipt.json --diagnostic-only
# Readable/consistent blocked evidence: exit 0.
"$PY" -m scripts.validate_phase_4j35_c1_receipt docs/evidence/phase_4j35_c1/receipt.json
# Acceptance remains exit 1.
```

The diagnostic CLI deliberately cannot turn runtime inventory completion into
scientific acceptance. It never falls back to external OPLS charges. Both
RESP charge generation and downstream GAFF2 preparation, minimization, NVT/NPT,
continuation, relocation and fresh-process MD remain **not executed**. Therefore
there are no new RESP-based coordinate, energy, checkpoint, RNG or barostat
comparison results. The retained native AM1-BCC −0.003 e failure remains blocked.

Tests include synthetic schema-only reviewed manifests (explicit software test
fixtures, not scientific evidence), genuine OPLS-reference rejection, re-signed
metadata mutations, altered intermediate/executable/environment assets, vector
mutations, relocation, fresh-process loading, nonzero CLI diagnostics and false
receipt claims. Actual test counts and exact commands are archived in the
additive receipt. Historical evidence hashes are checked read-only.

Required external assets to unblock: a declared versioned Psi4 executable/library
and environment; a pinned documented RESP fitter/exporter; actual calculations
and independently reviewed ESP/fit outputs for both retained exact final graphs.
Only then can fresh RESP-based GAFF2 dynamics establish the remaining acceptance
gates. No RESP success is claimed here.

`production_validated=False`; `simulation_readiness="not_established"`.
No equilibrium, density agreement, Tg, mechanics or production claims are made.
