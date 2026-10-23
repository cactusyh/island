# Phase 4J38-C1: exact-graph benchmark charge audit

This correction establishes the cause and failure boundaries of the three failed
J38 cells. It does **not** make those cells executable. The acceptance CLI returns
**1**; a faithful diagnostic receipt is separately readable with
`--diagnostic-only` (exit 0). No charge was normalized, redistributed, rounded,
subtracted or otherwise repaired. No force-field source, numerical tolerance,
geometry definition, historical serialized record or historical driver changed.

Base: `14b1b40ae43de22bf510731859e8e6cd337d9926`, verified after fetching origin
against `codex/phase-4j38-real-polymer-benchmark` and `git ls-remote`.
Tested additive implementation: `a07f36db4f56311e7ccffd6aa63a71900f7ded1b`.
The later evidence commit changes documentation/evidence only. Main was neither
modified nor merged. The benchmark inventory remains the same eight cells.

## Aromatic OPLS-AA: independently reproduced source incompleteness

The original native rejection was reproduced before implementation changes:
`OPLSAssignmentError: Native component neutrality failed: -0.23 e`.
The evidence includes the reproduction command, script, exact constructed fixture
and output; these are software failure diagnostics, not scientific validation.

The finite DP2 para-phenylene recipe produces two complete explicit-H biphenyl
chains, each C12H10 with formal charge zero. Their packed final graph is
`f7730225baafb0a2164cdcf3db2f8813d83d3b8d145923053b0fc4ffd193f052`.
The independent auditor does not import the ISLAND graph/charge importer. It:

- Reconstructs connected components from raw bonds and validates the serialized
  component and molecule maps, stable IDs and packing remaps.
- Checks local-to-packed element, formal-charge and repeat/end-group metadata.
  Every repeat has C6H5; each finite terminal attachment carbon has its explicit
  hydrogen. A separately written RDKit biphenyl reference matches the complete
  chemical graph, with no omitted/duplicated atoms.
- Builds a ParmEd topology directly from raw graph records, executes pinned
  upstream `Foyer.run_atomtyping`, and reads charge literals directly from the
  authoritative XML. It verifies the Foyer implementation hashes against the
  retained pin. The driver separately compares native stable-ID type assignments.
- Sums unmodified XML decimal literals using `Decimal` and independently uses
  `math.fsum`. It records every atom's source row, SMARTS definition, element,
  ordering, neighbors, chain/component, formal charge and charge units.

Pinned XML SHA-256:
`c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
The selected `opls_145` aromatic-carbon rule gives **−0.115 e** to all twelve
carbons; `opls_146` gives **+0.115 e** to ten hydrogens. In each chain two inter-ring
ipso carbons have no H: stable IDs 4/7 and 26/29. Thus:

`12 × (−0.115) + 10 × (+0.115) = −0.230 e` per chain.

The total system residual is **−0.460 e**. Decimal totals are exactly −0.230
per component, while compensated binary sums are −0.23. This is not floating-point
summation, serialization or importer rounding. Component/chain membership,
end-group completion and atom ordering agree. The pinned rule selection itself
is not charge-complete for this exact substituted aromatic final graph. The
result is `source_incomplete` at `charge_completeness`; it is not a claim that
all OPLS-AA aromatic systems fail. A separately reviewed complete source would
be necessary before preparing this cell. The zero-charge `opls_147` source row
is explicitly for naphthalene ring-fusion carbons, not biphenyl ipso carbons; it
was not substituted. Raw definitions are retained in `source-context.json`.
Parameter/backend/MD stages did not run.

## PE and aromatic GAFF2: exact RESP assets remain absent

The PE graph is
`1baea92607ac7d99dd85d04d55c22808c6bae9d6501824f7828088a2fca0b66b`.
The aromatic graph is the one above. Neither matches the reviewed retained
linear PEO graph
`1480f6da41f766852310f979bb8f820de11fe1832061b61157acd8577e52172c`
or the retained network graph
`01ec38d6b0177fe2f3f40e747f6c53315034be9b7902fd168b21e69a83e47407`.

The auditor extracts the declared J36/J37 raw archives into a new temporary
location, verifies archive/input/source/executable hashes against retained J37
anchors, and runs the independent raw-source review for both available PEO
fixtures. That review checks atom order, elements, optimized coordinates,
component formal charge, native binary64 charge bits, Psi4/RESP provenance and
finite-precision policy. It **does not** transfer those vectors to the new
polymers. The inventory explicitly limits its search to declared reviewed raw
assets; it makes no exhaustive filesystem-absence claim.

Both new GAFF2 cells are `asset_unavailable` at charge provenance. Required assets
are reviewed raw Psi4/RESP geometry, ESP, native charge vectors and executable
provenance bound to **each exact PE/aromatic graph**, stable atom ordering and
component charge. Psi4/RESP installation alone would not supply these artifacts.
No GAFF2 preparation or downstream MD ran for these cells. Native AM1-BCC remains
blocked with its historical diagnostic; external OPLS-AA charges remain a
separate origin and are never labeled RESP/Psi4.

## Additive matrix and persistence

Every cell now separately reports graph construction, packing, typing, charge
completeness, parameter completeness, executable backend, numerical reference,
equilibration and scientific validation. The five successful J38 cells are
explicitly marked `retained_j38`, with original row/graph/packing identities;
they are not new C1 simulations. Historical receipt validation rechecks their
workflow/analysis evidence, source hashes, relocation and continuation equality.
The failed cells have `new_dynamics_executed=False` and no trajectory claims.
No corrected charge/source gate passed, so **zero new benchmark MD workflows**
were launched. The required regression suites do execute their established short
integration runs; these are tests, not a new scientific campaign.

The companion matrix, fixture, audit, native rejection and RESP inventory survive
relocation and fresh-process reconstruction using declared sources. The validator
rebuilds exact frozen fixtures and independently recomputes the raw audits;
re-signing altered JSON cannot make it pass. Six fresh-process, consistently
re-signed receipt/manifest controls reject altered source charges, graph
membership, RESP graph substitution, a false passed charge gate, a false
acceptance claim and a missing stage. Their archived mutated records are
explicitly **software negative controls, not scientific evidence**. The receipt additionally checks all
artifact hashes, tested implementation blobs, JUnit counts/command results and
historical bytes. All eight cells remain distinct. The PCFF negative control and
native AM1-BCC blocked record are inherited unchanged through the bound J38
matrix, not rerouted into a newly supported domain.

## Reproduction

From the repository on the declared Linux server, with the retained environments
and source files present:

```bash
git fetch origin
git rev-parse origin/codex/phase-4j38-real-polymer-benchmark
git ls-remote origin refs/heads/codex/phase-4j38-real-polymer-benchmark
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c1 /tmp/j38-c1-new-empty-output
# Expected exit 1. The output directory must not already exist.
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c1_receipt /tmp/j38-c1-new-empty-output --bundle --diagnostic-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_c1_receipt docs/evidence/phase_4j38_c1/receipt.json --diagnostic-only
# Omit --diagnostic-only to enforce acceptance: expected exit 1.
```

The receipt's command ledger and archived runners give exact pytest file lists,
flags, environment variables, source/executable fingerprints and validation
commands. In particular the established test configuration retains
`ISLAND_J33_C1_RUN=/tmp/j33-c1-validation/runs`,
`ISLAND_J34_RUN=/tmp/j34-validation/runs` and
`ISLAND_J35_RUN=/tmp/j35-validation/runs`. No marker or required test was disabled.
The five execution environments, including package versions and Python executable
hashes, are recorded in `execution-environments.json`.

## Validation and limits

All runs tested `a07f36db4f56311e7ccffd6aa63a71900f7ded1b`:

| Gate | Actual outcome |
| --- | --- |
| Focused J38-C1 | 39 passed |
| J22–J38 regressions | 803 passed |
| Full supported suite | 2,771 passed, 10 established skips |
| Ruff / compileall / git diff --check | Passed |
| All five environment pip checks | Passed |
| Historical receipt validation and byte preservation | Passed with original blocked acceptance semantics |
| C1 benchmark acceptance | Exit 1: three unmet source/charge cells |

The 10 exact skip reasons, JUnit files and command ledgers are archived. No new
skip or test exclusion was introduced. Final actual counts and commands are
bound in the additive receipt and JUnit files. Historical J17–J38 evidence/manifest/.gitattributes inventory contains
**1,570 files**, checked against base Git blobs and current bytes. Standalone
historical validators are executed where supplied; earlier receipts also receive
JSON readability/hash checks and applicable regression coverage. Historical
blocked acceptances remain nonzero. J35's original source-hash validator runs in
its retained original implementation checkout, as required by that old contract.

Missing exact RESP assets and incomplete native aromatic source charges remain
unmet gates. Short retained trajectories lack sufficient sampling, replicas and
thermodynamic convergence evidence. No experimental agreement, equilibrium,
Tg, mechanical properties, universal force-field coverage, PCFF completion,
AM1-BCC support, universal msi2lmp parity or LAMMPS KSPACE validation is claimed.

`production_validated=False`

`simulation_readiness="not_established"`
