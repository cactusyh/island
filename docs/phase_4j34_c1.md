# Phase 4J34-C1: fresh AM1-BCC charge-path investigation

Base: `74ff640ccad8d5f0c1dd88e172f061b800489d84` on
`codex/phase-4j34-gaff2-periodic-npt-acceptance`. The corrective branch is
`codex/phase-4j34-c1-fresh-am1bcc-charge-conservation`. Main is untouched.
The additive receipt identifies the implementation actually tested. Its later
evidence commit is independently identified by Git; a receipt cannot contain
its own eventual commit hash without changing that hash.

## Outcome and boundary

Fresh linear AM1-BCC regeneration remains **blocked** in the declared
AmberTools 24.8 runtime. The original strict conservation check correctly
rejects native output. This correction adds an auditable explanation and a
fresh-only acceptance entry point; it does not repair, normalize, redistribute,
or alter charges, relax tolerances, change parameters, or rebuild AmberTools.

The frozen-input J34 periodic NPT result is retained as a separate historical
control. It proves the workflow for those prepared inputs, not reproducible
from-scratch AM1-BCC generation. No new minimization, NVT, NPT, continuation,
or network campaign is authorized by a failed fresh linear preparation. Short
integration tests in the regression suite are reported separately from that
unexecuted fresh scientific workflow.

`production_validated=False`; `simulation_readiness="not_established"`.
There is no equilibrium, converged density, experimental agreement, Tg,
mechanical-property, production-readiness, or universal force-field claim.

## Reproduction and charge trace

Before changing production code, the exact retained failure was replayed
through the original importer, and a new unmodified J34 CLI fresh preparation
was executed without `--linear-preparation`. Both rejected component charge
conservation. The new CLI also performs an independent fresh preparation.
Commands, exit codes, native outputs, hashes, and per-atom audits are archived
under `docs/evidence/phase_4j34_c1/`.

The complete 60-atom linear fixture has two neutral 30-atom C8O4H18 chains.
Three independent fresh attempts (original CLI, corrective CLI, and the final
committed implementation) all exited 1 with the same component-total failure.
The final tested implementation is
`aa813ed3d082d623fb5ad4b102a315265af6a731`.
For the retained failing second component (stable atom IDs 31–60):

| Stage | Charge total (e) |
| --- | ---: |
| Final graph formal charge / SQM `qmcharge` | 0 |
| SQM separately printed total | 0.000 |
| Sum of SQM printed atom charges | −0.003 |
| Native `ANTECHAMBER_AM1BCC_PRE.AC` | −0.003000 |
| Native `ANTECHAMBER_AM1BCC.AC` | −0.003000 |
| Native `typed.mol2` | −0.003000 |
| Parsed prmtop | approximately −0.002999997805 |

The audit regenerates the exact input from the final component and requires
byte equality. Input atom names encode stable IDs; element, name, graph and
stereochemical checks use the existing lineage/parser/importer. SQM sequential
atom numbers and elements are checked independently, as are fixed-width AC
names and prmtop name mappings. MOL2-to-prmtop charge serialization differs by
at most about 2.2×10⁻⁹ e for the retained failure. The original importer rejects
at `import/component_total_validation`, using the unchanged `1e-4 e` tolerance.
Failed native directories have no `ISLAND_SUCCESS` marker or published charge
record. Publication, the original engine, and Linux atomic publication are
unchanged.

The discrepancy first appears in the per-atom SQM text, before ISLAND's charge
parser. The installed SQM binary contains the fixed `F14.3` per-atom format;
the native output actually contains three decimal places. The separately
printed total is computed before this individual formatting. Native
antechamber averages equivalent atoms before PRE.AC (for example three methyl
hydrogens become their mean), so individual SQM/PRE values need not match.
That established native operation and the subsequent BCC corrections preserve
the already incorrect printed total. The diagnostic records every value and
does not apply either transformation itself.

Pinned upstream AmberClassic source excerpts corroborate this explanation:
`qm2_print_charges.F90` formats `scf_mchg` using `F14.3`; `rsqmcharge` in
`charge.c` reads that table; `bccharge` performs equivalent-atom averaging and
writes PRE.AC before invoking native BCC. These excerpts are explicitly an
upstream reference, **not an attestation that they are the exact source of
the installed conda build**. Installed executable hashes, embedded format
strings, package metadata, build recipe, and actual native artifacts establish
the observed runtime behavior independently.

The evidence rules out an ISLAND atom-order permutation, lost formal charge,
cancellation between disconnected components, significant prmtop parsing
error, and stale prepared inputs for these failures. Both fresh attempts use
the same verified J34 runtime/source identities; no runtime version change is
needed to reproduce the failure. Fresh SQM optimization outputs can differ
despite identical inputs, and a component whose rounded values happen to sum
to zero is not proof that this precision path conserves charge reliably.
The frozen successful control has a zero printed sum and is audited separately.

The unresolved external requirement is a separately pinned and verified native
charge-generation path that retains sufficient *actual* AM1 Mulliken precision
through BCC and satisfies the existing component tolerance. The discarded
digits cannot be reconstructed from these files. A native source/build change
would require explicit provenance, independent numerical verification and fresh
acceptance; none is silently substituted here.

## Diagnostic and acceptance contracts

`island.am1bcc-charge-audit.v1` binds runtime identity, exact input hash, native
file hashes, stable atoms/elements/formal charges, per-stage atom values and
totals, original importer failure stage, tolerance, and publication marker.
`validate_audit` recomputes the record from native files; consistently re-signing
a changed total, stage, tolerance, runtime, or file inventory cannot make it
validate. Missing/malformed/mismatched artifacts raise rather than being
classified as a native precision limitation.

`scripts.validate_phase_4j34_c1` only permits fresh preparation. It invokes the
existing strict J34 worker without frozen inputs and writes
`island.j34-c1-fresh-outcome.v1` plus original results/logs. Failure exits
nonzero. If preparation succeeds, the existing real minimization/NVT/NPT,
split/relocation/child continuation and comparison gates still apply. This
does not bypass the existing source, parameter, graph or charge checks.

The additive v1 C1 receipt represents the observed **blocked** outcome.
The default receipt validator returns nonzero because fresh acceptance is
unmet. `--diagnostic-only` returns zero only for a complete, internally
consistent blocked diagnostic; its output explicitly says acceptance is
unmet. It reconstructs native audits offline, checks fixture and result
bindings, artifact hashes, historical hashes and readiness flags, and rejects
contradictory dynamics claims. `--offline` skips present-runtime availability
checks, not archived identity checks. Historical schemas and receipt bytes
are unchanged. A future successful receipt format must explicitly validate
the full fresh workflow evidence rather than relabel this blocked record.

## Reproduction

Run from the repository root on the declared Linux server. The manifest retains
the exact Amber prefix, executable paths, native versions/build, libraries,
source hashes and environment fingerprints. Activation occurs only in a
scoped child; legacy tests keep their original environment.

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_c1 \
  --runtime docs/evidence/phase_4j34_c1/runtime.json \
  --runtime-identity cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab \
  --output /tmp/j34-c1-new-empty-run
# Expected nonzero while native fresh charge conservation fails.

../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_c1_receipt \
  docs/evidence/phase_4j34_c1/receipt.json --diagnostic-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_c1_receipt \
  docs/evidence/phase_4j34_c1/receipt.json
# Second command remains nonzero: this is not full GAFF2 acceptance.
```

Exact gate commands, environments, JUnit files/counts, skipped tests, native
preparation duration, current-runtime validation and historical read-only
reconstruction results are in the additive receipt and bound artifacts.
Software mutation controls and archived failure replays are labeled as such;
they do not count as fresh molecular dynamics or positive scientific evidence.

## Executed verification

| Group | Passed | Skipped |
| --- | ---: | ---: |
| J34-C1 focused, final implementation | 27 | 0 |
| J34 | 36 | 0 |
| J33-C1 / J33 | 50 / 7 | 0 |
| J32 | 53 | 0 |
| J31 and C1 | 111 | 0 |
| J26–J30 | 48 | 0 |
| J22/J23 | 91 | 0 |
| Full supported suite | 2391 | 10 |

The full run took 1040.95 seconds. The ten existing opt-in skips request native
Amber execution, additional archived Amber/PE references, or an explicit pinned
OPLS source; exact test names and reasons are in the JUnit and receipt. No tests
were disabled to obtain this result. Fresh native Amber execution was separately
performed three times and failed as described above. Ruff, compileall, all four
pip checks, diff checks, historical receipt checks and the 622-file historical
preservation check passed. C1 diagnostic validation passes; C1 fresh acceptance
validation intentionally remains nonzero. The implementation commit and each
test group's command record remain distinct from the later evidence commit.
