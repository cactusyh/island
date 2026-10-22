# Phase 4J34-C2: native AM1-BCC precision investigation

Fresh GAFF2 acceptance remains **blocked**. No authoritative higher-precision
native AM1-BCC vector was found for the declared AmberTools 24.8 environment.
No charges were normalized, redistributed, repaired or substituted. The original
component conservation tolerance remains `1e-4 e`.

Base: `03682038e067a3ce54132adb4beaba84626b8794`, verified against origin's
`codex/phase-4j34-c1-fresh-am1bcc-charge-conservation` branch. Work is on
`codex/phase-4j34-c2-native-charge-extraction`; main is not merged or modified.
The additive receipt names the implementation commit actually tested. Its
containing Git commit is the subsequent evidence commit (a receipt cannot
contain its own commit SHA); the final report records the remote branch SHA.

## Reproduction and source investigation

Before changing production code, the unchanged C1 strict importer replayed the
exact retained fresh native outputs from `verified-fresh.tar.gz`. This is a
read-only replay, **not another fresh QM generation or scientific campaign**.
Component IDs 31–60 have formal charge 0 e. The native SQM printed atom charges
sum to −0.003 e although its separately summed total prints 0.000 e. PRE.AC,
BCC.AC and MOL2 retain −0.003000 e; prmtop decodes to approximately
−0.002999997805 e. The original importer rejects this at
`import/component_total_validation`. No charge record is published.

The official source tarball comes from
<https://ambermd.org/downloads/AmberTools24_rc5.tar.bz2>, SHA-256
`52fb4fb3370a89b7ce738a2dc3e513c2fc1943fde4b4381846d9e75cc48d840f`.
That is the pinned conda recipe's source hash. The package applies official
AmberTools 24 updates 1–8; their complete path inventories and checksums are
retained. None changes the relevant SQM/antechamber source. Exact native source
bytes and recipe files are archived without editing their whitespace.

The [study](evidence/phase_4j34_c2/study.json) binds the official download,
manual excerpts, source review, native probes, original native outputs, runtime
and final graph. Source and independently executed native observations agree:

- `qm2_calc_mulliken` calculates charges into `_REAL_ scf_mchg` using the density
  matrix; `dprec.fh` defaults to double precision. These internal arrays are not
  a documented SQM charge-extraction API. The SQM static library is part of the
  Sander implementation, rather than a documented standalone charge accessor.
- Both final and intermediate AM1 charge output call `qm2_print_charges`, which
  writes atom values with `F14.3` and the independently summed total with
  `F12.3`. SQM's supported CLI exposes input/output files, overwrite and help.
- The Amber24 manual's SQM options document `printcharges` and `verbosity`,
  not an AM1 charge precision setting. Six real native AM1 water single-point
  probes (`maxcyc=0`, `printcharges=1`, verbosity 0–5) all exit 0 and retain
  three-decimal atom output. These are **capability probes**, not acceptance MD.
- Antechamber's `rsqmcharge` reads that printed table with `%lf`. Its in-memory
  `double` array consequently starts with already rounded data. Equivalent-atom
  averaging and BCC preserve the total residual. `wcharge` writes `%10.6lf`;
  native `-c wc -cf charges.qout` conversion reproduces −0.003000 e.
- AC, MOL2, QOUT and Amber prmtop can store more decimal places, but this charge
  path has lost the upstream digits before those files are written. Padding
  zeros or reformatting values is not higher-precision extraction.
- DFTB-specific/debug charge output, external Fireball output, RESP and MOPAC
  are different methods or runtimes. They cannot establish the required native
  AM1-BCC vector. The documented Sander Python API does not expose the required
  AM1 charge vector; its SQM print options are suppressed.

`AMBER SQM VERSION 19` is the engine banner, not the package release. The
runtime manifest separately binds AmberTools 24.8, build, binary/wrapper hashes,
libraries, source files and executable resolution. This conclusion is limited
to that pinned runtime and the inspected official source/manual; it does not
claim that every AmberTools version lacks such an interface.

## Diagnostic contract

`native_precision.py` implements an additive, bounded investigation contract,
not a substitute charge engine. It never returns an authoritative vector, even
for the other component whose rounded total happens to pass. It checks the
reviewed study identity, final graph, component case, stable atom IDs/names and
elements, runtime/executable/source identities, input and every native
intermediate hash, original strict importer outcome, and unchanged tolerance.
Any stale, reordered, truncated, rounded, zero-padded or mutated intermediate
is rejected. A re-signed diagnostic cannot claim publication or dynamics.

The CLI accepts no frozen charge input and preserves `outcome.json` before
returning nonzero. `--diagnostic-only` returns 0 solely for validation of the
blocked diagnostic; it does not satisfy acceptance. Reconstruction after a
path change and in a fresh Python process is tested **for the diagnostic**.
It is not claimed as a fresh GAFF2 MD continuation.

Fresh preparation, minimization, NVT, NPT, split/uninterrupted comparison,
checkpoint relocation and fresh-process MD continuation remain unexecuted in
C2 because the prerequisite native charge gate is unmet. The historical
prepared-input J34 success is preserved separately, byte-for-byte.

The external requirement is an official, documented, independently verified
native AM1-BCC interface/export that retains sufficient actual charge precision
and stable atom association, in a pinned runtime. An upstream-supported change
would require a separately reviewed extraction implementation and fresh charge
and downstream acceptance. This phase does not patch SQM's print format, scrape
process memory, reconstruct lost digits or silently switch methods.

## Reproduction and verification

From the repository root on the declared Linux server:

```sh
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_c2 --output /tmp/j34-c2-new-empty-output
# Expected exit 1: diagnostic retained, fresh acceptance blocked.
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j34_c2_receipt docs/evidence/phase_4j34_c2/receipt.json --diagnostic-only
# Expected exit 0: valid blocked evidence. Omit --diagnostic-only for acceptance exit 1.
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_phase_4j34_c2.py
```

Use a previously nonexistent output directory. The evidence archive includes
exact native-probe commands, input/output hashes, official-source retrieval
scripts, environment versions, JUnit results, complete regression/check command
records and skip reasons. Native probes use scoped activation of the declared
runtime; historical regressions do not globally activate AmberTools.

J22–J34 regressions and the supported full suite include existing short
integration runs and retained workflow revalidation. They are not a new C2
scientific campaign and cannot close the blocked fresh-charge acceptance gate.
Software mutation fixtures and retained failure replay are explicitly negative
controls, never positive scientific evidence.

`production_validated=False` and `simulation_readiness="not_established"`.
This work establishes neither equilibrium nor experimental agreement, Tg,
mechanical properties, universal force-field coverage or production readiness.
