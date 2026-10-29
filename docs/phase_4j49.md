# Phase 4J49: isolate external OPLS PDMS electrostatic error

J49 reproduces the J48 failures and identifies a **generic electrostatic
backend error**, rather than a PDMS charge/type mapping error. The explicit
exact-`erfc` numerical control passes the unchanged energy, force and finite-
difference bounds for HMDS and TMS-capped DP2/DP4/DP6. It also passes bounded
minimization, 40-step NVT, relocation, separate-process continuation, **native
checkpoint byte equality, complete state-dump equality and full RNG equality**.

The **original author PPPM/default-table execution remains blocked**. The
acceptance command returns **1**, and the public `run_md` entry point still
rejects that profile before engine invocation or output creation. The numerical
control is an explicitly modified executable and solver, recorded separately;
it is not silently installed as a replacement author convention or a new
force-field profile. Its bounded software acceptance is 0. Scientific acceptance
remains 1. No original energy is adjusted or subtracted to pass an acceptance
bound.

`production_validated=False`, `simulation_readiness="not_established"`.
No NPT is run. Pressure, polymer long-range scientific qualification, cutoff
recovery and long NPT remain blocked. The 10 fs NVT checks establish short
algorithm reproducibility, with no equilibrium or experimental claim.

## Exact base and unchanged physical model

Fetched `origin/main` was `9e2e5df9cfebf119c0c857c841a401fd98fcb42a`, the
merge of J48. Branch: `codex/phase-4j49-opls-pdms-electrostatic-closure`.
All files from that exact base, including J17–J48 historical receipts, bundles,
serialized identities, source records and implementation, remain byte-preserved.

The canonical public family remains **`forcefield="opls"`**, with no profile or
charge-method argument. The resolver still selects
`island_external_opls_pdms_xiang2025_dft_local_v1`, an external literature-derived
OPLS-compatible local-research profile. Its internal identity, domain, exact
charge literals, sources and source-row references remain visible. This is
neither native universal OPLS-AA nor PCFF/COMPASS PDMS support.

```python
from island.polymer import (
    OPLSAssets, use_opls_assets, prepare_polymer_forcefield,
)

assets = OPLSAssets(
    "island_ffref/pcff/PDMS_benchmark-main.zip",
    "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4",
)
with use_opls_assets(assets):
    prepared = prepare_polymer_forcefield(system, forcefield="opls")
# prepared.payload retains selected_profile, source_identity, charge_origin,
# domain, assigned terms, backend status and scientific limitations.
```

`island.electrostatic_audit.bind_inputs` adds an immutable numerical audit binding
without redesigning the registry or changing any historical API. Integrity
reprepares from the actual graph and pinned archive, rejecting even re-signed
edits. Coordinates/box, atom IDs/order, source/charge/type identity, masses,
source model, exclusions, special pairs, units and backend conventions are bound.

## Root cause and raw comparisons

The paper is Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane:
Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal
Conductivity*, JPCB **129** (2025), 1864–1873, DOI
**10.1021/acs.jpcb.4c08471**. The local PDF was re-read offline; Methods and
Figures 7/8 distinguish DFT/NPA and Qeq. This profile retains the author ZIP's
DFT/NPA branch. The Results' Qeq/thermophysical claims are not transferred to it.
Local research authorization does not establish redistribution rights; the ZIP
and PDF are not committed.

Every run keeps the same final graph, explicit-H atom order, coordinates,
40 Å periodic cubic box, masses, charge literals, 1–2/1–3 exclusions, 0.5 Coulomb
and LJ 1–4 scaling, geometric sigma/epsilon mixing, CHARMM r-squared 9–11 Å LJ
switch, 11 Å Coulomb cutoff, dielectric 1 and conducting periodic electrostatics.
The numerical solver/table/erfc differences are recorded explicitly.

| Fixture | Author Coulomb error, kJ/mol | Exact-kernel control error, kJ/mol | Control max force error, kJ/(mol Å) |
| --- | ---: | ---: | ---: |
| HMDS | 0.0699658923 | 0.00000978363 | 0.00000190698 |
| DP2 | 0.1436606402 | 0.0000100741 | 0.00000262941 |
| DP4 | 0.4018122426 | 0.00000949840 | 0.00000219366 |
| DP6 | 0.5085865935 | 0.0000147126 | 0.00000388548 |

The inherited bounds are **3e-4 kJ/mol**, **3e-4 kJ/(mol Å)** and
**1e-4 kJ/(mol Å)** for FD. No bound, charge, coefficient or force-field family
was changed. All comparisons use raw engine results.

The control ladder separates:

1. Author PPPM `1e-4`, default pair table versus OpenMM Reference Ewald `1e-7`:
   **Ewald/PPPM convention** and **reciprocal tolerance**. Their real/reciprocal
   partitions use different alpha values and cannot individually be equated.
2. Disabling only the pair table: interpolation **precision** contributes, but
   PPPM mesh energy and force errors remain. This is not sufficient closure.
3. Converged PPPM (`1e-10`, order 7, table off): energy and FD still fail.
4. Matching Ewald alpha, converging reciprocal truncation (20 per axis), and
   disabling the table: the pinned pair style's Abramowitz–Stegun `erfc`
   polynomial still disagrees with the exact-erfc equation. Its energy and force
   are not precisely mutual derivatives, causing the remaining FD discrepancy.
5. In a **separate local source copy**, replace only the 3 CHARMM and 2 pure-
   Coulomb `erfc` polynomial evaluations with `std::erfc`. Preserve force
   equations, the Coulomb conversion constant and every physical input. Raw
   energies, all-atom forces and sampled FD now pass. The executable hash changes
   and is recorded; the original source checkout and pinned binary are intact.

A charged dimer and a dimer straddling the periodic boundary reproduce the error
without PDMS typing. Both pass with the explicit exact-erfc control. The remaining
small raw residual includes the **unit conversion/precision** difference between
OpenMM's `1389.35457644382` and LAMMPS's `332.06371*4.184` Coulomb constant.
Neither constant is modified and no predicted residual is subtracted.

Bond/angle/torsion/LJ agree. Improper terms are source-declared inapplicable;
there is no fabricated zero row. LJ14/Coulomb14 are independently reconstructed
from graph distances in the preserved J47 comparison. LAMMPS `ecoul` includes
exclusion corrections and is never mislabeled as standalone Coulomb14.

The combined LAMMPS pair style does not expose separate per-atom Coulomb and LJ
force dumps. Total forces are raw measurements; the Coulomb force diagnostic is
explicitly labeled as LAMMPS total minus independently verified bonded/LJ
analytical forces. Independent component energy differences check this
partition, rather than claiming a nonexistent native observable.

## Finite differences and bounded continuation

The standalone FD driver imports no source parser. Each perturbation changes one
coordinate, records positive and negative raw bonded/LJ/Coulomb/total energies,
and computes `-(Eplus-Eminus)/(2h)` at **2e-4, 1e-4 and 5e-5 Å**. It preserves
analytical forces, errors and explicit sign/unit-error controls. Each PDMS
fixture samples all axes of a representative H/O/Si/C (36 comparisons per
component per backend setting). Both charged controls sample every coordinate.
Force comparison itself includes every atom. No mean-error gate hides a failed
sample. Independent stdlib reconstruction checks the raw LAMMPS log energies,
component sums, differences and unchanged bounds outside the production parser.

Only after the exact-kernel static/FD and generic charged controls pass does the
separate numerical harness run CG minimization and 40 steps of NVT at 300 K,
0.25 fs, fixed box. Minimization is limited to 100 iterations and lowers energy;
its recorded iteration-limit stop is **not convergence or equilibrium**.
The thermostat is a deterministic single Nosé–Hoover (`tchain 1`, 100 fs damping).
Initial velocities are drawn once with PCG64; complete before/after states and
velocities are preserved. The deterministic integrator consumes no RNG.

Uninterrupted 40 steps are compared to 20 steps plus 20 after filesystem
relocation and fresh Python/LAMMPS processes. Entire native restart files,
including thermostat state, all-atom coordinate/velocity/force dumps, and complete
RNG records compare byte-for-byte. A preliminary default three-chain thermostat
experiment showed last-bit restart drift; the bounded control explicitly declares
a single thermostat and makes no guarantee for other thermostats or stochastic
integrators. No rounding, canonicalization or ignored checkpoint bytes is used.

The original profile remains blocked because no replacement numerical contract
has been registered as its source-author execution convention. The distinct
control's success is not universal PPPM closure or pressure validation. Six
virial components are compared with OpenMM affine energy derivatives at two
strain sizes. Raw virials and configurational pressure proxies are retained;
there is no inherited scientific virial/pressure acceptance bound.

## Reproduction

No download or runtime asset discovery occurs. Use new output directories.
The exact executable produced by a local build must be pinned independently;
compiler/build paths can change its bytes.

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J49_PY=../island-validation/phase4g1-env/bin/python

# Explicit diagnostic build from the supplied local checkout. Never edits it.
"$J49_PY" scripts/build_phase_4j49_kernel.py \
  --source ../lammps/lammps --output /tmp/j49-exact --jobs 8
sha256sum /tmp/j49-exact/build/lmp

"$J49_PY" scripts/phase_4j49.py \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --exact-kernel /tmp/j49-exact/build/lmp \
  --exact-kernel-sha256 <independent-sha256-from-above> \
  --output /tmp/j49-audit
# Expected exit 1: original author-profile acceptance remains blocked.
# summary.json separately records bounded_diagnostic_acceptance_gate=0.

"$J49_PY" scripts/phase_4j49_virial.py --directory /tmp/j49-audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4
python scripts/audit_phase_4j49_inputs.py --directory /tmp/j49-audit \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --output /tmp/j49-independent-inputs.json

"$J49_PY" -m pytest -q tests/test_phase_4j49.py
"$J49_PY" -m pytest -q tests/test_phase_4j47.py tests/test_phase_4j48.py
"$J49_PY" -m pytest -q
ruff check src tests scripts
"$J49_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$J49_PY" scripts/validate_phase_4j49_receipt.py docs/evidence/phase_4j49 --repo .
"$J49_PY" scripts/validate_phase_4j48_receipt.py docs/evidence/phase_4j48 --repo . \
  --sha256 28bf6886f518977c0945eb1437f091ce0970cfb8c519a773a511b96fe8aa28eb
```

The reproducibility manifest binds exact source hashes, both executable hashes,
explicit numerical overrides, all new implementation/evidence files, exact full-
suite node accounting and every historical base-tree file. Receipt validation
and original-profile scientific acceptance are separate results.

Verification completed: **52 J49 focused tests passed**. The complete collected
suite contains **4107 tests: 4046 passed and 61 existing skips**; the separate
J47/J48 regression run passed **127 tests**. Ruff, compileall and all five pip
checks passed. Raw JUnit/log records and exact node accounting are retained.
Original author-profile acceptance remains nonzero despite these software checks.
