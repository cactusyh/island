# Phase 4J37: independent RESP/Psi4 provenance and GAFF2 reproducibility

J37 independently reviews the retained linear (60 atoms, two components) and
network (63 atoms, one component) quantum/RESP evidence before allowing fresh
GAFF2 parameterization and finite-run periodic workflows. It starts from C2
commit `4e0457e0dafcd2c6e1338f64e8143ebcf7a83236` on a separate branch.

`phase_4j37_reference.py` uses only the Python standard library. It imports
neither the production charge importer nor the finite-precision policy. Its
reviewed archive/input/result/worker/runtime anchors are in
`phase_4j37_anchors.json`. The auditor reconstructs chemical graph identity and
connected components from bonds, checks atom IDs/elements/formal totals,
initial centered coordinates, optimized native coordinates, Psi4 printed
geometry witnesses, raw artifact hashes, runtime identity, and unchanged
native charge bits against the imported vector. It independently reconstructs
exact rational residuals, compensated sums, Decimal reconstructions and C2's
half-ULP plus gamma-n bound. Descriptive policy changes are rejected too.
The source vector is not normalized, rounded or redistributed.

The initial-coordinate centering cross-check allows 5e-14 angstrom for different
summation orders. The printed Psi4 geometry witness allows 5.1e-12 bohr for its
12-decimal rendering. Neither allowance applies to charge values: the raw NPY
and imported JSON charges must be bit-identical. The existing downstream Amber
PRMTOP serialization allowance remains 1e-14 elementary charge. These are
mapping/serialization checks, not changes to simulation or charge tolerances.

The native replay invokes the pinned external CDS RESP implementation on the
retained native ESP/inverse-distance inputs, independently of ISLAND. It also
reconstructs the inverse-distance matrix (1e-13 bohr inverse allowance) and
requires the fitted vector to reproduce all retained float64 bits. It does not
repeat geometry optimization or Psi4 SCF. Historical optimized HF/STO-3G
geometries and HF/6-31G* ESP are explicitly retained inputs. Hashes establish
integrity against reviewed anchors, not authenticity of an arbitrary external
calculation or experimental validation.

Only after both reviews, relocation/child reconstruction and native refits pass
does the driver invoke fresh GAFF2 preparation and the existing NPT continuation
implementation. The frozen protocol minimizes, executes 8 NVT and 32 NPT steps
at 0.5 fs, 300 K and 1 bar, and splits at NPT step 16. All retained positions,
velocities, boxes, energies, temperatures, densities, stages, checkpoints and
complete RNG evidence are compared by the existing continuation contract.
Adaptive barostat state is reconstructed by independent prefix replay before
checkpoint loading. Measured pressure remains unavailable; target pressure is
never substituted for a measurement. These short runs verify software behavior
and provide four density blocks; they do not establish polymer equilibration.

Independent topology paths and native Amber single-point checks are retained.
J37 additionally compares periodic Ewald energy, forces, and bond/angle/torsion/
nonbonded components at first/last retained frames through OpenMM's own
AmberPrmtopFile reader against ISLAND's periodic construction. Both use OpenMM
Reference. This is an independent construction path, not an independent engine.
Existing limits are 1e-6 kJ/mol for energies/components and 1e-5 kJ/mol/angstrom
for forces. No LAMMPS KSPACE/Ewald/PPPM execution is claimed.

Run from the repository with declared environments intact:

```sh
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j37 \
  --output /tmp/j37-fresh-empty-output
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_phase_4j37.py
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j37_receipt \
  docs/evidence/phase_4j37/receipt.json
```

The output path must not exist. Diagnostic artifacts survive nonzero outcomes.
`--audit-only` is an inventory mode and always returns nonzero acceptance status;
it cannot claim molecular-workflow success. Missing/unverifiable raw assets,
missing pinned environments, failed refits or workflows fail acceptance.
Relocated provenance records are reconstructed offline using the bundle plus
declared reviewed anchors; the runtime need not execute for offline loading.

RESP/Psi4, native AM1-BCC and external OPLS-AA remain separate origins. Native
AM1-BCC fresh acceptance stays blocked by its historical precision diagnostic.
No universal GAFF2/PCFF coverage, experimental agreement, thermodynamic
equilibrium, Tg, mechanical properties or scientific production readiness is
claimed. All historical evidence is preserved byte-for-byte.

`production_validated=False`, `simulation_readiness="not_established"`.

## Retained-input audit and fresh execution results

The tested implementation is
`a41de9bd4cc3c17b207a7d871433e8f2808a0571`. Both external native RESP fits
reproduce the authoritative charge vector bit-for-bit. Linear component absolute
residuals are `4.423544863740858e-17 e` and `8.066464163292153e-17 e`;
the system residual is `1.249000902703301e-16 e`. The connected network has
one component with absolute residual `9.985502008591496e-17 e`. These are
nonzero source-value residuals accepted under the independently reconstructed,
unchanged C2 numerical policy; no atom charge is repaired.

Both fresh GAFF2 workflows complete, show actual NPT box-volume change, retain
four density blocks, and pass split/uninterrupted comparison with zero numerical
tolerances and byte-exact checkpoint/RNG hashes. Independent periodic reference
checks pass for both fixtures. Maximum total-energy discrepancies are
`7.0273813435051125e-9 kJ/mol` (linear) and
`1.6464241525682155e-8 kJ/mol` (network); maximum force discrepancies are
`8.04098476692161e-9` and `7.62846497082137e-9 kJ/mol/angstrom`, respectively.
Bond, angle, torsion and nonbonded component comparisons pass their existing
limits. These values are fresh software cross-checks, not thermodynamic results.

The independent live-runtime probe verifies all 1,462 pinned quantum runtime
files. Five Python environments are recorded with package versions and Python
executable fingerprints: supported tests, Foyer/OpenMM, MDAnalysis, AmberTools,
and Psi4/RESP. The native environment pins Psi4 1.9.1 and CDS RESP conda package
1.0 (module version 0.8); these version labels are kept distinct. All five pip
checks pass. The historical inventory contains 1,352 files and matches the base
Git blobs byte-for-byte. Local `main` stays at
`3de000e05e0af3071d1a23a8d6abc807f57909ac`.

Exact commands, JUnit reports, per-module counts, unchanged optional/archive skip
reasons, manifests, raw provenance reviews, refit results and workflow bundles
are retained in [the J37 evidence directory](evidence/phase_4j37/receipt.json).
Preliminary runs began before receipt hardening and are excluded from final
acceptance; the complete final gates were rerun against the tested SHA above.
Historical exact-equality RESP and native AM1-BCC validators retain their
expected nonzero acceptance outcomes. The original J35 external-OPLS control
is checked using its unchanged historical source checkout rather than rewriting
its source-hash-bound receipt.

Final verification: **33 focused tests**, **696 J22–J36 regression tests**, and
**2,697 full-suite passes with 10 established skips**. No test failures or errors
occurred in the final gates. Ruff, compileall, all five pip checks, diff checks,
historical validators and the new receipt/manifest validator pass. Seven fresh
charge/provenance CLI mutations and seven consistently re-signed full-receipt
mutations fail as intended and preserve their diagnostic records. These are
explicitly archived software negative controls, not scientific evidence.
The frozen implementation SHA above is distinct from the later evidence/docs
commit; resolve the latter with `git log -1 --format=%H --
docs/evidence/phase_4j37/receipt.json`.
