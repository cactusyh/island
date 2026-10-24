# Phase 4J40-C1: corrective equilibration and cross-engine protocol

**Scientific acceptance is blocked. No new equilibration campaign was executed.**
This phase adds a fail-closed staged protocol, independently parses actual
LAMMPS inputs, reproduces the J40 failures, and audits network forces component
by component. It does not repair Ewald discrepancies or provide a validated
instantaneous pressure estimator. Longer integration cannot satisfy those two
prerequisites. A blocked result is deliberately nonzero.

- Base: `codex/phase-4j40-polymer-equilibration-cross-engine`, fetched and verified
  at `a7dfb9a5aaf8ef37d94604e7771cb2771572d82c`.
- Tested implementation: `fd591b7eff8c6f1a4ad0554b8da705b6e0706dfe`.
- Final evidence commit: the Git commit introducing the C1 receipt; kept distinct
  from the implementation commit to avoid a self-referential receipt hash.
- Historical evidence is read-only. No old source, numerical tolerance, force
  field, parameter file, serialized identity or receipt was changed.
- `production_validated=False`, `simulation_readiness="not_established"`.

## Reproduced failures

Before implementation, the retained J40 bundle was validated offline and its
failures recorded in a fresh `/tmp/j40-c1-validation/before.json`. A second,
independent scalar reduction of the raw frame files appears in
`independent-baseline.json`. It computes density from actual diagonal box volume
and total mass, uses the original quarter-history burn-in, and never substitutes
target temperature or pressure for observations.

| Cell | Replica temperatures (K) | Density ESS | Energy discrepancy (kJ/mol) | Force discrepancy (kJ/mol/Å) |
|---|---|---|---|---|
| GAFF2 PE, RESP/Psi4 | 160.737, 160.626 | 10.028, 3.242 | 0.042948679 | 0.041929003 |
| OPLS-compatible aromatic, RESP/Psi4 | 144.150, 164.702 | 4.590, 3.121 | 0.007970134 | 0.025599057 |
| Native OPLS PEO | 156.980, 168.018 | 4.063, 3.217 | 0.062370852 | 0.054798189 |
| Native OPLS cyclic network | 147.791, 152.372 | 2.942, 7.805 | 0.128373826 | 0.036190044 |

Every replica has only 192 retained NPT observations after burn-in, below the
unchanged 512-observation floor. Temperature means fail the 5% target bound.
Pressure observations are `null`. All four original Ewald comparisons fail the
unchanged J22/J40 absolute energy and force bounds. The original first/last atom
network finite difference fails at `h=1e-6 nm`, with a maximum error about
`0.041686684 kJ/(mol Å)`.

The fresh single points use the identical retained KSPACE executable and
regenerated input files. They are new executions, not copied numerical results.
They preserve the discrepancies; see the full precision results and LAMMPS logs
in the C1 archive. Bond, angle and LJ component agreement do not excuse the
failed electrostatic component or total/force gates.

## Explicit protocol and stopping policy

`CorrectiveProtocol` is immutable and versioned. Its scheduler calls a real
engine adapter only after independent reference, finite-difference and reliable
pressure prerequisites pass. There is no default pressure estimator or synthetic
trajectory callback.

Declared, **unexecuted** protocol:

1. Existing J23 force-checked minimization of the actual packed graph.
2. At 0.5 fs, heat for 5,000 steps at 100 K, then to step 10,000 at 200 K,
   then to step 20,000 at 300 K.
3. Continue NVT to step 40,000 (10 ps constant-temperature pre-equilibration).
4. NPT at 300 K, 1 bar, isotropic Monte Carlo barostat every 25 steps;
   evaluate cumulative NPT histories every 200,000 steps (100 ps).
5. Stop only when every diagnostic passes, or after ten windows (1 ns NPT).
   Exhausting the budget preserves insufficient sampling; it does not prove
   equilibrium. Three independently seeded replicas are required.
6. Retain frames/checkpoints every 100 steps. Burn the first quarter of each
   cumulative NPT history, with the boundary recorded in frames. No favorable
   tail is selected after inspecting results.

The unchanged J40 statistical policy requires four contiguous blocks, at least
512 retained frames, ESS at least 50, block length at least five integrated ACF
times, relative drift at most 0.02, measured temperature within 5% of target,
and pairwise replica agreement within twice the combined standard error.
Standard error is the larger of block and ACF estimates. Constant series cannot
establish independence. Diagnostic records include means, SE, ACF, ESS, rolling
means, half-history differences and burn-in indices. The scheduler validates
stage schedule, timestep, targets, platform, thermostat, barostat, coupling,
retention and cumulative step counts. It does not run in this acceptance because
the prerequisite gates remain blocked.

## Pressure limitation

The existing OpenMM state capture records actual kinetic temperature, potential
and kinetic energy, positions, velocities and box vectors. It does not expose
instantaneous pressure. A barostat's pressure parameter is a target. Computing
only `sum(r·F)` would omit important periodic/reciprocal/exclusion contributions;
no such estimator is published here. A reliable, independently checked periodic
virial/pressure implementation is still required. Consequently every cell has
`pressure_unavailable`, with no instantaneous or block-mean pressure invented.

## Component finite differences and geometry

The new audit evaluates every Cartesian degree of freedom of all 63 network
atoms, with independent central energy differences. It separates bonds, angles,
RB torsions, Coulomb and LJ plus scaled LJ 1–4. No improper terms occur in the
validated source assignment; that category is explicitly not present.

Displacements (nm): `1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9`.
The bound remains `1e-5 kJ/(mol Å)`; two adjacent passing sizes are required for
each component. The audit records worst atom ID, axis, analytical/numerical
force, sign convention, native engine version, constraints, exclusion count,
cutoff margin, periodic-image boundary margin and minimum angle sine.

On **packed** coordinates, bonds, Coulomb and LJ have converged intervals;
angles and RB torsions fail the interval criterion. The full-atom maximum
errors at `h=1e-6 nm` are approximately 0.00413019 for angles and 0.32491992
for torsions, larger than the old two-atom probe. The minimum angle sine is
0.00333814. The minimum cutoff margin is `5.54282e-5 nm`, so the largest
displacement can cross the unswitched cutoff. The minimum image-boundary margin
is `2.16131e-4 nm`; no constraint is present.

On the **retained minimized** geometry, every applicable component has a passing
interval at the same bound; the minimum angle sine is 0.89447349. This supports
geometry conditioning/truncation and cancellation as the finite-difference
issue, rather than demonstrating a force-sign or unit error. The packed network
failure remains a failure. Neither geometry was silently substituted for the
other, and the minimized audit is not a new MD run.

## Engine input parity and remaining discrepancy

The independently parsed actual LAMMPS data verifies atom ordering, coordinates,
box, masses, charges, LJ coefficients, interaction memberships, units, cutoff,
switch/tail settings, 1–4 factors, mixing and Ewald request. Data and input files
must also match the original hash exactly, binding all remaining coefficients.
Only the 17-digit unit conversion of LJ epsilon allows a 2e-15 relative
serialization error; this does not change any energy or force acceptance bound.

The actual KSPACE executable hash is
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
Its declared build source is LAMMPS commit
`36f7c81ced8dd456430b13647e48f0c7d7bd8faa`.
Ewald was executed; PPPM was not. The common requested Ewald accuracy of `1e-5`
is an algorithmic accuracy request, not a guarantee of absolute cross-engine
agreement at `1e-5 kJ/mol` or `kJ/(mol Å)`. J40's tighter-accuracy single points
also failed some bounds. Matching input semantics is therefore distinct from
numerical reference agreement. Cross-engine forces also retain the original relative tolerance of `2e-10`. This phase preserves the electrostatic
mismatch rather than rescaling constants or loosening tolerances.

## Continuation and reconstruction

The historical uninterrupted/split histories remain exact. C1 validates their
nested graph/source/assignment/frame/checkpoint/RNG identities, copies the full
bundle to a different directory, and validates it in a fresh interpreter.
These are explicitly historical reconstruction checks, not new integrations.
There is no new equilibration trajectory or new continuation claim.

Native aromatic OPLS remains source-incompatible. PCFF gaps and native AM1-BCC
precision limitations remain unchanged. RESP/Psi4 and native OPLS charges retain
their distinct origins; no external charge vector was substituted.

## Reproduction and verification

Use the five declared Python environments inventoried in `environments.json`.
Set `PYTHONPATH=.:src`, `OPENMM_CPU_THREADS=1`, `OMP_NUM_THREADS=1`.
From a fresh directory, extract the historical archive using automatic XZ
recognition, retaining its `run/` and `accuracy-study/` subdirectories:

```sh
../island-validation/phase4e2-env/bin/python -m tarfile -e docs/evidence/phase_4j40/artifacts.tar.xz /tmp/j40-c1-input
../island-validation/phase4e2-env/bin/python -m scripts.validate_phase_4j40_c1 \
  --source /tmp/j40-c1-input/run --output /tmp/j40-c1-new \
  --lammps ../island-validation/phase4j22/lmp-kspace
# Expected exit 1: scientific acceptance is blocked.
../island-validation/phase4e2-env/bin/python -m scripts.validate_phase_4j40_c1_receipt \
  docs/evidence/phase_4j40_c1/receipt.json --diagnostic-only
# Expected exit 0 for valid blocked evidence; omit flag for acceptance exit 1.
```

Exact commands, implementation SHA, environment variables, logs, JUnit counts,
source/executable/artifact hashes and historical byte hashes are in the additive
receipt and validation archive. The supported suite contains short integration
regressions; it must not be described as execution of the proposed longer
protocol. Missing physical sampling, a reliable pressure estimator and failed
cross-engine agreement remain explicit scientific limitations.

No thermodynamic equilibrium, experimental agreement, Tg, mechanical properties,
universal force-field coverage, AM1-BCC support or production readiness is claimed.

An initial history helper treated the J40 archive symlink as a regular Git blob,
which caused a false hash alarm. The corrected helper verifies both the link
text against Git and the linked archive bytes. All 2,408 historical entries are
unchanged. Preliminary and interrupted test executions are retained separately
and are excluded from passing validation counts.

Final validation against the tested implementation:

| Gate | Result |
|---|---|
| J40-C1 focused | 37 passed |
| J22–J40 regressions (including C1) | 996 passed |
| Full supported suite | 2,964 passed, 10 skipped |
| Ruff, compileall, five environment pip checks | Passed |
| Historical entries, including archive symlink | 2,408 unchanged |
| Historical receipt integrity and expected blocked-gate semantics | Preserved |

The ten full-suite skips are recorded individually in JUnit. They concern
explicitly opt-in AmberTools QM/reference, analysis-reference and OPLS source
integration gates; they are not counted as executed scientific acceptance.
Earlier interrupted suites and preliminary results are archived separately.
The final receipt validator distinguishes integrity success (diagnostic mode)
from the intentionally nonzero scientific acceptance gate.
