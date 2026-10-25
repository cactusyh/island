# Phase 4J41 — PCFF silicon/phosphorus source tranche

J41 qualifies **finite neutral silane only** through an explicit opt-in source
profile. Methylsilane and a mixed C/H/O/Si siloxane fail complete Class-II source
coverage. Neutral phosphate ester and phosphonate remain blocked at chemical
source-type interpretation. The tranche acceptance gate is **1**, not success.
`production_validated=False`; `simulation_readiness="not_established"`.

The base is merged C14, `bed595291d52ce15c44ee59dcce6b681b6a8e91f`.
The tested implementation is `c4faced420a064d07c90ff6f15a9d1ea6d9a854d`, on
`codex/phase-4j41-pcff-silicon-phosphorus`. Historical implementation behavior
and evidence are unchanged. The receipt and manifest bind the exact commands,
source files, numerical artifacts, environments and historical byte hashes.

## Source inventory and interpretation

The independent ASCII scanner preserves original line numbers, versions,
namespace, ordering and text. It inventories every numeric source section,
selects direct Si/P labels and described companion labels, retains their
ordinary/automatic equivalence tokens, and includes all wildcard candidates.
Generic equivalence tokens intentionally overinclude rows; catalog membership
is never chemical or execution authorization. Empty/unsupported families are
not interpreted as zero interactions. The original FRC remains in the bundle.

Primary element labels are `si`, `sio`, `sz`, `p`, `p=`. Described companions
include `hsi`, `hp`, `hos`, `osi`, `osh`, `oss`, and `oas`. Inventory includes
atom definitions, both equivalence sections, increments, all ordinary bonded
families, 9–6 nonbonded parameters, bond–bond, bond–angle, angle–angle,
bond–bond-1–3, end/middle-bond–torsion, angle–torsion and angle–angle–torsion.
The empty torsion–torsion section remains visible in the preserved source;
no torsion–torsion support is inferred. Charge records are bond increments,
not a separately supplied atomic charge vector.

The independent audit compares **1,942** retained raw records against the native
catalog and executes **16** Si/P-focused queries against the compiled pinned C
leaf matcher. A separate retained matcher audit tests wildcard ambiguity and
ordinary equivalences. Automatic-equivalence queries are diagnostic tests of
leaf matching, not evidence that the complete converter supports that route.
Native guarded selection rejects ambiguous candidates rather than endorsing
file-order selection. No historical source resolver is changed.

Distribution SHA-256:

| File | SHA-256 |
|---|---|
| `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `pcff.rlb` | `82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36` |
| `pcff_templates.dat` | `dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093` |

The FRC bytes match ISLAND's existing source pin. The actual local LAMMPS
checkout is `36f7c81ced8dd456430b13647e48f0c7d7bd8faa`; this is recorded separately
from the historical FRC/implementation pin. All compiled converter C/header
hashes, compiler version, flags and source text are retained.

## Bounded profiles and failure stages

All fixtures use explicit hydrogens, seed 2026, finite coordinates and formal
molecular charge zero. A declared `0.015*sin(arange(3*N)+0.3)` Å perturbation
makes the numerical checks nontrivial. Stable atom IDs, canonical sorted-ID
ordering, coordinates, component membership, graph, source, typing, charge,
parameter and model identities are persisted. Identity validation reconstructs
the assessment against the exact supplied system and pinned source; neither
provided charges nor RESP/Psi4 charges enter these profiles.

| Profile suffix (`island_pcff_…_v1`) | Fixture / SMILES | Atoms | Result |
|---|---|---:|---|
| `silane_source` | `[SiH4]` | 5 | Source-complete; finite numerical/minimization/continuation checks pass |
| `organosilicon_source` | `C[SiH3]` | 8 | Source-incomplete: 9 bond–bond-1–3 and 9 end-bond–torsion requests missing |
| `siloxane_source` | `C[Si](C)(C)O[Si](C)(C)C` | 27 | Source-incomplete: 378 coupling requests missing |
| `phosphate_ester_audit` | `COP(=O)(OC)OC` | 17 | Chemical source-type interpretation unresolved; no charges/model published |
| `phosphonate_audit` | `CP(=O)(OC)OC` | 16 | Chemical source-type interpretation unresolved; no charges/model published |

For the siloxane, the 378 requests are 78 angle–angle requests and 60 each of
bond–bond-1–3, end-bond–torsion, middle-bond–torsion, angle–torsion and
angle–angle–torsion. Exact sites, supplied type tuples, equivalence paths,
searched source candidates and missing dependencies are retained for every
interaction. Native increment charges are component-neutral for both mixed-Si
fixtures, but neutrality cannot compensate for missing coupling rows.

The independently compiled converter exits **17** without `-ignore` for both:
`Unable to find endbonddihedral data for hc c3 si hsi` (methylsilane) and
`… hc c3 sio c3` (siloxane). J41's strict guarded-source policy also rejects
missing bond–bond-1–3 rows; the historical non-`cp` converter-zero convention is
not used to satisfy J41 completeness. Explicit numerical zeros in real source
rows remain distinguishable from missing parameters.

FRC line 163 describes `p` as general phosphorus with four connections; line
164 describes `p=` as five-connected phosphazene phosphorus. The companion
`pcff_templates.dat` line 1562 also uses `p=` for a four-neighbor oxo template.
The RLB/templates are preserved as source evidence, not executed as a typing
engine. There is no independently established exact chemical mapping resolving
this conflict for the retained oxo-P fixtures. Neither label is guessed;
converter assembly is not attempted with invented P/O types. This does not
assert that the distribution contains no phosphorus parameters. It records why
those rows do not yet authorize these exact final graphs.

Unsupported silicon valence, phosphorus oxidation/ambiguous states, charged,
ionic, metallic, radical and isotopic domains remain explicitly blocked. There
is no unrelated atom-type fallback or inferred missing zero. Incomplete models
cannot be returned by `prepare_pcff_sip`.

## Executed finite Class-II checks

`msi2lmp` is compiled from the local C sources using `gcc -O2 -Wall -Wextra …
-lm`; its observed version is **3.9.11 / 6 Sep 2024**. Its SHA-256 is
`59f3718d4e537d584011315b76d0bf80ca105b43468f87b384cfa4ce19d76123`.
It consumes independently supplied atom types and echoes CAR/MDF charges;
**it performs neither chemical perception nor native charge assignment**.

The LAMMPS executable SHA-256 is
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
J41 executes `bond_style class2`, `angle_style class2`, `dihedral_style class2`,
`improper_style class2`, `pair_style lj/class2/coul/cut`, sixth-power mixing,
and `special_bonds lj/coul 0 0 1`. These are finite, nonperiodic evaluations:
the declared reference cutoff exceeds every pair separation. The available
binary has KSPACE, but J41 does not execute or qualify PCFF Ewald/PPPM.
OpenMM is **8.6.1.dev-b399af4**, Reference platform, double precision.

For silane, independent graph walks verify 4 bonds, 6 angles, no proper
torsions and 4 neighbor triples. The model retains all 12 applicable
angle–angle couplings. All executable bond/angle/coupling coefficients are
compared with the converter data and independent raw-row reconstruction.
The converter's degree-four zero Wilson placeholders are not mistaken for
missing source terms: Wilson bending is inapplicable there, while nonzero
angle–angle terms are retained. Existing AA equilibrium-role auditing is used;
all silane angle equilibria coincide, so the emitted values do not change.

The source increments give Si `+0.504 e`, each H `−0.126 e`, and exact decimal
component total zero; float residual also zero, tolerance `1e-12 e`.
No individual charge is modified. Input coordinates are restored at full
precision for evaluation after converter text formatting/recentering; the
original CAR, MDF, converter data, full-precision coordinates and mappings are
all retained.

Initial maximum component/total energy difference is **2.40e-14 kJ/mol**;
maximum force difference is **1.11e-12 kJ/(mol Å)**. Unchanged bounds are
`atol=1e-5`, `rtol=2e-10`. Central finite differences cover every Cartesian
coordinate and compiled component at `1e-3`, `1e-4`, `1e-5 Å`; the final maximum
error is **1.59e-7 kJ/(mol Å)** against the `1e-5` bound. Absent components are
reported as absent/zero, not missing-row repairs. There are no constraints,
periodic images or neighbor-list discontinuities in this finite test.

Minimization passes the existing `fmax <= 0.1 kJ/(mol Å)` criterion, reaching
`0.013605`. NVE and Langevin BAOAB each run 200 steps at 0.1 fs (**20 fs**), with
100+100 split continuations. Retained frames are compared against independently
executed LAMMPS at steps 0, 100 and 200. Relocation and fresh-process
reconstruction pass. Final coordinates, velocities and complete RNG states
match the uninterrupted trajectory exactly. Checkpoint bytes match between
in-process continuation and relocated fresh-process continuation with identical
segmented lineage. Whole-run and segmented checkpoints intentionally retain
different lineage/evaluation counters; their entire files are not claimed
byte-identical. No serialization tolerance is substituted for the exact
checkpoint comparison that is reported.

Silane has no active nonbonded pairs under its shortest-path exclusions. Its
LJ/Coulomb source coefficients are audited, but this fixture does not validate
intermolecular Si nonbonded physics. Twenty femtoseconds demonstrates software
execution and continuation only; it establishes no sampling, equilibrium,
pressure, density or experimental agreement. No incomplete Si/P fixture runs MD.

## Reproduction and verification

The additive APIs are `assess_pcff_sip(system, source, profile=...)`,
`prepare_pcff_sip(system, assessment)` and
`load_pcff_sip_assessment(path, source, system=system)`. They do not change the
default PCFF profile registry or grant general polymer/periodic preparation.

The executable study command is retained verbatim in `execution-command.json`:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --converter /tmp/j41-validation/converter/msi2lmp \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j41-validation/final-study
```

The expected exit is **1**, with a complete report recording four source-blocked
profiles. Unexpected numerical or runtime failures are not accepted as source
rejections. `execution.tar.xz` contains the raw distribution, systems,
assessments, converter source/binary/build log, matcher, native model,
reference outputs, finite-difference sweeps, minima, trajectories, checkpoints,
environment and independent source audit. `validation.tar.xz` contains the
frozen implementation, full commands, JUnit reports, logs, package inventories
and historical hashes. Extracting the execution archive to a new directory
allows construction-free verification:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41.py \
  --reconstruct /path/to/extracted/study
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_receipt.py docs/evidence/phase_4j41
```

Final verification: **26 focused tests**, **1,156 J22–J40 regressions**, and
**3,150 full-suite passes with 10 optional skips**. Ruff, compileall and all five
pip checks passed. The receipt/manifest validator and extracted-bundle
fresh-process reconstruction passed; all **2,498** historical documentation
and evidence byte hashes are preserved. Exact command identities and JUnit
results are retained. Mutation controls cover valence/chemical state, missing
couplings, source hash, atom
ordering, coordinates, graph, charge values, component totals, special-pair
weights, source rows and incomplete interaction inventories. No universal PCFF,
COMPASS, msi2lmp, Si/P, charged-material or long-range parity is claimed.
AM1-BCC support and scientific production readiness are not claimed.
