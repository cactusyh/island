# Phase 4J47: external PDMS static LAMMPS execution and strict acceptance

J47 connects the existing source-bound PDMS assignment to an actual pinned
LAMMPS process. HMDS, trimethylsilyl-capped DP2, DP4 and the supported longer DP6
produce deterministic data/input, checksummed static results and reconstructed
workflows after relocation into a fresh Python process.

**The end-to-end MD acceptance remains blocked.** The author-authored PPPM
settings fail the inherited energy and force thresholds. The tighter historical
numerical control explains the remaining erfc/constant error but does not satisfy
exact energy parity. No LAMMPS minimization, NVT, NPT, MD trajectory, restart or RNG
continuation is authorized by J47. Each required unexecuted stage is explicit in
the case report. The full acceptance command returns **1** for this classified
blocked outcome, without a Python traceback.

`production_validated=False`, `simulation_readiness="not_established"`.
Polymer Ewald/PME/PPPM scientific closure remains unvalidated, backend pressure
unavailable, cutoff recovery blocked and long NPT blocked. There is no claim of
production readiness, equilibrium, experimental agreement, universal PDMS or
universal msi2lmp coverage.

## Exact base, profile and source

After `git fetch origin`, the exact merged J46 `origin/main` was
`9e7d9765d0248e14320b4214db8831771f06e84e`. This branch starts there. Its
prerequisite profile from snapshot
`8257a2ef25042494d0e4541b8da117cd8e97072a` remains unchanged:

`island_external_opls_pdms_xiang2025_dft_local_v1`.

It is an **external, literature-derived OPLS-compatible local research model**.
It supplies neither native universal OPLS-AA nor PCFF/COMPASS PDMS parameters.
The public preparation consumes a final graph built from SMILES/PSMILES. No
CAR/MDF input, converter or network parameter discovery is needed.

The inspected paper is Zhirui Xiang, Chao Gao, Teng Long, Lifeng Ding, Tianhang
Zhou and Zhenghao Wu, *Force-Field Benchmark for Polydimethylsiloxane: Density,
Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*,
J. Phys. Chem. B **129** (2025), 1864–1873,
DOI **10.1021/acs.jpcb.4c08471**. Its methods distinguish Qeq from DFT charges;
the Results use Qeq. This existing profile uses the separate author DFT/NPA
serialized input. J47 does not transfer the paper's thermophysical claims to
that branch. The local PDF was read with offline `pypdf`; no parameter was
extracted anew or replaced.

| Explicit asset | SHA-256 |
| --- | --- |
| `island_ffref/pcff/PDMS_benchmark-main.zip` | `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4` |
| `island_ffref/pcff/jp4c08471.pdf` | `04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc` |
| `../island-validation/phase4j22/lmp-kspace` | `6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191` |
| `../lammps/lammps/src/KSPACE/pair_lj_charmm_coul_long.cpp` | `8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5` |
| `../lammps/lammps/src/ewald_const.h` | `3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c` |

All five existing archive companion hashes, 17 source rows and original raw-line
citations remain in `source-preparation.json`. Execution records capture the
10 Sep 2025 development build/revision, complete help/package/compile flags,
compiler, FFT double precision, platform, units and executable hash. The executable
path in help is replaced by `<pinned-lammps>` so relocation cannot change build
identity. The hash still binds the actual executable bytes. The source archive
and PDF are user-supplied research assets with redistribution permission not
established; this branch does not redistribute them.

## Public API and integrity boundary

The additive `island.forcefields.pdms_lammps` namespace preserves all historical
preparation APIs and records. It accepts the requested named profile explicitly:

```python
import numpy as np
from island.builders import build_linear_polymer
from island.core import SimulationBox
from island.forcefields.pdms_lammps import (
    ForceFieldRequest, OPLSOptions, prepare_forcefield,
)
from island.pdms_lammps import LAMMPSOptions, run_static, reconstruct, run_md

system = build_linear_polymer(
    "[*:1]O[Si](C)(C)[*:2]", dp=4,
    head_end_group="[*][Si](C)(C)C", tail_end_group="[*]C",
    add_hydrogens=True, random_seed=2026,
)
system.box = SimulationBox(40, 40, 40)
xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
system.coordinates.translate(20 - xyz.mean(axis=0))
options = OPLSOptions(
    profile="island_external_opls_pdms_xiang2025_dft_local_v1",
    archive_path="island_ffref/pcff/PDMS_benchmark-main.zip",
    expected_sha256="8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4",
    e_oxygen_site_ids=(system.metadata["polymer"]["head_site_id"],),
)
request = ForceFieldRequest("opls", options)
prepared = prepare_forcefield(system, request)
report = prepared.payload
# report includes profile/source/graph identities, exact ordered charge literals,
# charge origin, assigned/missing interactions, conventions, backend availability,
# separate stage flags and scientific limitations.
engine = LAMMPSOptions(
    "../island-validation/phase4j22/lmp-kspace",
    "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191",
)
result = run_static(system, prepared, "/tmp/pdms-static", engine=engine)
restored_system, restored_prepared, observed = reconstruct(
    "/tmp/pdms-static", manifest_sha256=result["manifest_sha256"],
    request=request, engine=engine,
)
# Raises a classified PDMSOPLSError before engine invocation or output creation:
# run_md(system, prepared, "/tmp/blocked-md", engine=engine, steps=100)
```

`prepared_static_only` establishes a complete existing parameter assignment for
static export. It does not establish numerical or workflow MD acceptance. A
blocked source/typing/convention record carries its full diagnostic. Integrity
recomputes preparation against the current exact graph, order, component IDs,
source and charge literals, and additionally binds coordinates and box.

The workflow stores the system, immutable source preparation, public assignment,
engine build record, deterministic data/coefficients/input, raw process logs,
force/coordinate dump, static report and checksum manifest. Reconstruction
requires an independent manifest SHA-256 plus explicit source and executable
paths/hashes. It verifies membership and safe paths, rehashes everything,
reprepares from the actual source, regenerates input bytes and re-parses raw
output. Re-signing edited JSON cannot grant execution authority. Static evidence
contains no checkpoint or trajectory. The reported pressure is a configurational
virial diagnostic only.

## Discrepancy diagnosis

Both engines receive identical coordinates, atom order, types, masses, charge
literals, box, graph-distance exclusions, LJ/Coulomb 1–4 scale 0.5, geometric
sigma/epsilon mixing, CHARMM r-squared 9–11 Å LJ switching and 11 Å real-space
Coulomb cutoff. No charge normalization/redistribution, parameter adjustment,
zero filling, electrostatics disabling or tolerance change occurs.

The current profile declares LAMMPS **PPPM `1e-4`**; the existing OpenMM research
implementation declares **Ewald `1e-7`**. These numerical settings are explicitly
unequal. J47 does not pretend that equal coordinates establish solver equality.
A source-bound contract for identical reciprocal solver/splitting/mesh treatment
between these two backends is missing. The author run uses the pinned engine's
default pair table and no overrides. The separate historical control retains
`table 0`, PPPM `1e-10`, order 7, solely to diagnose algorithmic error.

| Fixture | Atoms | Author-settings total energy error, kJ/mol | Maximum force error, kJ/(mol Å) |
| --- | ---: | ---: | ---: |
| HMDS | 27 | 0.0699658923 | 0.0839530434 |
| Capped DP2 | 37 | 0.1436606402 | 0.2561358095 |
| Capped DP4 | 57 | 0.4018122426 | 0.4421014032 |
| Capped DP6 | 77 | 0.5085865935 | 0.3382671888 |

The inherited energy and force bounds are both `3e-4` in the corresponding units.
Bond, angle, OPLS torsion and LJ agree near floating-point precision. Native
Coulomb and total energies, forces and source-settings LAMMPS finite differences
fail. Cause: **long-range solver**; the diagnostic control further identifies
**precision** (erfc polynomial) and **unit conversion** (Coulomb constant).
The independent exact-equation reconstruction and the unchanged control retain
the raw error; no corrected engine output is used for acceptance. The control's
remaining exact energy errors are approximately 0.00154/0.00235/0.00349 kJ/mol
for HMDS/DP2/DP4, still above `3e-4`.

Real-space and reciprocal/self/background energies are recorded separately.
Their partitions depend on the different Ewald splitting parameters and are
classified as non-comparable partitions; only the total Coulomb sum is invariant.
LJ14 is measured separately in OpenMM and independently reconstructed; Coulomb14
and special-pair IDs are reconstructed explicitly. LAMMPS `ecoul` already includes
exclusion corrections and is never mislabeled as a standalone Coulomb14 energy.
Improper interactions are absent in this declared domain; no zero improper row
is introduced.

Global LAMMPS configurational virial is compared with six affine-box/coordinate
energy derivatives from OpenMM at two strain steps and an independent exact-Ewald
isotropic trace. Virial remains **unresolved/unqualified**: there is no inherited
acceptance threshold or identical reciprocal solver contract. No new tolerance
is invented. Pressure proxy and NPT remain unqualified.

Each fixture also has strict OpenMM finite differences, source-row independent
energy/force reconstruction and actual source-settings LAMMPS finite differences
at `1e-4` and `5e-5` Å. The latter failures are retained. The downstream
minimization, short NVT, fresh-process MD continuation and exact checkpoint/full
RNG comparison remain blocked behind these static gates. Historical OpenMM
minimization, 342 DP4 checks and 1000-step research NVT evidence are byte-preserved,
without being substituted for LAMMPS acceptance.

## Reproduction and verification

Run from the repository root with the supplied assets present. No download occurs.
Use a new output directory on each audit invocation:

```bash
export PYTHONPATH=.:src
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PDMS_PY=../island-validation/phase4g1-env/bin/python

"$PDMS_PY" scripts/phase_4j47.py \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --lammps-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --lammps-source ../lammps/lammps/src \
  --paper island_ffref/pcff/jp4c08471.pdf \
  --paper-sha256 04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc \
  --output /tmp/island-j47-reproduction
# Expected exit: 1; summary.json documents the classified blocked matrix.

# Independent standard-library reconstruction, outside all production parsers:
python scripts/audit_phase_4j47_inputs.py \
  --directory /tmp/island-j47-reproduction \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --output /tmp/j47-independent-inputs.json

"$PDMS_PY" -m pytest -q tests/test_phase_4j47.py
"$PDMS_PY" -m pytest -q
ruff check src tests scripts
"$PDMS_PY" -m compileall -q src tests scripts
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
../island-validation/phase4j32-analysis-env/bin/python -m pip check
../island-validation/phase4j36-native-env/bin/python -m pip check
../miniforge3/envs/md_rg/bin/python -m pip check
git diff --check
"$PDMS_PY" scripts/validate_phase_4j47_receipt.py \
  docs/evidence/phase_4j47 --repo .
"$PDMS_PY" scripts/validate_phase_4j46_receipt.py \
  docs/evidence/phase_4j46 --repo . \
  --sha256 136226931e81fe5be450ace9617548f185ba164f151ae8f9c5bc2307f1a45776
"$PDMS_PY" scripts/validate_pdms_public.py \
  docs/evidence/pdms_public_v1 --repo . \
  --sha256 4cad00e0639404c0bf32c35d0bb2dbfd9bc032a5b141a51bc786d1844dcd67c3
```

The receipt/manifest validator's success means the evidence, blocked gates,
complete suite coverage and historical bytes are intact. It does not turn the
nonzero end-to-end acceptance into a success. The manifest protects every tracked
entry from the exact merged J46 base, including all J17–J46 receipts, bundles,
identities and symlink target bytes. Only additive J47 files are committed.

Focused tests cover exact profile/charge/source reporting, deterministic export,
actual LAMMPS static execution, negative source/ZIP/parameter/charge/order/box/
1–4/mixing/cutoff/reciprocal/profile-family controls, forged/re-signed workflow
records, pre-engine MD rejection, executable relocation and fresh-process
reconstruction. Full-suite JUnit coverage and all check outcomes are recorded in
`docs/evidence/phase_4j47/receipt.json`.
