# Phase 4J39 — aromatic OPLS charge-source acceptance

The aromatic benchmark now has an executable **OPLS-compatible parameters with
RESP/Psi4 charges** path. This is option C, not native OPLS-AA charge completion.
The pinned native OPLS charge vector remains **source_incompatible**, with
−0.230 e per neutral chain. Native acceptance returns nonzero. No source charge
was changed, normalized, redistributed, or selectively mixed with RESP charges.

Branch: `codex/phase-4j39-authoritative-aromatic-opls-charges`.

Base: fetched `origin/main` at
`78a46d46bd66921f527fc08082903d151f0e4f5e`, the merged J38-C3 implementation.
Its tree (`a651b52d2f49219f646186aaf8f260ddc61c8e5b`) equals the reviewed C3 tip
`b129d79e6ff51dee5df892aa07ef8fe463c47a9b`. Implementation actually tested:
`1f23a4a8548a959b6f09038a19bf2d6ce5f7c49d`. The later evidence commit contains
only documentation and additive J39 artifacts. Main was not modified.

## Independent native-source reproduction

The pre-change reproduction is retained in `reproduction/` inside the evidence
archive. A fresh raw-graph audit uses pinned Foyer typing, raw XML decimal
literals, independently enumerated components, and a separately constructed
RDKit biphenyl reference. Production preparation raises:

```text
OPLSAssignmentError: Native component neutrality failed: -0.23 e
```

The frozen benchmark is `[*:1]c1ccc([*:2])cc1`, DP2, two explicit-H C12H10
chains (44 atoms), seeds 2038/2039. Its graph identity is
`f7730225baafb0a2164cdcf3db2f8813d83d3b8d145923053b0fc4ffd193f052`;
packing identity is
`4aadf3eee05388b2eca1c47cac2ec353c963e4389136f7251d368bb1785fb66c`.
The audit records every stable ID, original ordering, aromatic bond, repeat and
end-group metadata, remapping, molecule membership, formal charge, selected
source type, definition, location and original charge literal.

For **each** component:

| Inventory | Source type | Original charge | Contribution |
|---|---|---|---|
| 12 C | opls_145 | −0.115 e | −1.380 e |
| 10 H | opls_146 | +0.115 e | +1.150 e |
| Complete neutral C12H10 | | | **−0.230 e** |

System total: −0.460 e, formal total zero. Decimal arithmetic and compensated
summation agree. Inter-ring ipso carbon IDs 4/7 and 26/29 have no attached H;
all other aromatic carbons have H. Both C6H5 repeat inventories, terminal H,
valence, complete biphenyl topology and component membership are correct. There
are no omitted/duplicated atoms or lost formal charges. This is not a rounding,
serialization, chain-end remapping or intended ionic-repeat convention. The
pinned source's generic aromatic charge pair is not charge-complete for these
biaryl substitution sites. Its zero-charge `opls_147` fused/naphthalene `R2`
rule does not apply to those carbons and was not substituted. No counterions or
compensating-background model were added.

Authoritative XML SHA-256:
`c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
The source points to DOI `10.1021/ja9621760`; its hydrogen definition is explicitly
described as “benzene H.” The carbon SMARTS matches degree-three six-membered
ring carbons without requiring an attached H, so it also selects both ipso sites.
The conclusion concerns this pinned library profile, not every published OPLS
parameterization.

The source revision, pinned Foyer implementation and all original rows remain
unchanged. No newly established native exact-graph neutral OPLS source is claimed.

## Explicit alternative and immutable binding

J38-C3 generated reviewed native-memory RESP/Psi4 artifacts for **this exact
44-atom graph**. J39 independently re-audits those raw inputs, optimized
coordinates, atom/element order, component formal charges, native NPY bits,
intermediates, executable/runtime identities and SHA-256 anchors before MD.
This is a validated-source import, **not new quantum generation in J39**, and
not a transfer from benzene, PE DP3, PEO, or another graph. An additional official
CDS RESP refit from the native ESP matrices reproduces both aromatic component
arrays bit-for-bit in the pinned Psi4 runtime; this reexecutes fitting, not SCF
or geometry optimization. Its command, output and hashes are retained. The raw evidence is
included in the new bundle and retains its original hashes.

Accepted RESP identity:
`0fcbcbcab2d8cb73a89f943ae4bf5b0daade08a480eb24d232e982d8d7315bf3`.
Every charge is bit-identical to the authoritative native binary64 vector.
The complete RESP vector supplies all 44 charges; native OPLS charges supply
none. Native AM1-BCC remains blocked by the retained precision limitation.
External OPLS provided charges remain a separate historical origin.

The J36-C2 policy is unchanged: exact binary Fraction residual decisions,
compensated `math.fsum` and `Decimal.from_float` cross-checks; allowance
`sum(0.5*ulp(q_i)) + gamma_n*sum(abs(q_i))`, where
`gamma_n=n*2^-53/(1-n*2^-53)`. This measures finite-precision consistency,
not charge-model accuracy. Component absolute residuals are
2.2551405187698492e−16 and 1.0061396160665481e−16 e, below respective bounds
6.656880071025696e−15 and 6.581383438004177e−15 e. System residual is
−3.2612801348363973e−16 e, below 2.6057374458176912e−14 e.
No charge values were repaired and no tolerances were relaxed.

`island.opls-compatible-resp.v1` resolves complete bonded/LJ rows independently
of native charge assignment. It preserves original OPLS geometric LJ mixing,
1–4 Coulomb/LJ factors 0.5, exclusions and source parameter values. The pinned
XML has no explicit improper section; its existing RB-proper convention is
preserved. The upstream warning about unassigned impropers is retained in the
logs, and no improper term is fabricated or estimated. This verifies the
existing source profile, not universal interaction coverage. It binds the
original XML, validated pinned typing, exact final graph, complete finite RESP
record and policy. The existing periodic snapshot contains the additive record;
its source identity binds both source XML identity and RESP vector identity.
An offline validation hook rebuilds the numerical assignment and nested records
from these retained sources. Historical J23/J32 serialized records are unchanged.
Only this reviewed graph/vector profile is supported; this is not a general
OPLS charge override or universal force-field coverage.

The matrix retains every original C3 cell unchanged and adds explicit aromatic
paths: native OPLS (`source_incompatible`, no dynamics); separately provided
non-RESP vector (`asset_unavailable`, no source declared); and OPLS-compatible
RESP/Psi4 (software pass). The native/global gate stays nonzero even when the
alternative software gate passes.

## Actual periodic workflow and independent checks

The final run executes minimization, 64 NVT steps and 256 NPT steps, timestep
0.5 fs, target 300 K / 1 bar, Reference platform, seeds 2033/2034/2035, and
retains 321 frames. Backend cutoff is 0.5 nm, periodic Ewald tolerance 1e−5,
with geometric LJ and unchanged source conventions. The legacy NPT record's
0.9 nm default is not the executor's cutoff; the backend's 0.5 nm setting is
recorded and used. No historical setting was changed to conceal this distinction.

A separate split workflow stops at NPT step 128, relocates, and resumes in a
fresh Python process using only that bundle. The unchanged J33-C1 architecture
replays the independent NPT prefix to reconstruct OpenMM's adaptive barostat
proposal state before checkpoint restoration. All retained positions,
velocities, boxes, energies, measured kinetic temperatures, derived volumes and
densities, stages, workflow/source/assignment identities, checkpoint bytes and
opaque complete RNG states compare exactly, with zero comparison tolerances.
Pressure observations are **unavailable**; 1 bar is a target, never a measured
pressure. Actual box volumes change. No trajectory or thermodynamic values are
copied from the uninterrupted workflow into split integration.

The independent reference uses pinned Foyer.apply → ParmEd parameter selection
and geometric mixing, applies the full independently audited RESP vector,
reconstructs Coulomb exceptions from independently traversed bond distances,
and evaluates the first/last actual frames. Maximum energy error is
5.684341886080802e−14 kJ/mol; maximum force error is
8.526512829121202e−14 kJ/(mol Å), both below the existing 1e−5 thresholds.
Bond, angle, RB and nonbonded components agree. This is independent parameter
construction using the **same OpenMM Reference engine**, not an independent
engine or a LAMMPS KSPACE/Ewald/PPPM validation.

A second fresh MDAnalysis subprocess recomputes the full observable payload after
relocation and reproduces it exactly. Fresh MDAnalysis reconstruction records chain Rg/Re/Re-Rg, orientation/nematic
order, valid carbon RDF, contacts, topology/connectivity and density diagnostics,
with selections, units and PBC rules inherited unchanged from J38. The trajectory
NPZ includes positions, velocities, box vectors and IDs and is checked against
every retained frame. Four contiguous 64-frame NPT density blocks have mean
0.29758994769749514 g/cm³. Density always uses actual mass/volume. Measured
mean NPT kinetic temperature is about 170.884 K despite the 300 K target:
this short run is **insufficient_sampling**, not equilibrated polymer evidence.
No density or temperature target is substituted for a measurement.

## Reproduction and validation

Use the declared environments adjacent to the checkout. Exact argv, environments,
tested SHA, exit codes, elapsed times, JUnit and raw logs are archived under
`docs/evidence/phase_4j39/validation`. The environment inventory fingerprints
all five Python executables and package sets. The supported, Foyer and analysis
environments use Python 3.11.16 / OpenMM 8.6.1; Foyer is 1.2.0, ParmEd 4.3.1,
and MDAnalysis 2.10.0. The AmberTools control environment uses Python 3.10.19
/ OpenMM 8.4.0, and the pinned quantum environment uses Python 3.11.14.
The environments are supported tests, Foyer/OpenMM,
MDAnalysis, AmberTools, and Psi4/RESP. The imported quantum protocol retains
Psi4 1.9.1, CDS RESP Conda package 1.0/module 0.8 (Python distribution
metadata reports 0.0.0; all three identities are retained), optimization HF/STO-3G (GAU_LOOSE,
100 iterations), ESP HF/6-31G*, and the original native fitting settings.

From the repository root, with new empty output names:

```bash
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
mkdir /tmp/j39-source-replay
tar -xzf docs/evidence/phase_4j38_c3/fresh-source.tar.gz -C /tmp/j39-source-replay
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j39 \
  --reviewed /tmp/j39-source-replay/raw --output /tmp/j39-replay --software-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j39_receipt \
  --bundle /tmp/j39-replay --software-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j39_receipt \
  docs/evidence/phase_4j39/receipt.json --software-only
```

`--software-only` returns zero only for valid completed alternative software
acceptance. Omitting it returns 1 for the preserved native/global block; malformed
evidence returns 2. Missing source files preserve `asset_unavailable` and the
blocking stage with no downstream integration. All tests use the established
integration configuration; no required tests or markers were disabled. The
suite executes its existing short integration runs in addition to the final
acceptance run. Two earlier exploratory workflow validations also executed short
MD before the final receipt freeze; they are not used as final acceptance proof.
No long scientific campaign is claimed.

An earlier validation run was superseded after tightening receipt schema/model
assertions: its full/regression processes were terminated (−15), and its partial
logs are archived as exploratory evidence, not passing suite results. The final
suites run against the final tested implementation above.

Final results for implementation `1f23a4a`:

| Gate | Actual result |
|---|---|
| J39 focused | 43 passed |
| J22–J38-C3 regressions | 886 passed |
| Full supported suite | 2,897 passed, 10 skipped (23m 04s) |
| Ruff / compileall / diff check | passed |
| pip check, all five environments | passed |
| Raw source audit / official RESP refit | passed, unchanged charge bits |
| Relocation / fresh analysis and continuation | passed, exact payloads/checkpoints/RNG |
| Historical receipt diagnostics and byte preservation | passed; original blocked acceptance exits remain 1 |

The ten established optional skips are: one opt-in AmberTools QM integration,
seven archived Amber reference cases needing `ISLAND_AMBERTOOLS_REFERENCE_MANIFEST`,
one retained analysis reference needing `ISLAND_ANALYSIS_REFERENCE_ROOT`, and one
OPLS permutation check needing `ISLAND_OPLSAA_XML` in its prescribed environment.
Their exact reasons and test IDs are in the JUnit file. These remain unexecuted
gates, not positive scientific evidence.

The additive receipt records exact counts and all command results. Optional
skips are reported separately and never counted as passed validation. Historical
J17–J38-C3 documents/evidence and `.gitattributes` (2,078 tracked files) are
byte-identical to the actual main base. Unrelated untracked notebook checkpoints
and `radonpy.pdf` were preserved.

## Limits

Native aromatic OPLS remains source-incompatible; a separate non-RESP provided
source is unavailable. Finite-time software stability and exact continuation do
not establish thermodynamic equilibrium, experimental density agreement,
experimental validation, Tg or mechanical properties. There is no claim of
universal OPLS-AA/GAFF2 coverage, PCFF completion, AM1-BCC support, universal
msi2lmp parity or LAMMPS KSPACE validation.

```text
production_validated=False
simulation_readiness="not_established"
```
