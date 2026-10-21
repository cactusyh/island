# Phase 4J32: real periodic multichain workflow acceptance

J32 constructs chemically complete, explicit-hydrogen PEO oligomers, packs their
actual coordinates with J21, assigns native force fields to the packed final
graph, and executes J22/J23 minimization and short NVT dynamics. It extends the
analysis layer with a relocatable bundle that contains the actual workflow,
checkpoints, source assets and frame contents. This is software acceptance and
physical diagnostics, not equilibrium material-property validation.

The verified base is `d212fbaac72583e5a91f1519166c17a65e6e7218`, the merged J31-C1
evidence tip. Its tested implementation was
`14abbd11548112a15bc3332f6c7177eb73d48a91`. The implementation and evidence commits
for this phase are deliberately separate; the additive receipt records the
implementation commit actually exercised. Main is not merged by this task.

## Frozen chemistry and scope

`scripts/phase_4j32_real_fixture.py` constructs three chains from
`[*:1]CCO[*:2]`, with builder seeds 2032–2034 and explicit terminal hydrogens.
The linear fixture uses DP4 and has 90 atoms, three independent molecules and no crosslinks.
The optional network fixture uses three DP3 chains and removes one H at each terminal reactive C/O site,
uses the unchanged J20 planner to create three C/O links, verifies closed-shell
valence and uniqueness, and embeds the completed 63-atom macrocycle with RDKit
ETKDGv3. This is a connected cyclic multichain network with retained chain
ancestry; it is not evidence of gelation or a bulk percolating polymer.

J21 packs the three independent chains, or the whole completed network as one
rigid unit, into orthorhombic fully periodic boxes. Setup density is 0.4 g/cm³.
The linear seed is 2032, minimum interunit distance 1.8 Å, box length
13.044539881426486 Å. The network seed is 650 and box length
11.806895015924196 Å. Its interunit threshold is zero because only one connected
unit is packed. All topology, ancestry, mapping and packing identities are
retained and packing is replayed during bundle loading.

The network placement seed was selected to keep bonded input coordinates inside
the native Amber finite-molecule input domain. Earlier exploratory placements
were rejected for a bond crossing the input cell. No bond-length guard or
source typing rule was relaxed. Earlier prototype/copied short trajectories are
not used as J32 acceptance evidence.

Preflight DP4 and asymmetric network candidates were also rejected because
unmodified native AM1-BCC output did not meet the existing total-charge check.
A bounded set of embedding/packing candidates was evaluated before freezing
DP3 with embedding seed 2032 and packing seed 650. Selection was based on native
source completeness and existing validation, never a thermodynamic target.
The selected raw charge total is approximately −8.99e-6 e, inside the native
Amber tolerance 1e-4 e. No charge renormalization or tolerance relaxation was
performed. Candidate logs and source artifacts are retained as preflight
records, separate from acceptance and software injection controls. The J32
external-process time budget is 1800 s; it does not change numerical criteria.

## Force-field acceptance contract

| Family / fixture | Required outcome | Source and independent verification |
| --- | --- | --- |
| OPLS-AA / linear | Real workflow completes | Pinned Foyer XML, native upstream force inventory and numerical-content verification, independent interaction-path enumeration |
| OPLS-AA / network | Real workflow completes | Same pinned source, fresh typing/charges/assignment for this final graph |
| GAFF2 / network | Real workflow completes | Native AmberTools AM1-BCC, exact source/executable hashes, no parmchk2 estimates, OpenMM AmberPrmtopFile versus ParmEd/ISLAND energy/force check |
| PCFF / network | Verified typing rejection, no dynamics | Pinned PCFF source, exact cyclic-component diagnostic over all 63 stable IDs, expected PCFFError at bridge_data |

GAFF2 independent disconnected multichain input is outside the existing finite
Amber adapter domain. GAFF is unexecuted in this tranche; GAFF2 is the selected
Amber family. PCFF's frozen J31-C1 contract remains unchanged. J32 separately
checks the same completed cyclic graph used by its successful network cases.
Missing sources, dependencies, unrelated errors or unexpected PCFF success fail
acceptance. No universal force-field coverage or msi2lmp parity is claimed.

Each native preparation occurs on this packed final graph. Independent adjacency
paths verify every bond, angle and proper torsion, and atom types, charges and
nonbonded coverage are checked before dynamics. OPLS geometric LJ mixing and
0.5/0.5 1–4 scales are preserved; GAFF2 retains source-derived Lorentz–Berthelot
mixing and 1–4 scales. The periodic OpenMM backend uses its existing Ewald policy,
5 Å cutoff and 1e-5 tolerance. This does not validate LAMMPS KSPACE/Ewald/PPPM.

## Dynamics and interpretation

Each successful case executes minimization, then 1000 NVT steps at 300 K,
0.5 fs/step (0.5 ps), LangevinMiddle on OpenMM Reference, seeds 2032/2033.
The explicit J32 minimization configuration is 100 kJ/mol/nm and at most 10000
iterations. Preliminary candidates failed the stricter 10 kJ/mol/nm diagnostic;
those failures are not successful acceptance results. J23 numerical checks and
historical tolerances remain unchanged. The actual achieved force and independent
context checks are archived in minimization.json.

A separate 500-step run is persisted, relocated into a continuation directory,
and resumed for 500 steps by a fresh interpreter. The existing J23 comparison
checks the complete uninterrupted/split frame histories, checkpoints and RNG
state. No trajectory is synthesized by copying coordinates. Actual initial and
final boxes, kinetic and potential energies, temperatures, checkpoints and
engine details remain in the bundle.

NPT is explicitly unexecuted: this small, short fixture cannot establish
interpretable density relaxation. NVT is never labeled NPT. The NVT segment
exercises an equilibration workflow but does not establish equilibration.

## Analysis and persistence

`island.condensed-workflow-analysis.v1` is necessary because the earlier descriptor
records do not archive this complete real workflow and its frame-content
bindings. MDAnalysis 2.10.0 runs in a dedicated environment so the original
validation/Foyer dependencies remain unchanged.

Frames 201–1000 provide 800 retained samples, sixteen contiguous 50-frame blocks,
rolling means and provisional block SEM. Correlation independence is not
established. Stationarity and uncertainty are explicitly `insufficient_sampling`.
Measured temperature is calculated from actual kinetic energy with 3N degrees
of freedom (unconstrained J23 model, no COM remover), never from the target label.
Density is derived from mass and the actual periodic volume; constant NVT density
is not a convergence result.

Each observable records atom IDs, ancestry, PBC policy, weighting, aggregation,
units, frame range and block method. Outputs include chain mass-weighted Rg,
ancestry endpoint distance and orientation, nematic order, C–C RDF (including
intramolecular pairs), temporal minimum-image MSD without a diffusion fit,
inter-ancestry contact fraction, contact clusters, chemical components and
crosslink count per Å³. System Rg is blocked for disconnected linear chains
because their relative periodic images are not unique. For the connected
network it requires every cycle closure to be consistent. J31's deterministic
image assignment is reused; independent NumPy checks verify MDA mass moments
and periodic distances. Bond-length and minimum-distance sanity checks are
recorded. Tg and mechanical observables are explicitly blocked.

Loading verifies the exact input inventory/hashes, original complete chains,
J20 plan, J21 packing, graph/mapping, source, assignment, workflow, continuation,
checkpoint and frame identities, then recomputes observables. Both same-process
relocation and a fresh-process reconstruction run against a copied bundle.
Archives preserve byte-exact source output, including upstream fixed-column
padding, without changing historical publication attributes.

## Commands, tests and evidence

From the repository on the established Linux host:

```sh
export PYTHONPATH=src:.
export PYTHONDONTWRITEBYTECODE=1
export OPENMM_CPU_THREADS=1
export OMP_NUM_THREADS=1
../island-validation/phase4e2-env/bin/python scripts/validate_phase_4j32.py --output /fresh/nonexistent/j32-run
../island-validation/phase4e2-env/bin/python -m pytest tests/test_phase_4j32_condensed.py tests/test_phase_4j32_bundle.py
../island-validation/phase4j32-analysis-env/bin/python scripts/phase_4j32_receipt.py
```

The run driver selects the pinned Foyer environment for OPLS and the established
AmberTools installation for GAFF2. All failures retain result.json and logs and
exit nonzero. Its exact required matrix is shared with receipt validation.
Software injections exist only in test subprocesses; there is no production
bypass switch. Mutated/re-signed files are software negative controls, never
scientific evidence.

The additive receipt is `docs/evidence/phase_4j32/receipt.json`; the nested location
preserves J31's frozen historical top-level receipt inventory. It archives exact
commands, environments, test counts/JUnit, artifact hashes, outcomes and read-only
historical revalidation. The archive contains actual workflow and analysis
bundles, native source artifacts and reconstruction logs. Historical J17–J31
receipts, scientific artifacts, parameter sources and Linux atomic publication
code are unchanged. Applicable full-suite optional skips are listed in JUnit;
required tests are not disabled.

Always: `production_validated=False`,
`simulation_readiness="not_established"`. There are no AI models, experimental
regressions, Tg claims, mechanical-property claims or production-readiness claims.

The dedicated analysis environment can be recreated with the supported Python's
`venv --system-site-packages` and `python -m pip install MDAnalysis==2.10.0`;
the receipt pins the actual resolved dependency versions and records its
separate pip check. The original supported and Foyer environments are unchanged.
`scripts/validate_phase_4j32_gates.py DIRECTORY GATE` records every validation
command, its start/end implementation SHA, exit code and environment. Supported
gates include acceptance, j32, j31_c1, j26_j30, j22_j23, full_supported, ruff,
compileall, pip_check, foyer_pip_check, analysis_pip_check, git_diff_check,
historical_receipt, retained_evidence and j32_receipt.
