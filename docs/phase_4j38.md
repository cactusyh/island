# Phase 4J38 — real-polymer benchmark and cross-force-field reproducibility

J38 extends the software-method benchmark inventory. It does not establish
thermodynamic equilibrium or scientific production readiness. Its acceptance CLI
returns nonzero when the declared matrix or scientific gates are unmet, while
preserving successful execution and explicit failure evidence separately.

The base is origin/main `fb4dc8bed7a4e69e42c3186a22a5c4e5da39906b`, fetched and
confirmed by `git ls-remote origin refs/heads/main`. The reviewed J37 evidence tip
`2ebcb5277ebfa3cd52672a6c3321827f26b971d6` is an ancestor. Main is not modified.

## Frozen benchmark inventory

| Key | PSMILES | DP | Complete chains | Construction |
| --- | --- | --- | --- | --- |
| pe | `[*:1]CC[*:2]` | 4 | 2 | Fresh explicit-H polyethylene oligomers |
| peo | `[*:1]CCO[*:2]` | 4 | 2 | Rebuilt retained J34/J37 recipe |
| aromatic | `[*:1]c1ccc([*:2])cc1` | 2 | 2 | Fresh poly(para-phenylene) oligomers |
| network | `[*:1]CCO[*:2]` | 3 | 3 ancestors | Rebuilt J32/J37 controlled terminal C/O closure |

These are finite, small oligomer acceptance systems. RDKit completes terminal
hydrogen valence and explicit hydrogens; no implicit-H or disconnected input
chain is accepted. The network removes the specified terminal H atoms before
J20 crosslinking and re-embeds the completed connected graph before rigid-unit
periodic packing. Its three original molecule ancestries remain explicit. It
is a controlled end-linked cycle: one connected component, cycle rank one,
and no heavy-atom branch points. It is not evidence of a branched gel or
macroscopic gel behavior.

Every case has its own new output directory, PSMILES/DP/end-group/stereochemistry
and seed metadata, complete input/packed chemical graphs, stable atom mapping,
packing plan, box, and derived setup density. Retained PEO/network recipes are
rebuilt deliberately with the exact graph identities to which their independent
RESP/Psi4 artifacts apply. No historical trajectory or measurement is copied as
new dynamics.

## Execution and source boundaries

Each of the four benchmarks is attempted with OPLS-AA and GAFF2. Preparation
reassigns the exact packed final graph. No missing interaction is invented and
no native charge is repaired. Unsupported cells retain stage-specific failures.
New PE/aromatic GAFF2 graphs do not have independently reviewed raw RESP/Psi4
artifacts; their charge/source gate is blocked before parameterization or MD.
Those graphs cannot borrow the retained PEO charge vector.

The retained GAFF2 cells replay J37 independent raw-artifact, mapping, component,
geometry, finite-precision and executable provenance checks and actual native
RESP fitting from retained ESP. They then execute fresh GAFF2 preparation and
MD. This is not a new Psi4 geometry optimization or SCF campaign. RESP/Psi4,
external OPLS-AA provided charges, and native AM1-BCC remain distinct origins.
Native AmberTools AM1-BCC remains blocked with its historical precision
diagnostic. The PCFF network control requires the exact source-verified native
typing rejection with no dynamics; unrelated exceptions do not satisfy it.

OPLS-AA runs minimize, 64 NVT steps, and 256 NPT steps at 0.5 fs, 300 K and
1 bar on OpenMM Reference, retaining every step. GAFF2 reuses the frozen J37
8-NVT/32-NPT protocol. Each family has separately integrated uninterrupted and
split workflows, checkpoint relocation, fresh-process continuation, and exact
comparisons of every retained field and opaque checkpoint/RNG bytes. Existing
adaptive barostat restoration independently replays the split prefix before
loading the boundary checkpoint. Historical checkpoint schemas are unchanged. The inherited J33 `cutoff_nm`
configuration field is declarative and is not consumed by that executor; the
effective 0.5 nm cutoff and Ewald policy are taken from the identity-bound
periodic backend and serialized system. The source inventory records both; no
0.9 nm execution is claimed.

Cross-family comparisons bind identical PEO/network graphs and packing, but
charges, parameter models, protocol lengths and trajectories differ. Seeds
are retained explicitly in each configuration.
Reproducibility is assessed within each family; numerical equality or material
property agreement between force fields is not an acceptance criterion.

## Diagnostics and strict records

`island.j38-matrix.v1` and `island.j38-reproducibility.v1` separate six gates:
software execution, source completeness, numerical reference agreement,
finite-time stability, thermodynamic equilibration, and scientific validation.
A readable all-blocked diagnostic record is never accepted as benchmark success.
Only the first four can contribute to the separately named software gate.
Scientific gates remain `insufficient_sampling` and `blocked`, respectively.

The finite-time software sanity criteria require four computable density
blocks, finite density/volume/energies/temperature, 0 < measured T < 2000 K,
and NPT maximum/minimum density < 2. These bounded checks do not prove a stable
material or equilibration. Density is mass divided by actual frame volume.
The existing NPT schema provides block means, standard deviation, standard
error, rolling means and first-half/second-half differences. There is only one
short trajectory per case; independent replicas, long-time stationarity and
thermodynamic equilibration are unestablished. Pressure targets are not measured
pressure: the latter remains unavailable in this backend record.

The new MDAnalysis record binds frame hashes, all atom IDs, molecular ancestry,
endpoint selections, mass weighting, aggregation rules, RDF bin edges, units
and exact NPT frame indices. It calculates chain Rg (mass weighted), Re and
Re/Rg where endpoints exist, end-to-end unit axes, nematic order, carbon-pair
RDF (including intramolecular pairs), nonbonded contact fraction, minimum
nonbonded distance and maximum bonded distance. Orthorhombic graph BFS
unwrapping checks every cycle closure (absolute closure tolerance 1e-7
angstrom) and rejects half-box ambiguity (1e-8 in fractional-box units). Rg is
also checked with independent float64 mass moments; MDAnalysis PBC distances
are checked against independent minimum-image distances. Network endpoints are
explicitly blocked; disconnected system Rg is blocked because relative component
images are not uniquely defined. Network connectivity, component sizes, cycle
rank, strands and crosslink count come from final topology; crosslink density
is crosslinks per all final-graph atoms. Portable NPZ trajectory export contains
actual coordinates, atom IDs and per-frame boxes.

The periodic numerical references use independent native construction:
Foyer.apply/ParmEd geometric mixing for OPLS-AA, and AmberPrmtopFile for GAFF2.
First/last actual frames compare energies, forces and force components. Both
paths use OpenMM Reference; this is not independent-engine validation or
LAMMPS KSPACE/Ewald/PPPM execution. Existing numerical tolerances are retained.

The receipt validator validates the manifest, exact case inventory, JUnit
counts and tested commit; verifies historical byte hashes; extracts the new
archive into a fresh location; rechecks raw RESP provenance, all nested workflow
and checkpoint identities, continuation equality, trajectory export, density,
and independently reconstructed MDAnalysis/reference results. Consistently
re-signing JSON does not bypass these artifact/semantic checks. Diagnostic-only
validation exits zero for a faithful incomplete report; acceptance validation
still exits nonzero.

## Reproduction

From the repository, with declared sibling environments and sources retained:

```sh
export PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38 --output /tmp/j38-new-empty-run
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_receipt docs/evidence/phase_4j38/receipt.json --diagnostic-only
../island-validation/phase4g1-env/bin/python -m scripts.validate_phase_4j38_receipt docs/evidence/phase_4j38/receipt.json
```

The output path must not exist. Both acceptance commands are expected to return
nonzero for the retained incomplete matrix. Exact validation commands,
environments, counts and artifact hashes are recorded in the additive evidence
manifest/receipt. No historical evidence is rewritten.

`production_validated=False`, `simulation_readiness="not_established"`.
No AI prediction, experimental agreement, Tg, mechanical-property, universal
GAFF2/PCFF coverage, universal msi2lmp parity, or production-equilibrium claim is
made.

## Recorded outcome

The tested implementation is
`1f21e6d00a0b661eafd2354726ed68c52085878d`.

| Benchmark | OPLS-AA software gates | RESP/Psi4 GAFF2 software gates |
| --- | --- | --- |
| PE | Passed | Blocked: graph-specific reviewed raw quantum/charge artifacts absent |
| PEO | Passed | Passed |
| Aromatic | Failed before MD: source component charge −0.23 e | Blocked: graph-specific reviewed raw quantum/charge artifacts absent |
| Controlled network | Passed | Passed |

The independent aromatic XML audit finds −0.230 e for **each** neutral
22-atom chain using exact Decimal source reconstruction and compensated float
summation. No charge is normalized or redistributed. The raw source entries,
selected types, atom IDs, components and source hash are retained in
`aromatic-source-audit.json`.

The two GAFF2 case directories are archival copies of the **current executed**
J37 driver outputs within this fresh J38 run; they are not additional
simulations or replicas. The benchmark driver executes three successful
OPLS-AA cases and two successful GAFF2 cases. Its GAFF2 work includes fresh
parameter preparation, native RESP refitting from retained ESP, independently
integrated uninterrupted/split histories and child continuation. Historical
trajectories are not used as new executions.

All five completed cases have independent periodic reference agreement,
changed NPT box volume, four density blocks, byte-identical checkpoints and
exact opaque RNG-state equality between uninterrupted and resumed runs. The
PCFF network control passes its expected typing rejection without MD. The
native AM1-BCC historical acceptance command remains nonzero. All scientific
validation and thermodynamic-equilibration gates remain unmet; the overall
benchmark and full software matrix acceptance commands return **1**.

Validation against the tested implementation:

- Focused J38: **35 passed**.
- J22–J37 regressions: **729 passed**.
- Full supported suite: **2,732 passed, 10 established skips**, with the skip
  inventory unchanged from J37.
- Ruff, compileall, all five environment pip checks and diff checks passed.
- Historical receipts were validated read-only; their original blocked
  acceptance exits are retained. The source-frozen J35 receipt uses its
  unchanged detached historical implementation checkout.
- **1,463** historical evidence/manifest/attribute files match the base Git
  blobs byte-for-byte. Force-field implementation and parameter sources are
  unchanged.

The MD environments use Python 3.11.16, OpenMM 8.6.1 Reference and MDAnalysis
2.10.0 in the analysis environment. The declared AmberTools 24.8 and pinned
Psi4 1.9.1/RESP source environments retain their separate versions and source
fingerprints in the evidence manifests. The full suite repeats existing short
integration runs and native ESP refits; it is not a new scientific simulation
or quantum-optimization campaign.

Original and relocated receipt validation passed in fresh processes. Six
additional consistently re-signed receipt/manifest controls were rejected:
false scientific success, wrong charge origin, mismatched case identity,
changed test counts, changed tested-source hash, and changed retained frame
content. The last control also recomputes the outer archive hash; nested frame
validation still rejects it. These are software negative controls, not
scientific evidence. Their inputs, commands, diagnostics and reconstruction
script are archived in `artifacts/publication-controls.tar.gz`.
