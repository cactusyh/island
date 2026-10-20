# Phase 4J31 — bounded periodic network graph analysis

Base (reviewed J30-C1): `aa3014648ddda4b032d6b07cdea381792330f99d`.
Original J30 implementation: `55de228de9247017ba3d706d2f6c7bdda63f6087`.
Branch: `codex/phase-4j31-crosslinked-periodic-network-analysis`.
Origin was fetched and both complete SHAs were confirmed with `git ls-remote`
before implementation. Main was not merged.

## Construction and scientific boundary

The fixtures build two explicit-H DP2 PEO chains with distinct original molecule
labels. Selected hydrogens are explicitly removed before J20. Its existing
`ReactiveSiteRule`, `generate_crosslink_network` and `apply_crosslink_plan`
perform site selection, capacity/valence validation, bond creation, local
interaction rebuilding and final-graph generation. The branched fixture has one
C/O crosslink. The cyclic fixture has two links, including a C/O link. A zero-link
control preserves two disconnected, deliberately incomplete pre-reaction chains;
it is an analytical topology control and is not parameterized or simulated.

J20 deliberately rejects periodic-image planning, so chemistry precedes packing.
For closed-shell branched/crosslinked fixtures, ETKDGv3 independently generates
coordinates from the final graph. This is a coordinate-generation operation
outside the crosslink module. J21 packs the resulting network as one rigid
packing unit while preserving the two original molecule memberships. It does
not pretend disconnected-chain random packing can maintain already formed bonds.

`src/island/graph/network.py` derives immutable-analysis facts from the final
chemical graph. `src/island/analysis/network.py` owns selections and records.
No original J17–J30 receipt, graph identity, hash or force-field behavior is
changed. Graph-only controls explicitly retain blocked assignment evidence and
cannot be interpreted as force-field-assigned MD trajectories.

## Selection and records

`island.network-analysis-selection.v1` binds the complete final graph, exact J20
plan and J21 packing plan, original molecule membership, connected-component
membership, ordered stable IDs, mass/mapping identity, source/family, typing,
charge and parameter-assignment identities, workflow/checkpoint/trajectory
identities, every frame's coordinates/box/step/time content hash, geometry policy,
network policy, cutoff and tolerances. It retains the pre-crosslink system/graph
and declared post-crosslink coordinates. Reconstruction replays J20 and J21 and
requires the resulting graph to equal the supplied final graph. No typing or
parameterization is performed by analysis or reconstruction.

Canonical JSON records own their payloads. A selection rejects changed graphs,
plans, masses/mappings, family/source, typing/charge/parameter assignment, workflow,
checkpoint, trajectory or frame contents. A trusted expected identity is needed
to reject wholesale replacement by another internally valid history. Hashes are
integrity checks, not cryptographic attestations of physical validity.

`island.network-analysis-record.v1` binds the selection, identity fields, all
observables and readiness flags. `island.network-analysis-bundle.v1` additionally
stores the trajectory and independently recomputes observables when loading.
Bundle publication uses a sibling staging directory, fsync and existing Linux
`renameat2(RENAME_NOREPLACE)` exclusive atomic publication. Failures leave no
published partial bundle. Individual JSON records use exclusive atomic file
publication. Same-process relocation and fresh-process reconstruction are tested.

## Periodic geometry

Only positive orthorhombic, fully periodic cells are accepted. Each frame gets a
deterministic breadth-first spanning forest rooted at the minimum stable ID of
each component, with sorted neighbors. Bond images use the minimum-image rule.
Half-box ties within the declared tolerance (default `1e-6` Å) are ambiguous and
explicitly rejected. Every non-tree bond is checked against its minimum-image
bond vector **without taking the closure residual modulo the box**.

A nonzero integer lattice closure vector indicates a winding cycle. The strict
`unwrap_network` API rejects it as an inconsistent finite embedding. Analysis
reports nonzero winding axes as the explicitly defined graph-based periodic
spanning/percolation indicator and blocks finite-network Rg. It does not infer a
macroscopic percolation threshold from these small fixtures. Ambiguous images
block both geometric metrics and winding classification. Unsupported observables
return explicit blocked diagnostics. All retained frames are checked.

## Observable definitions

- Components, cluster sizes, gel fraction and cycle rank use all declared atoms.
  The gel component means the largest atom-count component, with minimum stable
  ID breaking ties. Gel fraction is its atom count divided by the total; it is
  not an experimental gel measurement.
- All-atom degree distribution is distinct from branch statistics on the induced
  heavy-atom skeleton. Hydrogen atoms are never counted as polymer branches or
  dangling chain ends.
- Strand lengths are maximal edge-disjoint heavy-bond paths between non-degree-2
  sites or crosslink endpoints. Terminal strands and parallel paths are retained;
  a pure degree-2 ring is counted once. Dangling-end fraction is heavy degree-one
  sites divided by all heavy sites. Cycle rank is `E - V + components`.
- Crosslink count comes from final-graph crosslink provenance, including
  intramolecular links; it is not inferred from different molecule labels.
  Density is crosslink count divided by each frame's periodic volume, in Å^-3.
- Contact fraction uses all unordered heavy-atom pairs except directly bonded
  pairs, including intrachain and interchain pairs. Its declared denominator is
  the number of eligible pairs; its numerator is the number at minimum-image
  distance no greater than the explicit cutoff. Missing cutoffs/no eligible
  pairs are blocked. Cutoffs must be below half the shortest box length.
- Network Rg is mass-weighted over the largest connected component and only
  returned when every cycle in every component has a consistent finite image
  assignment. It is not an Rg of all disconnected clusters under arbitrary images.

## Force-field source-completeness matrix

All execution begins by rerunning the existing native preparation on the exact
final chemical view of the packed graph; the periodic adapter then binds it to
the final periodic graph. Historical pre-crosslink assignments are never reused.

| Family | Branched C/O fixture | Cyclic crosslinked fixture | Independent checks |
|---|---|---|---|
| PCFF | blocked: incomplete automatic typing | blocked: incomplete automatic typing | pinned FRC checked; no dynamics, missing rows never filled |
| OPLS-AA | complete for this fixture | complete for this fixture | pinned Foyer/OpenMM force-content and source-row verification; exhaustive final-graph interaction inventory |
| GAFF | complete for this fixture | complete for this fixture | real AmberTools AM1-BCC; no parmchk2 estimates; exhaustive paths; OpenMM AmberPrmtopFile energy/force reference |
| GAFF2 | complete for this fixture | complete for this fixture | same independent checks with GAFF2 source |

The J31 PCFF path additionally rejects any source-missing interaction even if an
older converter model could explicitly supply a zero under its historical
convention. That existing model is unchanged. Estimated Amber parameters are
also rejected. Failures in executable workflows remain failures, never relabeled
as unsupported-source outcomes. No universal coverage or msi2lmp parity is claimed.

Each executable family/fixture runs four short OpenMM Reference NVT steps, then
compares against two steps, relocation, and fresh-process continuation. Equality
includes every retained frame and complete checkpoint bytes. These tests establish
software and numerical behavior only, not equilibration or production suitability.

## Reproduction

Use the supported `../island-validation/phase4e2-env/bin/python`. OPLS preparation
uses the explicitly selected pinned `../island-validation/phase4g1-env/bin/python`.
Sources are the pinned LAMMPS `pcff.frc` and Foyer `oplsaa.xml`; real AmberTools is
selected through `--amberhome`. Artifacts must have a new, empty destination.

```sh
PYTHONPATH=src:. python -m scripts.phase_4j31_network_fixture \
  --kind crosslinked --output /tmp/new-j31-graph-control
PYTHONPATH=src:. python -m scripts.validate_phase_4j31_network \
  --output /tmp/new-j31-source-acceptance \
  --amberhome ../miniforge3/envs/md_rg
PYTHONPATH=src:. python -m pytest -q tests/test_phase_4j31_network.py \
  tests/test_phase_4j31_execution.py tests/test_phase_4j31_receipt.py
PYTHONPATH=src:. python scripts/validate_phase_4j31_receipt.py \
  docs/evidence/phase_4j31.json docs/paper/reproducibility_manifest.json
```

The receipt records actual test counts, commands, JUnit/log hashes, source matrix,
fixture identities and independent references. Validation verifies current source
hashes, all retained acceptance files, all J17–J30 receipts against the exact base,
fixture reconstruction, source outcomes, manifest consistency and readiness.

## Limits

No triclinic geometry, automatic reaction chemistry, large-network scaling claim,
production equilibrium, experimental agreement, Tg or mechanical claim is made.
The supplied fixtures are deliberately small. J23 execution is Reference only;
LAMMPS KSPACE/Ewald/PPPM runtime acceptance remains blocked. No AI prediction or
experimental-property prediction is added.

```text
production_validated=False
simulation_readiness="not_established"
```
