# Phase 4J28 — final-graph-bound all-atom analysis

Phase 4J28 adds a versioned `AnalysisSelection` record and connects it to the
retained J23 workflow loader, observable calculation, checksummed analysis
records, relocation, and the CLI.  Selections bind stable atom IDs, molecule
membership, all-atom Rg atoms, an ordered heavy-atom backbone, explicit
endpoints, mass weighting, chain aggregation, and the orthorhombic PBC policy.

The supported software tranche is bounded, explicit-hydrogen, acyclic linear
molecules in orthorhombic fully periodic or nonperiodic cells.  Backbone
unwrapping uses authoritative bonds; hydrogens and pendant atoms are placed
relative to their bonded neighbour and are included in mass-weighted all-atom
Rg.  Re and orientation use the declared endpoints and backbone path.  The
implementation rejects missing bond evidence, ambiguous images, branched or
crosslinked periodic graphs, and unsupported triclinic geometry rather than
silently producing values.

The J25 analysis layer remains a software scaffold.  J26 TrajectoryStore v2
provides transactional persistence and opt-in export; the
`export_periodic_md_workflow_to_trajectory_store()` path exports a validated
completed J23 workflow and is not an online writer integrated into every MD
step.  No production-scale trajectory campaign, physical equilibrium claim,
experimental agreement, Tg/mechanical result, or predictive model is made.

## Reproduction

```text
python -m pytest -q tests/test_phase_4j28_selection.py
../island-validation/phase4e2-env/bin/python scripts/analyze_periodic_trajectory.py \
  ../island-validation/phase4j23/acceptance/test_four_family_minimization_0/min \
  --observables network_connectivity --output /tmp/island-j28-analysis.json
```

The first command exercises explicit-H selection, nonperiodic endpoint
semantics, authoritative-bond rejection, selection identity persistence, and
nested-selection tamper rejection.  The second command runs the real retained
J23 loader and CLI path; it is a software integration check, not physical
equilibrium validation.

Readiness remains:

```text
production_validated=False
simulation_readiness="not_established"
```

Full PCFF coverage, universal msi2lmp parity, LAMMPS KSPACE/Ewald/PPPM runtime
validation, branched/crosslinked/triclinic geometry, and scientific
equilibrium remain blocked.
