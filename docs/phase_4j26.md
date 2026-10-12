# Phase 4J26: trajectory integrity and periodic geometry corrections

Base: `5a4e0bcb53c2ed72aa9bc4a07902b7f894aee7e`  
Branch: `codex/phase-4j26-trajectory-integrity-pbc-corrections`

J25 analysis was a software scaffold. J26 corrects its wrapped-chain geometry,
nested identity binding and trajectory-store transaction behavior before any
physical production campaign. Linear chains are unwrapped along explicit bonded
paths with an orthorhombic minimum-image displacement. Branched and crosslinked
periodic unwrapping is rejected. Re/Rg is reported per chain and aggregated as
the equal-weight mean of per-frame chain ratios.

Analysis records now bind graph, packing, force-field, assignment, workflow,
checkpoint, mapping, trajectory and source identities. Trusted expected bindings
or a verified trajectory are required to detect nested tampering. External
MDAnalysis input requires an explicit stable atom-ID order and source hashes.

The opt-in trajectory store uses schema v2, a single-writer lock, staged restart
publication, manifest-declared restart history, checkpoint interval enforcement,
chunk validation and transactional relocation. Valid v1 stores remain readable.
The J23 workflow schema and APIs are unchanged; the additive
`export_periodic_md_workflow_to_trajectory_store()` function writes a validated
completed J23 workflow and its checkpoints into the v2 store. It is an opt-in
export path, not an online writer integrated into every MD step. No production
scale trajectory campaign has been performed.

The supported-environment gate is separate from the base environment. Scientific
production readiness remains blocked:

```text
production_validated=False
simulation_readiness="not_established"
```

Full PCFF coverage, universal msi2lmp parity, LAMMPS KSPACE/Ewald/PPPM runtime
agreement, equilibrium physical validation, experimental comparison, Tg scans
and mechanical tests remain outside this phase.
