# Phase 4J29 — real all-atom observable acceptance

Phase 4J29 exercises the public J23 workflow path with the retained
`test_four_family_exact_relocat0/full` fixture. The fixture contains seven
accepted frames for a 20 atom explicit-hydrogen polyethylene DP3 chain, with
authoritative final-graph bonds, stable IDs, molecule membership, an ordered
heavy-atom backbone, explicit endpoints, complete masses, an orthorhombic
periodic cell, and PCFF assignment provenance.

The real command-line path calculates radius of gyration, end-to-end distance,
Re/Rg, gyration tensor, chain anisotropy, chain orientation, nematic order,
density, and network connectivity. The resulting `AnalysisRecord` includes the
selection, graph, mapping, trajectory, workflow, checkpoint, and source
identities. Relocation and a separate Python process reproduce the same record
identity and observable payload without reconstructing chemistry or loading a
force-field backend. Same-process relocation is tested separately from the
fresh-process reconstruction test.

The evidence receipt uses `validated_commit` for the exact 40-character commit
whose implementation was executed (`b21ecd0705564b3cc9586ae4ad8038a2f2b47c2a`).
The later receipt and test-only commits remain part of the corrective branch;
the terminology avoids calling the validated implementation commit the current
branch tip.

An independent reference implementation unwraps the first frame from the
declared backbone bonds, places selected hydrogens through their authoritative
bonds, and computes the mass-weighted Rg and endpoint Re directly from stable
IDs. It gives Rg `2.1744452728578727` Å, Re `5.512538802850571` Å, and Re/Rg
`2.535147180598132` for the first frame. The production analysis reports
frame-averaged values over the seven-frame trajectory, so these values are a
control for the trajectory mean.

The acceptance evidence includes software contract verification, independent
analytical numerical control, same-process relocation, and actual
separate-process reconstruction. It does not establish production trajectory quality, physical
equilibrium, density accuracy, Tg, mechanical properties, or experimental
agreement. Unsupported branched/crosslinked periodic graphs, triclinic cells,
ambiguous images, missing bonds, changed atom order, stale checkpoints, and
nonfinite or zero-Rg inputs remain explicit rejection cases.

## Reproduction commands

```text
../island-validation/phase4e2-env/bin/python scripts/analyze_periodic_trajectory.py \
  ../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full \
  --observables radius_of_gyration end_to_end_distance re_over_rg gyration_tensor \
    chain_anisotropy chain_orientation nematic_order density network_connectivity \
  --output /tmp/island-j29-cli.json --block-size 2

../island-validation/phase4e2-env/bin/python scripts/validate_phase_4j29_real_workflow.py \
  ../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full \
  --record /tmp/island-j29-record.json \
  --relocated-record /tmp/island-j29-record-relocated.json --block-size 2

../island-validation/phase4e2-env/bin/python scripts/validate_phase_4j29_receipt.py \
  docs/evidence/phase_4j29.json docs/paper/reproducibility_manifest.json
```

Readiness remains:

```text
production_validated=False
simulation_readiness="not_established"
```
