# Phase 4J30 — multi-chain periodic all-atom analysis

J30 extends graph-bound analysis to a bounded two-chain explicit-hydrogen
fixture derived from the retained J23 PCFF workflow
`test_four_family_exact_relocat0/full`. The derived fixture retains seven real
workflow frames and the orthorhombic 30 Å box, duplicates a supported linear
molecule through the J21 packing contract, and assigns interleaved stable IDs.
Chain 0 has 20 atoms and six heavy-atom backbone sites; chain 1 has 13 atoms
and four heavy-atom backbone sites. Both chains have authoritative bonds,
explicit endpoints, complete masses, and separate molecule membership.

Rg, Re, Re/Rg, gyration tensor, anisotropy, orientation, nematic order, density,
and connectivity are calculated per chain and frame. Equal-chain aggregation
is the default. Explicit `mass_weighted` aggregation uses declared all-atom
chain masses and changes the selection identity. Integer-box translations of an
individual chain preserve its internal observables. The real CLI accepts
`--derived-multichain` and runs from the retained J23 workflow directory.

Independent first-frame controls give:

| chain | atoms | mass (Da) | Rg (Å) | Re (Å) | Re/Rg |
|---|---:|---:|---:|---:|---:|
| chain-0 | 20 | 86.178 | 2.1744452728578727 | 5.512538802850571 | 2.535147180598132 |
| chain-1 | 13 | 57.116 | 1.4798506695630267 | 3.2753975309780095 | 2.2133297624855457 |

The acceptance evidence covers software contracts, independent numerical
controls, same-process relocation, and fresh-process reconstruction. It does
not establish production trajectory quality, physical equilibrium, density
accuracy, Tg, mechanical properties, or experimental agreement. Rings,
branches, disconnected graphs, crosslinks, triclinic cells, missing bonds,
ambiguous images, stale mappings, and unsupported aggregation policies remain
explicit rejection cases.

## Reproduction

```text
../island-validation/phase4e2-env/bin/python scripts/analyze_periodic_trajectory.py \
  ../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full \
  --derived-multichain \
  --observables radius_of_gyration end_to_end_distance re_over_rg gyration_tensor \
    chain_anisotropy chain_orientation nematic_order density network_connectivity \
  --output /tmp/island-j30-cli.json --block-size 2

../island-validation/phase4e2-env/bin/python -m pytest -q \
  tests/test_phase_4j30_multichain.py tests/test_phase_4j29_real_acceptance.py \
  tests/test_phase_4j28_selection.py
```

Readiness remains:

```text
production_validated=False
simulation_readiness="not_established"
```
