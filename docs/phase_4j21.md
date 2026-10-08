# Phase 4J21: force-field-neutral multichain packing

J21 adds deterministic rigid-body packing for complete molecular units. A unit
is a `(label, MolecularSystem, FinalChemicalGraph)` tuple. The packer merges
topology and coordinates, remaps atom IDs into one stable sequence, prefixes
molecule labels with the unit label, and creates an orthorhombic
`SimulationBox`. It does not infer chemistry or import a force-field backend.

## Public records and API

`island.packing` provides `PeriodicPackingConfig`, `PeriodicPackingPlan`,
`PeriodicPackingResult`, `PeriodicPackingTransformation`,
`pack_multichain_periodic`, `apply_periodic_packing_plan`,
`save_periodic_packing_plan`, and `load_periodic_packing_plan`.

The config and plan are checksum-protected JSON records. A config chooses one
of target density in g/cm³ or explicit `(Lx, Ly, Lz)` lengths in Å, periodic
boundary flags, seed, `random_uniform` or `none` rotation, minimum interunit
distance, attempt budget, and density tolerance. Target-density volume uses
the total site mass and Avogadro's constant. Explicit box lengths are retained
exactly. The box shape is currently fixed to `orthorhombic`.

The plan binds ordered unit labels, every input graph identity, a coordinate and
topology system checksum, deterministic source-to-output site IDs, molecule
membership remapping, box and boundary settings, placement policy, accepted
rotation matrices and translations, every accepted coordinate, rejected-attempt
count, provenance, evidence, and output system/graph identities. Its identity
is recomputed from the serialized payload. The graph history contains the
existing final-graph transformation schema with operation
`periodic_multichain_packing`; its parameters contain the additive
`island_periodic_packing_transformation_v1` record and plan identity.

## Placement and reconstruction

Units are validated against their exact final graph before placement. Current
input units must be nonperiodic; a periodic box is created by this stage. Each
unit is centered, rotated by a deterministic local RNG, translated into the
box, wrapped into the box, and accepted only when every atom is at least the
configured minimum distance from atoms in earlier units under the selected
boundary convention. Intraunit distances are unchanged. A failed placement
consumes an attempt and is recorded. No packing or topology operation removes
atoms or bonds.

`apply_periodic_packing_plan` reconstructs the output from the serialized
placements. When units are supplied to `load_periodic_packing_plan` or the
apply function, it recomputes unit identities, remappings, rigid rotations,
translations, coordinates, and minimum distances. Re-signed coordinate or
rotation edits therefore reject even when their outer JSON checksum is valid.
The plan can be loaded in a separate process without rebuilding PSMILES
chains, rerunning crosslink construction, automatic typing, or importing
RDKit, OpenMM, Foyer, ParmEd, SciPy, or a force-field module.

Crosslink rows from each input system are remapped into output site IDs and
retained in output provenance. Input repeat provenance is retained under the
packed unit records. Original molecule membership labels are retained with a
unit prefix while connected components are recomputed for the merged graph.
The authoritative output `FinalChemicalGraph` includes the new periodic box,
membership mappings, crosslink provenance, and packing transformation.

## Example

```python
from island.graph import build_psmiles_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic

system0, graph0, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
system1, graph1, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=10)
config = PeriodicPackingConfig(
    target_density=0.90,
    periodic=(True, True, True),
    seed=2026,
    rotation_mode="random_uniform",
    minimum_interunit_distance=1.5,
    maximum_attempts=10000,
    density_tolerance=1e-8,
)
result = pack_multichain_periodic(
    [("chain-0", system0, graph0), ("chain-1", system1, graph1)],
    config=config,
)
```

After packing, the selected PCFF, OPLS-AA, GAFF, or GAFF2 preparation path can
consume `result.graph` and independently type and assign charges. J21 does not
perform force evaluation, periodic energy calculations, packing relaxation,
parameter lookup, MD, or scientific validation of density or polymer
properties. Readiness remains `production_validated=False` and
`simulation_readiness="not_established"`.

## Corrective revision

The corrective revision keeps the v1 plan payload and identity for the original
explicit-box-only path. Configurations with a target density use
`island_periodic_packing_plan_v2` and record `target_density`,
`calculated_density`, and `density_tolerance`. A config may provide either
value or both. When both are provided, the explicit box lengths determine the
box and the calculated mass density must be within the declared tolerance of
the target. Density tolerance is part of the config identity and the v2 plan
identity.

`load_periodic_packing_plan(path, units=None, expected_identity=None)` validates
the optional trusted identity after loading. Plan publication validates the
complete record, writes and fsyncs a temporary file in the destination
directory, and publishes it atomically. Existing destinations are rejected and
temporary files are removed after publication failures.

Periodic final graphs are rejected by the shared preparation gate before PCFF,
OPLS-AA, GAFF, or GAFF2 typing and parameterization. The diagnostic identifies
the selected family and states that periodic final graphs are unsupported by
the current backend. This correction does not add periodic force evaluation or
change any force-field record.
