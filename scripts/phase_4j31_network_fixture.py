"""Reproduce J20 -> J21 network fixtures, topology oracles and source audits.

Graph-only frames are explicit analytical controls, never assigned MD states.
Run with the supported environment. Sources and artifact directories are explicit.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from island.analysis import (
    NetworkAnalysisSelection,
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    calculate_network_observables,
    save_network_analysis_bundle,
)
from island.analysis.network import _digest
from island.builders import build_linear_polymer
from island.core import Coordinates, MolecularSystem, Topology
from island.crosslinking import ReactiveSiteRule, generate_crosslink_network
from island.graph import final_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic

ROOT = Path(__file__).resolve().parents[2]
FRC = ROOT / "lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"
XML = ROOT / "island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml"
RULE = ReactiveSiteRule(
    "J31 controlled neutral C/O",
    "j31_reactive",
    True,
    ("C", "O"),
    1,
    1.0,
    "Explicit H removal before J20",
    ("docs/phase_4j31.md",),
)


@dataclass
class NetworkFixture:
    kind: str
    input_system: object
    input_graph: object
    crosslinked: object
    packed: object


def build_network_fixture(kind="crosslinked", *, seed=2026):
    """Two DP2 PEO chains; zero, one C/O link, or two C/O+C/C links.

    J20 cannot plan on periodic graphs. Crosslink the nonperiodic molecules first,
    then pack the complete network as one rigid unit via J21, preserving both
    original molecule labels. Coordinates are deliberately unequilibrated.
    """
    if kind not in {"uncrosslinked", "branched", "crosslinked"}:
        raise ValueError("Unknown J31 fixture")
    topology, positions = Topology(), {}
    for chain_index in range(2):
        chain = build_linear_polymer(
            "[*:1]CCO[*:2]", dp=2, random_seed=seed + chain_index
        )
        # PEO heavy path is C1-C2-O3-C4-C5-O6. Link internal C2 to
        # terminal O6, and optionally C4 to C1, across the two molecules.
        targets = (2, 4) if chain_index == 0 else (6, 1)
        if kind != "crosslinked":
            targets = targets[:1]
        for atom_id in targets:
            hydrogen = next(
                i
                for i in sorted(chain.topology.neighbors(atom_id))
                if chain.topology.sites[i].element == "H"
            )
            chain.topology.remove_site(hydrogen)
            chain.coordinates._positions.pop(hydrogen)
            chain.topology.sites[atom_id].metadata["j31_reactive"] = True
        center = np.mean(
            [chain.coordinates.get(i) for i in chain.topology.sites], axis=0
        )
        remap = {i: i * 10 + chain_index for i in chain.topology.sites}
        for atom in chain.topology.sites.values():
            copied = deepcopy(atom)
            copied.id = remap[atom.id]
            copied.metadata["chain_id"] = f"chain-{chain_index}"
            topology.add_site(copied)
            positions[copied.id] = (
                chain.coordinates.get(atom.id)
                - center
                + np.array([chain_index * 5.0, 0.0, 0.0])
            )
        for bond in chain.topology.bonds.values():
            topology.add_bond(remap[bond.site1], remap[bond.site2], order=bond.order)
    system = MolecularSystem(
        topology,
        Coordinates(positions),
        metadata={
            "polymer": {"fixture": kind, "chains": 2, "repeat_units_per_chain": 2}
        },
    )
    graph = final_graph(system)
    linked = generate_crosslink_network(
        system,
        graph,
        RULE,
        target_crosslinks={"uncrosslinked": 0, "branched": 1, "crosslinked": 2}[kind],
        seed=seed,
        provenance="J31 explicit C/O two-chain fixture",
        evidence=("docs/phase_4j31.md",),
    )
    if kind != "uncrosslinked":
        # Coordinate generation is separate from the force-field-neutral J20
        # chemistry operation. Embed only the final closed-shell graph.
        from rdkit.Chem import AllChem

        from island.chemistry.rdkit_graph import system_to_rdkit_graph

        converted = system_to_rdkit_graph(linked.system)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(converted.mol, params) != 0:
            raise ValueError("Final network ETKDG embedding failed")
        conformer = converted.mol.GetConformer()
        linked.system.coordinates = Coordinates(
            {
                i: list(conformer.GetAtomPosition(index))
                for i, index in converted.site_id_to_rdkit_index.items()
            }
        )
        linked.graph.validate_integrity(linked.system)
    packed = pack_multichain_periodic(
        [("two-chain-network", linked.system, linked.graph)],
        config=PeriodicPackingConfig(
            box_lengths=(40.0, 40.0, 40.0), seed=seed, minimum_interunit_distance=0
        ),
        provenance="J31 rigid network J21 packing",
        evidence=(linked.plan.identity,),
    )
    return NetworkFixture(kind, system, graph, linked, packed)


def graph_control_trajectory(fixture, *, frames=2):
    """Unassigned static geometry control with explicit blocked evidence hashes."""
    packed = fixture.packed
    graph = packed.graph
    ids = tuple(sorted(packed.system.topology.sites))
    blocked = {
        "status": "blocked",
        "reason": "Analytical graph-only fixture; no force-field assignment or dynamics",
        "final_graph_identity": graph.identity,
    }
    evidence = {
        name: _digest({**blocked, "stage": name})
        for name in (
            "source",
            "typing",
            "charge",
            "assignment",
            "workflow",
            "checkpoint",
        )
    }
    xyz = np.array([packed.system.coordinates.get(i) for i in ids])
    metadata = TrajectoryMetadata(
        graph.identity,
        packed.plan.identity,
        "unassigned_graph_control",
        evidence["source"],
        evidence["typing"],
        evidence["charge"],
        evidence["workflow"],
        evidence["checkpoint"],
        ids,
        tuple(packed.system.topology.sites[i].mass for i in ids),
        {i: n for n, i in enumerate(ids)},
        assignment_identity=evidence["assignment"],
        source_metadata={
            "crosslink_plan_identity": fixture.crosslinked.plan.identity,
            "assignment_status": blocked,
        },
    )
    trajectory = VerifiedTrajectory(
        metadata,
        tuple(
            TrajectoryFrame(
                xyz.copy(), np.diag(packed.system.box.lengths), step, step * 0.001
            )
            for step in range(frames)
        ),
        "j31:analytical-static-control",
        {"geometry": _digest(xyz.tolist())},
        tuple(tuple(b["sites"]) for b in graph.payload["bonds"]),
        molecule_membership={
            int(i): v for i, v in graph.payload["molecule_membership"].items()
        },
    )
    selection = select_network(trajectory, fixture)
    return trajectory, selection


def select_network(trajectory, fixture, *, contact_cutoff_angstrom=4.0):
    return NetworkAnalysisSelection.from_trajectory(
        trajectory,
        graph=fixture.packed.graph,
        crosslink_plan=fixture.crosslinked.plan,
        packing_plan=fixture.packed.plan,
        input_system=fixture.input_system,
        input_graph=fixture.input_graph,
        crosslinked_coordinates=fixture.crosslinked.system.to_dict()["coordinates"],
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        provenance={
            "fixture": fixture.kind,
            "input_graph_identity": fixture.input_graph.identity,
            "crosslinked_graph_identity": fixture.crosslinked.graph.identity,
        },
    )


def independent_topology(graph):
    """Union-find components; enumerate simple junction paths independently.

    This oracle does not call the production component/strand implementation.
    """
    p = graph.payload
    ids = [r["id"] for r in p["sites"]]
    parent = {i: i for i in ids}

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    adjacency = {i: set() for i in ids}
    for row in p["bonds"]:
        a, b = row["sites"]
        parent[root(a)] = root(b)
        adjacency[a].add(b)
        adjacency[b].add(a)
    sizes = Counter(root(i) for i in ids)
    heavy_ids = {row["id"] for row in p["sites"] if row["element"] != "H"}
    heavy = {i: adjacency[i] & heavy_ids for i in heavy_ids}
    junctions = {i for i in heavy if len(heavy[i]) != 2} | {
        i for row in p["crosslink_provenance"] for i in row["sites"]
    }
    paths = set()

    def enumerate_paths(path):
        current = path[-1]
        if len(path) > 1 and current in junctions:
            paths.add(min(tuple(path), tuple(reversed(path))))
            return
        for neighbour in heavy[current]:
            if neighbour not in path:
                enumerate_paths(path + [neighbour])
            elif neighbour == path[0] and len(path) > 2:
                cycle = path + [neighbour]
                paths.add(min(tuple(cycle), tuple(reversed(cycle))))

    for start in junctions:
        enumerate_paths([start])
    return {
        "component_count": len(sizes),
        "component_sizes": sorted(sizes.values(), reverse=True),
        "degree_counts": {
            str(k): v for k, v in sorted(Counter(map(len, adjacency.values())).items())
        },
        "cycle_rank": len(p["bonds"]) - len(ids) + len(sizes),
        "crosslink_count": len(p["crosslink_provenance"]),
        "crosslink_density": len(p["crosslink_provenance"])
        / np.prod(p["periodic"]["lengths"]),
        "strand_lengths": sorted(len(path) - 1 for path in paths),
    }


def prepare_network_backend(
    fixture, family, *, frc=FRC, xml=XML, amberhome=None, work_root=None
):
    """Re-run the selected native backend on the exact final chemical view.

    No historical assignment is accepted. Native failures remain explicit. PCFF
    source-missing rows are rejected even if an older converter convention could
    supply a model zero. Estimated Amber parameters are likewise not accepted.
    """
    from island.forcefields import (
        AmberToolsOptions,
        ForceFieldRequest,
        OPLSOptions,
        PCFFOptions,
        prepare_forcefield,
    )
    from island.periodic import (
        PeriodicForceFieldConfig,
        periodic_chemical_view,
        prepare_periodic_forcefield,
    )

    system = periodic_chemical_view(fixture.packed.system)
    graph = final_graph(
        system,
        molecule_membership={
            int(i): v
            for i, v in fixture.packed.graph.payload["molecule_membership"].items()
        },
    )
    if family == "pcff":
        options = PCFFOptions(frc, (0, 0, 1), (0, 0, 1))
    elif family == "oplsaa":
        options = OPLSOptions(xml)
    else:
        options = AmberToolsOptions(
            family,
            "am1bcc",
            amberhome=amberhome,
            work_root=work_root,
            retain_success_artifacts=True,
        )
    prepared = prepare_forcefield(
        system, ForceFieldRequest(family, options, final_graph=graph)
    )
    native = prepared.native_result
    if family == "pcff":
        missing = [
            row
            for row in native.assignment.payload["assignments"]
            if row["status"] not in {"assigned", "not_applicable"}
        ]
        if missing:
            raise ValueError(
                f"PCFF final graph has {len(missing)} source-missing interactions; no converter-zero fallback accepted"
            )
        source = native.payload["source"]
        mixing, scales = "sixthpower_9_6", (1, 1)
    elif family == "oplsaa":
        source = {
            **prepared.source.identity,
            "sha256": prepared.source.identity["xml_sha256"],
        }
        mixing, scales = "geometric_12_6", (0.5, 0.5)
    else:
        if native.record["parmchk2_estimation_lines"]:
            raise ValueError(
                "Estimated Amber parameters are not independently verified"
            )
        imported = native.imported_result
        source = {
            "family": family.upper(),
            "sha256": imported.source_sha256,
            "source": imported.source,
        }
        mixing = "lorentz_berthelot_12_6"
        scales = (
            imported.nonbonded_policy.lj_scale_14,
            imported.nonbonded_policy.coulomb_scale_14,
        )
    config = PeriodicForceFieldConfig(
        family="OPLS-AA" if family == "oplsaa" else family.upper(),
        source=source,
        box_lengths=fixture.packed.system.box.lengths,
        nonbonded_cutoff=10.0,
        lj_mixing_rule=mixing,
        one_four_scaling=scales,
        electrostatics_method="ewald",
        pme_tolerance=1e-5,
        provenance="J31 reassignment from crosslinked final graph",
        evidence=(fixture.crosslinked.plan.identity, fixture.packed.plan.identity),
    )
    return prepare_periodic_forcefield(
        fixture.packed.system, fixture.packed.graph, prepared, config
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--kind",
        choices=("uncrosslinked", "branched", "crosslinked"),
        default="crosslinked",
    )
    args = parser.parse_args()
    fixture = build_network_fixture(args.kind)
    trajectory, selection = graph_control_trajectory(fixture)
    save_network_analysis_bundle(trajectory, selection, args.output)
    print(
        json.dumps(
            {
                "fixture": args.kind,
                "reference": independent_topology(fixture.packed.graph),
                "observables": calculate_network_observables(trajectory, selection),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
