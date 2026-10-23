"""MDAnalysis observables from actual variable-box NPT frames, never targets."""

import argparse
import json
from pathlib import Path

import numpy as np

from island.analysis.npt import FLAGS, _digest, _signed, _write_json, frame_measurements
from island.analysis.npt_continuation import load_continuation
from island.graph import FinalChemicalGraph
from island.graph.final import pack
from island.graph.network import network_topology
from scripts.phase_4j38_contract import require


def unwrap(xyz, lengths, adjacency):
    """Deterministic BFS and every-edge closure check; no ambiguous half boxes."""
    out = np.empty_like(xyz)
    seen = set()
    for root in range(len(xyz)):
        if root in seen:
            continue
        seen.add(root)
        out[root] = xyz[root]
        queue = [root]
        for a in queue:
            for b in sorted(adjacency[a]):
                raw = (xyz[b] - xyz[a]) / lengths
                require(
                    not np.any(
                        np.isclose(np.abs(raw - np.rint(raw)), 0.5, atol=1e-8, rtol=0)
                    ),
                    "Ambiguous half-box bond",
                )
                delta = (raw - np.rint(raw)) * lengths
                candidate = out[a] + delta
                if b in seen:
                    require(
                        np.allclose(out[b], candidate, atol=1e-7, rtol=0),
                        "Inconsistent periodic cycle; finite moments blocked",
                    )
                else:
                    out[b] = candidate
                    seen.add(b)
                    queue.append(b)
    return out


def analyze(case):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    case = Path(case)
    workflow, _restart = load_continuation(case / "uninterrupted")
    backend = json.loads((case / "uninterrupted/backend.json").read_text())
    evidence = json.loads((case / "fixture.json").read_text())
    graph = backend["final_graph"]
    ids = workflow["manifest"]["atom_ids"]
    mapping = {a: i for i, a in enumerate(ids)}
    sites = {s["id"]: s for s in graph["sites"]}
    require(
        graph["graph_identity"] == evidence["final_graph"]["graph_identity"],
        "Analysis graph/fixture mismatch",
    )
    adjacency = {i: set() for i in range(len(ids))}
    for bond in graph["bonds"]:
        a, b = [mapping[x] for x in bond["sites"]]
        adjacency[a].add(b)
        adjacency[b].add(a)
    masses = np.array([sites[i]["mass"] for i in ids])
    membership = graph["molecule_membership"]
    chains = [
        [mapping[a] for a in ids if membership[str(a)] == member]
        for member in sorted(set(membership.values()))
    ]
    carbons = [mapping[i] for i in ids if sites[i]["element"] == "C"]
    u = mda.Universe.empty(len(ids), trajectory=True)
    u.add_TopologyAttr("masses", masses)
    u.add_TopologyAttr("ids", ids)
    raw_frames = workflow["frames"]
    frames = frame_measurements(raw_frames, workflow["manifest"]["total_mass_dalton"])
    selected = [i for i, f in enumerate(frames) if f["stage"] == "npt"]
    require(len(selected) >= 32, "Insufficient NPT frames")
    rmax = min(min(np.diag(frames[i]["box_vectors_nm"])) * 10 / 2 for i in selected)
    edges = np.linspace(0, rmax, 25)
    result = []
    system_rgs = []
    mask = np.triu(np.ones((len(ids), len(ids)), dtype=bool), 1)
    for a, neighbors in adjacency.items():
        for b in neighbors:
            mask[min(a, b), max(a, b)] = False
    for i in selected:
        f = frames[i]
        xyz = np.array(f["positions_nm"]) * 10
        lengths = np.diag(f["box_vectors_nm"]) * 10
        points = unwrap(xyz, lengths, adjacency)
        u.atoms.positions = points
        system_rgs.append(float(u.atoms.radius_of_gyration()))
        rgs = [float(u.atoms[c].radius_of_gyration()) for c in chains]
        independent = [
            float(
                np.sqrt(
                    np.average(
                        np.sum(
                            (
                                points[c]
                                - np.average(points[c], axis=0, weights=masses[c])
                            )
                            ** 2,
                            axis=1,
                        ),
                        weights=masses[c],
                    )
                )
            )
            for c in chains
        ]
        require(
            np.allclose(rgs, independent, atol=1e-6, rtol=1e-6),
            "Independent Rg mismatch",
        )
        distances = distance_array(xyz, xyz, box=[*lengths, 90, 90, 90])
        delta = xyz[:, None] - xyz[None, :]
        delta -= np.rint(delta / lengths) * lengths
        require(
            np.allclose(distances, np.linalg.norm(delta, axis=2), atol=5e-6, rtol=1e-6),
            "Independent PBC distance mismatch",
        )
        pair = distances[np.ix_(carbons, carbons)][np.triu_indices(len(carbons), 1)]
        shell = 4 * np.pi / 3 * np.diff(edges**3)
        rdf = np.histogram(pair, bins=edges)[0] / (len(pair) * shell / np.prod(lengths))
        row = {
            "frame_index": i,
            "chain_rg_angstrom": rgs,
            "equal_chain_mean_rg_angstrom": float(np.mean(rgs)),
            "contact_fraction": float(np.mean(distances[mask] < 4)),
            "minimum_nonbonded_distance_angstrom": float(distances[mask].min()),
            "maximum_bond_length_angstrom": float(
                max(
                    distances[a, b]
                    for a, neighbors in adjacency.items()
                    for b in neighbors
                )
            ),
            "rdf_cc": rdf.tolist(),
            "density_g_cm3": f["density_g_cm3"],
            "volume_nm3": f["volume_nm3"],
            "measured_temperature_kelvin": f["measured_temperature_kelvin"],
        }
        if evidence["endpoints"]:
            axes = np.array(
                [
                    points[mapping[b]] - points[mapping[a]]
                    for a, b in evidence["endpoints"]
                ]
            )
            re = np.linalg.norm(axes, axis=1)
            require(np.all(re > 0), "Undefined chain orientation")
            axes /= re[:, None]
            nematic = np.linalg.eigvalsh(
                np.mean([1.5 * np.outer(a, a) - 0.5 * np.eye(3) for a in axes], axis=0)
            )[-1]
            row.update(
                chain_re_angstrom=re.tolist(),
                chain_re_over_rg=(re / rgs).tolist(),
                orientation=axes.tolist(),
                nematic_order=float(nematic),
            )
        result.append(row)
    topology = network_topology(FinalChemicalGraph(pack(graph)))
    return _signed(
        dict(
            schema="island.j38-observables.v1",
            frame_content_identity=_digest(raw_frames),
            workflow_identity=workflow["manifest"]["identity"],
            final_graph_identity=graph["graph_identity"],
            packing_identity=evidence["packing_plan"]["identity"],
            selection={
                "atom_ids": ids,
                "molecule_membership": membership,
                "endpoints": evidence["endpoints"],
                "rdf_atom_ids": [ids[i] for i in carbons],
                "rdf_policy": "all unique carbon pairs including intramolecular pairs; per-frame shell/actual-volume normalization",
                "rdf_edges_angstrom": edges.tolist(),
                "pbc_policy": "deterministic BFS minimum-image bonds; all cycle closures verified",
                "mass_weighting": "atomic mass Rg; equal-chain aggregation",
                "contact_policy": "all unique nonbonded pairs <4 angstrom",
                "frame_indices": selected,
                "units": {
                    "Rg": "angstrom",
                    "Re": "angstrom",
                    "Re/Rg": "dimensionless",
                    "RDF": "dimensionless",
                    "orientation": "dimensionless unit vector",
                    "nematic_order": "dimensionless",
                    "contact_fraction": "dimensionless",
                    "density": "g/cm^3",
                    "temperature": "kelvin",
                    "volume": "nm^3",
                },
                "statistics": "four contiguous NPT blocks; no independent-replica uncertainty",
            },
            endpoint_observables_status="ok"
            if evidence["endpoints"]
            else "blocked: connected network has no unique chain end-to-end observable",
            system_rg_status="blocked: disconnected relative component images are not unique"
            if len(topology["components"]) > 1
            else "ok",
            system_rg_angstrom=system_rgs if len(topology["components"]) == 1 else None,
            topology=topology,
            crosslink_density={
                "value": len(topology["crosslink_bonds"]) / len(ids),
                "denominator": "all final-graph atoms",
                "units": "crosslinks/atom",
            },
            npt_statistics=workflow["analysis"],
            frames=result,
            equilibrium_status="insufficient_sampling",
            **FLAGS,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = analyze(args.case)
        code = 0
    except Exception as error:  # noqa: BLE001 -- preserve diagnostic, never invent observables
        result = {"status": "blocked", "diagnostic": str(error), **FLAGS}
        code = 1
    _write_json(args.output, result)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
