"""Strict, relocatable analysis of real J23 multichain trajectories.

Unlike descriptor-only records, this bundle retains the actual workflow and
source assets. Loading checks every archived input and recomputes observables.
The schema declares short-run diagnostics and cannot establish equilibrium.
"""

import json
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.analysis.mdanalysis import AnalysisSelection, load_verified_trajectory
from island.analysis.network import _body, _digest, _image_geometry, _signed
from island.exceptions import AnalysisError
from island.periodic_md import compare_periodic_md_workflows, load_periodic_md_workflow

SCHEMA = "island.condensed-workflow-analysis.v1"
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
STEPS = 1000
DISCARD = 200
BLOCK = 50


def require(ok, message):
    if not ok:
        raise AnalysisError(message)


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def files(root):
    return {
        str(p.relative_to(root)): file_hash(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name != ".writer.lock" and p != root / "analysis.json"
    }


def statistics(values):
    values = np.asarray(values, dtype=float)
    require(
        len(values) >= 2 * BLOCK and np.isfinite(values).all(),
        "Insufficient trajectory frames for block statistics",
    )
    blocks = (
        values[: len(values) // BLOCK * BLOCK]
        .reshape((-1, BLOCK) + values.shape[1:])
        .mean(axis=1)
    )
    return {
        "mean": np.mean(values, axis=0).tolist(),
        "block_means": blocks.tolist(),
        "block_standard_error": (
            blocks.std(axis=0, ddof=1) / np.sqrt(len(blocks))
        ).tolist(),
        "blocks": len(blocks),
        "block_size_frames": BLOCK,
        "sampling_status": "insufficient_sampling",
        "uncertainty_method": "SEM of contiguous nonoverlapping block means; correlation independence not established",
    }


def validate_fixture(evidence, backend):
    from island.crosslinking import CrosslinkNetworkPlan, apply_crosslink_plan
    from island.graph import FinalChemicalGraph
    from island.graph.final import _system_from_payload, pack
    from island.packing import PeriodicPackingPlan, apply_periodic_packing_plan
    from scripts.phase_4j32_real_fixture import setup_metadata, validate_complete_chain

    chains = [_system_from_payload(s) for s in evidence["complete_chains"]]
    require(len(chains) == evidence["chain_count"] >= 2, "Wrong chain count")
    for chain in chains:
        validate_complete_chain(chain)
    units = [
        (
            r["label"],
            _system_from_payload(r["system"]),
            FinalChemicalGraph(pack(r["graph"])),
        )
        for r in evidence["packing_units"]
    ]
    for _, system, graph in units:
        validate_complete_chain(system)
        graph.validate_integrity(system)
    packed = apply_periodic_packing_plan(
        units, PeriodicPackingPlan(pack(evidence["packing_plan"]))
    )
    require(
        json.loads(json.dumps(packed.system.to_dict())) == evidence["packed_system"]
        and packed.graph.payload == evidence["final_graph"],
        "Stale packed graph or mapping",
    )
    require(
        backend["system"]["coordinates"] == evidence["packed_system"]["coordinates"]
        and backend["final_graph"] == evidence["final_graph"],
        "Packed final graph/backend mismatch",
    )
    require(
        setup_metadata(packed) == evidence["setup"],
        "Invalid minimum-distance/setup result",
    )
    if evidence["kind"] == "linear":
        require(
            evidence["setup"]["component_count"] == len(chains)
            and evidence["setup"]["minimum_interchain_distance_angstrom"]
            >= evidence["setup"]["minimum_interunit_distance_angstrom"]
            > 0,
            "Invalid independent-chain packing distance",
        )
        require(evidence["crosslink_plan"] is None, "Unexpected linear crosslink plan")
    else:
        original = _system_from_payload(evidence["crosslink_input"])
        graph = FinalChemicalGraph(pack(evidence["crosslink_input_graph"]))
        _, final, _ = apply_crosslink_plan(
            original, graph, CrosslinkNetworkPlan(pack(evidence["crosslink_plan"]))
        )
        require(final.payload == units[0][2].payload, "Stale crosslink plan or graph")
        require(
            evidence["setup"]["component_count"] == 1,
            "Network must be packed as one connected unit",
        )
    return packed


def validate_execution(execution, workflow):
    """Bind execution claims to the archived engine state, not a copied label."""
    manifest = workflow["manifest"]
    expected = {
        "family",
        "comparison",
        "config",
        "engine",
        "initial_box_vectors_nm",
        "final_box_vectors_nm",
        "final_checkpoint_sha256",
        "rng_evidence",
        "npt",
        "equilibrium_status",
        "thermodynamic_origin",
        "short_nvt_duration_ps",
        *FLAGS,
    }
    require(set(execution) == expected, "Unrecognized physical or execution claim")
    for key in ("config", "engine", "initial_box_vectors_nm"):
        require(execution[key] == manifest[key], "Execution/manifest mismatch: " + key)
    last = workflow["frames"][-1]
    require(
        execution["final_box_vectors_nm"] == last["box_vectors_nm"]
        and execution["final_checkpoint_sha256"] == last["checkpoint_sha256"],
        "Stale checkpoint or box claim",
    )
    require(
        execution["short_nvt_duration_ps"]
        == manifest["config"]["steps"] * manifest["config"]["timestep_fs"] / 1000,
        "Inherited thermodynamic duration",
    )


def validate_assignment(assignment, backend):
    """Require the independent source report to describe this backend exactly."""
    data = backend["numerical_data"]
    require(
        assignment["atom_types_and_charges"] == data["sites"],
        "Atom typing or charge mutation",
    )
    require(
        assignment["periodic_policy"] == backend["config"],
        "Wrong force-field family or periodic policy",
    )
    for key, term in (
        ("bonds", "bonds"),
        ("angles", "angles"),
        ("proper_paths", "torsions"),
    ):
        require(
            assignment[key] == len(data[term]), "Independent topology count mismatch"
        )
    require(
        abs(assignment["charge_total_e"] - sum(s["charge"] for s in data["sites"]))
        < 1e-12,
        "Charge total mutation",
    )
    require(
        bool(assignment["independent_source_verification"]),
        "Missing independent source verification",
    )


def trajectory_from_bundle(root):
    workflow = load_periodic_md_workflow(root / "workflow").payload
    manifest = workflow["manifest"]
    backend = manifest["backend"]
    config = manifest["config"]
    require(
        config["mode"] == "nvt"
        and config["pressure_bar"] is None
        and manifest["barostat"] is None,
        "Attempted NPT claim when NPT did not execute",
    )
    require(
        config["steps"] == STEPS
        and len(workflow["frames"]) == STEPS + 1
        and workflow["frames"][-1]["step"] == STEPS,
        "Insufficient trajectory frames or incomplete execution",
    )
    require(config["minimize"] is True, "Minimization required")
    fixture = json.loads((root / "fixture.json").read_text())
    validate_fixture(fixture, backend)
    execution = json.loads((root / "execution.json").read_text())
    validate_execution(execution, workflow)
    assignment = json.loads((root / "assignment.json").read_text())
    validate_assignment(assignment, backend)
    if config["family"] == "GAFF2":
        for relative, digest in assignment["source_executables"].items():
            if relative.startswith("dat/"):
                require(
                    file_hash(root / "source-assets/amber-parameters" / relative)
                    == digest,
                    "Amber source hash mutation",
                )
    require(
        execution["npt"]
        == {
            "status": "unexecuted",
            "reason": "Small finite oligomer acceptance run; density relaxation is not interpreted from this short sample",
        },
        "Contradictory NPT claim",
    )
    require(
        execution["equilibrium_status"] == "not_established"
        and execution["thermodynamic_origin"]
        == "J23 frame kinetic/potential energies and box; no inherited measurements",
        "Copied or inherited thermodynamic metadata",
    )
    require(
        all(
            execution.get(k) == v and type(execution.get(k)) is type(v)
            for k, v in FLAGS.items()
        ),
        "Readiness promotion",
    )
    require(execution["family"] == config["family"], "Wrong force-field family")
    comparison = compare_periodic_md_workflows(root / "workflow", root / "continuation")
    require(
        comparison == execution["comparison"] and comparison["passed"],
        "Workflow comparison/continuation failure",
    )
    source = backend["config"]["source"]["sha256"]
    source_path = (
        root
        / "source-assets"
        / ("oplsaa.xml" if config["family"] == "OPLS-AA" else "amber/result.prmtop")
    )
    require(file_hash(source_path) == source, "Source hash mutation")
    trajectory = load_verified_trajectory(root / "workflow")
    selection = AnalysisSelection.from_final_graph(
        backend["final_graph"],
        atom_mapping_identity=trajectory.mapping_identity,
        provenance={
            "origin": "authoritative J23 final graph",
            "ancestry": "original molecule IDs; crosslinks retained in full graph",
        },
    )
    require(
        len(selection.chains) == fixture["chain_count"], "Wrong molecule membership"
    )
    # Absolute paths are transport locations, never scientific identities.
    trajectory = replace(trajectory, source_path="workflow", selection=selection)
    return workflow, trajectory, fixture, execution


def analyze(root):
    """MDAnalysis mass moments and distances over unwrapped genuine MD frames."""
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    root = Path(root)
    workflow, tr, fixture, execution = trajectory_from_bundle(root)
    backend = workflow["manifest"]["backend"]
    graph = backend["final_graph"]
    ids = tr.metadata.atom_ids
    mapping = tr.metadata.atom_mapping
    masses = np.array(tr.metadata.masses_dalton)
    adj = {i: set() for i in ids}
    for a, b in tr.topology_bonds:
        adj[a].add(b)
        adj[b].add(a)
    u = mda.Universe.empty(len(ids), trajectory=True)
    u.add_TopologyAttr("masses", masses)
    u.add_TopologyAttr("ids", ids)
    chains = tr.selection.chains
    ci = [[mapping[i] for i in c["rg_atom_ids"]] for c in chains]
    carbons = [mapping[s["id"]] for s in graph["sites"] if s["element"] == "C"]
    xyz = np.array([f.coordinates_angstrom for f in tr.frames])
    lengths = np.array([np.diag(f.box_vectors_angstrom) for f in tr.frames])
    require(
        all(
            np.array_equal(
                f.box_vectors_angstrom, np.diag(np.diag(f.box_vectors_angstrom))
            )
            for f in tr.frames
        ),
        "Only orthorhombic boxes supported",
    )
    require(
        np.all(lengths > 0) and np.all(lengths == lengths[0]),
        "Missing box or unexplained density drift",
    )
    delta = np.diff(xyz, axis=0)
    delta -= np.rint(delta / lengths[1:, None, :]) * lengths[1:, None, :]
    require(
        np.max(np.abs(delta) / lengths[1:, None, :]) < 0.25,
        "Ambiguous temporal image assignment",
    )
    temporal = np.concatenate([xyz[:1], xyz[:1] + np.cumsum(delta, axis=0)], axis=0)
    selected = range(DISCARD + 1, len(tr.frames))
    rg = []
    re = []
    orient = []
    nematic = []
    system_rg = []
    rdf = []
    contacts = []
    clusters = []
    minima = []
    bond_max = []
    rmax = float(lengths.min() / 2)
    edges = np.linspace(0, rmax, 33)
    pair_mask = np.triu(np.ones((len(ids), len(ids)), dtype=bool), 1)
    membership = np.array([tr.molecule_membership[i] for i in ids])
    cross_mask = pair_mask & (membership[:, None] != membership[None, :])
    for a, b in tr.topology_bonds:
        cross_mask[min(mapping[a], mapping[b]), max(mapping[a], mapping[b])] = False
    full_consistent = True
    for fi in selected:
        frame = tr.frames[fi]
        geometry = _image_geometry(frame, mapping, adj, 1e-8)
        full_consistent &= geometry["consistent"]
        # Whole graph unwrapping also verifies every closing bond.
        require(
            geometry["consistent"],
            "Inconsistent periodic cycle; finite moments blocked",
        )
        points = np.array([geometry["positions"][i] for i in ids])
        u.atoms.positions = points
        per_rg = [float(u.atoms[index].radius_of_gyration()) for index in ci]
        for index, value in zip(ci, per_rg):
            center = np.average(points[index], axis=0, weights=masses[index])
            reference = np.sqrt(
                np.average(
                    np.sum((points[index] - center) ** 2, axis=1), weights=masses[index]
                )
            )
            require(
                np.isclose(value, reference, rtol=1e-6, atol=1e-6),
                "Independent mass moment mismatch",
            )
        axes = np.array(
            [
                points[mapping[c["endpoint_ids"][1]]]
                - points[mapping[c["endpoint_ids"][0]]]
                for c in chains
            ]
        )
        ends = np.linalg.norm(axes, axis=1)
        require(np.all(ends > 0), "Zero chain orientation")
        axes /= ends[:, None]
        rg.append(per_rg)
        re.append(ends)
        orient.append(axes)
        nematic.append(
            np.linalg.eigvalsh(
                sum(1.5 * np.outer(a, a) - 0.5 * np.eye(3) for a in axes) / len(chains)
            )[-1]
        )
        if fixture["setup"]["component_count"] == 1:
            system_rg.append(float(u.atoms.radius_of_gyration()))
        distances = distance_array(xyz[fi], xyz[fi], box=[*lengths[fi], 90, 90, 90])
        # Independent float64 minimum-image check validates MDA geometry.
        diff = xyz[fi, :, None, :] - xyz[fi, None, :, :]
        diff -= np.rint(diff / lengths[fi]) * lengths[fi]
        require(
            np.allclose(distances, np.linalg.norm(diff, axis=2), atol=5e-6, rtol=1e-6),
            "MDAnalysis periodic distance mismatch",
        )
        cc = distances[np.ix_(carbons, carbons)][np.triu_indices(len(carbons), 1)]
        hist = np.histogram(cc, bins=edges)[0]
        shell = 4 * np.pi / 3 * np.diff(edges**3)
        pairs = len(carbons) * (len(carbons) - 1) / 2
        rdf.append(hist / (pairs * shell / frame.volume))
        contacts.append(float(np.mean(distances[cross_mask] < 4.0)))
        minima.append(float(distances[cross_mask].min()))
        bond_max.append(
            max(distances[mapping[a], mapping[b]] for a, b in tr.topology_bonds)
        )
        connected = {i: {i} for i in range(len(ci))}
        for i in range(len(ci)):
            for j in range(i):
                if np.any(distances[np.ix_(ci[i], ci[j])] < 4.0):
                    merged = connected[i] | connected[j]
                    for k in merged:
                        connected[k] = merged
        clusters.append(len({tuple(sorted(v)) for v in connected.values()}))
    retained = np.array(list(selected))
    frame_range = {
        "first_step": int(retained[0]),
        "last_step": int(retained[-1]),
        "stride": 1,
        "discarded_initial_steps": DISCARD,
        "count": len(retained),
    }
    policy = {
        "atom_ids": list(ids),
        "molecule_membership": {str(k): v for k, v in tr.molecule_membership.items()},
        "selection_identity": tr.selection.identity,
        "frame_range": frame_range,
        "pbc_unwrapping": "J31 sorted BFS and every-cycle closure; MSD temporal minimum image < quarter box per step",
        "mass_weighting": "atomic mass for Rg/density; equal atoms for MSD/RDF/contact",
        "chain_aggregation": "equal chain weights",
        "block_statistics": {
            "block_size": BLOCK,
            "method": "contiguous block means and provisional standard error; no independent-sample claim",
        },
    }

    def observable(values, units, **extra):
        return {
            "status": "ok",
            "units": units,
            "policy": policy,
            **statistics(values),
            **extra,
        }

    def blocked(reason, units):
        return {"status": "blocked", "units": units, "policy": policy, "reason": reason}

    density = np.array(
        [tr.metadata.mass_total_dalton / f.volume * 1.66053906660 for f in tr.frames]
    )[retained]
    displacement = temporal[retained] - temporal[retained[0]]
    msd = np.mean(np.sum(displacement**2, axis=2), axis=1)
    energies = np.array([f["potential_energy_kj_mol"] for f in workflow["frames"]])[
        retained
    ]
    kinetic = np.array([f["kinetic_energy_kj_mol"] for f in workflow["frames"]])[
        retained
    ]
    temperatures = 2 * kinetic / (3 * len(ids) * 0.00831446261815324)
    observables = {
        "density": observable(
            density, "g/cm^3", meaning="constant NVT box, not convergence evidence"
        ),
        "radius_of_gyration": observable(
            np.mean(rg, axis=1), "angstrom", per_chain=statistics(rg)
        ),
        "system_radius_of_gyration": observable(system_rg, "angstrom")
        if system_rg
        else blocked(
            "Independent components have no unique relative periodic images", "angstrom"
        ),
        "end_to_end_distance": observable(
            np.mean(re, axis=1),
            "angstrom",
            per_chain=statistics(re),
            meaning="original chain ancestry endpoints; not a linear network interpretation",
        ),
        "chain_orientation": observable(orient, "unit vector"),
        "nematic_order": observable(nematic, "dimensionless"),
        "rdf": observable(
            rdf,
            "dimensionless",
            atom_ids=[ids[i] for i in carbons],
            edges_angstrom=edges.tolist(),
            normalization="all unique C-C pairs including intramolecular; ideal gas spherical shell count",
        ),
        "msd": observable(
            msd,
            "angstrom^2",
            series=msd.tolist(),
            definition="all atoms, equal weight, displacement from first retained frame; no diffusion fit",
        ),
        "contact_fraction": observable(
            contacts,
            "fraction",
            cutoff_angstrom=4.0,
            denominator="inter-ancestry atom pairs excluding direct bonds",
        ),
        "contact_clusters": observable(
            clusters, "count", definition="chain graph with any atom-pair distance <4 A"
        ),
        "network_connectivity": observable(
            [fixture["setup"]["component_count"]] * len(retained),
            "count",
            component_sizes=[len(c) for c in graph["components"]],
        ),
        "crosslink_density": observable(
            [
                fixture["setup"]["crosslink_count"] / f.volume
                for f in tr.frames[DISCARD + 1 :]
            ],
            "crosslinks/angstrom^3",
            denominator="instantaneous periodic volume",
            crosslink_count=fixture["setup"]["crosslink_count"],
        ),
        "glass_transition": blocked(
            "Short single-temperature acceptance run; no Tg inference", "K"
        ),
        "mechanical_properties": blocked(
            "No deformation protocol or stress validation", "not_applicable"
        ),
    }
    require(np.max(bond_max) < 2.2, "Bond-length sanity failure")
    require(min(minima) > 0.5, "Nonbonded minimum-distance sanity failure")

    def diagnostic(series):
        return {
            **statistics(series),
            "rolling_mean": np.convolve(
                series, np.ones(BLOCK) / BLOCK, mode="valid"
            ).tolist(),
            "last_minus_first_block": float(
                np.mean(series[-BLOCK:]) - np.mean(series[:BLOCK])
            ),
        }

    diagnostics = {
        "stationarity": "insufficient_sampling",
        "temperature": diagnostic(temperatures),
        "temperature_definition": "2 * actual J23 kinetic energy / (3N R); unconstrained model, no COM remover",
        "temperature_units": "K",
        "energy": diagnostic(energies),
        "energy_units": "kJ/mol",
        "density": diagnostic(density),
        "density_units": "g/cm^3",
        "minimum_nonbonded_interchain_distance_angstrom": min(minima),
        "minimum_distance_violations_below_0_5_angstrom": int(
            sum(x < 0.5 for x in minima)
        ),
        "maximum_bond_length_angstrom": float(max(bond_max)),
        "equilibrium_status": "not_established",
    }
    bindings = {
        k: tr.metadata.to_dict()[k]
        for k in (
            "final_graph_identity",
            "packing_identity",
            "force_field_family",
            "force_field_source_identity",
            "typed_graph_identity",
            "charge_identity",
            "assignment_identity",
            "workflow_identity",
            "checkpoint_identity",
        )
    }
    bindings.update(
        periodic_backend_identity=backend["identity"],
        atom_mapping_identity=tr.mapping_identity,
        selection_identity=tr.selection.identity,
        trajectory_identity=tr.identity,
        frame_content_identity=_digest(workflow["frames"]),
    )
    return {
        "schema": SCHEMA,
        "bindings": bindings,
        "selection": tr.selection.to_dict(),
        "observables": observables,
        "diagnostics": diagnostics,
        "frame_range": frame_range,
        "engine": workflow["manifest"]["engine"],
        "analysis_engine": {"name": "MDAnalysis", "version": mda.__version__},
        "execution": execution,
        "setup": fixture["setup"],
        **FLAGS,
    }


def save_bundle_record(root):
    root = Path(root)
    require(not (root / "analysis.json").exists(), "Analysis destination exists")
    result = analyze(root)
    result["input_files"] = files(root)
    result = _signed(result)
    import os
    import tempfile

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", dir=root, prefix=".analysis-", delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(
                json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, root / "analysis.json")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return result


def load_bundle(root, *, expected_identity=None):
    root = Path(root)
    record = json.loads((root / "analysis.json").read_text())
    require(
        record.get("identity") == _digest(_body(record)), "Analysis checksum mismatch"
    )
    if expected_identity is not None:
        require(record["identity"] == expected_identity, "Unexpected analysis identity")
    require(
        record["input_files"] == files(root), "Bundle input hashes/inventory changed"
    )
    result = analyze(root)
    require(
        result
        == {k: v for k, v in record.items() if k not in ("identity", "input_files")},
        "Analysis recomputation/bindings mismatch",
    )
    return record
