"""Native final-graph preparation and independent checks; no parameter fallback."""

import shutil
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np

from island.forcefields import (
    AmberToolsOptions,
    ForceFieldRequest,
    OPLSOptions,
    PCFFOptions,
    prepare_forcefield,
)
from island.graph import final_graph
from island.periodic import (
    PeriodicForceFieldConfig,
    periodic_chemical_view,
    prepare_periodic_forcefield,
)
from scripts.phase_4j31_contract import OPLS_SHA, PCFF_SHA
from scripts.phase_4j32_real_fixture import AMBER, FRC, XML, require


def checksum(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def prepare(fixture, family, root):
    system = periodic_chemical_view(fixture.packed.system)
    graph = final_graph(
        system,
        molecule_membership={
            int(k): v
            for k, v in fixture.packed.graph.payload["molecule_membership"].items()
        },
    )
    if family == "oplsaa":
        require(checksum(XML) == OPLS_SHA, "Pinned OPLS source changed")
        options = OPLSOptions(XML)
    elif family == "pcff":
        require(checksum(FRC) == PCFF_SHA, "Pinned PCFF source changed")
        options = PCFFOptions(FRC, (0, 0, 1), (0, 0, 1))
    else:
        options = AmberToolsOptions(
            family,
            "am1bcc",
            amberhome=AMBER,
            work_root=root,
            retain_success_artifacts=True,
            max_atoms=100,
            timeout_seconds=1800,
        )
    prepared = prepare_forcefield(
        system, ForceFieldRequest(family, options, final_graph=graph)
    )
    native = prepared.native_result
    if family == "oplsaa":
        source = {**prepared.source.identity, "sha256": checksum(XML)}
        mixing, scales = "geometric_12_6", (0.5, 0.5)
    elif family == "pcff":
        raise ValueError(
            "Unexpected PCFF success: frozen cyclic negative control requires typing rejection"
        )
    else:
        require(
            native.record["parmchk2_estimation_lines"] == [],
            "Estimated Amber interactions rejected",
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
        nonbonded_cutoff=5.0,
        lj_mixing_rule=mixing,
        one_four_scaling=scales,
        electrostatics_method="ewald",
        pme_tolerance=1e-5,
        provenance="J32 native reassignment of this complete packed graph",
        evidence=(fixture.packed.graph.identity, fixture.packed.plan.identity),
    )
    backend = prepare_periodic_forcefield(
        fixture.packed.system, fixture.packed.graph, prepared, config
    )
    assets = root / "source-assets"
    assets.mkdir()
    if family == "oplsaa":
        shutil.copyfile(XML, assets / "oplsaa.xml")
    else:
        shutil.copytree(Path(native.record["artifact_dir"]), assets / "amber")
    return backend


def independent_check(backend):
    p = backend.payload
    data = p["numerical_data"]
    graph = p["final_graph"]
    adj = {a["id"]: set() for a in graph["sites"]}
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        adj[a].add(b)
        adj[b].add(a)
    canon = lambda x: min(tuple(x), tuple(reversed(x)))
    bonds = {canon(b["sites"]) for b in graph["bonds"]}
    angles = {canon((a, b, c)) for b in adj for a, c in combinations(sorted(adj[b]), 2)}
    torsions = {
        canon((a, b, c, d))
        for b, c in bonds
        for a in adj[b] - {c}
        for d in adj[c] - {b}
        if a != d
    }
    for key, expected in [("bonds", bonds), ("angles", angles), ("torsions", torsions)]:
        actual = [canon(row["sites"]) for row in data[key]]
        require(
            set(actual) == expected and len(actual) == len(expected),
            "Independent topology mismatch: " + key,
        )
    native = backend._prepared.native_result
    charge_tolerance = (
        1e-6
        if p["config"]["family"] == "OPLS-AA"
        else native.imported_result.charge_result.tolerance
    )
    charges = [s["charge"] for s in data["sites"]]
    require(
        np.isfinite(charges).all() and abs(sum(charges)) <= charge_tolerance,
        "Nonneutral or nonfinite charges",
    )
    require(
        len(data["sites"]) == len(adj)
        and all(
            s["atom_type"] and s["sigma"] > 0 and s["epsilon"] >= 0
            for s in data["sites"]
        ),
        "Missing nonbonded terms",
    )
    native = backend._prepared.native_result
    report = {
        "bonds": len(bonds),
        "angles": len(angles),
        "proper_paths": len(torsions),
        "charge_total_e": sum(charges),
        "charge_total_tolerance_e": charge_tolerance,
        "atom_types_and_charges": data["sites"],
        "periodic_policy": p["config"],
        "topology_reference": "independent adjacency path enumeration",
    }
    if p["config"]["family"] == "OPLS-AA":
        report["independent_source_verification"] = native.payload[
            "upstream_verification"
        ]
    else:
        from island.forcefields import create_evaluator
        from scripts.validate_singlepoint_references import reference

        prepared = backend._prepared
        xyz = {i: prepared._system.coordinates.get(i).tolist() for i in adj}
        energy, forces, _ = reference(
            Path(native.record["artifact_dir"]) / "result.prmtop",
            dict(native.imported_result.mapping),
            xyz,
        )
        measured = create_evaluator(prepared._system, prepared).evaluate()
        error = max(
            float(np.max(np.abs(np.array(forces[i]) - measured.forces[i]))) for i in adj
        )
        require(
            np.isclose(energy, measured.potential_energy, rtol=1e-10, atol=1e-6)
            and error <= 1e-5,
            "Independent Amber energy/force mismatch",
        )
        report.update(
            energy_error_kj_mol=abs(energy - measured.potential_energy),
            max_force_error_kj_mol_angstrom=error,
            independent_source_verification="OpenMM AmberPrmtopFile vs ParmEd/ISLAND",
            parmchk2_estimation_lines=native.record["parmchk2_estimation_lines"],
        )
    return report
