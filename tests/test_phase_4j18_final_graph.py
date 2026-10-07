"""J18 final graph contracts; no scientific crosslinked/periodic claim."""

from copy import deepcopy
from math import inf, nan

import pytest

from island.exceptions import ValidationError
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    assign_typed_pcff_charges,
    bind_pcff_types,
    load_pcff_source,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
from island.graph import (
    FinalChemicalGraph,
    FinalGraphChargeRecord,
    build_psmiles_graph,
    complete_end_groups,
    construct_periodic_box,
    create_crosslinks,
    force_neutral_charges,
    load_final_graph_bundle,
    save_final_graph_bundle,
)
from island.graph.final import pack


def test_dp3_dp10_deterministic_final_graph():
    a = build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)
    b = build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)
    c = build_psmiles_graph("[*:1]CC[*:2]", dp=10, random_seed=2026)
    assert a[1].identity == b[1].identity
    assert len(a[1].payload["sites"]) == 20
    assert len(c[1].payload["sites"]) == 62
    assert (
        a[1].payload["history_identity"] == a[2].payload["identity"]
        or a[2].payload["output_graph_identity"] == a[1].identity
    )
    assert a[1].payload["repeat_provenance"]["degree_of_polymerization"] == 3


def test_transformations_are_immutable_and_history_bound():
    system, graph, _generation = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    completed_system, completed, completion = complete_end_groups(
        system, bonds=(), seed=4
    )
    assert completed.identity == graph.identity
    assert completed.payload["history_identity"] != graph.payload["history_identity"]
    assert completion.payload["input_graph_identity"] == graph.identity
    tampered = completed.payload
    tampered["transformations"][0]["parameters"]["dp"] = 99
    with pytest.raises(ValidationError, match="history identity"):
        FinalChemicalGraph(pack(tampered)).validate_integrity()
    assert system.to_dict() != completed_system.to_dict()
    assert completed_system.metadata["final_graph_transformations"]


def _crosslink_fixture():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    # Explicit end completion removes one terminal H at each target before the
    # caller-authorized crosslink bond is added; no implicit valence repair occurs.
    carbons = [i for i, a in system.topology.sites.items() if a.element == "C"]
    targets = (carbons[0], carbons[-1])
    for c in targets:
        h = next(
            i
            for i in system.topology.neighbors(c)
            if system.topology.sites[i].element == "H"
        )
        system.topology.remove_site(h)
        system.coordinates._positions.pop(h)
    return system, graph, targets


def test_crosslink_and_periodic_records():
    system, _, targets = _crosslink_fixture()
    system, graph, _crosslink = create_crosslinks(system, [targets], seed=2026)
    assert graph.payload["crosslink_provenance"]
    _periodic_system, periodic, box = construct_periodic_box(
        system, (40, 40, 40), seed=2026
    )
    assert periodic.periodic and periodic.payload["periodic"]["boundary"] == ["p"] * 3
    assert box.payload["output_graph_identity"] == periodic.identity
    with pytest.raises(ValidationError, match="Duplicate crosslink"):
        create_crosslinks(system, [targets, targets])
    with pytest.raises(ValidationError, match="existing bond"):
        create_crosslinks(system, [tuple(next(iter(system.topology.bonds)))])
    bad = deepcopy(periodic.payload)
    bad["periodic"]["lengths"][0] = 0
    with pytest.raises(ValidationError):
        FinalChemicalGraph(pack(bad)).validate_integrity()


def test_invalid_valence_and_graph_ids_are_rejected():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    bad = graph.payload
    bad["sites"][0]["id"] = 999
    with pytest.raises(ValidationError):
        FinalChemicalGraph(pack(bad)).validate_integrity()
    bad = graph.payload
    bad["bonds"].append({"sites": [1, 999], "order": 1, "aromatic": False})
    bad["graph_identity"] = "0" * 64
    with pytest.raises(ValidationError):
        FinalChemicalGraph(pack(bad)).validate_integrity()
    with pytest.raises(ValidationError, match="valence"):
        create_crosslinks(system, [(1, 4)])


def test_neutral_charge_contract_exact_and_nonfinite():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    values = {s["id"]: 0.0 for s in graph.payload["sites"]}
    source = {"family": "test", "sha256": "a" * 64}
    q = force_neutral_charges(
        graph,
        values,
        force_field="PCFF",
        source=source,
        method="provided",
        origin="provided",
        provenance="software fixture",
        evidence=["fixture"],
    )
    q.validate_integrity(graph, system=system)
    assert q.payload["total_charge"] == 0
    for value in (nan, inf):
        bad = dict(values)
        bad[1] = value
        with pytest.raises(ValidationError, match="finite"):
            force_neutral_charges(
                graph,
                bad,
                force_field="PCFF",
                source=source,
                method="provided",
                origin="provided",
                provenance="x",
                evidence=["x"],
            )
    bad = q.payload
    bad["charges"].pop("1")
    with pytest.raises(ValidationError):
        FinalGraphChargeRecord(pack(bad)).validate_integrity(graph)


def test_periodic_bundle_relocation_and_tampering(tmp_path):
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    system, graph, _ = construct_periodic_box(system, (25, 26, 27), seed=77)
    target = tmp_path / "bundle"
    save_final_graph_bundle(system, graph, target)
    moved = tmp_path / "moved"
    target.rename(moved)
    loaded_system, loaded_graph, _charges = load_final_graph_bundle(moved)
    assert loaded_graph.identity == graph.identity and loaded_system.box.lengths == (
        25,
        26,
        27,
    )
    raw = (moved / "final-graph.json").read_text()
    data = graph.payload
    data["molecule_membership"]["1"] = "changed"
    (moved / "final-graph.json").write_text(pack(data))
    with pytest.raises(Exception, match="checksum"):
        load_final_graph_bundle(moved)
    (moved / "final-graph.json").write_text(raw)
    data = graph.payload
    data["periodic"]["lengths"][0] = 99
    (moved / "final-graph.json").write_text(pack(data))
    import json

    manifest = json.loads((moved / "manifest.json").read_text())
    manifest["files"]["final-graph.json"]["sha256"] = __import__(
        "island.workflows.storage", fromlist=["checksum"]
    ).checksum((moved / "final-graph.json").read_bytes())
    (moved / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValidationError, match="identity|checksum"):
        load_final_graph_bundle(moved)


def test_periodic_forcefield_rejection_is_explicit():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    system, graph, _ = construct_periodic_box(system, (25, 25, 25))
    with pytest.raises(PreparedForceFieldError, match="does not support periodic"):
        prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                PCFFOptions("missing.frc", (0, 0, 1), (0, 0, 1)),
                final_graph=graph,
            ),
        )


def test_pcff_j17_external_records_bind_to_nonperiodic_graph():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    source = load_pcff_source("../lammps/lammps/tools/msi2lmp/frc_files/pcff.frc")
    labels = {
        i: {"C": "c", "H": "h"}[a.element] for i, a in system.topology.sites.items()
    }
    typed = bind_pcff_types(
        system,
        source,
        labels,
        provenance="J18 final graph fixture",
        evidence_references=["phase_4j18 fixture"],
    )
    charges = assign_typed_pcff_charges(
        system, typed, resolution_policy=COMPATIBILITY_POLICY
    )
    assert charges.complete and typed.payload["graph_identity"]
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                "../lammps/lammps/tools/msi2lmp/frc_files/pcff.frc",
                (0, 0, 1),
                (0, 0, 1),
                resolution_policy=COMPATIBILITY_POLICY,
                typed_graph=typed,
                graph_charges=charges,
            ),
            final_graph=graph,
        ),
    )
    assert prepared.metadata["production_validated"] is False
