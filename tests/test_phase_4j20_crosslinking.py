"""J20 explicit network generation and integrity controls."""

import json
import subprocess
import sys
from copy import deepcopy
from hashlib import sha256
from math import inf, nan

import pytest

from island.core import AtomSite, Coordinates, MolecularSystem, Topology
from island.crosslinking import (
    CrosslinkNetworkPlan,
    ReactiveSiteRule,
    apply_crosslink_plan,
    generate_crosslink_network,
    identify_reactive_sites,
    load_crosslink_plan,
    plan_crosslinks,
    save_crosslink_plan,
)
from island.exceptions import ValidationError
from island.graph import (
    build_psmiles_graph,
    construct_periodic_box,
    final_graph,
    load_final_graph_bundle,
    save_final_graph_bundle,
)
from island.graph.final import pack


def _fixture(count=2, *, reverse=False):
    """Disconnected ethylene fragments, each carbon marked by the caller."""
    topology = Topology()
    positions = {}
    entries = []
    for fragment in range(count):
        a, b = fragment * 2 + 1, fragment * 2 + 2
        entries.extend([(a, fragment), (b, fragment)])
    for atom_id, fragment in reversed(entries) if reverse else entries:
        topology.add_site(
            AtomSite(
                atom_id,
                "C",
                12.0,
                metadata={"fragment": fragment, "reactive": True},
                element="C",
                atomic_number=6,
            )
        )
        positions[atom_id] = (float(atom_id), 0.0, 0.0)
    for fragment in range(count):
        topology.add_bond(fragment * 2 + 1, fragment * 2 + 2, order=2)
    system = MolecularSystem(topology, Coordinates(positions))
    return system, final_graph(system)


def _rule():
    return ReactiveSiteRule(
        "declared carbon",
        "reactive",
        True,
        ("C",),
        1,
        1.0,
        "caller-marked fixture",
        ("J20 fixture",),
    )


def _polymer_fragments():
    """Two independent DP3 polyethylene chains with explicit terminal H removal."""
    topology = Topology()
    positions = {}
    for fragment in range(2):
        chain, _, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)
        carbon = next(
            i for i, atom in chain.topology.sites.items() if atom.element == "C"
        )
        hydrogen = next(
            i
            for i in chain.topology.neighbors(carbon)
            if chain.topology.sites[i].element == "H"
        )
        chain.topology.remove_site(hydrogen)
        chain.coordinates._positions.pop(hydrogen)
        chain.topology.sites[carbon].metadata["reactive"] = True
        offset = fragment * 100
        for atom in chain.topology.sites.values():
            copied = deepcopy(atom)
            copied.id += offset
            copied.metadata["chain_id"] = str(fragment)
            topology.add_site(copied)
            x, y, z = chain.coordinates.get(atom.id)
            positions[copied.id] = (x + fragment * 30, y, z)
        for bond in chain.topology.bonds.values():
            topology.add_bond(
                bond.site1 + offset,
                bond.site2 + offset,
                order=bond.order,
                aromatic=bond.aromatic,
            )
    system = MolecularSystem(
        topology,
        Coordinates(positions),
        metadata={"polymer": {"fixture": "two independent DP3 chains"}},
    )
    return system, final_graph(system)


def _plan(graph, rule, *, count=1, seed=2026, **kwargs):
    return plan_crosslinks(
        graph,
        identify_reactive_sites(graph, rule),
        target_crosslinks=count,
        seed=seed,
        provenance="caller-selected network",
        evidence=("J20 fixture",),
        **kwargs,
    )


def test_one_crosslink_copies_system_and_binds_history():
    system, graph = _fixture()
    original = deepcopy(system.to_dict())
    rule = _rule()
    result = generate_crosslink_network(
        system,
        graph,
        rule,
        target_crosslinks=1,
        provenance="caller-selected network",
        evidence=("J20 fixture",),
    )
    assert len(result.reactive_sites) == 4
    assert len(result.plan.payload["candidate_pairs"]) == 4
    assert len(result.plan.payload["selected_crosslink_bonds"]) == 1
    assert len(result.graph.payload["components"]) == 1
    assert len(result.system.topology.bonds) == len(system.topology.bonds) + 1
    assert (
        result.transformation.payload["parameters"]["plan_identity"]
        == result.plan.identity
    )
    assert result.plan.identity in str(result.graph.payload["transformations"])
    assert system.to_dict() == original
    assert all(
        (result.system.coordinates.get(i) == system.coordinates.get(i)).all()
        for i in system.topology.sites
    )
    assert result.system.topology.angles and result.system.topology.dihedrals


def test_two_disconnected_polymer_fragments_crosslink():
    system, graph = _polymer_fragments()
    original = deepcopy(system.to_dict())
    result = generate_crosslink_network(
        system,
        graph,
        _rule(),
        target_crosslinks=1,
        seed=2026,
        provenance="caller-declared DP3 fragment link",
        evidence=("explicit terminal H removal",),
    )
    assert len(result.reactive_sites) == 2
    assert len(result.plan.payload["candidate_pairs"]) == 1
    assert len(result.graph.payload["components"]) == 1
    assert result.system.number_of_bonds == system.number_of_bonds + 1
    assert system.to_dict() == original


def test_multiple_deterministic_capacity_and_insertion_order():
    system, graph = _fixture(3)
    _, graph_reverse = _fixture(3, reverse=True)
    assert graph.identity == graph_reverse.identity
    rule = _rule()
    plan = _plan(graph, rule, count=2)
    assert plan.identity == _plan(graph_reverse, rule, count=2).identity
    assert plan.identity == _plan(graph, rule, count=2).identity
    assert (
        len({i for pair in plan.payload["selected_crosslink_bonds"] for i in pair}) == 4
    )
    assert (
        len(apply_crosslink_plan(system, graph, plan)[1].payload["transformations"])
        == 1
    )
    assert (
        plan.payload["selected_crosslink_bonds"]
        != _plan(graph, rule, count=2, seed=2027).payload["selected_crosslink_bonds"]
    )
    moved = system.copy()
    moved.coordinates.translate((10, 0, 0))
    moved_graph = final_graph(moved)
    assert moved_graph.identity == graph.identity
    assert _plan(moved_graph, rule, count=2).identity == plan.identity


def test_psmiles_dp3_requires_explicit_hydrogen_removal():
    system, original, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3)
    carbons = [i for i, atom in system.topology.sites.items() if atom.element == "C"]
    targets = (carbons[0], carbons[-1])
    for carbon in targets:
        system.topology.sites[carbon].metadata["reactive"] = True
    full = final_graph(system, transformations=original.payload["transformations"])
    with pytest.raises(ValidationError, match="capacity"):
        identify_reactive_sites(full, _rule())
    for carbon in targets:
        hydrogen = next(
            i
            for i in system.topology.neighbors(carbon)
            if system.topology.sites[i].element == "H"
        )
        system.topology.remove_site(hydrogen)
        system.coordinates._positions.pop(hydrogen)
    graph = final_graph(system, transformations=original.payload["transformations"])
    result = generate_crosslink_network(
        system,
        graph,
        _rule(),
        target_crosslinks=1,
        require_intercomponent=False,
        provenance="explicit intrachain fixture",
        evidence=("terminal H removed before planning",),
    )
    assert result.system.number_of_sites == system.number_of_sites
    assert result.system.number_of_bonds == system.number_of_bonds + 1
    assert result.graph.identity != original.identity


def test_invalid_inputs_and_periodic_rejection():
    system, graph = _fixture()
    rule = _rule()
    sites = identify_reactive_sites(graph, rule)
    for kwargs in (
        {"target_crosslinks": -1},
        {"target_crosslinks": True},
        {"target_crosslinks": 3},
        {"target_crosslinks": 1, "bond_order": nan},
        {"target_crosslinks": 1, "bond_order": inf},
    ):
        with pytest.raises(ValidationError):
            plan_crosslinks(
                graph, sites, provenance="fixture", evidence=("fixture",), **kwargs
            )
    with pytest.raises(ValidationError):
        plan_crosslinks(
            graph, (), target_crosslinks=1, provenance="fixture", evidence=("fixture",)
        )
    with pytest.raises(ValidationError):
        plan_crosslinks(
            graph,
            sites + (sites[0],),
            target_crosslinks=1,
            provenance="fixture",
            evidence=("fixture",),
        )
    with pytest.raises(ValidationError):
        plan_crosslinks(
            graph,
            sites[:-1],
            target_crosslinks=1,
            provenance="fixture",
            evidence=("fixture",),
        )
    periodic_system, periodic_graph, _ = construct_periodic_box(system, (20, 20, 20))
    with pytest.raises(ValidationError, match="Periodic"):
        identify_reactive_sites(periodic_graph, rule)
    assert periodic_system.box is not None
    bad_rule = ReactiveSiteRule("wrong", "reactive", True, ("N",), 1, 1, "fixture", ())
    with pytest.raises(ValidationError, match="element"):
        identify_reactive_sites(graph, bad_rule)


@pytest.mark.parametrize(
    "field,change",
    [
        (
            "candidate_pairs",
            lambda p: p["candidate_pairs"].append(p["candidate_pairs"][0]),
        ),
        ("candidate_pairs", lambda p: p["candidate_pairs"].append([1, 1])),
        ("candidate_pairs", lambda p: p["candidate_pairs"].append([1, 999])),
        ("candidate_pairs", lambda p: p["candidate_pairs"].append([1, 2])),
        (
            "selected_crosslink_bonds",
            lambda p: p["selected_crosslink_bonds"].append(
                p["selected_crosslink_bonds"][0]
            ),
        ),
        (
            "selected_crosslink_bonds",
            lambda p: p["selected_crosslink_bonds"].append([1, 1]),
        ),
        (
            "selected_crosslink_bonds",
            lambda p: p["selected_crosslink_bonds"].append([1, 999]),
        ),
        ("random_seed", lambda p: p.update(random_seed=777)),
        ("rule_identity", lambda p: p.update(rule_identity="0" * 64)),
        ("provenance", lambda p: p.update(provenance="altered")),
        ("evidence", lambda p: p.update(evidence=["altered"])),
    ],
)
def test_plan_tampering_rejected(field, change):
    _, graph = _fixture()
    plan = _plan(graph, _rule())
    payload = plan.payload
    change(payload)
    assert field in payload
    with pytest.raises(ValidationError):
        CrosslinkNetworkPlan(pack(payload)).validate_integrity(graph)


def test_resigned_candidate_tampering_still_rejects_graph_contradiction():
    from island.crosslinking import _pair_identity

    _, graph = _fixture()
    plan = _plan(graph, _rule())
    for altered in ([1, 2], [1, 1], [1, 999], [1, 3]):
        payload = plan.payload
        payload["candidate_pairs"].append(altered)
        payload["candidate_pairs"].sort()
        payload["candidate_pair_identities"] = [
            _pair_identity(
                graph.identity, payload["rule_identity"], pair, payload["bond_order"]
            )
            for pair in payload["candidate_pairs"]
        ]
        payload["identity"] = sha256(
            pack({k: v for k, v in payload.items() if k != "identity"}).encode()
        ).hexdigest()
        with pytest.raises(ValidationError):
            CrosslinkNetworkPlan(pack(payload)).validate_integrity(graph)


def test_stale_graph_missing_coordinates_and_source_metadata():
    system, graph = _fixture()
    plan = _plan(graph, _rule())
    changed = system.copy()
    changed.topology.sites[1].formal_charge = 1
    with pytest.raises(ValidationError, match="mismatch"):
        apply_crosslink_plan(changed, graph, plan)
    missing = system.copy()
    missing.coordinates._positions.pop(1)
    with pytest.raises(ValidationError, match="missing coordinates"):
        apply_crosslink_plan(missing, graph, plan)
    nonfinite = system.copy()
    nonfinite.coordinates._positions[1][0] = nan
    with pytest.raises(ValidationError, match="Nonfinite coordinates"):
        apply_crosslink_plan(nonfinite, graph, plan)
    with pytest.raises(ValidationError, match="force-field"):
        ReactiveSiteRule("bad", "force_field_type", "c3", ("C",), 1, 1, "fixture", ())


def test_plan_persistence_relocation_and_child_process(tmp_path):
    system, graph = _fixture()
    plan = _plan(graph, _rule())
    path = save_crosslink_plan(plan, tmp_path / "plan.json")
    bundle = save_final_graph_bundle(system, graph, tmp_path / "bundle")
    relocated = tmp_path / "moved"
    bundle.rename(relocated)
    _, loaded_graph, _ = load_final_graph_bundle(relocated)
    assert (
        load_crosslink_plan(
            path, loaded_graph, expected_identity=plan.identity
        ).identity
        == plan.identity
    )
    output_system, output_graph, _ = apply_crosslink_plan(system, graph, plan)
    output_bundle = save_final_graph_bundle(
        output_system, output_graph, tmp_path / "output_bundle"
    )
    moved_output = tmp_path / "moved_output"
    output_bundle.rename(moved_output)
    _, loaded_output, _ = load_final_graph_bundle(moved_output)
    assert loaded_output.identity == output_graph.identity
    assert plan.identity in str(loaded_output.payload["transformations"])
    script = "import sys; from island.graph import load_final_graph_bundle; from island.crosslinking import load_crosslink_plan; g=load_final_graph_bundle(sys.argv[1])[1]; print(load_crosslink_plan(sys.argv[2],g).identity)"
    child = subprocess.run(
        [sys.executable, "-c", script, str(relocated), str(path)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert child.stdout.strip() == plan.identity
    tampered = json.loads(path.read_text())
    tampered["payload"]["random_seed"] = 99
    path.write_text(json.dumps(tampered))
    with pytest.raises(ValidationError):
        load_crosslink_plan(path, loaded_graph)
    resigned = plan.payload
    resigned["provenance"] = "re-signed but untrusted"
    resigned["identity"] = sha256(
        pack({k: v for k, v in resigned.items() if k != "identity"}).encode()
    ).hexdigest()
    path.write_text(pack(resigned))
    with pytest.raises(ValidationError, match="trusted reference"):
        load_crosslink_plan(path, loaded_graph, expected_identity=plan.identity)


def test_crosslink_import_keeps_scientific_backends_unloaded():
    script = (
        "import island.crosslinking,sys; "
        "forbidden=('rdkit','openmm','scipy','foyer','parmed',"
        "'island.forcefields.pcff','island.forcefields.oplsaa',"
        "'island.forcefields.ambertools'); "
        "print([name for name in forbidden if name in sys.modules])"
    )
    child = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    )
    assert child.stdout.strip() == "[]"


def test_stale_plan_rejected_after_graph_change():
    system, graph = _fixture()
    plan = _plan(graph, _rule())
    changed = system.copy()
    changed.topology.sites[1].formal_charge = 1
    newer = final_graph(changed)
    with pytest.raises(ValidationError, match="Stale input graph"):
        plan.validate_integrity(newer)
