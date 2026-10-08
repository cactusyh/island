"""Corrective J20 acceptance: cumulative usage, membership and bounded chemistry."""

from collections import Counter
from copy import deepcopy
from dataclasses import replace

import pytest
from test_phase_4j20_crosslinking import _fixture, _rule

from island.core import AtomSite, Coordinates, MolecularSystem, Topology
from island.crosslinking import (
    CrosslinkNetworkPlan,
    ReactiveSite,
    apply_crosslink_plan,
    generate_crosslink_network,
    identify_reactive_sites,
    load_crosslink_plan,
    plan_crosslinks,
    save_crosslink_plan,
)
from island.crosslinking import _digest as digest
from island.exceptions import ValidationError
from island.graph import (
    FinalChemicalGraph,
    final_graph,
    load_final_graph_bundle,
    save_final_graph_bundle,
)
from island.graph.final import _chemical_payload, pack


def plan(graph, rule=None, seed=2026, count=1, **kwargs):
    rule = rule or _rule()
    return plan_crosslinks(
        graph,
        identify_reactive_sites(graph, rule),
        target_crosslinks=count,
        seed=seed,
        provenance="corrective software fixture",
        evidence=("J20 correction",),
        **kwargs,
    )


def resign_graph(payload):
    for tr in payload["transformations"]:
        tr["identity"] = digest({k: v for k, v in tr.items() if k != "identity"})
    payload["graph_identity"] = digest(_chemical_payload(payload))
    payload["history_identity"] = digest(payload["transformations"])
    return FinalChemicalGraph(pack(payload))


def resign_plan(payload):
    payload["identity"] = digest({k: v for k, v in payload.items() if k != "identity"})
    return CrosslinkNetworkPlan(pack(payload))


def molecule(element, hydrogens, charge=0, metadata=None):
    t = Topology()
    for i in (1, 10):
        t.add_site(
            AtomSite(
                i,
                element,
                14,
                element=element,
                formal_charge=charge,
                metadata={"reactive": True, **(metadata or {})},
            )
        )
        for h in range(i + 1, i + hydrogens + 1):
            t.add_site(AtomSite(h, "H", 1, element="H"))
            t.add_bond(i, h, order=1)
    system = MolecularSystem(t, Coordinates({i: (i, 0, 0) for i in t.sites}))
    return system, final_graph(system)


def test_three_rounds_reload_enforce_cumulative_capacity_and_preserve_inputs(tmp_path):
    system, graph = _fixture(4)
    usage = Counter()
    snapshots = []
    for round_id, seed in enumerate((2026, 0, 1)):
        original = deepcopy(system.to_dict())
        raw = graph.json_text
        p = plan(graph, seed=seed)
        if round_id == 0:
            assert p.payload["selected_crosslink_bonds"] == [[1, 7]]
        else:
            assert not any(
                i in (1, 7)
                for pair in p.payload["selected_crosslink_bonds"]
                for i in pair
            )
        path = save_crosslink_plan(p, tmp_path / f"plan-{round_id}.json")
        bundle = save_final_graph_bundle(system, graph, tmp_path / f"bundle-{round_id}")
        moved = tmp_path / f"relocated-{round_id}"
        bundle.rename(moved)
        loaded_system, loaded_graph, _ = load_final_graph_bundle(moved)
        loaded_plan = load_crosslink_plan(
            path, loaded_graph, expected_identity=p.identity
        )
        snapshots.append((system, original, graph, raw))
        system, graph, tr = apply_crosslink_plan(
            loaded_system, loaded_graph, loaded_plan
        )
        assert tr.payload["parameters"]["plan_identity"] == p.identity
        usage.update(i for pair in p.payload["selected_crosslink_bonds"] for i in pair)
        assert max(usage.values()) == 1
        observed = identify_reactive_sites(graph, _rule())
        for site in observed:
            assert site.lifetime_capacity == 1
            assert site.consumed_crosslinks == usage[site.site_id]
            assert (
                site.remaining_capacity
                == site.available_capacity
                == 1 - usage[site.site_id]
            )
        for old_system, old_data, old_graph, old_raw in snapshots:
            assert old_system.to_dict() == old_data
            assert old_graph.json_text == old_raw
    before = deepcopy(system.to_dict())
    with pytest.raises(ValidationError, match="exceeds feasible"):
        generate_crosslink_network(
            system,
            graph,
            _rule(),
            target_crosslinks=2,
            provenance="failure",
            evidence=(),
        )
    assert system.to_dict() == before


def test_generate_two_rounds_and_lifetime_change_rejection():
    system, graph = _fixture(4)
    first = generate_crosslink_network(
        system, graph, _rule(), target_crosslinks=1, provenance="fixture", evidence=()
    )
    original = deepcopy(first.system.to_dict())
    second = generate_crosslink_network(
        first.system,
        first.graph,
        _rule(),
        target_crosslinks=1,
        seed=0,
        provenance="fixture",
        evidence=(),
    )
    assert second.plan.payload["selected_crosslink_bonds"] != [[3, 7]]
    assert first.system.to_dict() == original
    with pytest.raises(ValidationError, match="Lifetime capacity changed"):
        identify_reactive_sites(
            first.graph, replace(_rule(), maximum_crosslinks_per_site=2)
        )


def test_optional_batch_limit_does_not_replace_lifetime_limit():
    system, graph = _fixture(4)
    rule = replace(_rule(), maximum_crosslinks_per_site=2)
    p = plan(graph, rule=rule, count=3, maximum_crosslinks_per_site=None)
    usage = Counter(i for pair in p.payload["selected_crosslink_bonds"] for i in pair)
    assert max(usage.values()) <= 2
    _, output, _ = apply_crosslink_plan(system, graph, p)
    for site in identify_reactive_sites(output, rule):
        assert site.consumed_crosslinks == usage[site.site_id]
        assert site.remaining_capacity == 2 - usage[site.site_id]


@pytest.mark.parametrize("custom", [False, True])
def test_original_molecule_membership_survives_every_graph_construction_and_relocation(
    tmp_path, custom, monkeypatch
):
    from island import crosslinking

    system, _ = _fixture()
    assert all("chain_id" not in a.metadata for a in system.topology.sites.values())
    labels = {i: "original-A" if i < 3 else "original-B" for i in system.topology.sites}
    graph = final_graph(system, molecule_membership=labels if custom else None)
    original = graph.payload["molecule_membership"]
    calls = []
    original_final_graph = crosslinking.final_graph

    def check_membership(*args, **kwargs):
        assert kwargs["molecule_membership"] == {
            int(i): label for i, label in original.items()
        }
        calls.append(kwargs)
        return original_final_graph(*args, **kwargs)

    monkeypatch.setattr(crosslinking, "final_graph", check_membership)
    output_system, output, tr = apply_crosslink_plan(system, graph, plan(graph))
    assert len(calls) == 3
    assert output.payload["molecule_membership"] == original
    assert len(output.payload["components"]) == 1
    assert len(set(output.payload["component_membership"].values())) == 1
    assert tr.payload["output_graph_identity"] == output.identity
    bundle = save_final_graph_bundle(output_system, output, tmp_path / "bundle")
    moved = tmp_path / "relocated"
    bundle.rename(moved)
    loaded_system, loaded, _ = load_final_graph_bundle(moved)
    loaded.validate_integrity(loaded_system)
    assert loaded.payload["molecule_membership"] == original
    assert (
        loaded.payload["component_membership"] == output.payload["component_membership"]
    )


@pytest.mark.parametrize(
    "element,hs,charge,metadata,reason",
    [
        ("N", 3, 0, {}, "Unsupported reactive element/charge"),
        ("N", 4, 0, {}, "Invalid charge-aware valence"),
        ("N", 4, 1, {}, "Unsupported reactive element/charge"),
        ("C", 4, 0, {}, "capacity"),
        ("O", 2, 0, {}, "capacity"),
        ("C", 3, 1, {}, "Unsupported element/charge"),
        ("O", 1, -1, {}, "Unsupported element/charge"),
        ("C", 3, 0, {"aromatic": True}, "Aromatic"),
        ("C", 3, 0, {"radical_electrons": 1}, "radical"),
        ("C", 3, 0, {"radical_electrons": True}, "radical"),
        ("C", 3, 0, {"radical": True}, "radical"),
    ],
)
def test_bounded_chemical_controls_preserve_input(
    element, hs, charge, metadata, reason
):
    system, graph = molecule(element, hs, charge, metadata)
    rule = replace(_rule(), allowed_elements=(element,))
    before = deepcopy(system.to_dict())
    with pytest.raises(ValidationError, match=reason):
        identify_reactive_sites(graph, rule)
    with pytest.raises(ValidationError, match=reason):
        generate_crosslink_network(
            system, graph, rule, target_crosslinks=1, provenance="fixture", evidence=()
        )
    assert system.to_dict() == before


def test_proposed_output_state_validated_before_publication():
    # Directly forged N sites cannot make a neutral four-coordinate N legal.
    _, graph = molecule("N", 3)
    rule = replace(_rule(), allowed_elements=("N",))
    sites = tuple(
        ReactiveSite(
            i,
            rule.name,
            "N",
            digest(c),
            3,
            1,
            rule.identity,
            rule.provenance,
            rule.evidence,
            1,
            0,
            rule,
        )
        for c in graph.payload["components"]
        for i in c
        if i in (1, 10)
    )
    with pytest.raises(ValidationError, match="Unsupported reactive element/charge"):
        plan_crosslinks(
            graph, sites, target_crosslinks=1, provenance="fixture", evidence=()
        )


@pytest.mark.parametrize("hydrogens,charge", [(3, 0), (4, 1)])
def test_valid_nitrogen_spectators_keep_charge_and_bond_orders(hydrogens, charge):
    system, _ = _fixture()
    nitrogen_system, _ = molecule("N", hydrogens, charge)
    for i, atom in nitrogen_system.topology.sites.items():
        copied = deepcopy(atom)
        copied.id += 100
        copied.metadata.pop("reactive", None)
        system.topology.add_site(copied)
        system.coordinates.set(copied.id, nitrogen_system.coordinates.get(i))
    for bond in nitrogen_system.topology.bonds.values():
        system.topology.add_bond(bond.site1 + 100, bond.site2 + 100, order=bond.order)
    graph = final_graph(system)
    output, _, _ = apply_crosslink_plan(system, graph, plan(graph))
    for i in (101, 110):
        assert output.topology.sites[i].formal_charge == charge
        assert output.topology.degree(i) == hydrogens


@pytest.mark.parametrize(
    "defect",
    [
        "missing",
        "duplicate",
        "seed",
        "plan",
        "bond",
        "limit",
        "transform",
        "missing_limit",
        "partial_limit",
    ],
)
def test_resigned_provenance_must_reconcile_with_actual_bonds_and_history(defect):
    system, graph = _fixture(4)
    _, output, _ = apply_crosslink_plan(system, graph, plan(graph))
    p = output.payload
    if defect == "missing":
        p["crosslink_provenance"] = []
    elif defect == "duplicate":
        p["crosslink_provenance"] *= 2
    elif defect == "seed":
        p["crosslink_provenance"][0]["seed"] += 1
    elif defect == "plan":
        p["crosslink_provenance"][0]["plan_identity"] = "a" * 64
    elif defect == "bond":
        p["transformations"][-1]["parameters"]["bond_order"] = 2
    elif defect == "limit":
        p["transformations"][-1]["parameters"]["site_lifetime_capacities"]["1"] = 0
    elif defect == "missing_limit":
        p["transformations"][-1]["parameters"].pop("site_lifetime_capacities")
    elif defect == "partial_limit":
        p["transformations"][-1]["parameters"]["site_lifetime_capacities"].pop("1")
    else:
        p["transformations"] = []
    tampered = resign_graph(p)
    tampered.validate_integrity()
    with pytest.raises(ValidationError):
        identify_reactive_sites(tampered, _rule())


def test_resigned_plan_nested_capacity_and_policy_tampering_rejects_load_and_apply(
    tmp_path,
):
    system, graph = _fixture(4)
    system, graph, _ = apply_crosslink_plan(system, graph, plan(graph))
    p = plan(graph, seed=0)
    before = deepcopy(system.to_dict())
    for defect in ("capacity", "policy", "history"):
        payload = p.payload
        if defect == "capacity":
            for site in payload["reactive_sites"]:
                if site["site_id"] == 7:
                    site["consumed_crosslinks"] = 0
                    site["remaining_capacity"] = site["available_capacity"] = 1
            payload["selected_reactive_site_identities"] = [
                digest(row) for row in payload["reactive_sites"]
            ]
        elif defect == "policy":
            payload["structural_policy"] = "unbounded"
        else:
            payload["input_history_identity"] = "0" * 64
        bad = resign_plan(payload)
        path = tmp_path / f"{defect}.json"
        path.write_text(bad.json_text)
        with pytest.raises(ValidationError):
            bad.validate_integrity(graph)
        with pytest.raises(ValidationError):
            load_crosslink_plan(path, graph)
        with pytest.raises(ValidationError):
            apply_crosslink_plan(system, graph, bad)
        assert system.to_dict() == before
