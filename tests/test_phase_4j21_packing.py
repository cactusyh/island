"""J21 force-field-neutral multichain packing contracts."""

import subprocess
import sys
from copy import deepcopy
from math import isclose

import numpy as np
import pytest

from island.crosslinking import ReactiveSiteRule, generate_crosslink_network
from island.exceptions import ValidationError
from island.graph import build_psmiles_graph, final_graph
from island.graph.final import pack
from island.packing import (
    PeriodicPackingConfig,
    _digest,
    apply_periodic_packing_plan,
    load_periodic_packing_plan,
    pack_multichain_periodic,
    save_periodic_packing_plan,
)


def _chain(dp=1, seed=2026):
    return build_psmiles_graph("[*:1]CC[*:2]", dp=dp, random_seed=seed)


def _config(**kwargs):
    values = {
        "box_lengths": (40.0, 40.0, 40.0),
        "seed": 2026,
        "minimum_interunit_distance": 1.5,
        "maximum_attempts": 1000,
    }
    values.update(kwargs)
    return PeriodicPackingConfig(**values)


def _pack_two(**kwargs):
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    return pack_multichain_periodic(
        [("chain-0", s0, g0), ("chain-1", s1, g1)],
        config=_config(**kwargs),
        provenance="J21 test fixture",
        evidence=("controlled fixture",),
    )


def test_explicit_box_is_deterministic_and_remaps_sites_and_membership():
    first = _pack_two()
    second = _pack_two()
    assert first.plan.identity == second.plan.identity
    assert first.graph.identity == second.graph.identity
    assert first.transformation.identity == second.transformation.identity
    assert first.system.to_dict() == second.system.to_dict()
    payload = first.plan.payload
    assert payload["unit_labels"] == ["chain-0", "chain-1"]
    assert payload["site_id_remapping"] == [
        {str(i): i for i in range(1, 9)},
        {str(i): i + 8 for i in range(1, 9)},
    ]
    membership = first.graph.payload["molecule_membership"]
    assert set(membership.values()) == {"chain-0:A", "chain-1:A"}
    assert first.graph.payload["periodic"]["lengths"] == [40.0, 40.0, 40.0]
    assert first.graph.payload["periodic"]["boundary"] == ["p", "p", "p"]
    assert all(
        0 <= coordinate < 40
        for position in first.system.to_dict()["coordinates"].values()
        for coordinate in position
    )
    assert (
        first.graph.payload["transformations"][-1]["operation"]
        == "periodic_multichain_packing"
    )
    assert first.system.topology.angles


def test_target_density_constructs_orthorhombic_box():
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    result = pack_multichain_periodic(
        [("a", s0, g0), ("b", s1, g1)],
        config=PeriodicPackingConfig(
            target_density=0.9,
            periodic=(True, True, True),
            seed=7,
            rotation_mode="none",
            minimum_interunit_distance=0,
            maximum_attempts=5,
        ),
        provenance="density fixture",
        evidence=("controlled fixture",),
    )
    lengths = result.system.box.lengths
    assert lengths[0] == lengths[1] == lengths[2]
    mass = sum(site.mass for site in s0.topology.sites.values()) + sum(
        site.mass for site in s1.topology.sites.values()
    )
    density = mass / (6.02214076e23 * np.prod(lengths) * 1e-24)
    assert isclose(density, 0.9, rel_tol=0, abs_tol=1e-12)
    assert result.plan.payload["config"]["target_density"] == 0.9


def test_rejected_attempts_and_rotation_are_recorded():
    result = _pack_two(minimum_interunit_distance=1.5, maximum_attempts=1000)
    assert (
        result.plan.payload["rejected_attempt_count"]
        == _pack_two(minimum_interunit_distance=1.5).plan.payload[
            "rejected_attempt_count"
        ]
    )
    assert result.plan.payload["accepted_placements"][0]["rotation"]


def test_crosslinked_network_and_provenance_survive_packing():
    system, graph, _ = _chain(dp=3)
    carbon_ids = [i for i, atom in system.topology.sites.items() if atom.element == "C"]
    for carbon in (carbon_ids[0], carbon_ids[-1]):
        hydrogen = next(
            i
            for i in system.topology.neighbors(carbon)
            if system.topology.sites[i].element == "H"
        )
        system.topology.remove_site(hydrogen)
        system.coordinates._positions.pop(hydrogen)
        system.topology.sites[carbon].metadata["reactive"] = True
    graph = final_graph(system, transformations=graph.payload["transformations"])
    crosslinked = generate_crosslink_network(
        system,
        graph,
        ReactiveSiteRule(
            "C ends", "reactive", True, ("C",), 1, 1.0, "fixture", ("fixture",)
        ),
        target_crosslinks=1,
        require_intercomponent=False,
        provenance="fixture crosslink",
        evidence=("fixture",),
    )
    packed = pack_multichain_periodic(
        [("network", crosslinked.system, crosslinked.graph), ("free", *_chain()[:2])],
        config=_config(minimum_interunit_distance=0),
        provenance="mixed fixture",
        evidence=("crosslink preservation",),
    )
    assert packed.graph.payload["crosslink_provenance"]
    assert (
        packed.plan.payload["input_graph_identities"][0] == crosslinked.graph.identity
    )
    assert packed.graph.payload["molecule_membership"]


def test_save_load_relocate_and_separate_process_reconstruct(tmp_path):
    result = _pack_two()
    path = save_periodic_packing_plan(result.plan, tmp_path / "plan.json")
    moved = tmp_path / "relocated-plan.json"
    path.rename(moved)
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    units = [("chain-0", s0, g0), ("chain-1", s1, g1)]
    loaded = load_periodic_packing_plan(moved, units)
    reconstructed = apply_periodic_packing_plan(units, loaded)
    assert reconstructed.plan.identity == result.plan.identity
    assert reconstructed.graph.identity == result.graph.identity
    assert reconstructed.system.to_dict() == result.system.to_dict()
    script = "import sys; from island.packing import load_periodic_packing_plan; print(load_periodic_packing_plan(sys.argv[1]).identity)"
    child = subprocess.run(
        [sys.executable, "-c", script, str(moved)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert child.stdout.strip() == result.plan.identity


@pytest.mark.parametrize("defect", ["seed", "graph", "placement", "config", "outer"])
def test_tampered_plan_rejected(defect, tmp_path):
    result = _pack_two()
    payload = result.plan.payload
    if defect == "seed":
        payload["seed"] += 1
    elif defect == "graph":
        payload["input_graph_identities"][0] = "0" * 64
    elif defect == "placement":
        payload["accepted_placements"][0]["translation"][0] += 1
    elif defect == "config":
        payload["config"]["seed"] += 1
    else:
        payload["identity"] = "0" * 64
    path = tmp_path / f"{defect}.json"
    path.write_text(pack(payload))
    with pytest.raises(ValidationError):
        load_periodic_packing_plan(path)


def test_resigned_geometric_tampering_rejected_with_units(tmp_path):
    result = _pack_two()
    payload = result.plan.payload
    payload["accepted_placements"][0]["rotation"][0][0] = -1.0
    payload["identity"] = _digest(
        {key: value for key, value in payload.items() if key != "identity"}
    )
    path = tmp_path / "resigned.json"
    path.write_text(pack(payload))
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    with pytest.raises(ValidationError, match="rotation|coordinates"):
        load_periodic_packing_plan(path, [("chain-0", s0, g0), ("chain-1", s1, g1)])


def test_input_mutation_and_invalid_configs_are_rejected():
    result = _pack_two()
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    s1.topology.sites[1].formal_charge = 1
    with pytest.raises(ValidationError):
        result.plan.validate_integrity([("chain-0", s0, g0), ("chain-1", s1, g1)])
    for kwargs in (
        {"target_density": None, "box_lengths": None},
        {"box_lengths": (10, 10)},
        {"box_lengths": (10, -1, 10)},
        {"box_lengths": (10, 10, 10), "rotation_mode": "forcefield"},
    ):
        with pytest.raises(ValidationError):
            PeriodicPackingConfig(**kwargs)


def test_no_forcefield_backend_imported_by_packing():
    script = "import island.packing,sys; print([x for x in ('rdkit','openmm','parmed','scipy','foyer','island.forcefields.pcff','island.forcefields.oplsaa') if x in sys.modules])"
    child = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    )
    assert child.stdout.strip() == "[]"


def test_input_objects_unchanged_on_failure():
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    before = (
        deepcopy(s0.to_dict()),
        deepcopy(s1.to_dict()),
        g0.json_text,
        g1.json_text,
    )
    with pytest.raises(ValidationError):
        pack_multichain_periodic(
            [("a", s0, g0), ("b", s1, g1)],
            config=_config(box_lengths=(1, 1, 1), maximum_attempts=1),
        )
    assert (s0.to_dict(), s1.to_dict(), g0.json_text, g1.json_text) == before
