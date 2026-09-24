from dataclasses import replace
from math import sqrt

import pytest

from island import AtomSite, Topology
from island.builders import build_linear_polymer
from island.exceptions import NonbondedPolicyError
from island.forcefields import NonbondedPolicy
from island.forcefields.nonbonded import nonbonded_policy_signature
from island.forcefields.parameters import LennardJonesParameter

COMMON = {
    "source": "synthetic mixing test",
    "library_name": "synthetic",
    "library_version": "1",
}


def policy(mixing_rule: str = "lorentz_berthelot") -> NonbondedPolicy:
    return NonbondedPolicy(
        name="synthetic_policy",
        version="1",
        mixing_rule=mixing_rule,  # type: ignore[arg-type]
        lj_scale_12=0.0,
        coulomb_scale_12=0.0,
        lj_scale_13=0.0,
        coulomb_scale_13=0.0,
        lj_scale_14=0.5,
        coulomb_scale_14=0.833333333333,
        source="synthetic policy test",
    )


def lj(
    parameter_id: str, atom_type: str, epsilon: float, sigma: float
) -> LennardJonesParameter:
    return LennardJonesParameter(parameter_id, atom_type, epsilon, sigma, **COMMON)


def topology_with_edges(
    site_ids: tuple[int, ...], edges: tuple[tuple[int, int], ...]
) -> Topology:
    topology = Topology()
    for site_id in reversed(site_ids):
        topology.add_site(
            AtomSite(site_id, f"C{site_id}", 12.0, element="C", atomic_number=6)
        )
    for site1, site2 in reversed(edges):
        topology.add_bond(site2, site1, order=1)
    return topology


def test_lorentz_berthelot_and_geometric_mixing_with_unequal_values() -> None:
    first = lj("a", "A", 0.25, 0.30)
    second = lj("b", "B", 1.00, 0.50)
    lorentz = policy().mix_lj(first, second)
    geometric = policy("geometric").mix_lj(first, second)
    assert lorentz.epsilon == pytest.approx(0.5)
    assert lorentz.sigma == pytest.approx(0.4)
    assert geometric.epsilon == pytest.approx(0.5)
    assert geometric.sigma == pytest.approx(sqrt(0.15))
    assert lorentz.atom_types == ("A", "B")


@pytest.mark.parametrize("mixing_rule", ["lorentz_berthelot", "geometric"])
def test_zero_epsilon_is_preserved_by_mixing(mixing_rule: str) -> None:
    mixed = policy(mixing_rule).mix_lj(
        lj("zero", "A", 0.0, 0.3), lj("nonzero", "B", 1.0, 0.5)
    )
    assert mixed.epsilon == 0.0


def test_unsupported_policy_and_invalid_scales_are_rejected() -> None:
    with pytest.raises(NonbondedPolicyError, match="Unsupported"):
        policy("sixth_power")
    with pytest.raises(NonbondedPolicyError, match=r"\[0, 1\]"):
        replace(policy(), lj_scale_14=1.1)


def test_linear_chain_pair_classes_and_independent_scales() -> None:
    topology = topology_with_edges(
        (10, 20, 30, 40, 50),
        ((10, 20), (20, 30), (30, 40), (40, 50)),
    )
    selected = policy()
    assert selected.scaling_for_pair(topology, 10, 20).relationship == "1-2"
    assert selected.scaling_for_pair(topology, 10, 30).relationship == "1-3"
    one_four = selected.scaling_for_pair(topology, 10, 40)
    assert one_four.relationship == "1-4"
    assert one_four.lj_scale == 0.5
    assert one_four.coulomb_scale == pytest.approx(0.833333333333)
    full = selected.scaling_for_pair(topology, 10, 50)
    assert (full.relationship, full.lj_scale, full.coulomb_scale) == (
        "full",
        1.0,
        1.0,
    )


def test_branch_and_ring_use_unique_shortest_path_classification() -> None:
    branched = topology_with_edges((1, 2, 3, 4), ((1, 2), (1, 3), (1, 4)))
    assert policy().scaling_for_pair(branched, 2, 3).relationship == "1-3"

    square = topology_with_edges((1, 2, 3, 4), ((1, 2), (2, 3), (3, 4), (4, 1)))
    assert policy().scaling_for_pair(square, 1, 4).relationship == "1-2"
    assert policy().scaling_for_pair(square, 1, 3).relationship == "1-3"
    exceptions = policy().local_pair_scalings(square)
    assert len(exceptions) == 6
    assert len(exceptions) == len(set(exceptions))

    disconnected = topology_with_edges((10, 20), ())
    full = policy().scaling_for_pair(disconnected, 10, 20)
    assert (full.relationship, full.shortest_bond_distance) == ("full", None)


def test_noncontiguous_ids_and_insertion_order_are_irrelevant() -> None:
    first = topology_with_edges((10, 35, 90, 140), ((10, 35), (35, 90), (90, 140)))
    second = Topology()
    for site_id in (10, 35, 90, 140):
        second.add_site(AtomSite(site_id, "C", 12.0, element="C", atomic_number=6))
    for edge in ((10, 35), (35, 90), (90, 140)):
        second.add_bond(*edge, order=1)
    assert policy().local_pair_scalings(first) == policy().local_pair_scalings(second)


def test_policy_content_changes_signature() -> None:
    original = policy()
    changed = replace(original, coulomb_scale_14=0.75)
    assert nonbonded_policy_signature(original) != nonbonded_policy_signature(changed)


@pytest.mark.parametrize("dp", [50, 100])
def test_long_polyethylene_policy_uses_only_local_exceptions(dp: int) -> None:
    system = build_linear_polymer("[*]CC[*]", dp=dp, generate_3d=False)
    exceptions = policy().local_pair_scalings(system.topology)
    assert exceptions
    assert len(exceptions) < 30 * system.number_of_sites
    assert all(item.shortest_bond_distance in {1, 2, 3} for item in exceptions.values())
