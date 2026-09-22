from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.reactions import BreakBond, FormBond, apply_transformation


def make_fragments() -> MolecularSystem:
    topology = Topology()
    for site_id in (1, 2, 3, 4):
        topology.add_site(AtomSite(id=site_id, name="C", mass=12.0, element="C"))
    topology.add_bond(1, 2)
    topology.add_bond(3, 4)
    return MolecularSystem(
        topology, Coordinates({i: [float(i), 0, 0] for i in range(1, 5)})
    )


def test_form_bond_joins_fragments_and_rebuilds_interactions() -> None:
    system = make_fragments()
    apply_transformation(system, FormBond(2, 3))
    assert len(system.topology.connected_components()) == 1
    assert len(system.topology.angles) == 2
    assert len(system.topology.dihedrals) == 1


def test_break_bond_splits_fragments_and_rebuilds_interactions() -> None:
    system = make_fragments()
    apply_transformation(system, FormBond(2, 3))
    apply_transformation(system, BreakBond(2, 3))
    assert len(system.topology.connected_components()) == 2
    assert not system.topology.angles
