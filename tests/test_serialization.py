from island import (
    AtomSite,
    BeadSite,
    Coordinates,
    MolecularSystem,
    SimulationBox,
    Topology,
)


def test_atom_site_serialization_is_representation_complete() -> None:
    topology = Topology()
    topology.add_site(
        AtomSite(
            id=10,
            name="O1",
            mass=15.999,
            element="O",
            atomic_number=8,
            formal_charge=-1,
            metadata={"aromatic": False},
        )
    )
    topology.add_site(
        AtomSite(id=35, name="C1", mass=12.011, element="C", atomic_number=6)
    )
    topology.add_bond(10, 35, order=1.5, aromatic=True)
    system = MolecularSystem(
        topology,
        Coordinates({10: [1, 2, 3], 35: [4, 5, 6]}),
        box=SimulationBox(10, 20, 30, periodic=(True, False, True)),
        metadata={"source": "synthetic"},
    )

    serialized = system.to_dict()
    assert serialized["sites"][0] == {
        "id": 10,
        "name": "O1",
        "mass": 15.999,
        "metadata": {"aromatic": False},
        "kind": "AtomSite",
        "element": "O",
        "atomic_number": 8,
        "formal_charge": -1,
    }
    assert "charge" not in serialized["sites"][0]
    assert "type_name" not in serialized["sites"][0]
    assert serialized["bonds"] == [
        {"site1": 10, "site2": 35, "order": 1.5, "aromatic": True}
    ]
    assert serialized["coordinates"] == {10: [1.0, 2.0, 3.0], 35: [4.0, 5.0, 6.0]}
    assert serialized["box"] == {
        "lengths": (10, 20, 30),
        "periodic": (True, False, True),
    }
    assert serialized["metadata"] == {"source": "synthetic"}


def test_bead_site_serialization_preserves_mapping_and_bead_identity() -> None:
    topology = Topology()
    topology.add_site(
        BeadSite(
            id=20,
            name="Backbone",
            mass=28.0,
            bead_type="BB",
            mapped_atom_ids=(101, 104),
            metadata={"resolution": "coarse"},
        )
    )
    system = MolecularSystem(
        topology, Coordinates({20: [0, 0, 0]}), representation="coarse_grained"
    )
    assert system.to_dict()["sites"] == [
        {
            "id": 20,
            "name": "Backbone",
            "mass": 28.0,
            "metadata": {"resolution": "coarse"},
            "kind": "BeadSite",
            "bead_type": "BB",
            "mapped_atom_ids": [101, 104],
        }
    ]
