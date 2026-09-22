import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.forcefields import ParameterizedSystem


def test_chemical_site_has_formal_charge_but_no_partial_charge_field() -> None:
    atom = AtomSite(
        id=10, name="N", mass=14.007, element="N", atomic_number=7, formal_charge=1
    )
    assert atom.formal_charge == 1
    assert not hasattr(atom, "charge")
    assert not hasattr(atom, "type_name")


def test_formal_charge_must_be_an_integer() -> None:
    with pytest.raises(ValueError, match="formal charge"):
        AtomSite(
            id=10,
            name="N",
            mass=14.007,
            element="N",
            atomic_number=7,
            formal_charge=0.5,
        )


def test_future_partial_charge_lives_in_parameterized_system() -> None:
    topology = Topology()
    atom = AtomSite(
        id=35, name="O", mass=15.999, element="O", atomic_number=8, formal_charge=-1
    )
    topology.add_site(atom)
    system = MolecularSystem(topology, Coordinates({35: [0, 0, 0]}))
    parameterized = ParameterizedSystem(
        system,
        backend_name="synthetic",
        site_assignments={
            35: {
                "atom_type": "future-type",
                "partial_charge": -0.75,
                "nonbonded": {"sigma": 1.0},
            }
        },
    )
    assert parameterized.site_assignments[35]["partial_charge"] == -0.75
    assert atom.formal_charge == -1
    assert not hasattr(atom, "charge")
