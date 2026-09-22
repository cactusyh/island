import numpy as np
import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology, ValidationError


def make_system() -> MolecularSystem:
    topology = Topology()
    topology.add_site(AtomSite(id=100, name="C", mass=12.0, element="C"))
    topology.add_site(AtomSite(id=300, name="H", mass=1.0, element="H"))
    topology.add_bond(100, 300)
    return MolecularSystem(topology, Coordinates({100: [0, 0, 0], 300: [1, 0, 0]}))


def test_coordinates_are_id_keyed_and_translatable() -> None:
    system = make_system()
    system.coordinates.translate([1, 2, 3], site_ids=[300])
    assert np.allclose(system.coordinates.get(300), [2, 2, 3])
    assert np.allclose(system.coordinates.get(100), [0, 0, 0])


def test_coordinate_validation_checks_missing_and_extra_ids() -> None:
    topology = Topology()
    topology.add_site(AtomSite(id=7, name="H", mass=1, element="H"))
    with pytest.raises(ValidationError, match="missing"):
        Coordinates().validate(topology)
    with pytest.raises(ValidationError, match="unknown"):
        Coordinates({7: [0, 0, 0], 8: [0, 0, 0]}).validate(topology)


def test_system_validation_and_deep_copy() -> None:
    system = make_system()
    system.validate()
    copied = system.copy()
    copied.coordinates.set(100, [9, 9, 9])
    copied.topology.get_site(100).metadata["changed"] = True
    assert np.allclose(system.coordinates.get(100), [0, 0, 0])
    assert "changed" not in system.topology.get_site(100).metadata
