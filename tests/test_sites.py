import pytest

from island import AtomSite, BeadSite, Site


def test_site_has_stable_id_and_metadata() -> None:
    site = Site(id=42, name="x", mass=1.0, metadata={"role": "end"})
    assert site.id == 42
    assert site.metadata["role"] == "end"


def test_site_id_must_be_integer() -> None:
    with pytest.raises(ValueError):
        Site(id=1.5, name="x", mass=1.0)


def test_atom_and_bead_specializations() -> None:
    atom = AtomSite(id=4, name="C", mass=12.011, element="C", atomic_number=6)
    bead = BeadSite(id=7, name="P", mass=44.0, bead_type="P1", mapped_atom_ids=[4, 5])
    assert atom.element == "C"
    assert bead.mapped_atom_ids == (4, 5)
