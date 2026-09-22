"""The central system aggregate."""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

from island.core.box import SimulationBox
from island.core.coordinates import Coordinates
from island.core.site import AtomSite, BeadSite, Site
from island.core.topology import Topology


@dataclass
class MolecularSystem:
    """A chemical topology plus coordinates and an optional simulation cell."""

    topology: Topology
    coordinates: Coordinates
    box: SimulationBox | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    representation: str = "atomistic"

    @property
    def number_of_sites(self) -> int:
        return len(self.topology.sites)

    @property
    def number_of_bonds(self) -> int:
        return len(self.topology.bonds)

    def validate(self) -> None:
        self.topology.validate()
        self.coordinates.validate(self.topology)

    def copy(self) -> "MolecularSystem":
        return deepcopy(self)

    def to_dict(self) -> dict[str, Any]:
        """Return a representation-complete, serialization-friendly mapping."""
        return {
            "representation": self.representation,
            "metadata": deepcopy(self.metadata),
            "box": None
            if self.box is None
            else {"lengths": self.box.lengths, "periodic": self.box.periodic},
            "sites": [_site_to_dict(site) for site in self.topology.sites.values()],
            "bonds": [
                {
                    "site1": bond.site1,
                    "site2": bond.site2,
                    "order": bond.order,
                    "aromatic": bond.aromatic,
                }
                for bond in self.topology.bonds.values()
            ],
            "coordinates": {
                site_id: self.coordinates.get(site_id).tolist()
                for site_id in self.topology.sites
            },
        }


def _site_to_dict(site: Site) -> dict[str, Any]:
    common = {
        "id": site.id,
        "name": site.name,
        "mass": site.mass,
        "metadata": deepcopy(site.metadata),
        "kind": type(site).__name__,
    }
    if isinstance(site, AtomSite):
        return {
            **common,
            "element": site.element,
            "atomic_number": site.atomic_number,
            "formal_charge": site.formal_charge,
        }
    if isinstance(site, BeadSite):
        return {
            **common,
            "bead_type": site.bead_type,
            "mapped_atom_ids": list(site.mapped_atom_ids),
        }
    return common
