"""The central system aggregate."""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

from island.core.box import SimulationBox
from island.core.coordinates import Coordinates
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
        return {
            "representation": self.representation,
            "metadata": deepcopy(self.metadata),
            "box": None
            if self.box is None
            else {"lengths": self.box.lengths, "periodic": self.box.periodic},
            "sites": [
                {
                    "id": s.id,
                    "name": s.name,
                    "mass": s.mass,
                    "charge": s.charge,
                    "type_name": s.type_name,
                    "metadata": deepcopy(s.metadata),
                    "kind": type(s).__name__,
                }
                for s in self.topology.sites.values()
            ],
            "bonds": [
                {
                    "site1": b.site1,
                    "site2": b.site2,
                    "order": b.order,
                    "aromatic": b.aromatic,
                }
                for b in self.topology.bonds.values()
            ],
            "coordinates": {
                site_id: self.coordinates.get(site_id).tolist()
                for site_id in self.topology.sites
            },
        }
