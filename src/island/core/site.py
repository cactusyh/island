"""Site models independent of coordinates and force fields."""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Site:
    """A stable-identity molecular site."""

    id: int
    name: str
    mass: float
    charge: float = 0.0
    type_name: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.id, int) or isinstance(self.id, bool):
            raise ValueError("Site id must be an integer")
        if self.mass < 0:
            raise ValueError("Site mass cannot be negative")


@dataclass
class AtomSite(Site):
    """An atomistic site."""

    element: str = ""
    atomic_number: int | None = None


@dataclass
class BeadSite(Site):
    """A coarse-grained site, optionally mapped to atomistic site IDs."""

    bead_type: str = ""
    mapped_atom_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        super().__post_init__()
        self.mapped_atom_ids = tuple(self.mapped_atom_ids)
