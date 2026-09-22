"""Site models independent of coordinates and force fields."""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Site:
    """A stable-identity molecular site."""

    id: int
    name: str
    mass: float
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.id, int) or isinstance(self.id, bool):
            raise ValueError("Site id must be an integer")
        if self.mass < 0:
            raise ValueError("Site mass cannot be negative")


@dataclass
class AtomSite(Site):
    """An atomistic site with chemical identity and integer formal charge."""

    element: str = ""
    atomic_number: int | None = None
    formal_charge: int = 0

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.formal_charge, int) or isinstance(
            self.formal_charge, bool
        ):
            raise ValueError("Atom formal charge must be an integer")


@dataclass
class BeadSite(Site):
    """A coarse-grained site with representation-level bead identity.

    ``bead_type`` identifies the coarse-grained site kind; it is not an
    atomistic force-field type or parameter assignment.
    """

    bead_type: str = ""
    mapped_atom_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        super().__post_init__()
        self.mapped_atom_ids = tuple(self.mapped_atom_ids)
