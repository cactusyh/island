"""Force-field-independent reaction descriptions."""

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class FormBond:
    site1: int
    site2: int
    order: float | None = None


@dataclass(frozen=True)
class BreakBond:
    site1: int
    site2: int


@dataclass(frozen=True)
class ChangeBondOrder:
    site1: int
    site2: int
    order: float | None


@dataclass
class ReactionTemplate:
    """A named sequence of topology-only transformations."""

    name: str
    transformations: tuple[FormBond | BreakBond | ChangeBondOrder, ...]
    metadata: dict[str, Any] = field(default_factory=dict)
