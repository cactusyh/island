"""Force-field-independent repeat-unit sequence representation."""

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from island.exceptions import PolymerBuildError


@dataclass(frozen=True)
class PolymerSequence:
    """An ordered sequence of repeat-unit identities and its generation metadata."""

    repeat_unit_types: tuple[str, ...]
    generation_metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        identities = tuple(self.repeat_unit_types)
        if not identities:
            raise PolymerBuildError("A polymer sequence cannot be empty")
        if any(
            not isinstance(identity, str) or not identity for identity in identities
        ):
            raise PolymerBuildError(
                "Every polymer sequence identity must be a non-empty string"
            )
        object.__setattr__(self, "repeat_unit_types", identities)
        object.__setattr__(self, "generation_metadata", dict(self.generation_metadata))

    @property
    def identities(self) -> tuple[str, ...]:
        return self.repeat_unit_types

    @property
    def number_of_repeat_units(self) -> int:
        return len(self.repeat_unit_types)

    @property
    def composition_counts(self) -> dict[str, int]:
        return dict(Counter(self.repeat_unit_types))
