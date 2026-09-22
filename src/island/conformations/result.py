"""Structured results from conformation generators."""

from dataclasses import dataclass, field
from typing import Any

from island.core import Coordinates, MolecularSystem


@dataclass(frozen=True)
class ConformationResult:
    """Generated coordinates together with reproducibility diagnostics."""

    coordinates: Coordinates
    method: str
    seed: int
    success: bool
    attempts: int
    rejected_trials: int
    rollback_count: int
    minimum_nonbonded_distance: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def apply_to(
        self, system: MolecularSystem, *, copy: bool = True
    ) -> MolecularSystem:
        """Apply coordinates to a copy by default, leaving topology untouched."""
        target = system.copy() if copy else system
        self.coordinates.validate(target.topology)
        target.coordinates = self.coordinates.copy()
        return target
