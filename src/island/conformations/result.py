"""Structured results from conformation generators."""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

from island.core import Coordinates, MolecularSystem
from island.exceptions import UnsupportedConformationError


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
        """Apply coordinates and provenance, returning a copy by default.

        All validation and provenance construction happen before mutation, so an
        in-place failure cannot leave a partially updated target.
        """
        self.coordinates.validate(system.topology)
        if not self.success:
            raise UnsupportedConformationError(
                "Cannot apply an unsuccessful conformation result"
            )
        provenance = self.coordinate_provenance()
        new_coordinates = self.coordinates.copy()
        target = system.copy() if copy else system
        target.coordinates = new_coordinates
        polymer = target.metadata.get("polymer")
        if isinstance(polymer, dict):
            polymer["coordinates"] = provenance["coordinate_source"]
            polymer["coordinate_generation"] = provenance
        else:
            target.metadata["coordinate_generation"] = provenance
        return target

    def coordinate_provenance(self) -> dict[str, Any]:
        """Return normalized coordinate provenance for builder and direct use."""
        coordinate_source = self.metadata.get("coordinate_source", self.method)
        return {
            "method": self.method,
            "coordinate_source": coordinate_source,
            "seed": self.seed,
            "success": self.success,
            "attempts": self.attempts,
            "rejected_trials": self.rejected_trials,
            "rollback_count": self.rollback_count,
            "minimum_nonbonded_distance": self.minimum_nonbonded_distance,
            **deepcopy(self.metadata),
        }
