"""General interfaces for coordinate-only conformation generation."""

from abc import ABC, abstractmethod

from island.conformations.result import ConformationResult
from island.core import MolecularSystem


class ConformationGenerator(ABC):
    """Generate coordinates without changing molecular chemistry or topology."""

    @abstractmethod
    def generate(self, system: MolecularSystem) -> ConformationResult:
        """Return new coordinates and diagnostics without mutating ``system``."""
