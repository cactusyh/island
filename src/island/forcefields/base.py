"""Force-field backend contracts."""

from abc import ABC, abstractmethod

from island.core.system import MolecularSystem
from island.forcefields.parameterized import ParameterizedSystem


class ForceFieldBackend(ABC):
    """Backend that assigns parameters without mutating chemical topology."""

    @abstractmethod
    def parameterize(self, system: MolecularSystem) -> ParameterizedSystem:
        """Return force-field assignments for ``system``."""
