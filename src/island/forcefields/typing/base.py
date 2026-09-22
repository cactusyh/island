"""Atom-typing engine contracts independent of chemistry adapters."""

from abc import ABC, abstractmethod

from island.core import MolecularSystem, Topology
from island.forcefields.typing.models import AtomTypingResult, AtomTypingRuleSet


class AtomTypingEngine(ABC):
    """Assign external atom types from an authoritative chemical graph."""

    @abstractmethod
    def type_topology(
        self,
        topology: Topology,
        ruleset: AtomTypingRuleSet,
        *,
        strict: bool = True,
    ) -> AtomTypingResult:
        """Type a coordinate-free atomistic topology."""

    def type_system(
        self,
        system: MolecularSystem,
        ruleset: AtomTypingRuleSet,
        *,
        strict: bool = True,
    ) -> AtomTypingResult:
        """Type a MolecularSystem while ignoring coordinates and their metadata."""
        if system.representation != ruleset.supported_representation:
            from island.exceptions import UnsupportedAtomTypingError

            raise UnsupportedAtomTypingError(
                f"Ruleset {ruleset.name!r} supports "
                f"{ruleset.supported_representation!r}, not {system.representation!r}"
            )
        return self.type_topology(system.topology, ruleset, strict=strict)
