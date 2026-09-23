"""Force-field abstractions and graph-based atom-typing interfaces."""

from typing import TYPE_CHECKING, Any

from island.forcefields.base import ForceFieldBackend
from island.forcefields.parameterized import ParameterizedSystem
from island.forcefields.parameters import (
    FamilyCoverage,
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    ParameterAssignmentDiagnostic,
    ParameterAssignmentEngine,
    ParameterAssignmentResult,
    ParameterLibrary,
    ParameterSelection,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
    derive_interaction_inventory,
    island_demo_parameters_v1,
)
from island.forcefields.typing import (
    AtomTypeAssignment,
    AtomTypingEngine,
    AtomTypingResult,
    AtomTypingRule,
    AtomTypingRuleSet,
    SiteTypingDiagnostic,
    island_demo_v1_ruleset,
)

if TYPE_CHECKING:
    from island.forcefields.typing import RDKitSmartsAtomTypingEngine

__all__ = [
    "AtomTypeAssignment",
    "AtomTypingEngine",
    "AtomTypingResult",
    "AtomTypingRule",
    "AtomTypingRuleSet",
    "FamilyCoverage",
    "ForceFieldBackend",
    "HarmonicAngleParameter",
    "HarmonicBondParameter",
    "LennardJonesParameter",
    "ParameterAssignmentDiagnostic",
    "ParameterAssignmentEngine",
    "ParameterAssignmentResult",
    "ParameterLibrary",
    "ParameterSelection",
    "ParameterizedSystem",
    "PeriodicTorsionTerm",
    "ProperTorsionParameter",
    "RDKitSmartsAtomTypingEngine",
    "SiteTypingDiagnostic",
    "derive_interaction_inventory",
    "island_demo_parameters_v1",
    "island_demo_v1_ruleset",
]


def __getattr__(name: str) -> Any:
    if name == "RDKitSmartsAtomTypingEngine":
        from island.forcefields.typing import RDKitSmartsAtomTypingEngine

        return RDKitSmartsAtomTypingEngine
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
