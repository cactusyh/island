"""Graph-based atom-typing interfaces with optional SMARTS backends."""

from typing import TYPE_CHECKING, Any

from island.forcefields.typing.base import AtomTypingEngine
from island.forcefields.typing.demo import island_demo_v1_ruleset
from island.forcefields.typing.models import (
    AtomTypeAssignment,
    AtomTypingResult,
    AtomTypingRule,
    AtomTypingRuleSet,
    SiteTypingDiagnostic,
)
from island.forcefields.typing.validation import validate_complete_typing_result

if TYPE_CHECKING:
    from island.forcefields.typing.smarts import RDKitSmartsAtomTypingEngine

__all__ = [
    "AtomTypeAssignment",
    "AtomTypingEngine",
    "AtomTypingResult",
    "AtomTypingRule",
    "AtomTypingRuleSet",
    "RDKitSmartsAtomTypingEngine",
    "SiteTypingDiagnostic",
    "island_demo_v1_ruleset",
    "validate_complete_typing_result",
]


def __getattr__(name: str) -> Any:
    if name == "RDKitSmartsAtomTypingEngine":
        from island.forcefields.typing.smarts import RDKitSmartsAtomTypingEngine

        return RDKitSmartsAtomTypingEngine
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
