"""Force-field abstractions and graph-based atom-typing interfaces."""

from typing import TYPE_CHECKING, Any

from island.forcefields.base import ForceFieldBackend
from island.forcefields.charges import (
    AtomTypeChargeEngine,
    AtomTypeChargeEntry,
    AtomTypeChargeTable,
    ChargeAssignment,
    ChargeAssignmentDiagnostic,
    ChargeAssignmentEngine,
    ChargeAssignmentResult,
    ProvidedChargeEngine,
    island_demo_charges_v1,
)
from island.forcefields.composition import compose_parameterized_system
from island.forcefields.nonbonded import (
    MixedLJParameters,
    NonbondedPolicy,
    PairScaling,
    nonbonded_policy_signature,
    shortest_bond_distance,
)
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
    "AtomTypeChargeEngine",
    "AtomTypeChargeEntry",
    "AtomTypeChargeTable",
    "AtomTypingEngine",
    "AtomTypingResult",
    "AtomTypingRule",
    "AtomTypingRuleSet",
    "ChargeAssignment",
    "ChargeAssignmentDiagnostic",
    "ChargeAssignmentEngine",
    "ChargeAssignmentResult",
    "FamilyCoverage",
    "ForceFieldBackend",
    "HarmonicAngleParameter",
    "HarmonicBondParameter",
    "LennardJonesParameter",
    "MixedLJParameters",
    "NonbondedPolicy",
    "PairScaling",
    "ParameterAssignmentDiagnostic",
    "ParameterAssignmentEngine",
    "ParameterAssignmentResult",
    "ParameterLibrary",
    "ParameterSelection",
    "ParameterizedSystem",
    "PeriodicTorsionTerm",
    "ProperTorsionParameter",
    "ProvidedChargeEngine",
    "RDKitSmartsAtomTypingEngine",
    "SiteTypingDiagnostic",
    "compose_parameterized_system",
    "derive_interaction_inventory",
    "island_demo_charges_v1",
    "island_demo_parameters_v1",
    "island_demo_v1_ruleset",
    "nonbonded_policy_signature",
    "shortest_bond_distance",
]


def __getattr__(name: str) -> Any:
    if name == "RDKitSmartsAtomTypingEngine":
        from island.forcefields.typing import RDKitSmartsAtomTypingEngine

        return RDKitSmartsAtomTypingEngine
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
