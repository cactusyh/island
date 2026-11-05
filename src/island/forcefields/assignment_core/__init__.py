"""Opt-in coverage analysis for PSMILES. This package cannot activate a profile."""

from .profiles import builtin_contracts, profile_contract
from .resolver import resolve, validate_class_i, validate_class_ii
from .schema import (
    AssignmentDecision,
    ChemicalGraph,
    CompletenessReport,
    EvidenceRecord,
    ForceFieldProfile,
    ParameterRecord,
    PolymerSpec,
    TypedGraph,
)
from .topology import normalize, rebuild
from .typing import type_graph

__all__ = [
    "AssignmentDecision",
    "ChemicalGraph",
    "CompletenessReport",
    "EvidenceRecord",
    "ForceFieldProfile",
    "ParameterRecord",
    "PolymerSpec",
    "TypedGraph",
    "builtin_contracts",
    "normalize",
    "profile_contract",
    "rebuild",
    "resolve",
    "type_graph",
    "validate_class_i",
    "validate_class_ii",
]
