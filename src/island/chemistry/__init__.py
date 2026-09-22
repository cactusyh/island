"""Optional chemistry adapters, imported lazily."""

from typing import Any

_RDKIT_ADAPTER_EXPORTS = {
    "RDKitToSystemResult",
    "SystemToRDKitResult",
    "from_rdkit",
    "from_smiles",
    "to_rdkit",
}
_RDKIT_GRAPH_EXPORTS = {
    "GraphToRDKitResult",
    "system_to_rdkit_graph",
    "topology_to_rdkit_graph",
}
_PSMILES_EXPORTS = {"AttachmentPoint", "RepeatUnit", "parse_psmiles"}
_STEREOCHEMISTRY_EXPORTS = {
    "RepeatUnitStereochemistry",
    "StereocenterInfo",
    "StereochemicalSequence",
    "generate_stereochemical_sequence",
    "inspect_repeat_unit_stereochemistry",
    "invert_repeat_unit_stereochemistry",
    "require_controllable_stereocenter",
}

__all__ = [
    "AttachmentPoint",
    "GraphToRDKitResult",
    "RDKitToSystemResult",
    "RepeatUnit",
    "RepeatUnitStereochemistry",
    "StereocenterInfo",
    "StereochemicalSequence",
    "SystemToRDKitResult",
    "from_rdkit",
    "from_smiles",
    "generate_stereochemical_sequence",
    "inspect_repeat_unit_stereochemistry",
    "invert_repeat_unit_stereochemistry",
    "parse_psmiles",
    "require_controllable_stereocenter",
    "system_to_rdkit_graph",
    "to_rdkit",
    "topology_to_rdkit_graph",
]


def __getattr__(name: str) -> Any:
    if name in _RDKIT_ADAPTER_EXPORTS:
        from island.chemistry import rdkit_adapter

        return getattr(rdkit_adapter, name)
    if name in _RDKIT_GRAPH_EXPORTS:
        from island.chemistry import rdkit_graph

        return getattr(rdkit_graph, name)
    if name in _PSMILES_EXPORTS:
        from island.chemistry import psmiles

        return getattr(psmiles, name)
    if name in _STEREOCHEMISTRY_EXPORTS:
        from island.chemistry import stereochemistry

        return getattr(stereochemistry, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
