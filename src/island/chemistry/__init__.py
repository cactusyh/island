"""Optional chemistry adapters, imported lazily."""

from typing import Any

_RDKIT_ADAPTER_EXPORTS = {
    "RDKitToSystemResult",
    "SystemToRDKitResult",
    "from_rdkit",
    "from_smiles",
    "to_rdkit",
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
    "to_rdkit",
]


def __getattr__(name: str) -> Any:
    if name in _RDKIT_ADAPTER_EXPORTS:
        from island.chemistry import rdkit_adapter

        return getattr(rdkit_adapter, name)
    if name in _PSMILES_EXPORTS:
        from island.chemistry import psmiles

        return getattr(psmiles, name)
    if name in _STEREOCHEMISTRY_EXPORTS:
        from island.chemistry import stereochemistry

        return getattr(stereochemistry, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
