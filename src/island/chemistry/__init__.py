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

__all__ = [
    "AttachmentPoint",
    "RDKitToSystemResult",
    "RepeatUnit",
    "SystemToRDKitResult",
    "from_rdkit",
    "from_smiles",
    "parse_psmiles",
    "to_rdkit",
]


def __getattr__(name: str) -> Any:
    if name in _RDKIT_ADAPTER_EXPORTS:
        from island.chemistry import rdkit_adapter

        return getattr(rdkit_adapter, name)
    if name in _PSMILES_EXPORTS:
        from island.chemistry import psmiles

        return getattr(psmiles, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
