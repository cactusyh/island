"""Optional chemistry adapters, imported lazily."""

from typing import Any

__all__ = [
    "RDKitToSystemResult",
    "SystemToRDKitResult",
    "from_rdkit",
    "from_smiles",
    "to_rdkit",
]


def __getattr__(name: str) -> Any:
    if name in __all__:
        from island.chemistry import rdkit_adapter

        return getattr(rdkit_adapter, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
