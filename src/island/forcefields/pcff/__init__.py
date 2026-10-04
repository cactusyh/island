"""PCFF source ingestion and explicit-type native increments; no energy model."""

from .charges import (
    PCFFChargeResult,
    PCFFTypingResult,
    assign_pcff_charges,
    assign_pcff_types,
    load_pcff_record,
    save_pcff_record,
)
from .source import PCFFSource, load_pcff_source

__all__ = [
    "PCFFChargeResult",
    "PCFFSource",
    "PCFFTypingResult",
    "assign_pcff_charges",
    "assign_pcff_types",
    "load_pcff_record",
    "load_pcff_source",
    "save_pcff_record",
]

from .automatic import (
    PCFFAutomaticChargeResult,
    PCFFAutomaticTypingResult,
    assign_automatic_pcff_charges,
    load_pcff_automatic_record,
    save_pcff_automatic_record,
    type_pcff_atoms,
)

__all__ += [
    "PCFFAutomaticChargeResult",
    "PCFFAutomaticTypingResult",
    "assign_automatic_pcff_charges",
    "load_pcff_automatic_record",
    "save_pcff_automatic_record",
    "type_pcff_atoms",
]
