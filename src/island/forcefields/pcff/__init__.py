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

from .class2 import (
    PCFFClass2Result,
    assign_pcff_parameters,
    inspect_pcff_class2,
    load_pcff_parameters,
    save_pcff_parameters,
)

__all__ += [
    "PCFFClass2Result",
    "assign_pcff_parameters",
    "inspect_pcff_class2",
    "load_pcff_parameters",
    "save_pcff_parameters",
]

from .model import (
    PCFFModelSpecification,
    define_pcff_model,
    load_pcff_model,
    save_pcff_model,
    special_pair_policy,
)

__all__ += [
    "PCFFModelSpecification",
    "define_pcff_model",
    "load_pcff_model",
    "save_pcff_model",
    "special_pair_policy",
]

from .catalog import inspect_pcff_full_source, resolve_pcff_source_record
from .expanded import assign_pcff_source_types

__all__ += [
    "assign_pcff_source_types",
    "inspect_pcff_full_source",
    "resolve_pcff_source_record",
]

from .domain_coverage import pcff_domain_coverage

__all__ += ["pcff_domain_coverage"]
