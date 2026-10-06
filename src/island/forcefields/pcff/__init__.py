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

from .operational import inspect_pcff_operational_support

__all__ += ["inspect_pcff_operational_support"]

from .resolution import (
    PCFFResolutionAssessment,
    assess_pcff_resolution,
    load_pcff_resolution_assessment,
    save_pcff_resolution_assessment,
)

__all__ += [
    "PCFFResolutionAssessment",
    "assess_pcff_resolution",
    "load_pcff_resolution_assessment",
    "save_pcff_resolution_assessment",
]

from .variants import (
    PCFFSourceComparison,
    PCFFSourceVariant,
    PCFFVariantPolicy,
    PCFFVariantResolution,
    PCFFVariantSelection,
    compare_pcff_sources,
    load_pcff_source_comparison,
    load_pcff_source_variant,
    load_pcff_variant_resolution,
    resolve_pcff_variant_record,
    save_pcff_source_comparison,
    save_pcff_variant_resolution,
    select_native_pcff_source,
)

__all__ += [
    "PCFFSourceComparison",
    "PCFFSourceVariant",
    "PCFFVariantPolicy",
    "PCFFVariantResolution",
    "PCFFVariantSelection",
    "compare_pcff_sources",
    "load_pcff_source_comparison",
    "load_pcff_source_variant",
    "load_pcff_variant_resolution",
    "resolve_pcff_variant_record",
    "save_pcff_source_comparison",
    "save_pcff_variant_resolution",
    "select_native_pcff_source",
]

from .variant_audit import (
    PCFFSourceVariantAudit,
    audit_pcff_source_variants,
    load_pcff_source_variant_audit,
    save_pcff_source_variant_audit,
)

__all__ += [
    "PCFFSourceVariantAudit",
    "audit_pcff_source_variants",
    "load_pcff_source_variant_audit",
    "save_pcff_source_variant_audit",
]

from .operational_profile import (
    PCFFOperationalProfile,
    PCFFOperationalSelection,
    load_pcff_operational_profile,
    save_pcff_operational_profile,
)

__all__ += [
    "PCFFOperationalProfile",
    "PCFFOperationalSelection",
    "load_pcff_operational_profile",
    "save_pcff_operational_profile",
]

from .registry import (
    inspect_pcff_profile,
    list_pcff_profiles,
    select_pcff_profile,
    validate_pcff_profile_system,
)

__all__ += [
    "inspect_pcff_profile",
    "list_pcff_profiles",
    "select_pcff_profile",
    "validate_pcff_profile_system",
]

from .polymer import (
    PCFFFinalGraphPreparation,
    edit_polymer_topology,
    prepare_pcff_polymer,
    reassign_pcff_polymer,
)

__all__ += [
    "PCFFFinalGraphPreparation",
    "edit_polymer_topology",
    "prepare_pcff_polymer",
    "reassign_pcff_polymer",
]

from .interop import export_msi2lmp_car_mdf, load_msi2lmp_car_mdf

__all__ += ["export_msi2lmp_car_mdf", "load_msi2lmp_car_mdf"]

from .registry import (
    PCFFProfileRegistration,
    load_pcff_profile_registration,
    save_pcff_profile_registration,
)

__all__ += [
    "PCFFProfileRegistration",
    "load_pcff_profile_registration",
    "save_pcff_profile_registration",
]

from .row_coverage import pcff_source_row_ledger
from .validation_cache import clear_pcff_validation_cache, pcff_validation_cache_info

__all__ += [
    "clear_pcff_validation_cache",
    "pcff_source_row_ledger",
    "pcff_validation_cache_info",
]
