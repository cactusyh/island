"""Empty scientific contracts, deliberately separate from all operational registries."""

from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    FAMILIES,
    ForceFieldProfile,
)

BASE_CONVENTIONS = (
    "energy_unit",
    "length_unit",
    "angular_unit",
    "atom_typing",
    "chemical_scope",
    "nonbonded_form",
    "sigma_definition",
    "epsilon_definition",
    "mixing_rule",
    "exclusions_12",
    "exclusions_13",
    "scaling_14_lj",
    "scaling_14_coulomb",
    "charge_origin",
    "charge_tolerance",
    "bond_polynomial",
    "angle_polynomial",
    "torsion_phase",
    "torsion_normalization",
    "improper_ordering",
    "stereochemistry",
    "electrostatics",
    "long_range",
    "cutoff",
    "tail",
    "finite_or_periodic_scope",
)
CLASS_II_CONVENTIONS = (
    "cross_term_ordering",
    "cross_term_dependencies",
    "cross_term_applicability",
)
BASE_TERMS = ("bond", "angle", "torsion", "improper", "nonbonded", "charge")


def required_conventions(profile):
    return BASE_CONVENTIONS + (
        CLASS_II_CONVENTIONS if profile.functional_class == "II" else ()
    )


def required_families(profile):
    return BASE_TERMS + (CROSS_TERMS if profile.functional_class == "II" else ())


def profile_contract(family):
    """Names specify obligations, never an assertion of available parameter coverage."""
    return ForceFieldProfile(
        "j82_contract_" + family.lower().replace("-", "_"),
        family,
        "II" if family in ("PCFF", "COMPASS") else "I",
    )


def builtin_contracts():
    return tuple(profile_contract(f) for f in FAMILIES)
