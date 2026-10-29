"""Explicit capabilities. Registry inspection never authorizes a final graph."""

from copy import deepcopy

from island.forcefields.external_pdms import LIMITS

STAGES = (
    "typing_complete",
    "charge_complete",
    "parameter_complete",
    "executable",
    "numerically_verified",
    "scientifically_validated",
)
PDMS = "island_external_opls_pdms_xiang2025_dft_local_v1"
OPLS = "island_oplsaa_foyer_v1"
USER_OPLS = "island_user_provided_opls_v1"
PCFF = "island_pcff_acyclic_cho_v1"
CLASS2 = (
    "quartic_bond",
    "quartic_angle",
    "torsion_3",
    "wilson_out_of_plane",
    "bond-bond",
    "bond-angle",
    "bond-bond_1_3",
    "angle-angle",
    "angle-torsion_3",
    "end_bond-torsion_3",
    "middle_bond-torsion_3",
    "angle-angle-torsion_1",
)


def capability_registry():
    """Return owned declarations, including legacy PCFF registrations as such.

    Null stage values mean a graph and explicit assets must be inspected. Backend
    declarations describe existing adapters, not validated periodic science.
    """
    entries = []

    def add(
        name,
        family,
        kind,
        elements,
        environments,
        forms,
        charges,
        mixing,
        scaling,
        electrostatics,
        backend,
        **extra,
    ):
        entries.append(
            {
                "profile_identity": name,
                "family": family,
                "source_kind": kind,
                "supported_elements": elements,
                "element_scope": "bounded declaration; environments still require typing"
                if elements
                else "source-dependent, not universal",
                "supported_atom_environments": environments,
                "bond_terms": forms[0],
                "angle_terms": forms[1],
                "torsions": forms[2],
                "impropers": forms[3],
                "class2_cross_terms": list(CLASS2[4:]) if family == "pcff" else [],
                "nonbonded_terms": "9-6" if family == "pcff" else "12-6",
                "charge_methods": charges,
                "charges": "exact source vector/increments; no normalization or redistribution",
                "mixing_rule": mixing,
                "scaling_14": scaling,
                "exclusions": "explicit source policy; never inferred from family name",
                "long_range_electrostatics": electrostatics,
                "backend_support": backend,
                "periodic_backend_support": "bounded external PDMS research Ewald only"
                if name == PDMS
                else "no new periodic preparation authorization",
                "stages": {
                    s: False if s == "scientifically_validated" else None
                    for s in STAGES
                },
                "scientific_verification_status": "not_established; historical bounded software checks only",
                "production_validated": False,
                "simulation_readiness": "not_established",
                "limitations": deepcopy(LIMITS),
                **extra,
            }
        )

    forms = (["quartic"], ["quartic"], ["three-term Class-II"], ["Wilson out-of-plane"])
    from island.forcefields.pcff.expanded import NUMBERS

    for name in (PCFF, *(f"island_pcff_source_graph_v{v}" for v in range(1, 7))):
        add(
            name,
            "pcff",
            "native",
            ["C", "H", "O"] if name == PCFF else sorted(NUMBERS),
            "bounded final-graph rules; source_graph_v1 inspects additional source environments",
            forms,
            ["native"],
            "sixth-power Class-II",
            "explicit caller LJ/Coulomb tuples",
            "finite NoCutoff",
            {
                "OpenMM": "bounded native Class-II evaluator",
                "LAMMPS": "historical export; no new MD runner",
            },
        )
    add(
        OPLS,
        "opls",
        "native",
        ["C", "H", "N", "O", "F", "S", "Cl", "Br", "I"],
        "pinned Foyer SMARTS; complete exact graph assignment required",
        (
            ["harmonic"],
            ["harmonic"],
            ["Ryckaert-Bellemans"],
            ["no explicit source section"],
        ),
        ["native"],
        "geometric_sigma_and_epsilon",
        {"lj": 0.5, "coulomb": 0.5},
        "finite NoCutoff",
        {"OpenMM": "finite native evaluator", "LAMMPS": False},
    )
    add(
        PDMS,
        "opls",
        "external_literature",
        ["Si", "C", "O", "H"],
        "neutral explicit-H trimethylsilyl-capped linear siloxanes; explicit terminal oxygen anchor",
        (
            ["harmonic"],
            ["harmonic"],
            ["OPLS four-term, half factor"],
            ["source-declared absent"],
        ),
        ["author_dft_npa"],
        "geometric_sigma_and_epsilon",
        {"lj": 0.5, "coulomb": 0.5},
        "OpenMM Ewald 1e-7; exported author LAMMPS PPPM 1e-4 is a different solver",
        {
            "OpenMM": "local research NVT only",
            "LAMMPS": "export and static comparison; MD workflow unqualified",
        },
    )
    add(
        USER_OPLS,
        "opls",
        "user_provided",
        [],
        "explicit graph-bound external OPLS manifest",
        (
            ["declared harmonic"],
            ["declared harmonic"],
            ["declared RB or OPLS"],
            ["declared policy"],
        ),
        ["provided"],
        "explicit manifest",
        "explicit manifest",
        "finite declared model",
        {"OpenMM": "existing external adapter only", "LAMMPS": False},
    )
    for family in ("gaff", "gaff2"):
        add(
            f"island_{family}_ambertools_v1",
            family,
            "native",
            [],
            "AmberTools exact typing and parameter checks; no Si/P coverage inferred",
            (
                ["Amber harmonic"],
                ["Amber harmonic"],
                ["Amber periodic"],
                ["Amber improper"],
            ),
            ["am1bcc", "provided", "resp_psi4", "external_provided"],
            "Amber source-declared",
            {"lj": "Amber source", "coulomb": "Amber source"},
            "existing Amber workflow gates",
            {"OpenMM": "existing Amber import", "LAMMPS": False},
            charge_provenance="RESP/Psi4 requires its existing independently verified charge record, not a renamed provided vector",
        )
    # Preserve legacy coverage declarations without widening strict J46 semantics.
    from island.forcefields.pcff.registry import list_pcff_profiles

    for old in list_pcff_profiles():
        add(
            old["name"],
            "pcff",
            "native" if old["state"] == "operational" else "comparison_only",
            [],
            old["definition"].get("chemical_rule"),
            forms,
            ["native"],
            "source-declared",
            "legacy explicit policy",
            "existing legacy gates",
            {"OpenMM": "existing legacy adapter", "LAMMPS": "legacy evidence only"},
            selection_adapter=False,
            legacy_registration=old,
            authorized_atom_labels=old["definition"].get("authorized_atom_labels", []),
            blocking_reason="legacy policy requires separate explicit review; J46 strict source never enables fallback",
        )
    return deepcopy(sorted(entries, key=lambda e: (e["family"], e["profile_identity"])))


def inspect_profile(name):
    matches = [p for p in capability_registry() if p["profile_identity"] == name]
    if len(matches) != 1:
        raise ValueError(f"Unknown explicit source profile: {name}")
    return matches[0]
