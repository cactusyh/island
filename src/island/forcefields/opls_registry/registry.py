"""Curated OPLS family declarations; registry membership is not MD authority."""

from copy import deepcopy

from island.forcefields.external_pdms import LIMITS
from island.forcefields.oplsaa.pdms_author import PINS
from island.forcefields.oplsaa.source import PIN
from island.forcefields.pdms_opls import DOI, PROFILE
from island.forcefields.selection.inventory import canonical, digest

REGISTRY = "island_opls_registry_v1"
NATIVE = "island_oplsaa_foyer_v1"
ARCHIVE_SHA256 = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"
PAPER_SHA256 = "04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc"


def _document():
    entries = []
    for profile, external in ((NATIVE, False), (PROFILE, True)):
        entry = {
            "family": "opls",
            "profile_identity": profile,
            "version": 1,
            "source_kind": "external_literature" if external else "native",
            "source_hashes": {
                "archive": ARCHIVE_SHA256,
                "paper": PAPER_SHA256,
                **PINS,
            }
            if external
            else {"oplsaa.xml": PIN["xml_sha256"]},
            "citation": DOI if external else PIN["revision"],
            "supported_elements": ["C", "H", "O", "Si"]
            if external
            else ["Br", "C", "Cl", "F", "H", "I", "N", "O", "S"],
            "supported_atom_environments": [
                "Si(C)3(O)",
                "Si(C)2(O)2",
                "O(Si)2",
                "C(Si)(H)3",
                "H(C)",
            ]
            if external
            else ["pinned Foyer SMARTS; exact assignment required"],
            "supported_graph_patterns": ["acyclic linear methyl siloxane; two TMS ends"]
            if external
            else ["exact source SMARTS and complete source rows"],
            "supported_end_groups": ["trimethylsilyl"]
            if external
            else ["only source-typed complete end groups"],
            "charge_method": "author_dft_npa" if external else "native",
            "charge_origin": "Figure 7 DFT/NPA; exact author serialized literals"
            if external
            else "pinned Foyer NonbondedForce atom charges",
            "charge_input": "source_bound_vector",
            "bond": "harmonic",
            "angle": "harmonic",
            "torsion": "OPLS four-term half factor"
            if external
            else "Ryckaert-Bellemans",
            "improper": "source-declared absent in domain"
            if external
            else "no source section; explicit impropers block",
            "nonbonded": "12-6 CHARMM r-squared switch 9-11 A" if external else "12-6",
            "lj_mixing": "geometric_sigma_epsilon",
            "coulomb_14": "0.5",
            "lj_14": "0.5",
            "exclusions": "graph-distance 1,2 excluded; distance 3 scaled 0.5",
            "long_range_convention": {
                "LAMMPS": "PPPM 1e-4; pinned default pair table",
                "OpenMM": "Ewald 1e-7 Reference; NOT identical numerical contract",
                "tail": False,
            }
            if external
            else {"OpenMM": "finite NoCutoff", "LAMMPS": None},
            "backend_support": {
                "OpenMM": "bounded research minimization/NVT"
                if external
                else "finite evaluator",
                "LAMMPS": "deterministic export and unqualified statics"
                if external
                else False,
                "MD_LAMMPS": False,
            },
            "numerical_verification_status": "OpenMM bounded historical checks; cross-engine electrostatics blocked"
            if external
            else "graph-dependent bounded native tests",
            "scientific_validation_status": "not_established",
            "priority": 100 if external else 10,
            "specificity": 100 if external else 10,
            "domain_restrictions": [
                "neutral explicit-H nonisotopic TMS-capped linear siloxanes",
                "HMDS unique O; longer chains require validated PSMILES head identity",
                "no Si-H, Si-OH, Si-O-C, cyclic, crosslinked or P extrapolation",
                "author DFT/NPA branch, not paper Results Qeq branch",
                "local research only; redistribution permission not established",
                "LAMMPS electrostatics/FD/virial/MD continuation unqualified",
            ]
            if external
            else ["no Si or P; complete source typing/charge/row assignment required"],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        entries.append({**entry, "identity": digest(entry)})
    value = {
        "schema": REGISTRY,
        "entries": sorted(entries, key=lambda r: r["profile_identity"]),
        "limitations": deepcopy(LIMITS),
        "policy": "exact domain, exact typing/charges, specificity, priority; no fallback after matching domain fails",
        "unregistered_profiles": "user manifests require explicit review; incomplete literature is not auto-registered",
    }
    return {**value, "identity": digest(value)}


_CANONICAL = canonical(_document())


def canonical_registry():
    """Owned JSON-compatible copy, independent of filesystem enumeration."""
    import json

    return json.loads(_CANONICAL)


def validate_registry(value):
    value = deepcopy(value)
    value["entries"].sort(key=lambda e: e["profile_identity"])
    if canonical(value) != _CANONICAL:
        raise ValueError("Unapproved or tampered canonical OPLS registry")
    return deepcopy(value)
