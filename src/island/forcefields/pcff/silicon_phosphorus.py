"""Opt-in Si/P source tranche; a source label is never execution authority.

Only finite neutral silane is executable in v1. Other profiles retain complete
source-query diagnostics without publishing an incomplete Class-II model.
Historical converter-derived zero policies are deliberately not used here.
"""

from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack

from .automatic import chemical_graph, type_pcff_atoms
from .charges import identity
from .class2 import assign_pcff_parameters
from .expanded import environments
from .fallbacks import MSI_POLICY
from .model import define_pcff_model, special_pair_policy
from .operational import inspect_pcff_operational_support
from .source import FLAGS, PIN, boundary, require

SCHEMA = "island_pcff_sip_source_assessment_v1"
PROFILES = {
    "island_pcff_silane_source_v1": "silane",
    "island_pcff_siloxane_source_v1": "siloxane",
    "island_pcff_organosilicon_source_v1": "organosilicon",
    "island_pcff_phosphate_ester_audit_v1": "phosphate_ester",
    "island_pcff_phosphonate_audit_v1": "phosphonate",
}
SPECIAL = special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))


def _assess(system, source, profile):
    require(profile in PROFILES, "Unregistered Si/P source profile")
    source.require_assignment()
    require(source.identity["sha256"] == PIN["sha256"], "Si/P source hash mismatch")
    system.validate()
    graph = chemical_graph(system)
    atoms, adjacency, _, env = environments(graph)
    elements = {a["element"] for a in atoms.values()}
    domain = PROFILES[profile]
    diagnostics = []
    if graph["representation"] != "atomistic" or graph["has_box"]:
        diagnostics.append(
            {"stage": "domain", "reason": "finite_nonperiodic_atomistic_only"}
        )
    if not elements <= {"C", "H", "O", "Si", "P"} or any(
        a["formal_charge"]
        or a["metadata"]["isotope"]
        or a["metadata"]["radical_electrons"]
        for a in atoms.values()
    ):
        diagnostics.append(
            {
                "stage": "domain",
                "reason": "charged_ionic_metallic_isotopic_or_radical_domain_unsupported",
            }
        )
    for i, a in atoms.items():
        if a["element"] == "Si" and env[i]["bond_orders"] != [1, 1, 1, 1]:
            diagnostics.append(
                {
                    "stage": "typing",
                    "sites": [i],
                    "reason": "unsupported_silicon_valence",
                }
            )
        if a["element"] == "P":
            diagnostics.append(
                {
                    "stage": "typing",
                    "sites": [i],
                    "reason": "phosphorus_oxidation_and_source_template_mapping_unresolved",
                    "evidence": "pcff.frc:163 p connections4; 164 p= connections5 phosphazene; pcff_templates.dat:1562 p= oxo template has four neighbors; no silent p/p= substitution",
                    "environment": env[i],
                }
            )
    if domain == "silane":
        valid = (
            len(atoms) == 5
            and sum(a["element"] == "Si" for a in atoms.values()) == 1
            and sum(a["element"] == "H" for a in atoms.values()) == 4
        )
        valid = valid and all(
            len(adjacency[i]) == (4 if a["element"] == "Si" else 1)
            for i, a in atoms.items()
        )
    elif domain == "siloxane":
        valid = (
            "Si" in elements and "O" in elements and elements <= {"C", "H", "O", "Si"}
        )
        valid = valid and any(
            a["element"] == "O" and env[i]["neighbor_elements"].get("Si") == 2
            for i, a in atoms.items()
        )
    elif domain == "organosilicon":
        valid = {"Si", "C"} <= elements <= {"C", "H", "Si"}
    else:
        valid = "P" in elements and "O" in elements and elements <= {"C", "H", "O", "P"}
    if not valid:
        diagnostics.append(
            {"stage": "domain", "reason": "outside_selected_source_profile"}
        )
    typing = inspection = None
    if not diagnostics:
        typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v1")
        if not typing.complete:
            diagnostics.extend(
                {"stage": "typing", **d} for d in typing.payload["diagnostics"]
            )
        else:
            inspection = inspect_pcff_operational_support(
                system, typing, resolution_policy=MSI_POLICY, special_pairs=SPECIAL
            )
            if not inspection["charges_complete"]:
                diagnostics.append(
                    {"stage": "charges", "reason": "native_component_charge_incomplete"}
                )
            for query in inspection["interaction_queries"]:
                if (
                    query["resolution"]["status"] not in ("assigned", "not_applicable")
                    or not query["dependencies_complete"]
                ):
                    diagnostics.append(
                        {
                            "stage": "parameters",
                            "reason": "required_source_row_or_dependency_unresolved",
                            "interaction": query,
                        }
                    )
            if not inspection["model_complete"]:
                diagnostics.extend(
                    {"stage": "model", **d} for d in inspection["model_diagnostics"]
                )
    return {
        "schema": SCHEMA,
        "profile": profile,
        "domain": domain,
        "source": source.identity,
        "graph": graph,
        "graph_identity": identity(graph),
        "atom_order": sorted(atoms),
        "coordinates": [
            [i, list(map(float, system.coordinates.get(i)))] for i in sorted(atoms)
        ],
        "typing_identity": typing.identity if typing else None,
        "typing": typing.payload if typing else None,
        "charge_origin": "PCFF/source_native_bond_increments",
        "inspection": inspection,
        "resolution_policy": MSI_POLICY,
        "special_pairs": SPECIAL,
        "complete": not diagnostics,
        "diagnostics": diagnostics,
        "publication_status": "source_complete"
        if not diagnostics
        else "source_incomplete",
        "periodic_md_status": "unsupported_in_this_finite_profile",
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFSiPAssessment:
    json_text: str
    source: object

    @boundary
    def validate_integrity(self, system):
        payload = unpack(self.json_text)
        require(
            pack(payload) == pack(_assess(system, self.source, payload["profile"])),
            "Contradictory Si/P source assessment or final graph/atom ordering",
        )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def assess_pcff_sip(system, source, *, profile):
    """Inspect explicit opt-in domains; incomplete results authorize no dynamics."""
    return PCFFSiPAssessment(pack(_assess(system, source, profile)), source)


@boundary
def prepare_pcff_sip(system, assessment):
    """Return a native finite model only after every applicable source row passes."""
    assessment.validate_integrity(system)
    p = assessment.payload
    require(p["complete"], "Incomplete Si/P source profile: " + str(p["diagnostics"]))
    from .automatic import PCFFAutomaticTypingResult, assign_automatic_pcff_charges

    typing = PCFFAutomaticTypingResult(pack(p["typing"]), assessment.source)
    charges = assign_automatic_pcff_charges(
        system, typing, resolution_policy=MSI_POLICY
    )
    assignment = assign_pcff_parameters(
        system, typing, charges, resolution_policy=MSI_POLICY
    )
    model = define_pcff_model(assignment, special_pairs=SPECIAL)
    require(
        assignment.payload["parameter_coverage_complete"], "Incomplete interaction set"
    )
    require(
        model.payload["model_definition_complete"]
        and all(
            t["origin"] == "source_row" and t["source_rows"]
            for t in model.payload["terms"]
        ),
        "Missing coupling cannot be assigned zero",
    )
    return model


@boundary
def load_pcff_sip_assessment(path, source, *, system):
    result = PCFFSiPAssessment(Path(path).read_text(), source)
    result.validate_integrity(system)
    return result
