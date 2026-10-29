"""Offline Figure-7 mapping candidate, with no execution or license authority.

The paper's CC-BY grant and the author's repository are separate assets. An
explicit terminal oxygen fixes the oriented E/M fragments for new graphs; this
is a reconstructed rule, not a recovered Moltemplate script. Old audits and the
finite external OPLS adapter remain unchanged.
"""

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.forcefields.external_pdms import LIMITS, digest, pdms_final_graph, require
from island.forcefields.oplsaa.pdms_author import DFT, AuthorOPLSAudit
from island.forcefields.oplsaa.pdms_conventions import FIGURE7_SHA, _partition
from island.forcefields.oplsaa.pdms_literature import OPLSPDMSLiteratureError
from island.forcefields.pdms_literature import INPUT_ERRORS, pinned_bytes

SCHEMA = "island_external_pdms_author_graph_mapping_v4"
PAPER = "pcff/jp4c08471.pdf"
PAPER_SHA = "04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc"


def source_rights(references):
    """Known pinned publication grant; it cannot license separate GitHub files."""
    pinned_bytes(Path(references) / PAPER, PAPER_SHA)
    return {
        "paper": {
            "source_file": PAPER,
            "sha256": PAPER_SHA,
            "doi": "10.1021/acs.jpcb.4c08471",
            "citation": "Xiang et al., J. Phys. Chem. B 129 (2025), 1864-1873",
            "page": 1,
            "printed_page": 1864,
            "raw_text": "This article is licensed under CC-BY 4.0",
            "metadata_stamp": "02/05/2025 10:35:59 (CCBY4_0)",
            "license": "CC-BY-4.0",
            "scope": "article contents; attribution and changes must be recorded",
        },
        "author_repository": {
            "license": None,
            "use_permission": "not_established",
            "redistribution_permission": "not_established",
            "README_usage_instructions_are_a_license": False,
            "paper_license_transferred_to_repository": False,
        },
        "supporting_information": {
            "separate_explicit_grant_in_local_file": False,
            "paper_grant_automatically_transferred": False,
        },
    }


def _document(system, author, references, e_oxygen_site_id):
    require(type(author) is AuthorOPLSAudit, "pinned_author_audit_required")
    author.validate_integrity()
    rights = source_rights(references)
    graph = pdms_final_graph(system)
    partition = _partition(system, graph, e_oxygen_site_id)
    hist = author.payload["source_graph"]["type_charge_histograms"]
    type_ids = {"H": "98", "O": "122", "Si": "868", "C": "871"}
    rows = []
    sites = {a["id"]: a for a in graph["sites"]}
    for r in partition["rows"]:
        element, end = r["element"], r["fragment"] == "E"
        # Select the existing literal, do not round, neutralize or regenerate it.
        printed = DFT[element][int(end)]
        literal = format(float(printed), ".16e")
        mass = author.payload["masses"][type_ids[element]]["original_mass"]
        require(
            Fraction(str(sites[r["site_id"]]["mass"])) == Fraction(mass),
            "source_mass_mismatch",
        )
        require(literal in hist[type_ids[element]], "author_charge_literal_missing")
        rows.append(
            {
                "site_id": r["site_id"],
                "element": element,
                "fragment": r["fragment"],
                "author_type_id": int(type_ids[element]),
                "exact_source_mass_dalton": mass,
                "printed_Figure7_DFT_charge": printed,
                "exact_author_serialized_charge": literal,
            }
        )
    charge = sum(
        (Fraction(r["exact_author_serialized_charge"]) for r in rows), Fraction()
    )
    printed_charge = sum(
        (Fraction(r["printed_Figure7_DFT_charge"]) for r in rows), Fraction()
    )
    require(printed_charge == 0, "printed_source_fragment_charge_inconsistent")
    payload = {
        "schema": SCHEMA,
        "author_audit_identity": author.payload["identity"],
        "source_hashes": author.payload["source_hashes"],
        "source_rights": rights,
        "graph": graph,
        "graph_identity": digest(graph),
        "atom_order": list(system.topology.sites),
        "mapping": rows,
        "E_oxygen_site_id": partition["e_oxygen_site_id"],
        "selection_origin": partition["selection_origin"],
        "terminal_silicon": partition["terminal_silicon"],
        "mapping_rule_provenance": {
            "paper": PAPER,
            "sha256": PAPER_SHA,
            "page": 7,
            "location": "Figure 7: E fragments -OSi(CH3)3 and -Si(CH3)3; M -Si(CH3)2O-",
            "image_sha256": FIGURE7_SHA,
            "reconstructed_from": "published fragments plus independently audited author graph",
            "original_Moltemplate_rule_recovered": False,
            "orientation_is_explicit_not_inferred_from_IDs": True,
        },
        "charge_policy": {
            "branch": "author_DFT_NPA_serialized; not paper_Qeq",
            "printed_decimal_sum": str(printed_charge),
            "exact_serialized_sum_fraction": str(charge),
            "normalization_or_redistribution": False,
            "printed_values_substituted_for_serialized_values": False,
            "finite_v1_neutrality_policy_overridden": False,
        },
        "status": "asset_unavailable",
        "missing_dependencies": [
            {
                "stage": "author_source_use_permission_missing",
                "asset": "pinned GitHub input closure",
            },
            {"stage": "physical_PDMS_backend_and_workflow_qualification_missing"},
        ],
        "execution_profile": None,
        "execution_authorized": False,
        "numerical_MD_run": False,
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    return {**payload, "identity": digest(payload)}


def _classified(error):
    if isinstance(error, OPLSPDMSLiteratureError):
        return error
    status = (
        "asset_unavailable"
        if isinstance(error, FileNotFoundError)
        else "typing_blocked"
    )
    return OPLSPDMSLiteratureError(
        status,
        "author_PDMS_graph_mapping_invalid",
        [{"stage": "author_PDMS_graph_mapping_invalid", "reason": str(error)}],
    )


@dataclass(frozen=True)
class PDMSAuthorGraphMapping:
    json_text: str
    author: AuthorOPLSAudit
    references: Path
    e_oxygen_site_id: int | None

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self, system):
        try:
            expected = _document(
                system, self.author, self.references, self.e_oxygen_site_id
            )
            require(pack(expected) == pack(self.payload), "stale_or_tampered_mapping")
        except OPLSPDMSLiteratureError:
            raise
        except INPUT_ERRORS as error:
            raise _classified(error) from error

    def require_executable(self, system):
        self.validate_integrity(system)
        error = OPLSPDMSLiteratureError(
            "asset_unavailable",
            "author_source_use_permission_missing",
            self.payload["missing_dependencies"],
        )
        error.result = self
        raise error


def map_author_pdms(system, author, *, references, e_oxygen_site_id=None):
    """Build only a candidate record; an explicit anchor never grants MD authority."""
    try:
        doc = _document(system, author, references, e_oxygen_site_id)
        return PDMSAuthorGraphMapping(
            pack(doc), author, Path(references), e_oxygen_site_id
        )
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        raise _classified(error) from error


def load_mapping(
    path, expected_sha256, system, author, *, references, e_oxygen_site_id=None
):
    try:
        raw = pinned_bytes(path, expected_sha256)
        result = PDMSAuthorGraphMapping(
            raw.decode(), author, Path(references), e_oxygen_site_id
        )
        result.validate_integrity(system)
        return result
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        raise _classified(error) from error
