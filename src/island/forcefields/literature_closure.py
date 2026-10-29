"""J45 graph-bound gap closure for one external OPLS-compatible PDMS table.

Candidate matches are literal source-token associations, never assignments.
This module registers no force field and cannot construct a numerical engine.
"""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack

from .external_pdms import LIMITS, digest, require
from .literature_audit import LiteratureSourceReview, fixture_requirements
from .pdms_literature import INPUT_ERRORS, pinned_bytes

SCHEMA = "island.j45-literature-gap-closure.v1"
CATALOG_SHA = "f8247ff78670dab209af0c0c9076c86182ee20a85e61c10440731bfab4cc318a"
SELECTION_SHA = "64ecdd14316fcc740210853b96917f4f674dfb8ec63d0ebc244ec8030643ede9"
STAGES = (
    "source_complete",
    "typing_complete",
    "charge_complete",
    "parameter_complete",
    "executable",
    "numerically_verified",
    "workflow_verified",
    "scientifically_validated",
)


@dataclass(frozen=True)
class DomainSelection:
    """Pinned J44 review plus a new, independently reconstructable assessment."""

    json_text: str
    catalog_text: str

    def __post_init__(self):
        require(
            sha256(self.json_text.encode()).hexdigest() == SELECTION_SHA,
            "selection_evidence_or_policy_mutated",
        )
        require(
            sha256(self.catalog_text.encode()).hexdigest() == CATALOG_SHA,
            "j44_catalog_mutated",
        )
        review = LiteratureSourceReview(self.catalog_text)
        p = self.payload
        require(p["j44_catalog_sha256"] == review.identity, "wrong_j44_dependency")
        selected = [
            r
            for r in review.payload["extracted_row_candidates"]
            if r["model_scope"] == "OPLS_Table_S4"
        ]
        require(p["candidate_rows"] == selected, "candidate_source_rows_changed")

    @property
    def payload(self):
        from island.dynamics._checkpoint_data import strict_load

        return strict_load(self.json_text)

    @property
    def catalog(self):
        return LiteratureSourceReview(self.catalog_text).payload

    def verify_sources(self, references):
        """Verify J44's known files, without rediscovery or a new bibliography."""
        diagnostics = []
        for row in self.catalog["files"]:
            path = Path(references) / row["filename"]
            try:
                pinned_bytes(path, row["sha256"])
            except INPUT_ERRORS as error:
                diagnostics.append(
                    {
                        "source_file": row["filename"],
                        "expected_sha256": row["sha256"],
                        "status": "source_incompatible"
                        if path.exists()
                        else "asset_unavailable",
                        "reason": str(error),
                    }
                )
        return diagnostics


def _candidate_index(selection):
    index = {}
    for row in selection.payload["candidate_rows"]:
        label = row["atom_or_interaction_labels"]
        family = row["raw_text"].split()[0]
        labels = label if family == "pair_coeff" else label[0].split("_")
        key = (family, tuple(x.lower() for x in labels))
        index[key] = row["evidence_checksum"]
    return index


def _assessment(required, selection, source_errors):
    graph = required["graph"]
    sites = {s["id"]: s for s in graph["sites"]}
    neighbors = {i: [] for i in sites}
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        neighbors[a].append(b)
        neighbors[b].append(a)
    local = {}
    unsupported = []
    for i, site in sites.items():
        element = site["element"]
        near = sorted(sites[n]["element"] for n in neighbors[i])
        valid = (
            (element == "Si" and near in (["C", "C", "C", "O"], ["C", "C", "O", "O"]))
            or (element == "C" and near == ["H", "H", "H", "Si"])
            or (element == "O" and near == ["Si", "Si"])
            or (element == "H" and near == ["C"])
        )
        valid = valid and not (
            site["formal_charge"]
            or site["metadata"]["radical_electrons"]
            or site["metadata"]["aromatic"]
            or site["metadata"]["isotope"]
        )
        valid = valid and site["atomic_number"] == {
            "H": 1,
            "C": 6,
            "O": 8,
            "Si": 14,
        }.get(element)
        valid = valid and all(
            b["order"] == 1 and not b["aromatic"]
            for b in graph["bonds"]
            if i in b["sites"]
        )
        local[str(i)] = {
            "element_token": element,
            "neighbor_elements": near,
            "candidate_environment_only": valid,
            "assigned_type": None,
        }
        if not valid:
            unsupported.append(
                {
                    "site": i,
                    "element": element,
                    "neighbors": near,
                    "reason": "no_exact_Table_S4_typing_charge_and_valence_evidence",
                }
            )
    index = _candidate_index(selection)
    matches, missing = {}, {}
    policies = selection.payload["required_unresolved_policies"]
    for key, request in required["requests"].items():
        kind = key.split(":")[0]
        path = request["sites"]
        tokens = tuple(sites[i]["element"].lower() for i in path)
        family = {
            "bonds": "bond_coeff",
            "angles": "angle_coeff",
            "rb": "dihedral_coeff",
            "nonbonded": "pair_coeff",
        }.get(kind)
        query = (family, (tokens[0], tokens[0]) if kind == "nonbonded" else tokens)
        reverse = (family, query[1][::-1])
        candidates = (
            [index[k] for k in dict.fromkeys((query, reverse)) if k in index]
            if family
            else []
        )
        matches[key] = candidates
        reasons = (
            ["source_units_and_conventions_unresolved"]
            if candidates
            else ["no_citable_complete_assignment"]
        )
        if kind == "policy":
            reasons = [
                policies.get(key.split(":", 1)[1], "exact_policy_evidence_missing")
            ]
        elif kind in ("atom_typing", "partial_charge", "mass"):
            reasons = [
                policies.get(
                    kind
                    if kind == "mass"
                    else "charge_origin"
                    if kind == "partial_charge"
                    else "chemical_domain"
                )
            ]
        elif kind == "improper_applicability":
            reasons = [policies["improper_applicability"]]
        elif kind == "special_pair":
            reasons = [policies["special_pairs"]]
        missing[key] = {
            **request,
            "element_labels": list(tokens),
            "candidate_evidence": candidates,
            "missing_dependencies": reasons,
            "accepted_assignment": None,
        }
    status = "typing_ambiguous" if unsupported else "parameter_incomplete"
    if source_errors:
        status = (
            "source_incompatible"
            if any(e["status"] == "source_incompatible" for e in source_errors)
            else "asset_unavailable"
        )
    return {
        "schema": SCHEMA,
        "domain": selection.payload["selected_domain"],
        "selection_sha256": SELECTION_SHA,
        "j44_catalog_sha256": CATALOG_SHA,
        "requirements": required,
        "local_environment_candidates": local,
        "unsupported_environments": unsupported,
        "literal_row_candidates": matches,
        "missing_rows": missing,
        "unresolved_policies": policies,
        "source_diagnostics": source_errors,
        "status": status,
        "stages": dict.fromkeys(STAGES, False),
        "accepted_rows": [],
        "parameter_assignment": None,
        "exact_charge_vector": None,
        "profile_identity": None,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "native_pcff_acceptance_gate": 1,
        "limitations": LIMITS,
    }


@dataclass(frozen=True)
class DomainClosure:
    """Immutable blocked evidence; a candidate match never grants execution."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(
            p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
            "closure_checksum_changed",
        )
        require(
            p["schema"] == SCHEMA
            and p["status"]
            in (
                "asset_unavailable",
                "parameter_incomplete",
                "typing_ambiguous",
                "source_incompatible",
            )
            and set(p["stages"]) == set(STAGES)
            and all(v is False for v in p["stages"].values())
            and p["selection_sha256"] == SELECTION_SHA
            and p["j44_catalog_sha256"] == CATALOG_SHA
            and p["profile_identity"] is None
            and p["parameter_assignment"] is None
            and p["exact_charge_vector"] is None
            and p["accepted_rows"]
            == p["execution_profiles"]
            == p["backend_executions"]
            == []
            and p["acceptance_gate"] == p["native_pcff_acceptance_gate"] == 1
            and p["limitations"] == LIMITS,
            "incomplete_domain_cannot_grant_execution",
        )

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self, system, selection, source_errors=()):
        expected = close_domain(system, selection, source_errors=source_errors)
        require(self == expected, "stale_graph_source_units_rows_charge_or_profile")

    def require_executable(self):
        raise DomainClosureBlockedError(self)


class DomainClosureBlockedError(ValueError):
    def __init__(self, result):
        self.result = result
        self.status = result.payload["status"]
        self.diagnostics = result.payload["missing_rows"]
        super().__init__(self.status + ": external PDMS literature domain incomplete")


def close_domain(system, selection, *, source_errors=()):
    required = fixture_requirements(system, "opls_aa")
    p = _assessment(required, selection, list(source_errors))
    p["identity"] = digest(p)
    return DomainClosure(pack(p))
