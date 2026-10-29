"""Source-bound Sun siloxane parameter records, separate from native PCFF.

The numerical supplement is usable for exact row inspection and graph-bound
completeness audits. It has no execution profile: missing hydrocarbon sources,
blank coefficients, conventions and coupling applicability remain blockers.
"""

from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.dynamics._checkpoint_data import strict_load
from island.exceptions import IslandError
from island.workflows.storage import child

from .external_pdms import LIMITS, digest
from .literature_audit import fixture_requirements
from .pdms_literature import INPUT_ERRORS, pinned_bytes, reconstruct_rows

SCHEMA = "island_sun_siloxane_parameter_supplement_v1"
ASSESSMENT_SCHEMA = "island_sun_siloxane_graph_assessment_v1"
# Filled from the independently reconstructed, versioned supplement.
CATALOG_IDENTITY = "22471d4ff4449b4147c15bd4207bbe4fe21018cfe2d483a12737998c6bad4832"


class SiloxaneParameterError(IslandError):
    """Classified source/integrity failure, without numerical execution."""

    def __init__(self, status, stage):
        self.status = status
        self.diagnostics = [{"stage": stage}]
        super().__init__(status + ": " + stage)


class SiloxaneParametersBlockedError(IslandError):
    def __init__(self, result):
        self.result = result
        self.status = result.payload["status"]
        self.diagnostics = result.payload["missing_dependencies"]
        super().__init__(self.status + ": literature model is not executable")


@dataclass(frozen=True)
class SiloxaneParameterSupplement:
    """Immutable audited numbers, not a PreparedForceField or C5 bundle."""

    json_text: str

    def __post_init__(self):
        try:
            p = unpack(self.json_text)
            valid = (
                digest(p) == CATALOG_IDENTITY
                and p["schema"] == SCHEMA
                and p["checksum"]
                == digest({k: v for k, v in p.items() if k != "checksum"})
                and p["execution_profiles"] == []
                and p["execution_authorized"] is False
                and p["limitations"] == LIMITS
            )
        except INPUT_ERRORS as error:
            raise SiloxaneParameterError(
                "source_incompatible", "malformed_parameter_supplement"
            ) from error
        if not valid:
            raise SiloxaneParameterError(
                "source_incompatible", "parameter_rows_or_authority_changed"
            )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return CATALOG_IDENTITY

    def verify_sources(self, references):
        for sid, source in self.payload["sources"].items():
            try:
                pinned_bytes(child(references, source["file"]), source["sha256"])
            except FileNotFoundError as error:
                raise SiloxaneParameterError(
                    "asset_unavailable", "source_missing:" + sid
                ) from error
            except INPUT_ERRORS as error:
                raise SiloxaneParameterError(
                    "source_incompatible", "source_hash_or_access_failed:" + sid
                ) from error

    def lookup(self, source, family, labels):
        """Literal row lookup. No implicit equivalence, reversal or wildcard."""
        if source not in self.payload["sources"]:
            raise SiloxaneParameterError("source_incompatible", "unknown_source")
        return tuple(
            r
            for r in self.payload["rows"]
            if r["source"] == source and r["family"] == family and r["labels"] == labels
        )


def load_siloxane_parameters(path, expected_sha256, *, references):
    """Read only an explicit catalog path/pin and its listed local evidence."""
    try:
        raw = pinned_bytes(path, expected_sha256)
        result = SiloxaneParameterSupplement(pack(strict_load(raw.decode("utf-8"))))
        result.verify_sources(references)
        p = result.payload
        evidence = p["numeric_evidence"]
        text = pinned_bytes(
            child(Path(path).parent, evidence["file"]), evidence["sha256"]
        )
        rows = {r["id"]: r["numeric_record"] for r in p["rows"]}
        reconstructed, errors = reconstruct_rows(
            {"rows": rows}, {evidence["file"]: text}
        )
        if errors or len(reconstructed) != len(rows):
            raise SiloxaneParameterError(
                "source_incompatible", "decimal_row_reconstruction_failed"
            )
        return result
    except SiloxaneParameterError:
        raise
    except FileNotFoundError as error:
        raise SiloxaneParameterError(
            "asset_unavailable", "catalog_or_numeric_evidence_missing"
        ) from error
    except INPUT_ERRORS as error:
        raise SiloxaneParameterError(
            "source_incompatible", "catalog_pin_or_evidence_invalid"
        ) from error


def _row_types(row):
    label = row["labels"]
    if "/" not in label:
        return label.split("-")
    left, right = [s.split("-") for s in label.split("/")]
    # Candidate association only; source angle-angle atom roles still need
    # verification before any numerical assignment is permitted.
    return left + [right[-1]]


def _matches(row, query):
    types = _row_types(row)
    if len(types) != len(query):
        return False
    reverse = (
        [query[3], query[1], query[2], query[0]]
        if row["family"] == "angle-angle"
        else query[::-1]
    )
    return any(
        all(a == "*" or a == b for a, b in zip(types, q)) for q in (query, reverse)
    )


def _candidates(rows, family, sites, atoms, neighbors):
    family = {"nonbonded": "nonbond(9-6)"}.get(family, family)
    native = {"Si": "si4", "C": "c4", "O": "o2", "H": "h1"}
    query = [native.get(atoms[i]["element"], "unsupported") for i in sites]
    if (
        family == "nonbond(9-6)"
        and query == ["si4"]
        and any(atoms[j]["element"] == "C" for j in neighbors[sites[0]])
    ):
        query = ["si4c"]
    comparison = [atoms[i]["element"] for i in sites]
    primary, inherited, wildcard, partial, comparison_only = [], [], [], [], []
    for row in rows:
        if row["family"] != family:
            continue
        labels = query if row["source"] == "sun1997" else comparison
        if not _matches(row, labels):
            continue
        if row["source"] == "sun1995" and family == "nonbond(9-6)":
            comparison_only.append(row["id"])
        elif row["missing_fields"]:
            partial.append(row["id"])
        elif "*" in row["labels"]:
            wildcard.append(row["id"])
        elif row["source"] == "sun1997":
            primary.append(row["id"])
        else:
            # 1995 nonbond parameters are intentionally comparison-only. No
            # source switch can replace 1997's reoptimized Si/O or base C/H.
            inherited.append(row["id"])
    return {
        "labels": query,
        "primary_source_rows": primary,
        "prior_silane_candidates": inherited,
        "partial_source_rows": partial,
        "wildcard_candidates": wildcard,
        "comparison_only_source_rows": comparison_only,
        "selected_row": None,
    }


@dataclass(frozen=True)
class SiloxaneParameterAssessment:
    json_text: str

    def __post_init__(self):
        try:
            p = unpack(self.json_text)
            valid = (
                p["schema"] == ASSESSMENT_SCHEMA
                and p["identity"]
                == digest({k: v for k, v in p.items() if k != "identity"})
                and p["status"] in ("parameter_incomplete", "typing_ambiguous")
                and p["profile_identity"] is None
                and p["parameter_assignment"] is None
                and p["exact_charge_vector"] is None
                and p["execution_authorized"] is False
                and p["acceptance_gate"] == 1
                and p["limitations"] == LIMITS
            )
        except INPUT_ERRORS as error:
            raise SiloxaneParameterError(
                "source_incompatible", "malformed_graph_assessment"
            ) from error
        if not valid:
            raise SiloxaneParameterError(
                "source_incompatible", "graph_assessment_authority_changed"
            )

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self, system, supplement, *, references):
        rebuilt = assess_siloxane_parameters(system, supplement, references=references)
        if self != rebuilt:
            raise SiloxaneParameterError(
                "source_incompatible", "stale_or_mutated_graph_assessment"
            )

    def require_executable(self):
        raise SiloxaneParametersBlockedError(self)


def assess_siloxane_parameters(system, supplement, *, references):
    """Associate real rows with final-graph requests; never fill or assign gaps."""
    supplement.verify_sources(references)
    catalog = supplement.payload
    rows = {r["id"]: r for r in catalog["rows"]}
    try:
        required = fixture_requirements(system, "class2")
        graph = required["graph"]
        atoms = {s["id"]: s for s in graph["sites"]}
        neighbors = {i: [] for i in atoms}
        for bond in graph["bonds"]:
            a, b = bond["sites"]
            neighbors[a].append(b)
            neighbors[b].append(a)
            required["requests"][f"bond_increment:{a},{b}"] = {
                "sites": [a, b],
                "dependencies": [],
            }
        unsupported = []
        typing_diagnostics = []
        for i, atom in atoms.items():
            element = atom["element"]
            near = [atoms[j]["element"] for j in neighbors[i]]
            valid = (
                (element == "Si" and len(near) == 4 and set(near) <= {"H", "C", "O"})
                or (
                    element == "O"
                    and len(near) == 2
                    and "Si" in near
                    and set(near) <= {"Si", "H"}
                )
                or (element == "C" and len(near) == 4 and set(near) <= {"H", "C", "Si"})
                or (element == "H" and len(near) == 1 and near[0] in {"Si", "C", "O"})
            )
            if (
                not valid
                or atom["formal_charge"]
                or atom["metadata"]["radical_electrons"]
            ):
                unsupported.append(i)
                typing_diagnostics.append(
                    {
                        "site": i,
                        "element": element,
                        "neighbors": near,
                        "formal_charge": atom["formal_charge"],
                        "radical_electrons": atom["metadata"]["radical_electrons"],
                        "reason": "neutral_saturated_siloxane_label_not_established_for_environment",
                    }
                )
        if (
            not any(a["element"] == "Si" for a in atoms.values())
            or not any(a["element"] == "O" for a in atoms.values())
            or any(b["order"] != 1 or b["aromatic"] for b in graph["bonds"])
        ):
            unsupported = sorted(atoms)
            typing_diagnostics.append(
                {
                    "sites": unsupported,
                    "reason": "requires_Si_and_O_with_explicit_H_and_only_single_nonaromatic_bonds",
                }
            )
    except INPUT_ERRORS as error:
        raise SiloxaneParameterError(
            "typing_ambiguous", "final_graph_invalid"
        ) from error
    associations, missing = {}, {}
    for key, request in required["requests"].items():
        family = key.split(":")[0]
        match = _candidates(rows.values(), family, request["sites"], atoms, neighbors)
        if match["primary_source_rows"]:
            state = "source_row_present_dependencies_unresolved"
        elif match["partial_source_rows"]:
            state = "partial_source_row"
        elif match["prior_silane_candidates"]:
            state = "prior_silane_row_requires_inheritance_verification"
        elif match["wildcard_candidates"]:
            state = "wildcard_precedence_unverified"
        elif match["comparison_only_source_rows"]:
            state = "comparison_only_source_role"
        else:
            state = "source_row_absent"
        reasons = {
            "source_row_present_dependencies_unresolved": "source_numbers_present; typing, equilibrium_dependencies_and_model_policies_not_qualified",
            "partial_source_row": "printed_coefficient_cells_blank; zero_inference_forbidden",
            "prior_silane_row_requires_inheritance_verification": "1995_row_present; exact_type_mapping_and_term_specific_1997_inheritance_not_verified",
            "wildcard_precedence_unverified": "literal_wildcard_row_present; applicability_and_precedence_not_verified",
            "comparison_only_source_role": "1995_nonbond_row_is_not_1997_base_or_reoptimized_nonbond_evidence; transfer_forbidden",
            "source_row_absent": "no_matching_row_in_the_two_papers; external_base_source_or_explicit_applicability_evidence_required",
        }
        associations[key] = {**request, **match, "state": state}
        missing[key] = {
            "state": state,
            "labels": match["labels"],
            "equilibrium_dependencies": request["dependencies"],
            "blank_fields": {
                rid: rows[rid]["missing_fields"] for rid in match["partial_source_rows"]
            },
            "reason": reasons[state],
        }
    result = {
        "schema": ASSESSMENT_SCHEMA,
        "supplement_identity": supplement.identity,
        "requirements": required,
        "unsupported_sites": unsupported,
        "typing_diagnostics": typing_diagnostics,
        "row_associations": associations,
        "missing_dependencies": missing,
        "unresolved_policies": catalog["unresolved_policies"],
        "status": "typing_ambiguous" if unsupported else "parameter_incomplete",
        "profile_identity": None,
        "parameter_assignment": None,
        "exact_charge_vector": None,
        "execution_authorized": False,
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    result["identity"] = digest(result)
    return SiloxaneParameterAssessment(pack(result))
