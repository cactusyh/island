"""J86 exact-graph qualification sets and family-scoped ranking, never activation.

The J85 per-profile qualifier is reused unchanged. Its historical aggregate and
cross-family ranking are deliberately not used by this workflow.
"""

import json
from dataclasses import dataclass, field
from decimal import Decimal

from island.forcefields.assignment_core.profiles import profile_contract
from island.forcefields.assignment_core.schema import identity, require, unique

from .qualification import BundleQualification, qualify

QUALIFIED = "authorized_bundle_qualified_runtime_blocked"
ALIASES = {
    "opls": "OPLS-AA",
    "opls-aa": "OPLS-AA",
    "gaff": "GAFF",
    "gaff2": "GAFF2",
    "pcff": "PCFF",
    "compass": "COMPASS",
}


def canonical_family(value):
    require(
        isinstance(value, str) and value.lower() in ALIASES, "unknown_requested_family"
    )
    return ALIASES[value.lower()]


@dataclass(frozen=True)
class QualificationRequest:
    family: str
    profile_id: str
    graph_identities: tuple[str, ...]

    def __post_init__(self):
        require(
            self.family == canonical_family(self.family), "canonical_family_required"
        )
        require(
            isinstance(self.profile_id, str) and self.profile_id.strip(),
            "profile_identity_required",
        )
        ids = tuple(self.graph_identities)
        require(ids, "declared_exact_graphs_required")
        require(
            all(
                isinstance(g, str)
                and len(g) == 64
                and set(g) <= set("0123456789abcdef")
                for g in ids
            ),
            "invalid_graph_identity",
        )
        unique(ids, "duplicate_declared_graph")
        object.__setattr__(self, "graph_identities", tuple(sorted(ids)))

    @property
    def key(self):
        return self.family, self.profile_id


@dataclass(frozen=True)
class ScopedQualification:
    request: QualificationRequest
    qualification: BundleQualification

    def __post_init__(self):
        q, r = self.qualification, self.request
        require(q.family == r.family, "cross_family_qualification_forbidden")
        require(
            q.bundle_identity is None or q.profile_id == r.profile_id,
            "declared_profile_identity_mismatch",
        )
        reports = json.loads(q.graph_reports_json)
        require(
            tuple(sorted(g["graph_identity"] for g in reports)) == r.graph_identities,
            "declared_graph_scope_mismatch",
        )
        require(
            not q.automatic_assignment_enabled
            and not q.production_validated
            and not q.runtime_registration
            and q.accepted_coefficients == 0
            and q.simulation_readiness == "not_established",
            "runtime_gate_mutation",
        )
        require(
            q.functional_class == profile_contract(r.family).functional_class,
            "functional_class_mismatch",
        )
        require(
            type(q.source_authorized) is bool, "source_authorization_must_be_boolean"
        )
        if q.source_authorized:
            require(
                q.bundle_identity is not None
                and q.source_hashes
                and q.source_status == "audited_bundle_candidate",
                "authorized_source_identity_required",
            )
        dimensions = {
            "atom_typing",
            "bond",
            "angle",
            "torsion",
            "improper",
            "charge",
            "nonbonded",
            "interaction_14",
            "class_ii_cross_terms",
            "topology",
        }
        complete = bool(reports) and all(
            g["decision"] == QUALIFIED
            and not g["blocking_reasons"]
            and set(g["coverage"]) == dimensions
            and all(
                c["required"] == c["covered"]
                and c["status"] in ("complete", "not_required_by_class_i_contract")
                for c in g["coverage"].values()
            )
            and g["core_report"]["status"] == "coverage_complete_not_activated"
            for g in reports
        )
        expected = (
            "source_bundle_required"
            if not q.source_authorized
            else QUALIFIED
            if complete and not q.blocking_reasons
            else "parameter_bundle_incomplete"
        )
        require(q.decision == expected, "inconsistent_per_profile_decision")


def qualify_request(captured, graphs, *, request):
    """Qualify one explicitly declared profile against exactly its declared graphs."""
    graphs = tuple(graphs)
    require(
        tuple(sorted(g.identity for g in graphs)) == request.graph_identities,
        "declared_graph_scope_mismatch",
    )
    return ScopedQualification(
        request, qualify(captured, graphs, family=request.family)
    )


@dataclass(frozen=True)
class QualificationSetReport:
    decision: str
    qualification_set: tuple[QualificationRequest, ...]
    profiles: tuple[ScopedQualification, ...]
    qualified_profiles: tuple[tuple[str, str], ...]
    missing_or_unauthorized_profiles: tuple[tuple[str, str], ...]
    incomplete_profiles: tuple[tuple[str, str], ...]
    scope_of_qualification: str = field(
        default="declared_exact_graphs_only", init=False
    )
    chemical_family_transfer: str = field(default="not_established", init=False)
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)
    accepted_coefficients: int = field(default=0, init=False)
    runtime_registration: bool = field(default=False, init=False)

    @property
    def replay_identity(self):
        return identity(self)


def aggregate_qualifications(results, *, qualification_set):
    """No omitted, extra or duplicate required slot may be hidden by another slot.

    A partial set means some, but not all, declared profiles have authorized
    sources; incomplete authorized profiles remain visible independently.
    """
    required, results = tuple(qualification_set), tuple(results)
    require(required, "declared_qualification_set_required")
    unique([r.key for r in required], "duplicate_required_profile")
    unique([r.request.key for r in results], "duplicate_profile_result")
    expected = {r.key: r for r in required}
    actual = {r.request.key: r for r in results}
    require(expected.keys() == actual.keys(), "qualification_set_coverage_mismatch")
    for key, request in expected.items():
        require(
            actual[key].request == request, "qualification_set_graph_scope_mismatch"
        )
        # Re-check even when constructed by a caller rather than qualify_request.
        ScopedQualification(request, actual[key].qualification)
    ordered = tuple(actual[k] for k in sorted(actual))
    missing = tuple(
        r.request.key for r in ordered if not r.qualification.source_authorized
    )
    incomplete = tuple(
        r.request.key
        for r in ordered
        if r.qualification.source_authorized and r.qualification.decision != QUALIFIED
    )
    complete = tuple(
        r.request.key for r in ordered if r.qualification.decision == QUALIFIED
    )
    decision = (
        "source_bundle_required"
        if len(missing) == len(required)
        else "partial_bundle_set"
        if missing
        else "parameter_bundle_incomplete"
        if incomplete
        else "family_scoped_qualification_ready_runtime_blocked"
    )
    return QualificationSetReport(
        decision,
        tuple(expected[k] for k in sorted(expected)),
        ordered,
        complete,
        missing,
        incomplete,
    )


def rank_qualifications(results, *, requested_family):
    """Rank only candidates for one explicit family and the same exact graph set.

    Mixed input is an error; callers must make separate requests for each family.
    Aliases such as forcefield='opls' mean OPLS-AA only, never another family.
    """
    family = canonical_family(requested_family)
    results = tuple(results)
    require(all(r.family == family for r in results), "cross_family_ranking_forbidden")
    unique(
        [(r.family, r.bundle_identity) for r in results], "duplicate_candidate_bundle"
    )
    scopes = {
        tuple(sorted(g["graph_identity"] for g in json.loads(r.graph_reports_json)))
        for r in results
    }
    require(len(scopes) <= 1, "candidate_graph_scopes_differ")

    def key(result):
        reports = json.loads(result.graph_reports_json)
        required = sum(c["required"] for r in reports for c in r["coverage"].values())
        covered = sum(c["covered"] for r in reports for c in r["coverage"].values())
        return (
            result.decision != QUALIFIED,
            -Decimal(covered) / max(1, required),
            len(result.blocking_reasons),
            result.profile_id,
            result.bundle_identity or "",
        )

    return tuple(sorted(results, key=key))
