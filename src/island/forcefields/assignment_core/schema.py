"""Immutable, coordinate-free assignment contracts. No runtime force-field registry."""

import json
from dataclasses import asdict, dataclass, field, is_dataclass
from hashlib import sha256
from math import isfinite

EVIDENCE_STATUSES = (
    "verified_stock",
    "verified_extension",
    "local_patch",
    "reported_claim",
    "derived_provisional",
    "fallback",
    "missing",
    "blocked",
)
FAMILIES = ("OPLS-AA", "GAFF", "GAFF2", "PCFF", "COMPASS")
CROSS_TERMS = (
    "bond_bond",
    "bond_angle",
    "angle_angle",
    "bond_bond_13",
    "end_bond_torsion",
    "middle_bond_torsion",
    "angle_torsion",
    "angle_angle_torsion",
    "torsion_torsion",
)
ARITY = {
    "bond": 2,
    "angle": 3,
    "torsion": 4,
    "improper": 4,
    "nonbonded": 1,
    "charge": 1,
    "bond_bond": 3,
    "bond_angle": 3,
    **{k: 4 for k in CROSS_TERMS[2:-1]},
    "torsion_torsion": 5,
}
READINESS = {
    "automatic_assignment_enabled": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "accepted_coefficients": 0,
}


def data(value):
    """Canonical JSON-compatible view; reject NaN and non-JSON metadata."""

    def encode(obj):
        if is_dataclass(obj):
            return asdict(obj)
        raise TypeError("non_JSON_schema_value:" + type(obj).__name__)

    return json.loads(json.dumps(value, default=encode, allow_nan=False))


def identity(value):
    return sha256(
        json.dumps(
            data(value), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def unique(values, reason):
    require(len(values) == len(set(values)), reason)


def freeze_sequences(instance):
    """Detach caller-owned lists; schema containers are recursively immutable."""
    from dataclasses import fields

    def freeze(value):
        require(not isinstance(value, dict), "mapping_not_allowed_in_immutable_schema")
        return (
            tuple(freeze(v) for v in value)
            if isinstance(value, (tuple, list))
            else value
        )

    for f in fields(instance):
        object.__setattr__(instance, f.name, freeze(getattr(instance, f.name)))


def claim_identity(kind, value):
    return identity((kind, value))


def is_hash(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(c in "0123456789abcdef" for c in value)
    )


@dataclass(frozen=True)
class Branch:
    """One single-dummy branch attached to a pre-edit heavy-atom site ID."""

    anchor_site: int
    psmiles: str


@dataclass(frozen=True)
class Crosslink:
    site1: int
    site2: int
    order: int = 1


@dataclass(frozen=True)
class PolymerSpec:
    repeat_units: tuple[tuple[str, str], ...]
    degree_of_polymerization: int
    sequence: tuple[str, ...] = ()
    connection_points: tuple[tuple[str, int, int], ...] = ()
    tacticity: str | None = None
    stereo_seed: int = 2026
    end_groups: tuple[str | None, str | None] = (None, None)
    branches: tuple[Branch, ...] = ()
    crosslinks: tuple[Crosslink, ...] = ()
    periodic_topology: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        freeze_sequences(self)
        require(
            type(self.degree_of_polymerization) is int
            and self.degree_of_polymerization > 0,
            "invalid_degree_of_polymerization",
        )
        require(bool(self.repeat_units), "empty_repeat_library")
        unique([n for n, _ in self.repeat_units], "duplicate_repeat_identity")
        require(all(n and p for n, p in self.repeat_units), "empty_repeat_definition")
        require(len(self.end_groups) == 2, "two_end_groups_required")
        require(
            not self.sequence or len(self.sequence) == self.degree_of_polymerization,
            "sequence_DP_mismatch",
        )
        require(
            self.sequence or len(self.repeat_units) == 1,
            "copolymer_requires_explicit_sequence",
        )
        require(
            set(self.sequence) <= {n for n, _ in self.repeat_units},
            "undefined_repeat_in_sequence",
        )
        unique([c[0] for c in self.connection_points], "duplicate_connection_identity")
        unique([k for k, _ in self.periodic_topology], "duplicate_periodic_metadata")
        require(
            self.tacticity in (None, "isotactic", "syndiotactic", "atactic"),
            "unsupported_tacticity",
        )
        unique([(b.anchor_site, b.psmiles) for b in self.branches], "duplicate_branch")
        unique(
            [tuple(sorted((c.site1, c.site2))) for c in self.crosslinks],
            "duplicate_crosslink",
        )


@dataclass(frozen=True)
class ChemicalGraph:
    """Final explicit-H graph, ordered site IDs and stereochemistry bound together."""

    spec: PolymerSpec
    molecule_json: str
    atom_order: tuple[int, ...]
    architecture: str
    blockers: tuple[str, ...] = ()
    normalizer_version: str = "island.j82.topology.v1"
    rdkit_version: str = ""

    def __post_init__(self):
        freeze_sequences(self)

    @property
    def identity(self):
        return identity(self)


@dataclass(frozen=True)
class EvidenceRecord:
    id: str
    status: str
    source_sha256: str
    citation: str
    location: str
    source_family: str
    synthetic: bool = False
    audited_claims: tuple[str, ...] = ()
    resolution: str = "all_atom"

    def __post_init__(self):
        freeze_sequences(self)
        unique(self.audited_claims, "duplicate_audited_claim")
        require(
            all(is_hash(c) for c in self.audited_claims), "invalid_audited_claim_hash"
        )
        require(type(self.synthetic) is bool, "evidence_kind_must_be_boolean")
        require(
            self.resolution in ("all_atom", "united_atom", "coarse_grained"),
            "unknown_evidence_resolution",
        )
        require(self.status in EVIDENCE_STATUSES, "unknown_evidence_status")
        require(
            bool(self.id and self.citation and self.location), "incomplete_provenance"
        )
        require(is_hash(self.source_sha256), "invalid_source_hash")
        require(self.source_family in FAMILIES, "unknown_evidence_family")


@dataclass(frozen=True)
class TypingRule:
    id: str
    atom_type: str
    smarts: str
    evidence_id: str
    priority: int = 0
    overrides: tuple[str, ...] = ()

    def __post_init__(self):
        freeze_sequences(self)
        require(
            bool(self.id and self.atom_type and self.smarts and self.evidence_id),
            "incomplete_typing_rule",
        )
        require(type(self.priority) is int, "noninteger_typing_priority")
        unique(self.overrides, "duplicate_rule_override")


@dataclass(frozen=True)
class Convention:
    name: str
    value: str
    evidence_id: str

    def __post_init__(self):
        require(
            all(
                isinstance(v, str) and v.strip()
                for v in (self.name, self.value, self.evidence_id)
            ),
            "incomplete_convention",
        )


@dataclass(frozen=True)
class TermForm:
    family: str
    equation: str
    fields: tuple[str, ...]
    units: tuple[str, ...]
    ordering: str
    equilibrium_fields: tuple[str, ...] = ()

    def __post_init__(self):
        freeze_sequences(self)
        require(self.family in ARITY, "unknown_parameter_family")
        require(
            bool(self.equation and self.ordering and self.fields), "incomplete_equation"
        )
        require(
            len(self.fields) == len(self.units) and all(self.units), "incomplete_units"
        )
        unique(self.fields, "duplicate_coefficient_identity")
        require(
            set(self.equilibrium_fields) <= set(self.fields),
            "unknown_equilibrium_field",
        )


@dataclass(frozen=True)
class ParameterRecord:
    id: str
    profile_family: str
    family: str
    atom_types: tuple[str, ...]
    environment: str
    form: TermForm
    values: tuple[str, ...]
    evidence_id: str
    convention_identity: str
    dependency_keys: tuple[str, ...] = ()
    applicability: str = "required"
    absence_reason: str = ""
    source_location: str = ""
    equilibrium_references: tuple[tuple[str, tuple[tuple[str, str], ...]], ...] = ()

    def __post_init__(self):
        freeze_sequences(self)
        require(
            bool(self.id and self.source_location), "parameter_source_location_missing"
        )
        unique(self.dependency_keys, "duplicate_dependency_key")
        require(
            self.family in ARITY and self.family == self.form.family,
            "parameter_family_mismatch",
        )
        require(
            len(self.atom_types) == ARITY[self.family] and all(self.atom_types),
            "parameter_arity_mismatch",
        )
        require(
            not any("*" in t for t in self.atom_types), "wildcard_parameter_forbidden"
        )
        require(
            is_hash(self.environment) and is_hash(self.convention_identity),
            "unbound_environment_or_conventions",
        )
        require(
            self.applicability in ("required", "not_applicable"),
            "invalid_applicability",
        )
        if self.applicability == "not_applicable":
            require(
                self.family in (*CROSS_TERMS, "improper")
                and not self.values
                and bool(self.absence_reason),
                "unevidenced_term_omission",
            )
        else:
            require(
                len(self.values) == len(self.form.fields),
                "missing_coefficient_no_zero_fill",
            )
            require(
                all(isinstance(v, str) and isfinite(float(v)) for v in self.values),
                "invalid_numeric_source_value",
            )


@dataclass(frozen=True)
class ForceFieldProfile:
    id: str
    family: str
    functional_class: str
    rules: tuple[TypingRule, ...] = ()
    forms: tuple[TermForm, ...] = ()
    conventions: tuple[Convention, ...] = ()
    evidence: tuple[EvidenceRecord, ...] = ()
    parameters: tuple[ParameterRecord, ...] = ()
    graph_scope: tuple[str, ...] = ()
    topology_support: tuple[str, ...] = ("finite_linear",)
    purpose: str = "contract_only"
    resolution: str = "all_atom"

    def __post_init__(self):
        freeze_sequences(self)
        require(self.family in FAMILIES, "unknown_forcefield_family")
        require(self.resolution == "all_atom", "atomistic_profile_required")
        require(
            all(e.resolution == self.resolution for e in self.evidence),
            "resolution_overlay_forbidden",
        )
        require(
            self.functional_class != "I"
            or not any(
                p.family in CROSS_TERMS for p in (*self.forms, *self.parameters)
            ),
            "class_i_cross_term_overlay",
        )
        require(
            self.functional_class
            == ("II" if self.family in ("PCFF", "COMPASS") else "I"),
            "forcefield_class_mismatch",
        )
        require(
            self.purpose
            in ("contract_only", "synthetic_control", "evidence_candidate"),
            "invalid_profile_purpose",
        )
        for rows, attr, label in (
            (self.rules, "id", "rule"),
            (self.forms, "family", "form"),
            (self.conventions, "name", "convention"),
            (self.evidence, "id", "evidence"),
            (self.parameters, "id", "parameter"),
        ):
            unique([getattr(r, attr) for r in rows], "duplicate_" + label + "_identity")
        require(
            all(p.profile_family == self.family for p in self.parameters),
            "cross_family_parameter_overlay",
        )
        require(
            all(e.source_family == self.family for e in self.evidence),
            "cross_family_evidence_overlay",
        )
        require(
            self.purpose != "contract_only" or not self.parameters,
            "numeric_parameters_in_contract_only_profile",
        )
        require(
            self.purpose != "synthetic_control"
            or all(e.synthetic for e in self.evidence),
            "scientific_evidence_in_synthetic_control",
        )
        require(
            self.purpose != "evidence_candidate"
            or not any(e.synthetic for e in self.evidence),
            "synthetic_evidence_in_scientific_candidate",
        )

    @property
    def identity(self):
        return identity(self)

    @property
    def convention_identity(self):
        return identity(self.conventions)


@dataclass(frozen=True)
class TypedGraph:
    graph_identity: str
    profile_identity: str
    assignments: tuple[tuple[int, str], ...]
    diagnostics: tuple[str, ...]
    rule_trace: tuple[tuple[int, tuple[str, ...]], ...]

    def __post_init__(self):
        freeze_sequences(self)

    def validate(self, graph, profile):
        require(
            self.graph_identity == graph.identity, "topology_changed_retyping_required"
        )
        require(
            self.profile_identity == profile.identity,
            "profile_changed_retyping_required",
        )


@dataclass(frozen=True)
class TermRequest:
    id: str
    family: str
    sites: tuple[int, ...]
    atom_types: tuple[str, ...]
    environment: str
    dependencies: tuple[str, ...] = ()

    @property
    def parameter_key(self):
        return identity((self.family, self.atom_types, self.environment))


@dataclass(frozen=True)
class MissingTerm:
    request_id: str
    family: str
    sites: tuple[int, ...]
    reason: str
    candidate_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class Coverage:
    required: int
    covered: int
    status: str


@dataclass(frozen=True)
class CompletenessReport:
    profile_id: str
    profile_identity: str
    functional_class: str
    graph_identity: str
    status: str
    coverage: tuple[tuple[str, Coverage], ...]
    missing_terms: tuple[MissingTerm, ...]
    blocking_reasons: tuple[str, ...]
    source_hashes: tuple[str, ...]
    matched_records: tuple[tuple[str, tuple[str, ...]], ...]
    replay_identity: str
    scientific_validation: str = "not_established"
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)


@dataclass(frozen=True)
class AssignmentDecision:
    status: str
    ranked_profiles: tuple[str, ...]
    reports: tuple[CompletenessReport, ...]
    blocking_reasons: tuple[str, ...]
    replay_identity: str
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)
    accepted_coefficients: int = field(default=0, init=False)
