"""Family-separated adapters into unchanged J82/J83 contracts. No conversion."""

import json
from dataclasses import dataclass, field, replace

from island.forcefields.assignment_core.profiles import (
    required_conventions,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    READINESS,
    Convention,
    EvidenceRecord,
    ForceFieldProfile,
    ParameterRecord,
    TermForm,
    TypingRule,
    claim_identity,
    data,
    identity,
    is_hash,
    require,
    unique,
)
from island.forcefields.coverage_registry import ForceFieldProfileRegistry

from .schema import (
    AUDIT,
    CAPTURE,
    ORIGINS,
    SOURCE,
    keys,
    loads,
    sha,
    text,
    timestamp,
    validate_manifest,
)

ADAPTERS = {
    f: f.lower().replace("-", "_") + ".canonical-json.v1"
    for f in ("OPLS-AA", "GAFF", "GAFF2", "PCFF", "COMPASS")
}


@dataclass(frozen=True)
class BundleIngestion:
    status: str
    bundle_identity: str | None
    adapter_identity: str | None
    profile: ForceFieldProfile | None
    trusted_evidence: tuple[EvidenceRecord, ...]
    source_hashes: tuple[str, ...]
    blocking_reasons: tuple[str, ...]
    evidence_graph_json: str
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)
    accepted_coefficients: int = field(default=0, init=False)
    runtime_registration: bool = field(default=False, init=False)

    @property
    def identity(self):
        return identity(self)

    def registry(self):
        return ForceFieldProfileRegistry(
            (self.profile,) if self.profile else (), self.trusted_evidence
        )


def _evidence(m, f):
    return EvidenceRecord(
        "source:" + m["family"] + ":" + f["sha256"],
        m["evidence_status"],
        f["sha256"],
        m["provenance"]["citation"] + "; release=" + m["version"],
        "bundle=" + m["bundle_id"] + "; file=" + f["path"],
        m["family"],
        m["synthetic"],
    )


def proposals(captured):
    """Untrusted review proposals, including exact mapped-claim hashes. Not an audit."""
    require(captured["schema"] == CAPTURE, "capture_schema")
    m = captured["manifest"]
    validate_manifest(m)
    files = {f["path"]: f for f in captured["files"]}
    require(
        set(files) == {f["path"] for f in m["files"]}
        and len(files) == len(captured["files"]),
        "capture_file_coverage",
    )
    evidence = tuple(
        _evidence(m, f) for f in sorted(m["files"], key=lambda f: f["path"])
    )
    e_by_path = {f["path"]: _evidence(m, f) for f in m["files"]}
    rows = []
    for f in sorted(m["files"], key=lambda f: f["path"]):
        cap = files[f["path"]]
        keys(cap, ("path", "sha256", "size", "raw_text"), "captured_file_schema")
        require(
            cap["sha256"] == f["sha256"] and cap["size"] == f["size"],
            "capture_source_hash_mismatch",
        )
        if f["role"] == "primary_artifact":
            require(
                cap["raw_text"] is None, "opaque_artifact_cannot_supply_parameter_rows"
            )
            continue
        text(cap["raw_text"], "source_record_bytes_missing")
        raw = cap["raw_text"].encode("utf-8")
        require(
            sha(raw) == f["sha256"] and len(raw) == f["size"],
            "captured_record_bytes_mutated",
        )
        src = loads(raw)
        keys(src, ("schema", "family", "version", "records"), "source_record_schema")
        require(
            src["schema"] == SOURCE
            and src["family"] == m["family"]
            and src["version"] == m["version"],
            "source_family_or_version_mismatch",
        )
        require(isinstance(src["records"], list), "records_require_list")
        for index, row in enumerate(src["records"]):
            keys(
                row,
                ("record_id", "kind", "location", "origin", "value"),
                "source_row_schema",
            )
            text(row["record_id"], "record_identity_missing")
            text(row["location"], "exact_source_location_missing")
            require(
                row["origin"] == "directly_reported",
                "inferred_analog_or_reconstructed_row_blocked",
            )
            rows.append(
                {
                    "file": f["path"],
                    "pointer": "/records/" + str(index),
                    "source_sha256": f["sha256"],
                    "record": row,
                    "evidence_id": e_by_path[f["path"]].id,
                }
            )
    unique([r["record"]["record_id"] for r in rows], "duplicate_source_record_identity")
    forms, conventions, rules, parameters, scopes = [], [], [], [], []
    for row in rows:
        r = row["record"]
        value, kind, eid = r["value"], r["kind"], row["evidence_id"]
        if kind == "form":
            keys(
                value,
                (
                    "family",
                    "equation",
                    "fields",
                    "units",
                    "ordering",
                    "equilibrium_fields",
                ),
                "functional_form_record_schema",
            )
            obj = TermForm(**value)
            forms.append(obj)
        elif kind == "convention":
            keys(value, ("name", "value"), "convention_record_schema")
            obj = Convention(**value, evidence_id=eid)
            conventions.append(obj)
        elif kind == "atom_type":
            keys(
                value,
                ("id", "atom_type", "smarts", "priority", "overrides"),
                "atom_type_record_schema",
            )
            obj = TypingRule(**value, evidence_id=eid)
            rules.append(obj)
        elif kind == "graph_scope":
            keys(value, ("graph_identity",), "graph_scope_record_schema")
            obj = value["graph_identity"]
            require(is_hash(obj), "graph_scope_hash_missing")
            scopes.append(obj)
        elif kind == "parameter":
            continue  # Convention hashes are computed only after all conventions exist.
        else:
            raise ValueError("unknown_source_record_kind:" + str(kind))
        row["mapped"] = obj
    p = ForceFieldProfile(
        m["profile_id"],
        m["family"],
        "II" if m["family"] in ("PCFF", "COMPASS") else "I",
        tuple(rules),
        tuple(forms),
        tuple(conventions),
        evidence,
        (),
        tuple(scopes),
        tuple(m["topology_support"]),
        "synthetic_control" if m["synthetic"] else "evidence_candidate",
    )
    for row in rows:
        if row["record"]["kind"] != "parameter":
            continue
        value = dict(row["record"]["value"])
        keys(
            value,
            (
                "id",
                "family",
                "atom_types",
                "environment",
                "form",
                "values",
                "dependency_keys",
                "applicability",
                "absence_reason",
                "equilibrium_references",
            ),
            "parameter_record_schema",
        )
        keys(
            value["form"],
            ("family", "equation", "fields", "units", "ordering", "equilibrium_fields"),
            "parameter_form_record_schema",
        )
        value["form"] = TermForm(**value["form"])
        obj = ParameterRecord(
            **value,
            profile_family=m["family"],
            evidence_id=row["evidence_id"],
            convention_identity=p.convention_identity,
            source_location=row["file"]
            + "#"
            + row["pointer"]
            + "; "
            + row["record"]["location"],
        )
        row["mapped"] = obj
        parameters.append(obj)
    p = replace(p, parameters=tuple(parameters))
    kind_map = {
        "atom_type": "typing_rule",
        "graph_scope": "graph_scope",
        "form": "form",
        "convention": "convention",
        "parameter": "parameter",
    }
    for row in rows:
        row["claim_sha256"] = claim_identity(
            kind_map[row["record"]["kind"]], row["mapped"]
        )
    return p, rows


def _conventions(m, p):
    values = {c.name: c.value for c in p.conventions}
    reasons = [
        "convention_missing:" + name
        for name in required_conventions(p)
        if name not in values
        or values[name].strip().lower() in ("unknown", "unresolved", "default")
    ]
    required = {
        **m["unit_convention"],
        "mixing_rule": m["mixing_rule"],
        "charge_origin": m["charge_convention"],
        "scaling_14_lj": m["policy_14"]["lj"],
        "scaling_14_coulomb": m["policy_14"]["coulomb"],
    }
    reasons += [
        "manifest_source_convention_mismatch:" + k
        for k, v in required.items()
        if values.get(k) != v
    ]
    expected = sorted(m["functional_form_convention"], key=lambda f: f["family"])
    if sorted(data(p.forms), key=lambda f: f["family"]) != expected:
        reasons.append("functional_form_units_or_ordering_mismatch")
    for family in required_families(p):
        if family not in {f.family for f in p.forms}:
            reasons.append("functional_form_missing:" + family)
    nonbonded = next((f for f in p.forms if f.family == "nonbonded"), None)
    if nonbonded and values.get("nonbonded_form") != nonbonded.equation:
        reasons.append("nonbonded_9_6_vs_12_6_or_equation_mismatch")
    return reasons


def _audit(m, captured, rows):
    a = captured["audit"]
    if a is None:
        return ["independent_source_audit_missing"]
    keys(
        a,
        (
            "schema",
            "manifest_identity",
            "reviewer",
            "review_timestamp",
            "authority_basis",
            "policy",
            "synthetic",
            "source_reviews",
        ),
        "audit_schema",
    )
    require(
        a["schema"] == AUDIT and a["manifest_identity"] == m["identity"],
        "audit_manifest_binding_mismatch",
    )
    for name in ("reviewer", "authority_basis"):
        text(a[name], "audit_authority_metadata_missing")
    timestamp(a["review_timestamp"])
    require(
        a["policy"] == "independent_source_review"
        and type(a["synthetic"]) is bool
        and a["synthetic"] == m["synthetic"],
        "audit_policy_or_scientific_kind_mismatch",
    )
    reviews = a["source_reviews"]
    require(isinstance(reviews, list), "audit_source_reviews_schema")
    require(
        len(reviews) == len(m["files"])
        and len({r["path"] for r in reviews}) == len(reviews),
        "audit_source_coverage",
    )
    for f in m["files"]:
        review = next((r for r in reviews if r["path"] == f["path"]), None)
        keys(
            review,
            ("path", "sha256", "provenance_verified", "license_verified", "claims"),
            "audit_source_review_schema",
        )
        require(review["sha256"] == f["sha256"], "audit_source_hash_mismatch")
        require(
            review["provenance_verified"] is True
            and review["license_verified"] is True,
            "source_authority_or_license_not_audited",
        )
        claims = review["claims"]
        for c in claims:
            keys(c, ("record_id", "claim_sha256"), "audit_claim_schema")
        expected = sorted(
            (r["record"]["record_id"], r["claim_sha256"])
            for r in rows
            if r["file"] == f["path"]
        )
        require(
            sorted((r["record_id"], r["claim_sha256"]) for r in claims) == expected,
            "audit_claim_set_mismatch",
        )
    return []


def adapt(captured):
    """Reconstruct only source-bound, independently audited candidates; never activate."""
    keys(
        captured,
        (
            "schema",
            "manifest",
            "manifest_raw",
            "manifest_sha256",
            "files",
            "audit",
            "audit_raw",
            "audit_sha256",
            "capture_timestamp",
        ),
        "capture_fields",
    )
    require(captured["schema"] == CAPTURE, "capture_schema")
    timestamp(captured["capture_timestamp"])
    m = captured["manifest"]
    if m is None:
        require(
            captured["manifest_sha256"] is None
            and captured["manifest_raw"] is None
            and not captured["files"]
            and captured["audit"] is None
            and captured["audit_raw"] is None
            and captured["audit_sha256"] is None,
            "source_claims_without_bundle",
        )
        return BundleIngestion(
            "source_bundle_missing",
            None,
            None,
            None,
            (),
            (),
            ("authorized_source_bundle_missing",),
            "{}",
        )
    require(is_hash(captured["manifest_sha256"]), "external_manifest_hash_missing")
    text(captured["manifest_raw"], "captured_manifest_bytes_missing")
    require(
        (captured["audit"] is None) == (captured["audit_sha256"] is None),
        "independent_audit_hash_missing",
    )
    require(
        sha(captured["manifest_raw"].encode()) == captured["manifest_sha256"]
        and loads(captured["manifest_raw"]) == m,
        "captured_manifest_bytes_mutated",
    )
    if captured["audit"] is not None:
        text(captured["audit_raw"], "captured_audit_bytes_missing")
        require(
            is_hash(captured["audit_sha256"])
            and sha(captured["audit_raw"].encode()) == captured["audit_sha256"]
            and loads(captured["audit_raw"]) == captured["audit"],
            "captured_audit_bytes_mutated",
        )
    else:
        require(captured["audit_raw"] is None, "unexpected_audit_bytes")
    p, rows = proposals(captured)
    reasons = []
    provenance = m["provenance"]
    if m["version"].strip().lower() in ("unknown", "unresolved", "default"):
        reasons.append("source_release_unresolved")
    if (
        provenance["origin"] not in ORIGINS
        or "polyjarvis"
        in (
            m["source_uri"]
            + provenance["upstream_identity"]
            + provenance["acquisition_method"]
        ).lower()
    ):
        reasons.append("unauthorized_analog_patch_or_source_origin")
    if any(
        provenance[k].strip().lower() in ("unknown", "unresolved", "default")
        for k in ("license", "citation", "upstream_identity")
    ):
        reasons.append("source_license_or_provenance_unresolved")
    if m["synthetic"] != (provenance["origin"] == "synthetic_control"):
        reasons.append("synthetic_source_cannot_be_promoted")
    if provenance["authorization"] != "authorized_for_local_readonly_ingestion":
        reasons.append("source_ingestion_not_authorized")
    if m["evidence_status"] not in (
        ("derived_provisional",)
        if m["synthetic"]
        else ("verified_stock", "verified_extension")
    ):
        reasons.append("source_evidence_status_unqualified")
    reasons += _conventions(m, p)
    reasons += _audit(m, captured, rows)
    evidence = []
    if not reasons:
        for e in p.evidence:
            claims = tuple(
                sorted(r["claim_sha256"] for r in rows if r["evidence_id"] == e.id)
            )
            evidence.append(replace(e, audited_claims=claims))
        p = replace(p, evidence=tuple(evidence))
    graph = {
        "bundle_identity": m["identity"],
        "manifest_sha256": captured["manifest_sha256"],
        "audit_sha256": captured["audit_sha256"],
        "source_files": m["files"],
        "records": [{k: data(v) for k, v in r.items() if k != "mapped"} for r in rows],
        "source_authority": "Independent audit supplied by caller; a hash or local copy alone proves no authenticity.",
        "convention_conversions": [],
        "capture_timestamp": captured["capture_timestamp"],
    }
    outcome = BundleIngestion(
        "source_parsed_audit_or_convention_blocked"
        if reasons
        else "audited_bundle_candidate",
        m["identity"],
        ADAPTERS[m["family"]],
        p if not reasons else None,
        tuple(evidence),
        tuple(sorted(f["sha256"] for f in m["files"])),
        tuple(sorted(set(reasons))),
        json.dumps(graph, sort_keys=True, separators=(",", ":")),
    )
    require(
        all(data(outcome)[k] == v for k, v in READINESS.items()), "readiness_mutation"
    )
    return outcome
