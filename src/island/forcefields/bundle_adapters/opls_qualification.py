"""J87: one explicitly supplied OPLS profile, exact fixture bytes, no activation.

Reuses J84 source capture/audit, J85 individual completeness and J86 scope gates.
No runtime registry, source discovery, scientific execution or numeric defaults.
"""

import json
from pathlib import Path

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.profiles import BASE_CONVENTIONS, BASE_TERMS
from island.forcefields.assignment_core.schema import data, is_hash, require, unique

from .intake import _read, capture
from .qualification import qualify
from .schema import keys, loads, relative_path, sha, text
from .scoped_qualification import QUALIFIED, QualificationRequest, ScopedQualification

REQUEST = "island.j87.opls-qualification-request.v1"
FIXTURE = "island.j87.exact-graph-fixture.v1"
CONVENTIONS = (*BASE_CONVENTIONS, "energy_convention")


def request_fields(request):
    keys(
        request,
        (
            "schema",
            "family",
            "profile_id",
            "graph_identities",
            "fixtures",
            "expected_version",
            "expected_conventions",
            "expected_forms",
        ),
        "opls_request_schema",
    )
    require(request["schema"] == REQUEST, "opls_request_schema")
    require(request["family"] == "OPLS-AA", "J87_OPLS_AA_only")
    # Reuse the immutable J86 declaration and its duplicate/identity checks.
    declaration = QualificationRequest(
        request["family"], request["profile_id"], request["graph_identities"]
    )
    require(
        isinstance(request["fixtures"], list) and request["fixtures"],
        "exact_fixture_files_required",
    )
    for f in request["fixtures"]:
        keys(
            f,
            ("id", "path", "sha256", "size", "graph_identity"),
            "graph_fixture_manifest_schema",
        )
        text(f["id"], "fixture_identity_required")
        relative_path(f["path"])
        require(
            is_hash(f["sha256"]) and type(f["size"]) is int and f["size"] > 0,
            "fixture_hash_or_size_required",
        )
    for key in ("id", "path", "sha256", "graph_identity"):
        unique([f[key] for f in request["fixtures"]], "duplicate_fixture_" + key)
    require(
        tuple(sorted(f["graph_identity"] for f in request["fixtures"]))
        == declaration.graph_identities,
        "fixture_graph_identity_coverage_mismatch",
    )
    require(
        request["expected_version"] is None
        or isinstance(request["expected_version"], str),
        "expected_version_schema",
    )
    require(
        isinstance(request["expected_conventions"], dict)
        and set(request["expected_conventions"]) <= set(CONVENTIONS),
        "expected_conventions_schema",
    )
    for value in request["expected_conventions"].values():
        text(value, "expected_convention_value_required")
    require(isinstance(request["expected_forms"], list), "expected_forms_schema")
    return declaration


def capture_request(path, expected_sha256):
    """Read only the explicit request and files named in its hash-bound manifest."""
    raw = _read(path, expected_sha256)
    request = loads(raw)
    request_fields(request)
    root = Path(path).parent.resolve()
    files = []
    for f in sorted(request["fixtures"], key=lambda f: f["path"]):
        candidate = root
        for part in relative_path(f["path"]).parts:
            candidate = candidate / part
            require(not candidate.is_symlink(), "fixture_path_symlink_forbidden")
        require(candidate.resolve().is_relative_to(root), "fixture_root_escape")
        source = _read(candidate, f["sha256"])
        require(len(source) == f["size"], "fixture_size_mutation")
        files.append({"path": f["path"], "raw_text": source.decode("utf-8")})
    return {
        "request_raw": raw.decode("utf-8"),
        "request_sha256": expected_sha256,
        "files": files,
    }


def validate_request(captured, normalize_fixture):
    """Pure replay: verify exact bytes, normalize final graphs, compare identities."""
    keys(
        captured, ("request_raw", "request_sha256", "files"), "captured_request_schema"
    )
    require(
        is_hash(captured["request_sha256"])
        and sha(captured["request_raw"].encode()) == captured["request_sha256"],
        "captured_request_hash_mutation",
    )
    request = loads(captured["request_raw"])
    declaration = request_fields(request)
    require(isinstance(captured["files"], list), "captured_fixtures_schema")
    unique([f["path"] for f in captured["files"]], "duplicate_captured_fixture")
    files = {f["path"]: f for f in captured["files"]}
    require(
        set(files) == {f["path"] for f in request["fixtures"]},
        "fixture_file_coverage_mismatch",
    )
    graphs = []
    for f in sorted(request["fixtures"], key=lambda f: f["id"]):
        cap = files[f["path"]]
        keys(cap, ("path", "raw_text"), "captured_fixture_schema")
        raw = cap["raw_text"].encode()
        require(
            sha(raw) == f["sha256"] and len(raw) == f["size"],
            "captured_fixture_hash_mutation",
        )
        fixture = loads(raw)
        keys(fixture, ("schema", "id", "spec"), "exact_graph_fixture_schema")
        require(
            fixture["schema"] == FIXTURE and fixture["id"] == f["id"],
            "fixture_identity_mutation",
        )
        graph = normalize(normalize_fixture(fixture["spec"]))
        require(graph.identity == f["graph_identity"], "exact_graph_scope_mismatch")
        graphs.append(graph)
    return declaration, request, tuple(graphs)


def capture_opls(
    manifest=None,
    manifest_sha256=None,
    audit=None,
    audit_sha256=None,
    *,
    capture_timestamp,
):
    """Reject another family before opening any of its parameter source files."""
    if manifest is not None:
        m = loads(_read(manifest, manifest_sha256))
        require(m["family"] == "OPLS-AA", "J87_OPLS_AA_only")
    return capture(
        manifest,
        manifest_sha256,
        audit,
        audit_sha256,
        capture_timestamp=capture_timestamp,
    )


def qualify_opls(request_capture, source_capture, normalize_fixture):
    """No IO/execution: return source-bound exact-graph coverage, never acceptance."""
    declaration, request, graphs = validate_request(request_capture, normalize_fixture)
    m = source_capture["manifest"]
    require(m is None or m["family"] == "OPLS-AA", "J87_OPLS_AA_only")
    require(
        m is None or m["profile_id"] == declaration.profile_id,
        "declared_profile_identity_mismatch",
    )
    q = qualify(source_capture, graphs, family="OPLS-AA")
    # A blocked J84 adapter deliberately supplies an empty coverage contract,
    # whose ID is not the requested source profile. Preserve that diagnostic;
    # it is never a substitute source. Only complete candidates enter J86's
    # strict completeness guard (which also requires the actual profile ID).
    if q.decision == QUALIFIED:
        ScopedQualification(declaration, q)
    evidence = json.loads(q.evidence_graph_json)
    records = evidence["adapter"].get("records", [])
    actual = {
        r["record"]["value"]["name"]: r["record"]["value"]["value"]
        for r in records
        if r["record"]["kind"] == "convention"
    }
    # Expected semantics are independently pinned in the request, not guessed from
    # a family name. No unit conversion or numerical reinterpretation is attempted.
    comparisons, blockers = {}, []
    if m is not None:
        expected_version = request["expected_version"]
        if not expected_version or expected_version.strip().lower() in (
            "unknown",
            "unresolved",
            "default",
        ):
            blockers.append("explicit_expected_opls_version_required")
        elif expected_version != m["version"]:
            blockers.append("opls_version_mismatch")
        for name in CONVENTIONS:
            expected = request["expected_conventions"].get(name)
            valid = expected is not None and expected.strip().lower() not in (
                "unknown",
                "unresolved",
                "default",
                "missing",
            )
            status = (
                "match"
                if valid and actual.get(name) == expected
                else "missing_expected_convention"
                if not valid
                else "mismatch"
            )
            comparisons[name] = {
                "expected": expected,
                "source": actual.get(name),
                "status": status,
            }
            if status != "match":
                blockers.append("request_source_convention_" + status + ":" + name)
        expected_forms = request["expected_forms"]
        form_families = [f["family"] for f in expected_forms]
        unique(form_families, "duplicate_expected_form")
        if set(form_families) != set(BASE_TERMS):
            blockers.append("complete_expected_class_i_forms_required")
        if sorted(expected_forms, key=lambda f: f["family"]) != sorted(
            m["functional_form_convention"], key=lambda f: f["family"]
        ):
            blockers.append("request_source_equation_units_or_ordering_mismatch")
    decision = (
        "family_scoped_qualification_ready_runtime_blocked"
        if q.decision == QUALIFIED
        else q.decision
    )
    if q.source_authorized and blockers:
        decision = "parameter_bundle_incomplete"
    locations = [
        {
            "source_file": r["file"],
            "source_sha256": r["source_sha256"],
            "pointer": r["pointer"],
            "record_id": r["record"]["record_id"],
            "kind": r["record"]["kind"],
            "location": r["record"]["location"],
            "mapped_parameter_id": r["record"]["value"].get("id")
            if r["record"]["kind"] == "parameter"
            else None,
            "claim_sha256": r["claim_sha256"],
        }
        for r in records
    ]
    by_parameter = {
        r["mapped_parameter_id"]: r for r in locations if r["kind"] == "parameter"
    }
    graph_reports = json.loads(q.graph_reports_json)
    for report in graph_reports:
        report["source_qualification_decision"] = report["decision"]
        report["decision"] = (
            "source_bundle_required"
            if not q.source_authorized
            else "parameter_bundle_incomplete"
            if blockers or report["blocking_reasons"]
            else "family_scoped_qualification_ready_runtime_blocked"
        )
        report["blocking_reasons"] = sorted(
            set(report["blocking_reasons"]) | set(blockers)
        )
        report["matched_source_records"] = [
            {"request_id": k, **by_parameter[v]}
            for k, record_ids in report["core_report"]["matched_records"]
            for v in record_ids
        ]
    return {
        "decision": decision,
        "request": data(declaration),
        "fixture_files": request["fixtures"],
        "scope_of_qualification": "declared_exact_graphs_only",
        "chemical_family_transfer": "not_established",
        "base_qualification": data(q),
        "source_evidence": evidence,
        "source_locations": locations,
        "graph_reports": graph_reports,
        "convention_comparison": {
            "performed": m is not None,
            "expected_version": request["expected_version"],
            "source_version": m["version"] if m else None,
            "fields": comparisons,
            "blocking_reasons": sorted(blockers),
        },
        "source_hashes": q.source_hashes,
        "blocking_reasons": sorted(set(q.blocking_reasons) | set(blockers)),
        "automatic_assignment_enabled": False,
        "runtime_registration": False,
        "accepted_coefficients": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
