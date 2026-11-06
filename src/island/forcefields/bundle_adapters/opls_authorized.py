"""J88: hardened, metadata-only, exact J87 OPLS qualification; no activation.

Raw source records are ephemeral. Replay must recapture all external bytes.
The historical J84--J87 implementations are deliberately unchanged.
"""

from pathlib import Path

from island.forcefields.assignment_core.profiles import BASE_TERMS
from island.forcefields.assignment_core.schema import TermForm, is_hash, require, unique

from . import opls_qualification as j87
from .intake import _read
from .schema import keys, loads, text, validate_manifest

FORM_FIELDS = (
    "family",
    "equation",
    "fields",
    "units",
    "ordering",
    "equilibrium_fields",
)

ALLOWED_FIXTURES = {
    "amide": {
        "sha256": "a71b3e48b71a6172a0498c05e51b2536af7983946c255895934a129d66670337",
        "size": 217,
        "graph_identity": "d7f155a11801e9426015c54b1ab531820530323980266346aa85678bdce83209",
    },
    "aromatic": {
        "sha256": "31590b249e331241e3307e33fb45ea5c79f84a21726b1971aad4d019c36aade1",
        "size": 226,
        "graph_identity": "8ca0c3a01313cb94122cc1a438570d032d5ad2bf3a2821c0478d46f8881c5604",
    },
    "carbonyl": {
        "sha256": "dd1c61da0e27e6b0f0170cc4af32b2401cedfebc68cd48ddae9848854bd1824d",
        "size": 223,
        "graph_identity": "976a09e143a8f646f52c4787c6b37c9187e11ccd6f6c7adcb30976859f9672c9",
    },
    "ester": {
        "sha256": "995f1855c690a88a4c549e6a7b2172628f30a6bcae52333bc61a3f5a7a7d3d59",
        "size": 267,
        "graph_identity": "8e4f86935cc373ee00c371b8a3da18452844e9bd4c69c5756739e2cb42b73c4d",
    },
    "ether": {
        "sha256": "aecd79b287cd4b4ed2cdba2ef6e75547eeb8b35fbf0b1980911a8d23b681d442",
        "size": 263,
        "graph_identity": "450370367f10cf10d442e0e7b7c93e856ce4b0860fd1996c2e542883d3102ff3",
    },
    "halogenated": {
        "sha256": "c8116a92def31f518592885a175d595523264e02680ac05988dd9f9c18975332",
        "size": 231,
        "graph_identity": "29180d84c97042bd55f57bb48b6431a2eb1dcc10e9fc4857eab8b2e856b4e600",
    },
    "hydrocarbon": {
        "sha256": "9e636f57cacbbcea72a31721a6f405ade97677d6c0fd702130d538a96aa0631f",
        "size": 224,
        "graph_identity": "30c53d0353935dc75f5592b48b8fc59f8c5d9763ae60d2d2fcce1507c9890cb3",
    },
    "nitrile": {
        "sha256": "20a176ea6aa536a2d76d73e2de231b196869a7cbe398135fc25886ad3529580a",
        "size": 221,
        "graph_identity": "85335c735681d8fa0537ca09ae3153d536f68c7246d5e69e982303884b97d7f7",
    },
}


class InputBlockers(ValueError):
    """A deterministic collection of independent input defects, without raw rows."""

    def __init__(self, blockers):
        self.blockers = tuple(sorted(set(blockers)))
        super().__init__(";".join(self.blockers))


def validate_forms(forms, *, complete=False):
    """Validate the entire structured collection before consumers index any form."""
    require(isinstance(forms, list), "expected_forms_schema")
    families = []
    for form in forms:
        keys(form, FORM_FIELDS, "expected_form_schema")
        for name in ("family", "equation", "ordering"):
            text(form[name], "expected_form_text:" + name)
        require(form["family"] in BASE_TERMS, "expected_form_family")
        for name in ("fields", "units", "equilibrium_fields"):
            require(isinstance(form[name], list), "expected_form_array:" + name)
            for value in form[name]:
                text(value, "expected_form_array_item:" + name)
            if name != "units":
                unique(form[name], "duplicate_expected_form_field:" + name)
        TermForm(**form)
        families.append(form["family"])
    unique(families, "duplicate_expected_form")
    if complete:
        require(
            set(families) == set(BASE_TERMS), "complete_expected_class_i_forms_required"
        )


def request_fields(request):
    # Validate containers before J87's sorting/set/dataclass consumers.
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
    text(request["profile_id"], "profile_id_required")
    require(
        isinstance(request["graph_identities"], list)
        and all(is_hash(v) for v in request["graph_identities"]),
        "graph_identities_schema",
    )
    validate_forms(request["expected_forms"])
    require(isinstance(request["fixtures"], list), "exact_fixture_files_required")
    for f in request["fixtures"]:
        keys(
            f,
            ("id", "path", "sha256", "size", "graph_identity"),
            "graph_fixture_manifest_schema",
        )
        text(f["id"], "fixture_identity_required")
        require(is_hash(f["graph_identity"]), "fixture_graph_identity_schema")
        require(f["id"] in ALLOWED_FIXTURES, "fixture_outside_immutable_J87_scope")
        expected = ALLOWED_FIXTURES[f["id"]]
        require(
            {k: f[k] for k in expected} == expected,
            "J87_fixture_byte_or_graph_identity_mutation",
        )
    return j87.request_fields(request)


def capture_request(path, expected_sha256):
    # In particular, no fixture or source table is opened on malformed forms.
    request = loads(_read(path, expected_sha256))
    request_fields(request)
    captured = j87.capture_request(path, expected_sha256)
    request_fields(loads(captured["request_raw"]))
    return captured


def capture_opls(
    manifest=None,
    manifest_sha256=None,
    audit=None,
    audit_sha256=None,
    *,
    capture_timestamp,
):
    if manifest is not None:
        m = loads(_read(manifest, manifest_sha256))
        validate_manifest(m)
        require(m["family"] == "OPLS-AA", "J88_OPLS_AA_only")
        validate_forms(m["functional_form_convention"])
    captured = j87.capture_opls(
        manifest,
        manifest_sha256,
        audit,
        audit_sha256,
        capture_timestamp=capture_timestamp,
    )
    # Capture verified every source byte before this additional metadata parsing.
    blockers = []
    for file in captured["files"]:
        if file["raw_text"] is None:
            continue
        source = loads(file["raw_text"])
        if not isinstance(source, dict) or not isinstance(source.get("records"), list):
            continue  # The inherited strict adapter reports the schema defect.
        for index, row in enumerate(source["records"]):
            if isinstance(row, dict) and (
                not isinstance(row.get("location"), str) or not row["location"].strip()
            ):
                blockers.append(
                    "source_location_blocked:"
                    + file["path"]
                    + "#/records/"
                    + str(index)
                )
    if blockers:
        raise InputBlockers(blockers)
    return captured


def qualify_metadata(request_capture, source_capture, normalize_fixture):
    """J87/J85/J86 qualification, projecting only audit metadata for persistence."""
    request_fields(loads(request_capture["request_raw"]))
    result = j87.qualify_opls(request_capture, source_capture, normalize_fixture)
    q = result["base_qualification"]
    # Never serialize q's raw evidence_graph_json or source parameter values.
    safe_q = {
        k: v
        for k, v in q.items()
        if k not in ("evidence_graph_json", "graph_reports_json")
    }
    return {
        k: v
        for k, v in result.items()
        if k not in ("base_qualification", "source_evidence")
    } | {
        "base_qualification": safe_q,
        "source_evidence": {
            "source_authority": result["source_evidence"]["source_authority"],
            "adapter_status": result["source_evidence"]["adapter_status"],
            "source_hashes": result["source_hashes"],
            "source_authorized": q["source_authorized"],
            "qualification_checks": result["source_evidence"]["qualification_checks"],
            "records": result["source_locations"],
            "raw_source_records_persisted": False,
        },
    }


def external_path(path):
    """Require explicit regular-file paths (including no symlink ancestors)."""
    if path is not None:
        path = Path(path).absolute()
        require(
            not any(p.is_symlink() for p in (path, *path.parents)),
            "external_path_symlink_forbidden",
        )
    return path
