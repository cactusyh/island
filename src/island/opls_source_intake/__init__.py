"""New-source intake, distinct from canonical assignment and scientific approval.

No engine, parameter converter or downloaded source is invoked by intake.
Untrusted package declarations cannot modify the immutable OPLS registry.
"""

import io
import json
import math
import re
import tarfile
import zipfile
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path, PurePosixPath

from island.exceptions import IslandError
from island.forcefields.opls_registry.api import require
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import canonical, digest, final_graph
from island.opls_conventions import signed
from island.opls_final_resolution import corpus_audit

from .schema import REQUIRED_CONVENTIONS, STATUSES, source_intake_schema

MANIFEST = "source-intake.json"
# No genuinely new source or independently qualified adapter has been supplied.
# In particular the J49 exact-erfc diagnostic is NOT an author qualification.
QUALIFIED_SOURCES = ()
_PACKAGE = ContextVar("island_opls_source_intake", default=None)


@dataclass(frozen=True)
class SourcePackage:
    root: object
    manifest_sha256: str


@contextmanager
def use_source_package(package):
    require(type(package) is SourcePackage, "source_incompatible", "package_required")
    token = _PACKAGE.set(package)
    try:
        yield
    finally:
        _PACKAGE.reset(token)


class IntakeProblem(ValueError):
    def __init__(self, status, code, **details):
        self.status = status
        self.diagnostic = {"code": code, **details}
        super().__init__(code)


def check(condition, code, *, conflict=False, **details):
    if not condition:
        raise IntakeProblem(
            "supplied_conflicting" if conflict else "supplied_incomplete",
            code,
            **details,
        )


def _json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            check(key not in result, "duplicate_JSON_key", conflict=True, key=key)
            result[key] = value
        return result

    return json.loads(
        raw,
        object_pairs_hook=pairs,
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)),
    )


def _shape(value, schema, path="manifest"):
    if "const" in schema:
        check(value == schema["const"], "schema_constant", field=path)
    if "enum" in schema:
        check(value in schema["enum"], "schema_enum", field=path)
    kind = schema.get("type")
    if kind:
        expected = {
            "object": dict,
            "array": list,
            "string": str,
            "integer": int,
            "number": (int, float),
        }[kind]
        check(
            type(value) in (expected if isinstance(expected, tuple) else (expected,)),
            "schema_type",
            field=path,
        )
    if kind == "object":
        check(
            set(value) == set(schema["properties"]),
            "schema_fields",
            field=path,
            missing=sorted(set(schema["properties"]) - set(value)),
            unexpected=sorted(set(value) - set(schema["properties"])),
        )
        for key, sub in schema["properties"].items():
            _shape(value[key], sub, path + "." + key)
    elif kind == "array":
        check(len(value) >= schema.get("minItems", 0), "schema_min_items", field=path)
        for n, item in enumerate(value):
            _shape(item, schema["items"], f"{path}[{n}]")
    elif kind == "string":
        check(
            bool(value.strip()) and len(value) >= schema.get("minLength", 0),
            "schema_empty_text",
            field=path,
        )
        if "pattern" in schema:
            check(
                re.fullmatch(schema["pattern"], value) is not None,
                "schema_pattern",
                field=path,
            )
    elif kind in ("number", "integer"):
        check(
            math.isfinite(value) and value >= schema.get("minimum", -math.inf),
            "schema_number",
            field=path,
        )


def _safe(name):
    p = PurePosixPath(name)
    check(
        bool(name)
        and not p.is_absolute()
        and ".." not in p.parts
        and "\\" not in name
        and str(p) == name,
        "unsafe_relative_path",
        conflict=True,
        path=name,
    )
    return name


def inventory(root):
    """Read every physical file; never follow links or execute source code."""
    root = Path(root)
    if not root.is_dir():
        raise IntakeProblem("supplied_unreadable", "package_directory_unreadable")
    check(not root.is_symlink(), "source_symlink", conflict=True)
    rows = []
    for path in sorted(root.rglob("*")):
        name = _safe(path.relative_to(root).as_posix())
        check(not path.is_symlink(), "source_symlink", conflict=True, path=name)
        if path.is_dir():
            continue
        check(path.is_file(), "nonregular_source_file", conflict=True, path=name)
        raw = path.read_bytes()
        rows.append({"path": name, "size": len(raw), "sha256": sha256(raw).hexdigest()})
    return rows


def archive_inventory(root, rows):
    """Hash archive members without extraction; reject unsafe/ambiguous archives."""
    result = []
    remaining = [256 * 1024 * 1024]

    def walk(raw, name, depth=0):
        check(depth <= 3, "archive_nesting_limit", path=name)
        entries = []
        if name.lower().endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                for info in z.infolist():
                    if info.is_dir():
                        continue
                    check(
                        (info.external_attr >> 16) & 0o170000 != 0o120000,
                        "archive_symlink",
                        conflict=True,
                        path=name,
                    )
                    check(
                        info.file_size <= remaining[0], "archive_size_limit", path=name
                    )
                    remaining[0] -= info.file_size
                    entries.append((info.filename, z.read(info)))
        elif name.lower().endswith((".tar", ".tar.gz", ".tgz", ".tar.xz")):
            with tarfile.open(fileobj=io.BytesIO(raw)) as t:
                for info in t:
                    if info.isdir():
                        continue
                    check(
                        info.isfile(),
                        "archive_nonregular_entry",
                        conflict=True,
                        path=name + "!" + info.name,
                    )
                    check(info.size <= remaining[0], "archive_size_limit", path=name)
                    remaining[0] -= info.size
                    entries.append((info.name, t.extractfile(info).read()))
        seen = set()
        for member, data in sorted(entries):
            _safe(member)
            check(
                member not in seen, "duplicate_archive_member", conflict=True, path=name
            )
            seen.add(member)
            full = name + "!" + member
            result.append(
                {"path": full, "size": len(data), "sha256": sha256(data).hexdigest()}
            )
            walk(data, full, depth + 1)

    for row in rows:
        if row["path"] != MANIFEST:
            walk((Path(root) / row["path"]).read_bytes(), row["path"])
    return sorted(result, key=lambda r: r["path"])


def _evidence(root, evidence, files, document):
    name = _safe(evidence["file"])
    check(
        name in files and evidence["sha256"] == files[name]["sha256"],
        "evidence_source_hash",
        conflict=True,
        path=name,
    )
    check(
        evidence["doi"] == document["doi"]
        and evidence["citation"] == document["citation"],
        "evidence_citation_mismatch",
        conflict=True,
        path=name,
    )
    # Binary PDFs require an inventoried text extraction. That text can be cited
    # here but cannot independently authorize a profile without adapter review.
    lines = (root / name).read_text(encoding="utf-8").splitlines()
    start = evidence["first_line"] - 1
    raw_lines = evidence["raw_text"].splitlines()
    check(
        lines[start : start + len(raw_lines)] == raw_lines,
        "raw_evidence_mismatch",
        conflict=True,
        path=name,
    )


def _record(root, row, files, document):
    check(
        row["units"].strip().lower()
        not in {"unknown", "unresolved", "default", "n/a", "unspecified"},
        "unresolved_units",
        key=row["key"],
    )
    payload = {k: v for k, v in row.items() if k != "evidence_checksum"}
    check(
        digest(payload) == row["evidence_checksum"],
        "row_evidence_checksum",
        conflict=True,
        key=row["key"],
    )
    for evidence in row["evidence"]:
        _evidence(root, evidence, files, document)
    recon = row["reconstruction"]
    for key in ("code", "result"):
        _evidence(root, recon[key], files, document)
    check(
        files[recon["code"]["file"]]["role"]
        in ("converter", "independent_reconstruction"),
        "reconstruction_code_required",
        key=row["key"],
    )
    check(
        recon["value"] == row["value"] and recon["units"] == row["units"],
        "independent_reconstruction_conflict",
        conflict=True,
        key=row["key"],
    )
    # Machine-readable reconstruction result must contain the exact value,
    # units and row ID. We do not execute untrusted reconstruction code.
    result = _json(recon["result"]["raw_text"])
    check(
        result == {"key": row["key"], "value": row["value"], "units": row["units"]},
        "reconstruction_result_mismatch",
        conflict=True,
        key=row["key"],
    )
    if row["resolution"] == "converter_resolved":
        check(
            any(files[e["file"]]["role"] == "converter" for e in row["evidence"]),
            "converter_rule_missing",
            key=row["key"],
        )


def _unique(rows, key):
    names = [row[key] for row in rows]
    check(
        len(names) == len(set(names)),
        "duplicate_conflicting_record",
        conflict=True,
        field=key,
    )
    return {row[key]: row for row in rows}


def _domains(root, document, files):
    from island.charge_references.records import unpack
    from island.forcefields.selection.inventory import requirements
    from island.pdms import restore_pdms_system

    _unique(document["domains"], "graph_identity")
    _unique(document["domains"], "name")
    for domain in document["domains"]:
        name = _safe(domain["system_file"])
        check(
            name in files and files[name]["role"] == "system", "domain_system_missing"
        )
        system = restore_pdms_system(unpack((root / name).read_text()))
        graph = final_graph(system)
        order = graph["atom_order"]
        check(
            digest(graph) == domain["graph_identity"], "graph_mutation", conflict=True
        )
        check(
            len(order) == len(domain["atom_types"]) == len(domain["charges"]),
            "typing_or_charge_incomplete",
        )
        check(
            digest(
                {
                    "graph_identity": domain["graph_identity"],
                    "atom_types": list(zip(order, domain["atom_types"])),
                }
            )
            == domain["typed_graph_identity"],
            "typed_graph_mutation",
            conflict=True,
        )
        rows = _unique(domain["rows"], "key")
        needed = {r["id"] for r in requirements(system, "opls")}
        check(
            set(rows) == needed,
            "interaction_rows_incomplete",
            missing=sorted(needed - set(rows)),
            unexpected=sorted(set(rows) - needed),
        )
        for row in domain["rows"]:
            _record(root, row, files, document)
            check(
                row["labels"] or row["key"].endswith(":"),
                "interaction_labels_missing",
                key=row["key"],
            )
        check(
            domain["charge_method"]
            == next(
                r["value"]
                for r in document["conventions"]
                if r["key"] == "charge_method"
            ),
            "charge_method_conflict",
            conflict=True,
        )
        for i, atom_type, charge in zip(order, domain["atom_types"], domain["charges"]):
            check(
                rows[f"atom_type:{i}"]["value"] == atom_type,
                "atom_type_conflict",
                conflict=True,
                site=i,
            )
            check(
                float(rows[f"charge:{i}"]["value"]) == charge
                and rows[f"charge:{i}"]["units"] == "elementary_charge",
                "charge_conflict",
                conflict=True,
                site=i,
            )
            check(
                float(rows[f"mass:{i}"]["value"]) == system.topology.sites[i].mass,
                "mass_conflict",
                conflict=True,
                site=i,
            )


def blocker_reopen_report(result):
    """Preserve each J52 result verbatim; an intake declaration cannot clear it."""
    audit = corpus_audit()
    return signed(
        {
            "schema": "island.opls-blocker-reopen.v1",
            "J52_source_audit_identity": audit["identity"],
            "intake_identity": result["identity"],
            "source_status": result["status"],
            "rows": [
                {
                    "setting": row["setting"],
                    "permanent_blocker": row,
                    "required_new_artifact": row["exact_missing_artifact"],
                    "reopen_status": "new_authority_and_independent_qualification_required",
                    "candidate_evidence": next(
                        (
                            r
                            for r in result["candidate_conventions"]
                            if r["key"] == row["setting"]
                        ),
                        None,
                    ),
                    "reopened": False,
                }
                for row in audit["final_resolutions"]
            ],
            "MD_authorized": False,
            "new_profiles": [],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def inspect_source(package=None):
    """Deterministic intake, with external manifest pin and no profile promotion."""
    status, diagnostics, rows, members, document = "not_supplied", [], [], [], None
    manifest_pin = None
    if package is not None:
        try:
            check(type(package) is SourcePackage, "package_required")
            root = Path(package.root)
            rows = inventory(root)
            raw = (root / MANIFEST).read_bytes()
            manifest_pin = sha256(raw).hexdigest()
            check(
                manifest_pin == package.manifest_sha256,
                "manifest_hash_mismatch",
                conflict=True,
            )
            document = _json(raw)
            _shape(document, source_intake_schema())
            check(
                document["proposed_profile"].startswith("island_external_opls_")
                and re.search(r"_v[1-9][0-9]*$", document["proposed_profile"])
                is not None
                and document["proposed_profile"] != PROFILE,
                "new_external_version_required",
                conflict=True,
            )
            files = _unique(document["files"], "path")
            for name in files:
                _safe(name)
            actual = [r for r in rows if r["path"] != MANIFEST]
            check(
                actual
                == sorted(
                    (
                        {k: r[k] for k in ("path", "size", "sha256")}
                        for r in files.values()
                    ),
                    key=lambda r: r["path"],
                ),
                "recursive_inventory_mismatch",
                conflict=True,
            )
            members = archive_inventory(root, rows)
            auth = document["authorization"]
            check(auth["status"] == "permitted", "source_permission_unavailable")
            _evidence(root, auth["evidence"], files, document)
            check(
                files[auth["evidence"]["file"]]["role"] == "permission",
                "permission_artifact_required",
            )
            coverage = _unique(auth["covered_files"], "file")
            required_coverage = {
                name: f["sha256"]
                for name, f in files.items()
                if f["role"] != "permission"
            }
            check(
                {name: f["sha256"] for name, f in coverage.items()}
                == required_coverage,
                "permission_scope_incomplete",
                conflict=True,
            )
            _evidence(root, document["version_evidence"], files, document)
            check(
                document["source_version"] in document["version_evidence"]["raw_text"],
                "version_evidence_mismatch",
                conflict=True,
            )
            conventions = _unique(document["conventions"], "key")
            check(
                set(conventions) == set(REQUIRED_CONVENTIONS),
                "execution_conventions_incomplete",
                missing=sorted(set(REQUIRED_CONVENTIONS) - set(conventions)),
                unexpected=sorted(set(conventions) - set(REQUIRED_CONVENTIONS)),
            )
            for row in document["conventions"]:
                _record(root, row, files, document)
            _domains(root, document, files)
            status = "source_complete_not_executed"
        except IntakeProblem as error:
            status, diagnostics = error.status, [error.diagnostic]
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
            zipfile.BadZipFile,
            tarfile.TarError,
        ) as error:
            status = "supplied_unreadable"
            diagnostics = [
                {
                    "code": "unreadable_source_package",
                    "error_type": type(error).__name__,
                }
            ]
        except (
            ValueError,
            KeyError,
            TypeError,
            AttributeError,
            IndexError,
            OverflowError,
            IslandError,
        ) as error:
            status = "supplied_incomplete"
            diagnostics = [
                {"code": "invalid_source_structure", "error_type": type(error).__name__}
            ]
    check(status in STATUSES, "invalid_source_status")
    return signed(
        {
            "schema": "island.opls-source-intake-result.v1",
            "family": "opls",
            "status": status,
            "manifest_sha256": manifest_pin,
            "inventory": rows,
            "archive_members": members,
            "diagnostics": diagnostics,
            "candidate_conventions": sorted(
                document.get("conventions", []), key=lambda r: str(r.get("key"))
            )
            if isinstance(document, dict) and status == "source_complete_not_executed"
            else [],
            "candidate_profile": document["proposed_profile"]
            if status == "source_complete_not_executed"
            else None,
            "candidate_domains": document["domains"]
            if status == "source_complete_not_executed"
            else [],
            "qualification_status": "no_new_independently_qualified_source_adapter",
            "automatically_selectable": False,
            "new_profiles": [],
            "execution_evidence": None,
            "independent_verification_evidence": None,
            "MD_authorized": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


@dataclass(frozen=True)
class IntakePreparation:
    json_text: str
    prior: object = field(repr=False)
    package: SourcePackage | None = field(default=None, repr=False)

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    @property
    def assignment(self):
        return self.prior.assignment

    def validate_integrity(self, system):
        value = self.prior.validate_integrity(system)
        require(
            self.json_text == canonical(_preparation(self.prior, self.package)),
            "source_incompatible",
            "intake_source_or_report_mutation",
        )
        return value

    def require_execution(self, system, action="static", *, research_opt_in=False):
        self.validate_integrity(system)
        require(
            action == "static",
            "parameter_incomplete",
            "new_authoritative_source_required",
            action=action,
            blockers=self.payload["blocker_reopen"]["rows"],
            source_status=self.payload["source_intake"]["status"],
        )
        require(
            self.package is None,
            "source_incompatible",
            "new_source_has_no_qualified_execution_adapter",
        )
        return self.prior.require_execution(
            system, action, research_opt_in=research_opt_in
        )


def _preparation(prior, package):
    result = inspect_source(package)
    assignment = prior.assignment.payload
    return signed(
        {
            "schema": "island.opls-source-intake-preparation.v1",
            **{
                k: v
                for k, v in prior.payload.items()
                if k not in ("schema", "identity")
            },
            "J52_resolution_identity": prior.identity,
            "assignment_binding": {
                "canonical_assignment_identity": prior.assignment.identity,
                "graph_identity": assignment["graph_identity"],
                "typed_graph_identity": digest(assignment["typed_graph"]),
                "charge_identity": digest(assignment["graph_charges"]),
                "source_identity": assignment["source_identity"],
            },
            "source_intake": result,
            "blocker_reopen": blocker_reopen_report(result),
            "new_source_selected": False,
            "selection_reason": "no_new_independently_qualified_authoritative_package",
        }
    )


def prepare_polymer_forcefield(system, forcefield="opls", **kwargs):
    """Family-only selection; new source installation is an administrator step."""
    from island.opls_final_resolution import prepare_polymer_forcefield as prepare

    package = _PACKAGE.get()
    prior = prepare(system, forcefield=forcefield, **kwargs)
    return IntakePreparation(canonical(_preparation(prior, package)), prior, package)


__all__ = [
    "STATUSES",
    "IntakePreparation",
    "SourcePackage",
    "blocker_reopen_report",
    "inspect_source",
    "inventory",
    "prepare_polymer_forcefield",
    "use_source_package",
]
