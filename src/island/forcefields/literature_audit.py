"""Local, data-only literature audits; no inferred chemical or MD authority.

J43 supplies pinned reads and independent decimal reconstruction. This layer
adds recursive inventories and explicit convention/conversion evidence. Native
sources and numerical adapters are deliberately outside this contract.
"""

import os
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.dynamics._checkpoint_data import strict_load

from .external_pdms import LIMITS, digest, require
from .oplsaa.parameters import inventory as opls_inventory
from .pcff.automatic import chemical_graph
from .pcff.class2 import FAMILIES, interaction_requests
from .pdms_literature import INPUT_ERRORS, pinned_bytes, reconstruct_rows

SCHEMA = "island_local_literature_audit_v1"
BUNDLE_SCHEMA = "island_external_literature_bundle_candidate_v1"
ROW_SCHEMA = "island_external_literature_evidence_row_v1"
GROUPS = (
    "pcff_compass_class2",
    "opls_aa",
    "pdms_siloxane",
    "other_extension",
    "explanatory_without_executable_rows",
    "incompatible_incomplete_or_unlicensed",
)
KINDS = (
    "atom_typing",
    "partial_charge",
    "bond",
    "angle",
    "torsion",
    "improper",
    "class2_cross_term",
    "nonbonded",
    "special_14",
    "exclusion",
    "mixing_rule",
)
STAGES = (
    "source_inventory_completed",
    "citation_audit_completed",
    "typing_complete",
    "charge_complete",
    "parameter_rows_complete",
    "executable_backend_complete",
    "independent_numerical_verification_complete",
    "scientific_validation_complete",
)


@dataclass(frozen=True)
class LiteratureAudit:
    """Immutable recursive inventory; missing bibliography is recorded as unknown."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(p["schema"] == SCHEMA, "wrong_audit_schema")
        require(
            p["execution_profiles"] == []
            and p["acceptance_gate"] == 1
            and p["limitations"] == LIMITS,
            "audit_cannot_grant_execution",
        )
        require(
            set(p["stages"]) == set(STAGES)
            and all(p["stages"][k] is False for k in STAGES[1:]),
            "inventory_is_not_qualification",
        )
        require(
            len({f["filename"] for f in p["files"]}) == len(p["files"]),
            "duplicate_inventory_path",
        )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return sha256(self.json_text.encode()).hexdigest()

    def validate_integrity(self, root):
        require(
            self == inventory_literature(root), "stale_or_mutated_literature_inventory"
        )


def file_record(path, name):
    """Stream every byte; never execute archives, PDFs or supplied source code."""
    metadata = {
        "title": None,
        "doi": None,
        "authors": None,
        "year": None,
        "journal": None,
    }
    base = {
        "filename": name,
        "bibliography": metadata,
        "license": "unknown",
        "access_status": "local_user_supplied",
        "category": GROUPS[-1],
        "numerical_parameter_rows": "not_independently_audited",
        "diagnostics": ["citation_license_and_row_evidence_required"],
    }
    if path.is_symlink():
        raw = os.readlink(path).encode()
        return {
            **base,
            "size": len(raw),
            "sha256": sha256(raw).hexdigest(),
            "hash_basis": "symlink_text",
            "file_type": "symlink_not_followed",
            "symlink_target": os.readlink(path),
        }
    h, size, head = sha256(), 0, b""
    with path.open("rb") as stream:
        before = os.fstat(stream.fileno())
        while chunk := stream.read(1024 * 1024):
            if not head:
                head = chunk[:4096]
            h.update(chunk)
            size += len(chunk)
        after = os.fstat(stream.fileno())
    require(
        (before.st_size, before.st_mtime_ns, before.st_ino)
        == (after.st_size, after.st_mtime_ns, after.st_ino)
        and size == before.st_size,
        "source_changed_during_inventory",
    )
    kind = "binary"
    for magic, value in (
        (b"%PDF", "pdf"),
        (b"PK\x03\x04", "zip_archive"),
        (b"\x1f\x8b", "gzip_archive"),
        (b"\xfd7zXZ", "xz_archive"),
    ):
        if head.startswith(magic):
            kind = value
            break
    else:
        try:
            head.decode("utf-8")
            kind = "utf8_text_candidate" if b"\0" not in head else "binary"
        except UnicodeError:
            pass
    return {
        **base,
        "size": size,
        "sha256": h.hexdigest(),
        "hash_basis": "file_bytes",
        "file_type": kind,
        "symlink_target": None,
    }


def inventory_literature(root):
    """Inventory an explicitly supplied directory, including hidden files.

    Neither filename hints nor readable numbers establish bibliography, licensing
    or parameter completeness. Unknown values remain unknown until independently
    evidenced. Symlinks are listed without following them outside this scope.
    """
    root = Path(root)
    files, errors = [], []

    def walk(folder):
        try:
            entries = sorted(folder.iterdir(), key=lambda p: p.name)
        except OSError as error:
            errors.append(
                {
                    "stage": "inventory_unreadable",
                    "path": str(folder.relative_to(root)),
                    "context": str(error),
                }
            )
            return
        for path in entries:
            name = path.relative_to(root).as_posix()
            try:
                if path.is_symlink() or path.is_file():
                    files.append(file_record(path, name))
                elif path.is_dir():
                    walk(path)
                else:
                    errors.append({"stage": "unsupported_file_kind", "path": name})
            except INPUT_ERRORS as error:
                errors.append(
                    {
                        "stage": "inventory_file_failed",
                        "path": name,
                        "context": str(error),
                    }
                )

    if root.is_dir() and not root.is_symlink():
        walk(root)
    else:
        errors.append({"stage": "source_directory_unavailable"})
    complete = not errors
    status = "asset_unavailable" if not files else "parameter_incomplete"
    return LiteratureAudit(
        pack(
            {
                "schema": SCHEMA,
                "files": files,
                "diagnostics": errors,
                "status": status,
                "stages": {
                    k: complete if i == 0 else False for i, k in enumerate(STAGES)
                },
                "groups": {
                    g: [f["filename"] for f in files if f["category"] == g]
                    for g in GROUPS
                },
                "execution_profiles": [],
                "acceptance_gate": 1,
                "limitations": LIMITS,
            }
        )
    )


def validate_evidence_row(record, *, root, audit):
    """Reconstruct a proposed row without promoting it to an executable profile.

    Only identity conversion is implemented. Other conversions require a new,
    independently verified version, never an inferred unit factor.
    """
    from island.workflows.storage import child

    audit.validate_integrity(root)
    require(
        set(record)
        == {
            "schema",
            "kind",
            "model_family",
            "functional_form",
            "source",
            "row",
            "conventions",
            "conversion",
            "evidence_checksum",
        },
        "row_fields_invalid",
    )
    require(record["schema"] == ROW_SCHEMA, "wrong_row_schema")
    require(
        record["evidence_checksum"]
        == digest({k: v for k, v in record.items() if k != "evidence_checksum"}),
        "evidence_checksum_changed",
    )
    require(record["kind"] in KINDS, "unsupported_row_kind")
    require(
        record["model_family"]
        in (
            "external:pcff",
            "external:compass",
            "external:opls_aa",
            "external:pdms_class2",
        ),
        "external_family_required",
    )
    row, source = record["row"], record["source"]
    allowed = {
        "atom_typing": {"atom_type"},
        "partial_charge": {"charge", "bond_increment"},
        "bond": {"quartic_bond", "opls_harmonic_bond"},
        "angle": {"quartic_angle", "opls_harmonic_angle"},
        "torsion": {"torsion_3", "opls_fourier", "opls_rb"},
        "improper": {"wilson_out_of_plane", "opls_improper"},
        "nonbonded": {"nonbond(9-6)", "opls_lj12_6"},
        "class2_cross_term": set(FAMILIES)
        - {
            "quartic_bond",
            "quartic_angle",
            "torsion_3",
            "wilson_out_of_plane",
            "nonbond(9-6)",
        },
        "special_14": {"special_14"},
        "exclusion": {"exclusions"},
        "mixing_rule": {"nonbonded_policy"},
    }
    require(row["family"] in allowed[record["kind"]], "kind_family_mismatch")
    require(
        type(row["types"]) is list
        and bool(row["types"])
        and all(type(t) is str and t.strip() and "*" not in t for t in row["types"]),
        "typing_ambiguous",
    )
    require(
        record["functional_form"] in ("class2", "opls_aa")
        and (record["functional_form"] == "opls_aa")
        == (record["model_family"] == "external:opls_aa"),
        "source_incompatible",
    )
    require(
        all(
            type(source[k]) is str and source[k].strip()
            for k in ("file", "sha256", "doi", "citation", "location")
        ),
        "row_citation_missing",
    )
    known = {f["filename"]: f for f in audit.payload["files"]}
    require(
        source["file"] in known
        and source["sha256"] == known[source["file"]]["sha256"]
        and known[source["file"]]["hash_basis"] == "file_bytes",
        "source_hash_mismatch",
    )
    conventions = record["conventions"]
    require(
        set(conventions)
        == {
            "mixing_rule",
            "one_four_scaling",
            "charge_method",
            "torsion_convention",
            "state_point",
        }
        and all(type(v) is str and v.strip() for v in conventions.values()),
        "conventions_incomplete",
    )
    conversion = record["conversion"]
    require(
        conversion
        == {
            "policy": "identity_v1",
            "original_values": row["values"],
            "converted_values": row["values"],
            "original_units": row["units"],
            "converted_units": row["units"],
            "code": "island.identity_decimal_v1",
            "version": 1,
        },
        "unsupported_or_mutated_conversion",
    )
    require(
        all(
            type(u) is str
            and u.strip()
            and "unresolved" not in u
            and not any(c in u for c in ("E", "L", "A"))
            for u in row["units"]
        ),
        "units_unresolved",
    )
    concrete_units = {
        "atom_type": ["dalton"],
        "charge": ["elementary_charge"],
        "bond_increment": ["elementary_charge"],
        "opls_harmonic_bond": ["nm", "kJ/(mol*nm^2)"],
        "opls_harmonic_angle": ["radian", "kJ/(mol*radian^2)"],
        "opls_fourier": ["kJ/mol"] * 4,
        "opls_rb": ["kJ/mol"] * 6,
        "opls_lj12_6": ["nm", "kJ/mol"],
    }
    if row["family"] in concrete_units:
        require(
            row["units"] == concrete_units[row["family"]], "units_or_shape_incompatible"
        )
    if record["functional_form"] == "class2":
        require(
            row["family"] in FAMILIES
            or record["kind"]
            in (
                "atom_typing",
                "partial_charge",
                "mixing_rule",
                "special_14",
                "exclusion",
            ),
            "unsupported_cross_term",
        )
    else:
        require(record["kind"] != "class2_cross_term", "unsupported_cross_term")
        require(
            conventions["mixing_rule"] == "geometric_sigma_and_epsilon"
            and conventions["one_four_scaling"] == "lj=0.5;coulomb=0.5"
            and conventions["torsion_convention"]
            in ("opls_fourier_v1", "source_rb_v1"),
            "opls_convention_unverified",
        )
    assets = {}
    for citation in row["citations"]:
        name = citation["file"]
        require(
            name in known and known[name]["hash_basis"] == "file_bytes",
            "uncatalogued_citation",
        )
        assets[name] = pinned_bytes(child(root, name), known[name]["sha256"])
    require(source["file"] in assets, "source_file_not_cited")
    labels = row["types"]
    primary_labels_verified = False
    for citation in row["citations"]:
        if citation["file"] != source["file"]:
            continue
        tokens = citation["raw"].split("!", 1)[0].split()
        for start in range(len(tokens) - len(labels) + 1):
            columns = set(range(start, start + len(labels)))
            if tokens[
                start : start + len(labels)
            ] == labels and not columns.intersection(citation["columns"]):
                primary_labels_verified = True
    require(primary_labels_verified, "atom_labels_not_reconstructed")
    extracted, errors = reconstruct_rows({"rows": {"candidate": row}}, assets)
    require(not errors, "source_row_not_reconstructed:" + str(errors))
    # These are verified tokens, not a licence, chemical typing or MD authority.
    return {
        "extracted": extracted["candidate"],
        "evidence_checksum": record["evidence_checksum"],
        "execution_authorized": False,
        "citation_metadata_status": "supplied_not_independently_qualified",
    }


def fixture_requirements(system, form):
    """Exact graph-bound requests, including unresolved source applicability."""
    system.validate()
    graph = chemical_graph(system)
    components = sorted(
        [sorted(c) for c in system.topology.connected_components()], key=min
    )
    graph["components"] = components
    graph["component_membership"] = {
        str(i): n for n, c in enumerate(components) for i in c
    }
    graph["molecule_membership"] = {
        str(i): str(system.topology.sites[i].metadata.get("chain_id", n))
        for n, c in enumerate(components)
        for i in c
    }
    ids = [s["id"] for s in graph["sites"]]
    requests = {}
    for i in ids:
        for family in ("atom_typing", "mass", "partial_charge", "nonbonded"):
            requests[f"{family}:{i}"] = {"sites": [i], "dependencies": []}
    if form == "class2":
        rows, _ = interaction_requests(graph, expanded=True)
        for family, sites, deps, reason in rows:
            if reason:
                continue
            key = family + ":" + ",".join(map(str, sites))
            requests[key] = {
                "sites": list(sites),
                "dependencies": [
                    f + ":" + ",".join(map(str, min(s, s[::-1]))) for f, s in deps
                ],
            }
        requests["policy:cross_term_applicability"] = {"sites": ids, "dependencies": []}
    else:
        require(form == "opls_aa", "unsupported_functional_form")
        inventory, pairs = opls_inventory(system)
        for family, paths in inventory.items():
            for sites in paths:
                requests[family + ":" + ",".join(map(str, sites))] = {
                    "sites": list(sites),
                    "dependencies": ["policy:torsion_convention"]
                    if family == "rb"
                    else [],
                }
        adjacent = {i: set() for i in ids}
        for b in graph["bonds"]:
            a, z = b["sites"]
            adjacent[a].add(z)
            adjacent[z].add(a)
        for i in ids:
            for outer in combinations(sorted(adjacent[i]), 3):
                requests[
                    "improper_applicability:" + ",".join(map(str, (i, *outer)))
                ] = {"sites": [i, *outer], "dependencies": []}
        for pair in pairs:
            requests["special_pair:" + ",".join(map(str, pair["sites"]))] = {
                "sites": list(pair["sites"]),
                "distance": pair["distance"],
                "dependencies": ["policy:special_pairs"],
            }
        requests["policy:torsion_convention"] = {"sites": ids, "dependencies": []}
    for p in ("mixing", "special_pairs", "charge_origin", "chemical_domain"):
        requests["policy:" + p] = {"sites": ids, "dependencies": []}
    return {
        "graph": graph,
        "graph_identity": digest(graph),
        "atom_order": ids,
        "requests": requests,
    }


@dataclass(frozen=True)
class LiteratureBundleAudit:
    """Checksummed blocked bundle report. Never a native or external MD model."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(
            p["schema"] == BUNDLE_SCHEMA
            and p["status"]
            in (
                "asset_unavailable",
                "parameter_incomplete",
                "typing_ambiguous",
                "source_incompatible",
            )
            and p["execution_profile"] is None
            and p["acceptance_gate"] == 1
            and p["limitations"] == LIMITS,
            "incomplete_bundle_cannot_authorize_md",
        )

    @property
    def payload(self):
        return unpack(self.json_text)

    def require_executable(self):
        raise LiteratureBlockedError(self)


class LiteratureBlockedError(ValueError):
    def __init__(self, result):
        self.result = result
        self.status = result.payload["status"]
        super().__init__(self.status + ": external literature profile not qualified")


def audit_parameter_bundle(
    system, *, root, audit, form, bundle=None, bundle_sha256=None
):
    """Validate explicit candidate rows and preserve every missing dependency.

    J44 has no qualified numerical adapter for new literature profiles. Even
    complete candidate rows cannot bypass the independent qualification gates.
    """
    required = fixture_requirements(system, form)
    errors, supplied, extracted = [], {}, {}
    status = (
        "asset_unavailable" if not audit.payload["files"] else "parameter_incomplete"
    )
    if bundle is not None:
        status = "parameter_incomplete"
        try:
            audit.validate_integrity(root)
            candidate = strict_load(pinned_bytes(bundle, bundle_sha256).decode())
            require(
                set(candidate)
                == {"schema", "graph_identity", "atom_order", "rows", "assignments"},
                "candidate_fields_invalid",
            )
            require(
                candidate["graph_identity"] == required["graph_identity"]
                and candidate["atom_order"] == required["atom_order"],
                "stale_graph_or_atom_order",
            )
            require(candidate["schema"] == BUNDLE_SCHEMA, "wrong_bundle_schema")
            families, selections, papers = set(), {}, set()
            for rid, record in candidate["rows"].items():
                try:
                    require(
                        record["functional_form"] == form,
                        "fixture_functional_form_mismatch",
                    )
                    extracted[rid] = validate_evidence_row(
                        record, root=root, audit=audit
                    )
                    families.add(record["model_family"])
                    papers.add(record["source"]["doi"])
                    key = (record["row"]["family"], tuple(record["row"]["types"]))
                    value = (record["row"]["values"], record["row"]["units"])
                    require(
                        key not in selections or selections[key] == value,
                        "conflicting_duplicate_row",
                    )
                    selections[key] = value
                except INPUT_ERRORS as error:
                    extracted.pop(rid, None)
                    errors.append(
                        {"stage": "row_rejected", "row": rid, "context": str(error)}
                    )
            require(len(families) <= 1, "cross_family_borrowing_forbidden")
            require(len(papers) <= 1, "cross_paper_borrowing_not_authorized")
            assignments = candidate["assignments"]
            require(
                set(assignments) <= set(required["requests"]),
                "unexpected_interaction_assignment",
            )
            supplied = {k: rid for k, rid in assignments.items() if rid in extracted}
            require(len(supplied) == len(assignments), "unreconstructed_row_assignment")
        except INPUT_ERRORS as error:
            errors.append({"stage": "candidate_rejected", "context": str(error)})
    missing = {k: v for k, v in required["requests"].items() if k not in supplied}
    dependencies = [
        {"consumer": k, "dependency": d}
        for k, v in required["requests"].items()
        for d in v["dependencies"]
        if k in missing or d in missing
    ]
    if not missing and not errors:
        status = "source_incompatible"
        errors.append(
            {"stage": "chemical_typing_and_numerical_qualification_unavailable"}
        )
    return LiteratureBundleAudit(
        pack(
            {
                "schema": BUNDLE_SCHEMA,
                "status": status,
                "requirements": required,
                "source_inventory_identity": audit.identity,
                "functional_form": form,
                "missing_rows": missing,
                "unresolved_dependencies": dependencies,
                "extracted_candidate_rows": extracted,
                "diagnostics": errors,
                "stages": audit.payload["stages"],
                "execution_profile": None,
                "acceptance_gate": 1,
                "limitations": LIMITS,
            }
        )
    )


REVIEW_SCHEMA = "island.j44-reviewed-literature.v1"


@dataclass(frozen=True)
class LiteratureSourceReview:
    """Immutable reviewed metadata and candidate tokens, without MD authority.

    Raw-file reconstruction is a separate explicit offline command. A checksummed
    transcription alone cannot establish citation accuracy or chemical validity.
    """

    json_text: str

    def __post_init__(self):
        import json
        from decimal import Decimal

        p = self.payload
        checksum = p.pop("checksum")
        encoded = json.dumps(
            p, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        )
        require(
            sha256(encoded.encode()).hexdigest() == checksum, "review_checksum_changed"
        )
        require(p["schema"] == REVIEW_SCHEMA, "review_schema_changed")
        require(
            p["accepted_parameter_rows"] == []
            and p["execution_profiles"] == []
            and p["acceptance_gate"] == 1
            and p["production_validated"] is False
            and p["simulation_readiness"] == "not_established",
            "review_cannot_authorize_execution",
        )
        files = {f["filename"]: f for f in p["files"]}
        require(len(files) == len(p["files"]), "duplicate_review_path")
        for f in files.values():
            require(f["category"] in GROUPS, "unknown_source_category")
            require(
                set(f["bibliography"]) == {"title", "doi", "authors", "year", "journal"}
                and f["diagnostics"]
                and f["license"]
                and f["access_status"],
                "review_metadata_missing",
            )
        for row in p["extracted_row_candidates"]:
            body = {k: v for k, v in row.items() if k != "evidence_checksum"}
            encoded = json.dumps(
                body, sort_keys=True, separators=(",", ":"), ensure_ascii=False
            )
            require(
                sha256(encoded.encode()).hexdigest() == row["evidence_checksum"],
                "candidate_evidence_changed",
            )
            f = files[row["source_file"]]
            require(
                row["source_sha256"] == f["sha256"]
                and row["doi"] == f["bibliography"]["doi"],
                "candidate_source_changed",
            )
            require(
                row["status"] == "parameter_incomplete"
                and row["converted_values"] is None,
                "unqualified_conversion_or_row",
            )
            tokens = row["raw_text"].split()
            require(
                row["atom_or_interaction_labels"]
                == tokens[1 : min(row["numeric_columns"])],
                "candidate_labels_changed",
            )
            require(
                [Decimal(tokens[i]) for i in row["numeric_columns"]]
                == [Decimal(v) for v in row["original_values"]],
                "candidate_decimal_mismatch",
            )

            token_row = {
                "family": "literature_token_candidate",
                "types": ["unqualified_transcription"],
                "values": row["original_values"],
                "units": ["unresolved; no parameter acceptance"]
                * len(row["original_values"]),
                "provenance": row["location"],
                "citations": [
                    {
                        "file": "excerpt",
                        "line": 1,
                        "raw": row["raw_text"],
                        "columns": row["numeric_columns"],
                    }
                ],
            }
            _, errors = reconstruct_rows(
                {"rows": {"token": token_row}}, {"excerpt": row["raw_text"].encode()}
            )
            require(not errors, "j43_token_reconstruction_failed")

    @property
    def payload(self):
        return strict_load(self.json_text)

    @property
    def identity(self):
        return sha256(self.json_text.encode()).hexdigest()

    def validate_inventory(self, audit):
        actual = {
            f["filename"]: (f["sha256"], f["size"]) for f in audit.payload["files"]
        }
        reviewed = {
            f["filename"]: (f["sha256"], f["size"]) for f in self.payload["files"]
        }
        require(actual == reviewed, "review_inventory_changed")
        return {
            "files": len(actual),
            "review_identity": self.identity,
            "candidate_rows": len(self.payload["extracted_row_candidates"]),
            "accepted_parameter_rows": 0,
            "execution_profiles": [],
        }
