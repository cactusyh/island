"""Read-only PDF/text/archive parsing. Source code and archives are never run."""

import hashlib
import io
import re
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

from island.pcff_pdms.extension import require

MAX_BYTES = 64 * 1024 * 1024
MAX_MEMBERS = 4096
MAX_EXPANDED = 128 * 1024 * 1024
FAMILIES = {
    "quartic_bond": "bond",
    "quartic_angle": "angle",
    "torsion_3": "torsion",
    "wilson_out_of_plane": "improper",
    "bond-bond": "bond-bond",
    "bond-bond_1_3": "bond-bond_1_3",
    "bond-angle": "bond-angle",
    "angle-angle": "angle-angle",
    "end_bond-torsion_3": "end-bond-torsion",
    "middle_bond-torsion_3": "middle-bond-torsion",
    "angle-torsion_3": "angle-torsion",
    "angle-angle-torsion_1": "angle-angle-torsion",
    "torsion-torsion": "torsion-torsion",
    "torsion-torsion_1": "torsion-torsion",
    "nonbond(9-6)": "nonbonded",
    "nonbond(12-6)": "nonbonded",
    "bond_increments": "charge",
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def local_bytes(root, path):
    require(
        isinstance(path, str) and path and "://" not in path,
        "explicit_local_path_required",
    )
    p = Path(path)
    require(".." not in p.parts, "path_traversal")
    p = p if p.is_absolute() else root / p
    require(not any(q.is_symlink() for q in (p, *p.parents)), "source_symlink_rejected")
    require(p.is_file(), "source_missing")
    require(p.stat().st_size <= MAX_BYTES, "source_size_limit")
    with p.open("rb") as stream:
        raw = stream.read(MAX_BYTES + 1)
    require(len(raw) <= MAX_BYTES, "source_size_limit")
    return raw


def safe_name(name):
    p = PurePosixPath(name)
    require(
        name and not p.is_absolute() and ".." not in p.parts and "\\" not in name,
        "unsafe_archive_member",
    )
    return name.rstrip("/")


def archive(raw):
    """Bounded member reads in memory, with no filesystem extraction."""
    files, seen, expanded = {}, set(), 0
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            require(len(z.infolist()) <= MAX_MEMBERS, "archive_member_limit")
            for info in z.infolist():
                name = safe_name(info.filename)
                require(name not in seen, "duplicate_archive_member")
                seen.add(name)
                require(
                    (info.external_attr >> 16) & 0o170000 != 0o120000,
                    "archive_link_rejected",
                )
                require(not info.flag_bits & 1, "encrypted_archive_rejected")
                if info.is_dir():
                    continue
                expanded += info.file_size
                require(
                    info.file_size <= MAX_BYTES and expanded <= MAX_EXPANDED,
                    "archive_expansion_limit",
                )
                files[name] = z.read(info)
    else:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as t:
            for i, info in enumerate(t):
                require(i < MAX_MEMBERS, "archive_member_limit")
                name = safe_name(info.name)
                require(name not in seen, "duplicate_archive_member")
                seen.add(name)
                require(
                    info.isfile() or info.isdir(), "archive_link_or_special_rejected"
                )
                if info.isdir():
                    continue
                expanded += info.size
                require(
                    info.size <= MAX_BYTES and expanded <= MAX_EXPANDED,
                    "archive_expansion_limit",
                )
                files[name] = t.extractfile(info).read(MAX_BYTES + 1)
    return files


def pdf(raw):
    require(raw.startswith(b"%PDF-"), "PDF_signature_mismatch")
    try:
        import pypdf
    except ImportError as e:
        raise ValueError("PDF_parser_unavailable: install pypdf==6.1.1") from e
    require(pypdf.__version__ == "6.1.1", "PDF_parser_version_mismatch")
    try:
        reader = pypdf.PdfReader(io.BytesIO(raw), strict=True)
        require(not reader.is_encrypted, "encrypted_PDF_rejected")
        require(len(reader.pages) <= 10000, "PDF_page_limit")
        pages = [p.extract_text() or "" for p in reader.pages]
    except pypdf.errors.PyPdfError as e:
        raise ValueError("PDF_unreadable:" + type(e).__name__) from e
    return pages, {"tool": "pypdf", "version": pypdf.__version__}


def frc(raw):
    text = raw.decode("utf-8", errors="strict")
    section, records = None, []
    for i, line in enumerate(text.splitlines(keepends=True), 1):
        if line.lstrip().startswith("#"):
            section = (
                line.lstrip()[1:].strip().split()[0]
                if line.lstrip()[1:].strip()
                else None
            )
        elif section in FAMILIES and line.strip() and not line.lstrip().startswith("!"):
            records.append(
                {
                    "line_start": i,
                    "line_end": i,
                    "section": section,
                    "family": FAMILIES[section],
                    "raw_record": line,
                }
            )
    return text, records


def extract(raw, artifact, candidates):
    """Extract only cited snippets into output, never full source documents."""
    kind = artifact["source_type"]
    members = archive(raw) if kind == "authorized_archive" else {None: raw}
    inventory = (
        [
            {"path": n, "sha256": sha(b), "file_size": len(b)}
            for n, b in sorted(members.items(), key=lambda x: str(x[0]))
        ]
        if kind == "authorized_archive"
        else []
    )
    facts, failures, documents = {}, {}, {}
    if not candidates and kind != "authorized_archive":
        if kind == "primary_pdf":
            pages, parser = pdf(raw)
            documents["None"] = {
                "parser": parser,
                "page_count": len(pages),
                "page_text_sha256": [sha(t.encode()) for t in pages],
            }
        else:
            text, records = frc(raw)
            documents["None"] = {
                "parser": {"tool": "island.j81-text", "version": "1"},
                "text_sha256": sha(text.encode()),
                "recognized_FRC_records": len(records),
            }
    for candidate in candidates:
        cid, loc = candidate["id"], candidate["source_location"]
        try:
            require(isinstance(loc, dict), "source_location_required")
            member = loc.get("archive_member") if kind == "authorized_archive" else None
            require(member in members, "archive_member_missing")
            data = members[member]
            require(
                kind != "authorized_archive" or loc.get("member_sha256") == sha(data),
                "archive_member_hash_mismatch",
            )
            is_pdf = kind == "primary_pdf" or (
                member is not None and member.lower().endswith(".pdf")
            )
            raw_text = candidate["raw_text"]
            require(
                bool(raw_text) and len(raw_text) <= 8192, "short_raw_record_required"
            )
            if is_pdf:
                pages, parser = pdf(data)
                require(
                    artifact["extraction_tool"] == parser["tool"]
                    and artifact["extraction_version"] == parser["version"],
                    "extraction_tool_version_mismatch",
                )
                page = loc.get("page")
                require(
                    type(page) is int and 1 <= page <= len(pages),
                    "PDF_page_out_of_range",
                )
                text = pages[page - 1]
                require(
                    loc.get("ocr_status") == "not_used"
                    and loc.get("ocr_confidence") is None,
                    "OCR_unresolved_no_silent_OCR",
                )
                require(
                    {"table", "row", "column", "figure", "appendix"} <= loc.keys(),
                    "PDF_cell_metadata_required",
                )
                require(
                    bool(loc["table"] or loc["figure"] or loc["appendix"]),
                    "PDF_exact_locator_required",
                )
                if loc["table"]:
                    require(
                        bool(loc["row"]) and bool(loc["column"]),
                        "PDF_table_cell_required",
                    )
                require(
                    loc.get("extracted_text_sha256") == sha(text.encode()),
                    "PDF_text_checksum_mismatch",
                )
                require(text.count(raw_text) == 1, "PDF_snippet_absent_or_ambiguous")
                documents[str(member)] = {
                    "parser": parser,
                    "page_count": len(pages),
                    "page_text_sha256": [sha(t.encode()) for t in pages],
                }
            else:
                text, records = frc(data)
                parser = {"tool": "island.j81-text", "version": "1"}
                require(
                    artifact["extraction_tool"] == parser["tool"]
                    and artifact["extraction_version"] == parser["version"],
                    "extraction_tool_version_mismatch",
                )
                start, end = loc.get("line_start"), loc.get("line_end")
                lines = text.splitlines(keepends=True)
                require(
                    type(start) is int
                    and type(end) is int
                    and 1 <= start <= end <= len(lines),
                    "line_range_invalid",
                )
                require(end - start < 30, "snippet_length_limit")
                require(
                    raw_text == "".join(lines[start - 1 : end]),
                    "source_line_or_raw_record_mutation",
                )
                require(
                    bool(loc.get("section")) and bool(loc.get("parameter_record_id")),
                    "section_and_record_id_required",
                )
                if kind in ("original_pcff_frc", "materials_studio_export") or (
                    member and member.lower().endswith(".frc")
                ):
                    matches = [
                        r
                        for r in records
                        if r["line_start"] == start
                        and r["line_end"] == end
                        and r["section"] == loc["section"]
                    ]
                    require(
                        len(matches) == 1
                        and matches[0]["family"] == candidate["parameter_family"],
                        "FRC_section_or_family_mismatch",
                    )
                documents[str(member)] = {
                    "parser": parser,
                    "text_sha256": sha(text.encode()),
                    "recognized_FRC_records": len(records),
                }
            tokens = raw_text.replace(",", " ").split()
            values = candidate["numerical_values"]
            require(
                all(v in tokens for v in values), "numerical_transcription_mismatch"
            )
            types = candidate["atom_type_sequence"]
            require(
                not types
                or any(tokens[i : i + len(types)] == types for i in range(len(tokens))),
                "ordered_source_types_not_in_raw_record",
            )
            require(
                all(
                    re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", v)
                    for v in values
                ),
                "finite_decimal_source_values_required",
            )
            facts[cid] = {
                "source_sha256": artifact["sha256"],
                "member_sha256": sha(data) if member else None,
                "source_location": loc,
                "raw_text": raw_text,
                "raw_text_sha256": sha(raw_text.encode()),
                "parser": parser,
            }
        except (
            ValueError,
            TypeError,
            KeyError,
            OSError,
            UnicodeError,
            zipfile.BadZipFile,
            tarfile.TarError,
        ) as e:
            failures[cid] = str(e)
    return {
        "facts": facts,
        "failures": failures,
        "documents": documents,
        "archive_members": inventory,
    }
