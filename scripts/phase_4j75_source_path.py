"""Fenced, read-only exact-root capture and offline PCFF evidence matching.

No scientific package is imported. Historical parsers read literals only; their
results are new J75 facts with explicit links to unchanged historical identities.
"""

import argparse
import json
import os
import re
import stat
import subprocess
from collections import Counter
from pathlib import Path

from scripts import phase_4j74_source_intake as old
from scripts.phase_4j73_source_intake import (
    IntakeError,
    digest,
    exact,
    pin,
    read,
    require,
    signed,
    verify_signed,
    write,
)

ROOT = Path("docs/evidence/phase_4j75")
AUTHORITATIVE = "/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/island_ffref/pcff"
BASE = "7536e143a563b8111d0f4a6c66c296e5456fd4fd"
MAIN = "c118ddc06a7bd2532bea93c5e8d6ce17727726ef"
BRANCH = "codex/phase-4j75-authoritative-pcff-source-path-correction"
# Set once after read-only capture, before scientific matching or verification.
FREEZE_SHA = "8919d63fe32ff9edfa508eb963ff5e90d3fb151fa21322f81892ee4bdd4c79d9"
FAMILIES = (
    "native_PCFF_CFF91",
    "native_published_COMPASS",
    "external_PDMS_benchmark",
    "external_Huang_PDMS",
    "external_Smith_Borodin_Smith_PDMS",
    "Sun_Rigby_polysiloxane",
    "analog_secondary_contextual",
)
PARENTS = (
    *old.PARENTS,
    "docs/evidence/phase_4j74/source-file-inventory.json",
    "docs/evidence/phase_4j74/row-resolution-matrix.json",
    "docs/evidence/phase_4j74/receipt.json",
    "docs/evidence/phase_4j74/frozen-inputs.json",
)
REPLAY_FILES = (
    "scripts/phase_4j75_source_path.py",
    "scripts/phase_4j74_source_intake.py",
    "scripts/phase_4j73_source_intake.py",
    "scripts/phase_4j73_scope.py",
    "scripts/phase_4j72_c1_fence.py",
)
EXTRA_PARENTS = (
    "docs/evidence/phase_4j74/source-conflicts.json",
    *REPLAY_FILES[1:],
)
SCOPE = {
    "J75_new_QM_attempts": 0,
    "J75_MD_steps": 0,
    "converter_attempts": 0,
    "fitting_executions": 0,
    "evaluator_executions": 0,
    "parameter_generation_executions": 0,
    "reserved_validation": 0,
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "parameter_completeness": "incomplete",
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "profile": None,
    "parameter_bundle": None,
    "converter_result": None,
    "static_closure_result": None,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "source_adjudication_excluded_quarantine": True,
    "author_confirmed_reported": 28,
    "author_confirmed_audited": 0,
    "original_author_artifact_status": "not_yet_supplied",
}
CAPTURE_FILES = (
    "source-path-report.json",
    "source-root-manifest.json",
    "source-file-inventory.json",
    "source-parsing-facts.json",
    "extracted-parameter-records.json",
    "historical-files.json",
)


def scope():
    old.scope()
    from scripts.phase_4j72_c1_fence import assert_installed

    assert_installed()


def source_family(name, historical=None):
    """Seven separate categories; historical model identity is also retained."""
    n = name.lower()
    h = historical if isinstance(historical, dict) else {}
    model = h.get(
        "source_force_field_family",
        h.get("force_field", historical if isinstance(historical, str) else ""),
    )
    citation = h.get("citation", h)
    if not isinstance(citation, dict):
        citation = h  # Reading ledgers also retain a rendered citation string.
    title = citation.get("title", "") or ""
    authors = citation.get("authors", "") or ""
    if "/huang/" in n or str(model).lower().startswith("huang"):
        return FAMILIES[3]
    if "pdms_benchmark" in n:
        return FAMILIES[2]
    if Path(n).name in {"pcff.frc", "cff91.frc", "pcff_templates.dat"}:
        return FAMILIES[0]
    if Path(n).name == "compass_published.frc" or model == "COMPASS":
        return FAMILIES[1]
    if "smith" in str(model).lower() and "pdms" in str(model).lower():
        return FAMILIES[4]
    if (
        model in {"PCFF", "PCFF_related"}
        and "polysiloxane" in title.lower()
        and ("sun" in authors.lower() or "rigby" in authors.lower())
    ):
        return FAMILIES[5]
    return FAMILIES[6]


def check_declared(value):
    # Compare spelling before Path normalization or realpath resolution.
    require(str(value) == AUTHORITATIVE, "incorrect_authoritative_absolute_path")


def path_report(value, repo):
    scope()  # Must precede even exists/stat/open on the source root.
    check_declared(value)
    p = Path(value)
    exists, directory = p.exists(), p.is_dir()
    links = [str(c) for c in (*reversed(p.parents), p) if c.is_symlink()]
    st = p.lstat() if p.exists() or p.is_symlink() else None
    mounts = []
    # Read kernel metadata only; never invoke mount, find, realpath or a shell.
    for line in Path("/proc/self/mountinfo").read_text().splitlines():
        fields = line.split()
        mountpoint = re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), fields[4])
        if p == Path(mountpoint) or Path(mountpoint) in p.parents:
            split = fields.index("-")
            mounts.append(
                {
                    "mountpoint": mountpoint,
                    "mount_root": fields[3],
                    "major_minor": fields[2],
                    "filesystem": fields[split + 1],
                    "optional_fields": fields[6:split],
                }
            )
    exact_mount = next((m for m in mounts if m["mountpoint"] == value), None)
    return signed(
        {
            "schema": "island.j75-source-path.v1",
            "requested_absolute_path": value,
            "path_exists": exists,
            "directory_status": directory,
            "resolved_realpath": str(p.resolve(strict=False)),
            "inside_git_checkout": p.is_relative_to(repo.resolve()),
            "is_symlink": p.is_symlink(),
            "symlink_components": links,
            "root_identity_usable": exists and directory and not links,
            "exact_mountpoint": exact_mount,
            "nearest_mount": max(mounts, key=lambda m: len(m["mountpoint"]))
            if mounts
            else None,
            "bind_mount_assessment": "exact mountpoint; mount-root metadata retained, bind provenance not independently proved"
            if exact_mount
            else "no exact mountpoint in proc mountinfo; no separate bind-mount evidence",
            "device": st.st_dev if st else None,
            "inode": st.st_ino if st else None,
            "permissions": oct(stat.S_IMODE(st.st_mode)) if st else None,
            "readable": os.access(p, os.R_OK | os.X_OK) if directory else False,
            "fence_installed_before_source_inspection": True,
            "requested_shell_semantics": {
                "test_-d": directory,
                "realpath": str(p.resolve(strict=False)),
                "find_-type_f": "source-root-manifest.json#/files; implemented in fenced Python without shell execution",
            },
            "fallback_used": False,
        }
    )


def tree_inventory(directory):
    """Read bytes without following any symlink, including directory aliases."""
    scope()
    p = Path(directory)
    require(not p.is_symlink(), "source_root_symlink_forbidden")
    require(p.is_dir(), "source_path_unavailable")
    files = []
    for folder, directories, names in os.walk(p, followlinks=False):
        directories.sort()
        for name in [*directories, *sorted(names)]:
            require(not (Path(folder) / name).is_symlink(), "source_symlink_forbidden")
        for name in sorted(names):
            f = Path(folder) / name
            st = f.stat()
            require(stat.S_ISREG(st.st_mode), "nonregular_source_file")
            require(st.st_size <= old.MAX_FILE, "source_size_limit")
            raw = f.read_bytes()
            require(len(raw) == st.st_size, "source_changed_while_hashing")
            files.append(
                {
                    "relative_path": f.relative_to(p).as_posix(),
                    "sha256": old.sha_bytes(raw),
                    "bytes": len(raw),
                    "device": st.st_dev,
                    "inode": st.st_ino,
                    "permissions": oct(stat.S_IMODE(st.st_mode)),
                    "readable": os.access(f, os.R_OK),
                }
            )
    return sorted(files, key=lambda f: f["relative_path"])


def root_manifest(report, files):
    return signed(
        {
            "schema": "island.j75-source-root.v1",
            "source_root": AUTHORITATIVE,
            "path_report_identity": report["identity"],
            "source_root_device": report["device"],
            "source_root_inode": report["inode"],
            "recursive_file_count": len(files),
            "total_bytes": sum(f["bytes"] for f in files),
            "files": files,
            "frozen_before_format_parsing": True,
            "no_symlinks_followed": True,
        }
    )


def new_record(record):
    data = {k: v for k, v in record.items() if k not in {"identity", "record_id"}}
    data.update(
        {
            "schema": "island.j75-parameter-fact.v1",
            "source_root": AUTHORITATIVE,
            "source_model_identity": record.get("historical_record_unchanged", {}).get(
                "source_force_field_family", record["source_family"]
            ),
            "source_family": source_family(
                record["source_path"],
                record.get("historical_record_unchanged", record["source_family"]),
            ),
            "extraction_parser_fact_identity": record["identity"],
            "exact_row_resolved": False,
            "convention_resolved": False,
            "executable_bundle_qualified": False,
        }
    )
    return signed({"record_id": "j75-fact-" + digest(data), **data})


def capture(repo, root, source):
    scope()
    require(not root.exists(), "do_not_overwrite_capture")
    report = path_report(source, repo)
    require(
        not root.resolve().is_relative_to(Path(source)), "evidence_output_inside_source"
    )
    root.mkdir(parents=True)
    write(root / "source-path-report.json", report)
    if not report["path_exists"]:
        # No directory substitute, no invented source facts.
        files = []
        decision = "source_path_unavailable"
    else:
        require(
            report["root_identity_usable"] and report["readable"],
            "source_root_identity_blocked",
        )
        files = tree_inventory(source)
        decision = "parameter_incomplete"
    write(root / "source-root-manifest.json", root_manifest(report, files))
    # Freeze before any source format is parsed; no J74 capture is rerun.
    parents = {n: read(repo / n) for n in PARENTS}
    historical = {r["sha256"]: r for r in parents[old.PARENTS[5]]["files"]}
    old_records = parents[old.PARENTS[4]]["records"]
    facts, sources, parsing, seen = [], [], {}, {}
    for entry in files:
        name, sha = entry["relative_path"], entry["sha256"]
        previous = historical.get(sha)
        model = old.family_of(name, previous)
        source_record = {
            **entry,
            "MIME_type": old.mimetypes.guess_type(name)[0]
            or "application/octet-stream",
            "citation": {
                k: previous.get(k) for k in ("title", "authors", "journal_year", "doi")
            }
            if previous
            else {"status": "local_copy_only_not_author_authenticated"},
            "version": "not_established",
            "license_status": "not_established_no_redistribution",
            "redistribution_status": "source asset excluded from Git; no grant inferred",
            "exact_derived_copied_unknown": "copied_origin_unverified",
            "source_family": model,
            "source_authority": "local_copy_not_author_primary",
            "native_source_identity_complete": False,
            "primary_artifact_supplied": False,
            "exact_graph_applicability": False,
            "parser_used": "not_parsed",
            "text_extraction_status": "not_parsed",
            "duplicate_content_of": seen.get(sha),
        }
        raw = (Path(source) / name).read_bytes()
        require(old.sha_bytes(raw) == sha, "source_hash_mutation")
        rows = []
        try:
            if name.endswith(".pdf"):
                meta = old.parse_pdf(raw)
                parsing[name] = meta
                source_record["parser_used"] = meta["parser"]
                source_record["text_extraction_status"] = meta["status"]
                if previous and sha not in seen:
                    for original in old_records:
                        if original["source_file_sha256"] == sha:
                            rows.append(
                                old.make_record(
                                    source_record,
                                    original["interaction_family"],
                                    original["ordered_sites"],
                                    original["raw_value"],
                                    original["location"],
                                    historical_J71_record_id=original["record_id"],
                                    historical_J71_identity=original["identity"],
                                    historical_record_unchanged=original,
                                    units={
                                        "status": "inherited_J71_facts_not_re_adjudicated",
                                        "raw": original["raw_unit"],
                                    },
                                    analytic_equation={
                                        "status": "inherited_J71_facts_not_re_adjudicated",
                                        "raw": original["scalar_energy_equation"],
                                    },
                                )
                            )
            elif name.endswith(".frc"):
                text = raw.decode("utf-8")
                source_record["version"] = old.version_lines(text)
                source_record["citation"].update(old.frc_provenance(text))
                if name == "pcff.frc" and sha == old.PCFF_SHA:
                    citation = next(
                        r["citation"]
                        for r in parents[old.PARENTS[3]]["rows"].values()
                        if r.get("source_file") == "pcff.frc"
                    )
                    source_record["citation"].update(citation)
                    source_record["native_source_identity_complete"] = True
                    source_record["source_authority"] = (
                        "byte_exact_copy_of_retained_LAMMPS_PCFF_release; not_new_author_primary"
                    )
                    source_record["exact_derived_copied_unknown"] = (
                        "exact_copy_of_pinned_distribution"
                    )
                rows, parsing[name] = old.parse_frc(text, source_record)
                source_record["parser_used"] = (
                    "J74 literal FRC parser v1 (hash-bound dependency)"
                )
                source_record["text_extraction_status"] = "decoded_full_UTF8"
            elif name.endswith(".zip"):
                rows, parsing[name] = old.parse_archive(raw, source_record)
                source_record["version"] = parsing[name]["archive_comment_revision"]
                source_record["parser_used"] = (
                    "Python zipfile + J74 literal topology/command parser v1"
                )
                source_record["text_extraction_status"] = "all_members_read_as_data"
            else:
                text = raw.decode("utf-8")
                parsing[name] = {
                    "lines": len(text.splitlines()),
                    "text_sha256": old.sha_bytes(text.encode()),
                    "not_executed": True,
                }
                source_record["parser_used"] = "Python UTF-8 reader"
                source_record["text_extraction_status"] = "decoded_full_UTF8"
        except (
            OSError,
            ValueError,
            ImportError,
            EOFError,
            old.zipfile.BadZipFile,
            RecursionError,
        ) as error:
            source_record["text_extraction_status"] = "parse_failed"
            source_record["parser_used"] = "attempted_format_parser"
            parsing[name] = {
                "status": "unreadable_no_OCR_or_substitution",
                "error_type": type(error).__name__,
                "error": str(error)[:400],
            }
        facts.extend(new_record(r) for r in rows)
        source_record["source_model_identity"] = model
        source_record["source_family"] = source_family(name, previous)
        source_record["source_root"] = AUTHORITATIVE
        source_record["upstream_repository_or_DOI"] = source_record["citation"].get(
            "repository", source_record["citation"].get("doi", "not_established")
        )
        sources.append(signed(source_record))
        seen.setdefault(sha, name)
        print(
            json.dumps(
                {
                    "path": name,
                    "status": source_record["text_extraction_status"],
                    "facts": len(rows),
                }
            ),
            flush=True,
        )
    if files:
        exact(tree_inventory(source), files, "source_tree_changed_during_capture")
    write(
        root / "source-file-inventory.json",
        signed(
            {
                "schema": "island.j75-source-inventory.v1",
                "source_root": AUTHORITATIVE,
                "files": sources,
                "source_assets_committed": False,
            }
        ),
    )
    write(
        root / "source-parsing-facts.json",
        signed(
            {
                "schema": "island.j75-parsing.v1",
                "files": parsing,
                "scientific_execution": False,
            }
        ),
    )
    write(
        root / "extracted-parameter-records.json",
        signed(
            {
                "schema": "island.j75-parameter-records.v1",
                "records": facts,
                "normalization_or_parameter_generation_performed": False,
            }
        ),
    )
    tree = subprocess.run(
        ["git", "ls-tree", "-r", "--full-tree", BASE],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    history = {}
    for line in tree.splitlines():
        meta, name = line.split("\t", 1)
        p = repo / name
        history[name] = (
            old.sha_bytes(str(p.readlink()).encode())
            if meta.startswith("120000")
            else pin(p)
        )
    write(root / "historical-files.json", history)
    freeze = signed(
        {
            "schema": "island.j75-capture-freeze.v1",
            "base_commit": BASE,
            "fetched_origin_main": MAIN,
            "branch": BRANCH,
            "decision": decision,
            "capture_files": {n: pin(root / n) for n in CAPTURE_FILES},
            "parents": {n: pin(repo / n) for n in PARENTS},
            "scope": SCOPE,
        }
    )
    write(root / "frozen-inputs.json", freeze)
    print(
        json.dumps(
            {
                "freeze_sha256": pin(root / "frozen-inputs.json"),
                "files": len(files),
                "facts": len(facts),
                "historical_paths": len(history),
            }
        )
    )


def inputs(repo, root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "fabricated_or_resigned_capture_freeze",
    )
    freeze = read(root / "frozen-inputs.json")
    verify_signed(freeze)
    exact(freeze["scope"], SCOPE, "scientific_scope_activation")
    require(
        freeze["base_commit"] == BASE
        and freeze["fetched_origin_main"] == MAIN
        and freeze["branch"] == BRANCH,
        "capture_base_identity",
    )
    require(
        set(freeze["parents"]) == set(PARENTS)
        and set(freeze["capture_files"]) == set(CAPTURE_FILES),
        "explicit_capture_inventory",
    )
    for n, sha in freeze["parents"].items():
        require(pin(repo / n) == sha, "historical_input_mutation:" + n)
    for n, sha in freeze["capture_files"].items():
        require(pin(root / n) == sha, "source_fact_hash_mutation:" + n)
    history = read(root / "historical-files.json")
    for n in EXTRA_PARENTS:
        require(pin(repo / n) == history[n], "historical_dependency_mutation:" + n)
    return freeze


def validate_inventory(document):
    verify_signed(document)
    check_declared(document["source_root"])
    seen = set()
    for s in document["files"]:
        verify_signed(s)
        require(s["relative_path"] not in seen, "duplicate_source_identity")
        seen.add(s["relative_path"])
        require(re.fullmatch(r"[0-9a-f]{64}", s["sha256"]), "missing_source_hash")
        require(
            s.get("license_status") and s.get("redistribution_status"),
            "missing_license_metadata",
        )
        require(
            s["source_family"] in FAMILIES and not s["primary_artifact_supplied"],
            "local_or_analog_authority_promotion",
        )
        if Path(s["relative_path"]).name == "compass_published.frc":
            require(
                s["source_family"] == FAMILIES[1], "COMPASS_PCFF_cross_family_promotion"
            )
        if s["native_source_identity_complete"]:
            require(
                s["relative_path"] == "pcff.frc"
                and s["sha256"] == old.PCFF_SHA
                and s["citation"]["revision"]
                == "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
                "native_release_identity",
            )


def compare_copies(current, previous, report, previous_directory):
    by_old = {s["relative_path"]: s for s in previous}
    by_new = {s["relative_path"]: s for s in current}
    require(
        len(by_old) == len(previous) and len(by_new) == len(current),
        "duplicate_copy_identity",
    )
    return signed(
        {
            "schema": "island.j75-checkout-comparison.v1",
            "authoritative_root": AUTHORITATIVE,
            "J74_recorded_checkout_root": previous_directory,
            "source_root_identity": {
                "authoritative_realpath": report["resolved_realpath"],
                "device": report["device"],
                "inode": report["inode"],
                "J74_recorded_root_same_path": previous_directory == AUTHORITATIVE,
                "independent_external_root": False,
                "interpretation": "same physical current root compared with frozen J74 file hashes; no historical inode claim",
            },
            "files": [
                {
                    "relative_path": n,
                    "authoritative_SHA256": by_new.get(n, {}).get("sha256"),
                    "J74_checkout_SHA256": by_old.get(n, {}).get("sha256"),
                    "authoritative_bytes": by_new.get(n, {}).get("bytes"),
                    "J74_checkout_bytes": by_old.get(n, {}).get("bytes"),
                    "relative_path_equality": n in by_new and n in by_old,
                    "size_equality": n in by_new
                    and n in by_old
                    and by_new[n]["bytes"] == by_old[n]["bytes"],
                    "byte_equality": n in by_new
                    and n in by_old
                    and by_new[n]["sha256"] == by_old[n]["sha256"],
                    "exact_duplicate": n in by_new
                    and n in by_old
                    and by_new[n]["sha256"] == by_old[n]["sha256"]
                    and by_new[n]["bytes"] == by_old[n]["bytes"],
                    "source_authority_inferred": False,
                }
                for n in sorted(set(by_new) | set(by_old))
            ],
            "hash_equality_proves_authorship": False,
        }
    )


def match(record, target):
    # Comparison fields are deliberately independent; unresolved is never True.
    converted = {**record, "source_family": record["source_model_identity"]}
    result = old.compatible(converted, target)
    native = result["source_family_match"] and record["source_family"] == FAMILIES[0]
    contract = target["required_evidence"]["equation_and_units"][
        "historical_unresolved_contract"
    ]
    # A self-declared complete record cannot close an unresolved target
    # contract. J75 never adjudicates or upgrades that historical authority.
    target_closed = (
        contract["equation_status"] == "source_complete"
        and contract["units_status"] == "source_complete"
    )
    native = native and target_closed
    result["target_convention_authority_resolved"] = target_closed
    result["convention_match"] = result["convention_match"] and target_closed
    result.update(
        {
            "phase_multiplicity_match": native
            and record["phase_or_multiplicity"].get("status")
            == "complete_source_contract",
            "orientation_order_match": native
            and record["orientation_order_convention"].get("status")
            == "complete_source_contract",
            "equilibrium_dependency_match": native
            and record["equilibrium_dependencies"].get("status")
            == "complete_source_contract",
            "exclusion_match": native
            and record["applicable_exclusions"].get("status")
            == "complete_source_contract",
            "charge_convention_match": native
            and record["charge_origin"].get("status") == "complete_source_contract",
            "LJ_and_mixing_match": native
            and record["LJ_equation_and_mixing"].get("status")
            == "complete_source_contract",
            "source_authority": record["source_authority"],
        }
    )
    return result


def candidate_set(candidates):
    require(
        len(candidates) == len({c["record_id"] for c in candidates}),
        "duplicate_candidate_identity",
    )


def derive(repo, root):
    freeze = inputs(repo, root)
    sources = read(root / "source-file-inventory.json")
    validate_inventory(sources)
    facts = read(root / "extracted-parameter-records.json")["records"]
    candidate_set(facts)
    for r in facts:
        verify_signed(r)
        require(
            not r["native_PCFF_eligible"]
            and not r["used_to_fill_missing_row"]
            and not r["exact_row_resolved"],
            "analog_promotion",
        )
    parsing = read(root / "source-parsing-facts.json")["files"]
    requests = read(repo / old.PARENTS[0])
    previous = read(repo / "docs/evidence/phase_4j74/row-resolution-matrix.json")
    require(
        len(requests["requests"]) == len(previous["instances"]) == 378,
        "incomplete_378_input_coverage",
    )
    previous_by_id = {r["id"]: r for r in previous["instances"]}
    historical_by_id = {
        r["historical_J71_record_id"]: r
        for r in facts
        if "historical_J71_record_id" in r
    }
    resolutions = []
    for t in requests["requests"]:
        candidates, queries = {}, []
        for s in sources["files"]:
            if s["source_family"] != FAMILIES[0] or not s["relative_path"].endswith(
                ".frc"
            ):
                continue
            equivalences = parsing[s["relative_path"]]["ordinary_equivalences"]
            role = "angle" if t["family"] == "angle-angle" else "torsion"
            mapped = [
                [equivalences.get(v, {}).get(role, v) for v in types]
                for types in t["symmetry_equivalent_types"]
            ]
            selected = [
                r
                for r in facts
                if r["source_path"] == s["relative_path"]
                and r["interaction_family"] == t["family"]
                and (
                    r["ordered_atom_types"] in mapped
                    or r["ordered_atom_types"] in t["symmetry_equivalent_types"]
                )
                and not r.get("wildcard")
                and not r.get("auto_equivalence")
            ]
            for r in selected:
                candidates[r["record_id"]] = r
            queries.append(
                {
                    "source_path": s["relative_path"],
                    "sha256": s["sha256"],
                    "role": role,
                    "mapped_source_types": mapped,
                    "matches": [r["record_id"] for r in selected],
                    "wildcard_or_auto_used": False,
                }
            )
        for a in t["unresolved_alternatives"]:
            if a["record_id"] in historical_by_id:
                r = historical_by_id[a["record_id"]]
                candidates[r["record_id"]] = r
        # External benchmark rows are family candidates only, not exact type analogs.
        for r in facts:
            if "::" in r["source_path"] and r["interaction_family"] == t["family"]:
                candidates[r["record_id"]] = r
        cs = [
            {
                "record_id": r["record_id"],
                "source_family": r["source_family"],
                "source_model_identity": r["source_model_identity"],
                "evidence_reference": {
                    "source_root": AUTHORITATIVE,
                    "path": r["source_path"],
                    "sha256": r["source_sha256"],
                    "container_sha256": r["container_sha256"],
                    "location": r["location"],
                },
                **match(r, t),
            }
            for r in sorted(candidates.values(), key=lambda r: r["record_id"])
        ]
        candidate_set(cs)
        native_count = sum(len(q["matches"]) for q in queries)
        require(
            native_count == 0, "unexpected_exact_source_match_requires_separate_audit"
        )
        state = old.STATUS_MAP[t["J72_status"]]
        resolutions.append(
            signed(
                {
                    "id": t["id"],
                    "J70_identity": t["J70_identity"],
                    "J72_identity": t["J72_identity"],
                    "J73_identity": t["identity"],
                    "J74_identity": previous_by_id[t["id"]]["identity"],
                    "ordered_atom_types": t["ordered_types"],
                    "site_ids": t["stable_sites"],
                    "local_graph_identity": requests["graph_identity"],
                    "local_graph_neighborhood": t["chemical_environment"],
                    "interaction_family": t["family"],
                    "candidate_source_records": cs,
                    "native_source_queries": queries,
                    "resolution_status": state,
                    "statuses": {s: s == state for s in old.STATES},
                    "author_confirmed_reported": "28 historical user attestations; no per-row primary evidence",
                    "author_confirmed_audited": False,
                    "blocking_reason": "Native Phi/AAT and complete source convention/dependency contract unresolved"
                    if state == "convention_blocked"
                    else "No exact native row found; external/analog records cannot supply native PCFF authority",
                    "required_evidence": t["required_evidence"],
                }
            )
        )
    counts = dict(Counter(r["resolution_status"] for r in resolutions))
    exact(
        counts,
        {"analog_only": 294, "convention_blocked": 60, "source_missing": 24},
        "378_statuses_changed",
    )
    ledger = signed(
        {
            "schema": "island.j75-artifact-ledger.v1",
            "source_files": sources["files"],
            "archive_members": {
                n: p["members"] for n, p in parsing.items() if "members" in p
            },
            "author_confirmed_reported": 28,
            "author_confirmed_audited": 0,
            "trusted_new_author_primary_artifacts": 0,
            "hash_proves_authority": False,
            "scope": SCOPE,
        }
    )
    old_inventory = read(repo / "docs/evidence/phase_4j74/source-file-inventory.json")
    families = [
        {
            "family": f,
            "extracted_facts": sum(r["source_family"] == f for r in facts),
            "exact_native_rows_supplied": 0,
            "exact_rows_unresolved": 378,
            "scope_of_counts": "exact J73 native PCFF requests; no completeness claim for a different external model",
            "unresolved_row_ids": [r["id"] for r in resolutions],
            "missing_Class_II_families": dict(
                Counter(r["interaction_family"] for r in resolutions)
            ),
            "convention_fields_unresolved": [
                "phase/multiplicity",
                "equilibrium dependencies",
                "charge origin",
                "exclusions",
                "LJ and mixing",
                "exact graph authority",
            ],
            "automatic_assignment": False,
        }
        for f in FAMILIES
    ]
    return {
        "checkout-comparison.json": compare_copies(
            sources["files"],
            old_inventory["files"],
            read(root / "source-path-report.json"),
            old_inventory["actual_directory"],
        ),
        "source-artifact-ledger.json": ledger,
        "row-resolution-matrix.json": signed(
            {
                "schema": "island.j75-row-resolution.v1",
                "instances": resolutions,
                "counts": counts,
                "scope": SCOPE,
            }
        ),
        "convention-resolution.json": signed(
            {
                "schema": "island.j75-conventions.v1",
                "exact_rows_resolved": 0,
                "convention_resolved_rows": 0,
                "executable_bundle_qualified": False,
                "rows": [
                    {
                        "id": r["id"],
                        "remaining_fields": r["required_evidence"][
                            "J72_remaining_missing_conventions"
                        ],
                        "complete_source_contract": False,
                    }
                    for r in resolutions
                ],
                "scope": SCOPE,
            }
        ),
        "source-conflicts.json": signed(
            {
                "schema": "island.j75-conflicts.v1",
                "historical_J72_matrix_sha256": pin(repo / old.PARENTS[6]),
                "equal_duplicate_equivalences": 16,
                "unresolved_unequal_COMPASS_assignments": 2,
                "new_adjudications": 0,
                "original_alternatives_preserved": True,
                "numeric_average_or_residual_selection": False,
                "distinct_native_release_alternatives": read(
                    repo / "docs/evidence/phase_4j74/source-conflicts.json"
                )["distinct_native_release_alternatives"],
                "scope": SCOPE,
            }
        ),
        "family-coverage-report.json": signed(
            {
                "schema": "island.j75-family-coverage.v1",
                "reports": families,
                "interaction_families": {
                    f: {
                        "extracted_facts": sum(
                            r["interaction_family"] == f for r in facts
                        ),
                        "absent_means": "not supplied; never zero-filled",
                    }
                    for f in sorted({v[0] for v in old.FAMILIES.values()})
                },
                "scope": SCOPE,
            }
        ),
        "license-report.json": signed(
            {
                "schema": "island.j75-license.v1",
                "source_assets_committed": False,
                "files": [
                    {
                        k: s[k]
                        for k in (
                            "relative_path",
                            "sha256",
                            "license_status",
                            "redistribution_status",
                        )
                    }
                    for s in sources["files"]
                ],
                "scope": SCOPE,
            }
        ),
        "capability-report.json": signed(
            {
                "schema": "island.j75-capability.v1",
                "decision": freeze["decision"],
                "unresolved_instances": 378,
                "resolved_native_rows": 0,
                "strict_converter": "not_run; prohibited in J75",
                "static_closure": "not_run; prohibited in J75",
                "execution_profile_created": False,
                "ordinary_pcff_behavior_changed": False,
                "scope": SCOPE,
            }
        ),
    }


def replay(repo, root, create=False):
    documents = derive(repo, root)
    for name, expected in documents.items():
        if create:
            require(not (root / name).exists(), "do_not_overwrite_derived_report")
            write(root / name, expected)
        exact(read(root / name), expected, "resigned_derived_evidence:" + name)
    return {
        "decision": read(root / "frozen-inputs.json")["decision"],
        "counts": documents["row-resolution-matrix.json"]["counts"],
        "resolved_native_rows": 0,
        "scope": SCOPE,
    }


def originals(root):
    scope()
    check_declared(read(root / "source-root-manifest.json")["source_root"])
    exact(
        tree_inventory(AUTHORITATIVE),
        read(root / "source-root-manifest.json")["files"],
        "authoritative_source_tree_mutation",
    )
    return {"source_bytes_preserved": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--directory", type=Path, default=ROOT)
    parser.add_argument("--source", default=AUTHORITATIVE)
    parser.add_argument("mode", choices=("capture", "build", "replay", "originals"))
    args = parser.parse_args()
    try:
        scope()
        result = (
            capture(args.repo, args.directory, args.source)
            if args.mode == "capture"
            else originals(args.directory)
            if args.mode == "originals"
            else replay(args.repo, args.directory, args.mode == "build")
        )
        print(json.dumps(result, indent=2))
        return 0
    except (IntakeError, OSError, RuntimeError) as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
