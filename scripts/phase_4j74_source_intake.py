"""Read-only local source facts and exact-row requests; no evaluator or converter.

PDF parsing is optional and lazy, used only during explicitly fenced capture.
Offline replay uses frozen extracted facts and the immutable J73 requests.
"""

import argparse
import hashlib
import io
import json
import mimetypes
import re
import subprocess
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

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
from scripts.phase_4j73_source_intake import (
    scope as inherited_scope,
)

ROOT = Path("docs/evidence/phase_4j74")
J73 = "dff35a0dfbca2c0517a9a7d5fe8e5032ac2d8375"
MAIN = "c118ddc06a7bd2532bea93c5e8d6ce17727726ef"
BASE = "a0fc6003dd637dc4920eab03f0b881dccd1dbc33"
BRANCH = "codex/phase-4j74-pcff-local-source-resolution"
DECLARED = "/island_ffref/pcff"
FREEZE_SHA = "76890c42892ae5bf84169c3a20072096da8024a22bdad66b810c5122eb01379a"
PARENTS = (
    "docs/evidence/phase_4j73/source-intake.json",
    "docs/evidence/phase_4j73/receipt.json",
    "docs/evidence/phase_4j73/artifact-contract.json",
    "docs/evidence/phase_4j69/source-row-inventory.json",
    "docs/evidence/phase_4j71/parameter-records.json",
    "docs/evidence/phase_4j71/reading-status.json",
    "docs/evidence/phase_4j72/source-conflict-matrix.json",
    "docs/evidence/phase_4j72_c1/quarantine-manifest.json",
    "docs/evidence/phase_4j72_c1/receipt.json",
)
SCOPE = {
    "J74_new_QM_attempts": 0,
    "J74_MD_steps": 0,
    "converter_attempts": 0,
    "fitting_executions": 0,
    "parameter_generation_executions": 0,
    "scientific_energy_force_pressure_evaluations": 0,
    "reserved_validation": 0,
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "profile": None,
    "parameter_bundle": None,
    "converter_result": None,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested; no authority promotion",
}
STATES = (
    "analog_only",
    "convention_blocked",
    "source_missing",
    "primary_artifact_supplied",
    "exact_row_resolved",
    "convention_resolved",
    "executable_bundle_qualified",
)
STATUS_MAP = {
    "analog_unqualified": "analog_only",
    "convention_unresolved": "convention_blocked",
    "source_missing": "source_missing",
}
FAMILIES = {
    "quartic_bond": ("bond", 2, "r0 K2 K3 K4"),
    "quadratic_bond": ("bond", 2, "r0 K2"),
    "quartic_angle": ("angle", 3, "theta0 K2 K3 K4"),
    "quadratic_angle": ("angle", 3, "theta0 K2"),
    "torsion_3": ("torsion", 4, "V1 Phi1 V2 Phi2 V3 Phi3"),
    "torsion_1": ("torsion", 4, "V Phi n"),
    "wilson_out_of_plane": ("improper", 4, "K chi0"),
    "bond-bond": ("bond-bond", 3, "Kbb"),
    "bond-bond_1_3": ("bond-bond_1_3", 4, "Kbb13"),
    "bond-angle": ("bond-angle", 3, "Kba_left Kba_right"),
    "angle-angle": ("angle-angle", 4, "Kaa"),
    "end_bond-torsion_3": ("end_bond-torsion_3", 4, "B1 B2 B3 C1 C2 C3"),
    "middle_bond-torsion_3": ("middle_bond-torsion_3", 4, "A1 A2 A3"),
    "angle-torsion_3": ("angle-torsion_3", 4, "D1 D2 D3 E1 E2 E3"),
    "angle-angle-torsion_1": ("angle-angle-torsion_1", 4, "Kaat"),
    "torsion-torsion_1": ("torsion-torsion", 5, "source_fields"),
    "nonbond(9-6)": ("nonbonded", 1, "rmin epsilon"),
    "bond_increments": ("bond_increments", 2, "delta_i delta_j"),
}
COEFF = {
    "bond_coeff": "bond",
    "angle_coeff": "angle",
    "dihedral_coeff": "torsion",
    "improper_coeff": "improper",
    "pair_coeff": "nonbonded",
}
SUBSTYLE = {
    ("angle_coeff", "bb"): "bond-bond",
    ("angle_coeff", "ba"): "bond-angle",
    ("dihedral_coeff", "mbt"): "middle_bond-torsion_3",
    ("dihedral_coeff", "ebt"): "end_bond-torsion_3",
    ("dihedral_coeff", "at"): "angle-torsion_3",
    ("dihedral_coeff", "aat"): "angle-angle-torsion_1",
    ("dihedral_coeff", "bb13"): "bond-bond_1_3",
    ("improper_coeff", "aa"): "angle-angle",
}
MAX_FILE = 256 * 1024**2
MAX_MEMBER = 64 * 1024**2
PCFF_SHA = "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"
_scope_installed = False


def scope():
    """Require C1 and reject named Python scientific drivers cooperatively."""
    global _scope_installed
    inherited_scope()
    if _scope_installed:
        return
    from scripts import phase_4j72_c1_fence as fence

    previous = subprocess.Popen

    def guarded(argv, *args, **kwargs):
        if isinstance(argv, (list, tuple)) and len(argv) > 1:
            word = str(argv[2] if argv[1] == "-m" and len(argv) > 2 else argv[1])
            if re.search(
                r"(?:^|[/._-])(?:evaluate_energy|evaluate_forces|evaluate_pressure|energy_evaluation|force_evaluation|pressure_evaluation|run_converter|optimize_coefficients)(?:$|[/._-])",
                word,
            ):
                fence.reject(
                    "subprocess.Popen",
                    list(map(str, argv)),
                    "J74_scientific_evaluator_or_converter_driver_forbidden",
                )
        return previous(argv, *args, **kwargs)

    subprocess.Popen = guarded
    _scope_installed = True


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def family_of(name, historical=None):
    n = name.lower()
    if "/huang/" in n:
        return "external_Huang_PDMS"
    if "pdms_benchmark" in n:
        if "/compass/" in n:
            return "external_Xiang_COMPASS_implementation"
        if "/opls/" in n:
            return "external_Xiang_OPLS_implementation"
        if "/ua/" in n:
            return "external_Xiang_united_atom_implementation"
        if "/dreiding/" in n:
            return "external_Xiang_Dreiding_implementation"
        return "external_Xiang_PDMS_benchmark"
    if n.endswith(("/pcff.frc", "/cff91.frc")) or n in {"pcff.frc", "cff91.frc"}:
        return "native_PCFF_CFF91_distribution_candidate"
    if n.endswith("compass_published.frc"):
        return "native_COMPASS_distribution_candidate"
    if historical:
        return historical.get("force_field", "contextual_literature")
    if "template" in n:
        return "native_named_typing_context_unverified"
    return "unclassified_contextual_source"


def inventory(directory):
    """Deterministic bytes-only inventory; never follow symlinks or execute."""
    scope()
    directory = Path(directory)
    if not directory.is_dir():
        return []
    result = []
    for path in sorted(directory.rglob("*")):
        require(not path.is_symlink(), "source_symlink_requires_explicit_provenance")
        if not path.is_file():
            continue
        require(path.stat().st_size <= MAX_FILE, "source_size_limit")
        raw = path.read_bytes()
        result.append(
            signed(
                {
                    "relative_path": path.relative_to(directory).as_posix(),
                    "bytes": len(raw),
                    "sha256": sha_bytes(raw),
                    "MIME_type": mimetypes.guess_type(path.name)[0]
                    or "application/octet-stream",
                }
            )
        )
    return result


def version_lines(text):
    return [
        {"line": n, "text": line.strip()}
        for n, line in enumerate(text.splitlines(), 1)
        if line.startswith("#version")
    ]


def frc_provenance(text):
    """Literal citation/license labels, never an inferred redistribution grant."""
    references, current = {}, None
    urls = set()
    notices = []
    for number, line in enumerate(text.splitlines(), 1):
        words = line.split()
        if line.startswith("#reference"):
            current = words[1] if len(words) > 1 else "unspecified"
            references.setdefault(current, [])
        elif line.startswith("#"):
            current = None
        elif current and line.strip():
            references[current].append({"line": number, "text": line.strip()[:300]})
        urls.update(re.findall(r"https?://[^\s<>]+", line))
        if any(
            w in line.lower()
            for w in ("copyright", "redistribut", "license", "warranty", "permission")
        ):
            notices.append({"line": number, "text": line.strip()[:300]})
    return {
        "literal_source_references": references,
        "literal_URLs": sorted(urls),
        "license_notices_not_a_grant": notices,
    }


def make_record(source, family, ordered, values, location, **extra):
    """Facts retain raw numeric lexemes; no normalization or coefficient use."""
    data = {
        "schema": "island.j74-extracted-fact.v1",
        "source_path": source["relative_path"],
        "source_sha256": source["sha256"],
        "container_sha256": source.get("container_sha256"),
        "citation": source["citation"],
        "version": source.get("version", "not_established"),
        "source_family": source["source_family"],
        "source_authority": source["source_authority"],
        "primary_artifact_supplied": False,
        "location": location,
        "interaction_family": family,
        "ordered_atom_types": ordered,
        "mapped_chemical_environment": "source types only; no native PCFF graph transfer inferred",
        "equilibrium_value": None,
        "coefficients_raw": values,
        "normalized_values": None,
        "units": {"status": "unresolved", "raw": []},
        "analytic_equation": {"status": "unresolved", "raw": []},
        "analytic_force_equation": "F_i = -partial E_source / partial x_i; not evaluated",
        "sign_convention": "source order and signs retained; unresolved unless source states it",
        "phase_or_multiplicity": {"status": "unresolved"},
        "orientation_order_convention": {"status": "unresolved"},
        "equilibrium_dependencies": {"status": "unresolved"},
        "applicable_exclusions": {"status": "unresolved"},
        "charge_origin": {"status": "unresolved"},
        "LJ_equation_and_mixing": {"status": "unresolved"},
        "exact_graph_applicability": False,
        "confidence": "literal source facts; no physical validation or authority by hash",
        "native_PCFF_eligible": False,
        "used_to_fill_missing_row": False,
        **extra,
    }
    record_id = "j74-fact-" + digest(data)
    return signed({"record_id": record_id, **data})


def parse_frc(text, source):
    """Read section rows, explicit equivalences and relevant numeric facts.

    Wildcards and auto-equivalences are retained only as noneligible context.
    All rows are parsed/count-bound; extraction is limited to HMDS component
    types or explicitly Si-typed environments, not every unrelated organic row.
    """
    sections, records, equivalences, silicon = {}, [], {}, set()
    section, namespace = "", ""
    row_counts = Counter()
    for number, line in enumerate(text.splitlines(), 1):
        words = line.split()
        if line.startswith("#") and words:
            section = words[0][1:]
            namespace = words[1] if len(words) > 1 else "not_stated"
            sections.setdefault(
                section + ":" + namespace, {"line": number, "header": []}
            )
            continue
        key = section + ":" + namespace
        if not words:
            continue
        if words[0].startswith(("!", ">")):
            if key in sections and len(sections[key]["header"]) < 28:
                sections[key]["header"].append({"line": number, "text": line.strip()})
            continue
        if not re.fullmatch(r"\d+(?:\.\d+)?", words[0]) or len(words) < 3:
            continue
        row_counts[key] += 1
        if section == "atom_types" and len(words) >= 5 and words[4] == "Si":
            silicon.add(words[2])
        if section == "equivalence" and len(words) >= 8:
            equivalences[words[2]] = {
                "line": number,
                "nonbonded": words[3],
                "bond": words[4],
                "angle": words[5],
                "torsion": words[6],
                "improper": words[7],
            }
        if section not in FAMILIES or len(words) < 3 + FAMILIES[section][1]:
            continue
        family, arity, labels = FAMILIES[section]
        ordered = words[2 : 2 + arity]
        values = words[2 + arity :]
        try:
            for value in values:
                float(value)
        except ValueError:
            continue
        if not (
            any(t in silicon or t in {"sio", "osi"} for t in ordered)
            or set(ordered) <= {"c3", "hc", "c", "h", "osi", "sio"}
        ):
            continue
        header = sections[key]["header"]
        record = make_record(
            source,
            family,
            ordered,
            values,
            {
                "line": number,
                "section": section,
                "namespace": namespace,
                "source_reference": words[1],
            },
            source_field_labels=labels.split(),
            row_version=words[0],
            source_section_header=header,
            source_section_row_count=None,
            wildcard=any("*" in t for t in ordered),
            auto_equivalence="auto" in namespace,
            analytic_equation={
                "status": "literal_header_only_not_independently_convention_closed",
                "raw": [
                    h
                    for h in header
                    if any(x in h["text"] for x in ("=", "Phi", "Energy", "energy"))
                ],
            },
            units={
                "status": "literal_header_only",
                "raw": [
                    h
                    for h in header
                    if any(
                        x in h["text"].lower()
                        for x in ("unit", "kcal", "degree", "radian", "angstrom")
                    )
                ],
            },
            equilibrium_value=values[0]
            if section
            in {"quartic_bond", "quartic_angle", "quadratic_bond", "quadratic_angle"}
            and values
            else None,
            phase_or_multiplicity={
                "status": "raw_header_and_torsion_phase_fields_only",
                "raw_values": values[1::2] if section == "torsion_3" else [],
            },
            native_PCFF_eligible=False,
        )
        records.append(record)
    return records, {
        "section_row_counts": dict(row_counts),
        "section_headers": sections,
        "ordinary_equivalences": equivalences,
        "auto_equivalences_used": False,
        "silicon_source_types": sorted(silicon),
        "extraction_filter": "all HMDS component-only rows and all explicitly silicon/osi rows; unrelated rows counted, not transferred",
    }


def data_topology(text):
    """Read static type/topology facts, never evaluate coordinates or forces."""
    section, atoms, masses = "", {}, {}
    counts, maps = Counter(), defaultdict(set)
    sections = {
        "Atoms",
        "Masses",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Velocities",
    }
    arity = {"Bonds": 2, "Angles": 3, "Dihedrals": 4, "Impropers": 4}
    for number, line in enumerate(text.splitlines(), 1):
        content, _, comment = line.partition("#")
        words = content.split()
        if not words:
            continue
        if words[0] in sections:
            section = words[0]
            continue
        if words[0][0].isalpha():
            section = ""
        if not words[0].isdigit():
            continue
        if section == "Masses" and len(words) >= 2:
            masses[words[0]] = {
                "raw_mass": words[1],
                "source_comment": comment.strip(),
                "line": number,
            }
        elif section == "Atoms" and len(words) >= 7:
            atoms[words[0]] = (words[2], words[3])
            counts[section] += 1
        elif section in arity and len(words) >= arity[section] + 2:
            ids = words[2 : arity[section] + 2]
            if all(i in atoms for i in ids):
                maps[(section, words[1])].add(tuple(atoms[i][0] for i in ids))
            counts[section] += 1
    charges = defaultdict(set)
    for atom_type, charge in atoms.values():
        charges[atom_type].add(charge)
    return {
        "counts": dict(counts),
        "mass_type_metadata": masses,
        "observed_atom_types": sorted({t for t, q in atoms.values()}),
        "observed_charge_values_by_type": {
            t: sorted(q) for t, q in sorted(charges.items())
        },
        "interaction_type_tuples": {
            s + ":" + t: [list(v) for v in sorted(vs)]
            for (s, t), vs in sorted(maps.items())
        },
        "coordinates_evaluated": False,
        "native_PCFF_graph_applicability": False,
    }


def parse_commands(text, source, topology, conventions):
    records, counts = [], Counter()
    commands = {
        "bond_coeff": "Bonds",
        "angle_coeff": "Angles",
        "dihedral_coeff": "Dihedrals",
        "improper_coeff": "Impropers",
    }
    for number, line in enumerate(text.splitlines(), 1):
        content, _, comment = line.partition("#")
        words = content.split()
        if not words or words[0] not in COEFF:
            continue
        command = words[0]
        counts[command] += 1
        if len(words) < 3:
            continue
        index = words[1]
        tuples = topology.get("interaction_type_tuples", {}).get(
            commands.get(command, "") + ":" + index, []
        )
        if command == "pair_coeff":
            if not set(words[1:3]) <= set(topology.get("observed_atom_types", [])):
                continue
            ordered, value_start = words[1:3], 3
        else:
            if not tuples:
                continue  # Counted as uninstantiated source rows, never fallback.
            ordered, value_start = tuples[0], 2
        family = COEFF[command]
        substyle = words[value_start] if len(words) > value_start else ""
        if (command, substyle) in SUBSTYLE:
            family = SUBSTYLE[(command, substyle)]
            value_start += 1
        values = words[value_start:]
        records.append(
            make_record(
                source,
                family,
                ordered,
                values,
                {
                    "line": number,
                    "command": command,
                    "type_id": index,
                    "comment": comment.strip()[:200],
                },
                source_topology_type_tuples=tuples,
                source_style_contract=conventions,
                units={
                    "status": "implementation_style_reference_only",
                    "raw": conventions.get("units", []),
                },
                analytic_equation={
                    "status": "requires_exact_source_style_equation_and_version",
                    "raw": conventions.get(command.replace("_coeff", "_style"), []),
                },
                applicable_exclusions={
                    "status": "implementation_declaration_only",
                    "raw": conventions.get("special_bonds", []),
                },
                LJ_equation_and_mixing={
                    "status": "implementation_declaration_only",
                    "pair_style": conventions.get("pair_style", []),
                    "pair_modify": conventions.get("pair_modify", []),
                },
                charge_origin={
                    "status": "observed_values_not_author_origin",
                    "observed_values": topology.get(
                        "observed_charge_values_by_type", {}
                    ),
                },
            )
        )
    return records, dict(counts)


def parse_archive(raw, source):
    records, members, texts, topologies = [], [], {}, {}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = [n for n in archive.namelist() if not n.endswith("/")]
        require(len(names) == len(set(names)), "duplicate_archive_source_identity")
        require(
            sum(x.file_size for x in archive.infolist()) <= MAX_FILE,
            "archive_expansion_limit",
        )
        revision = archive.comment.decode("utf-8", "replace")
        for name in sorted(names):
            require(
                not Path(name).is_absolute() and ".." not in Path(name).parts,
                "archive_member_path",
            )
            info = archive.getinfo(name)
            require(info.file_size <= MAX_MEMBER, "archive_member_size_limit")
            data = archive.read(name)
            try:
                text = data.decode("utf-8")
                texts[name] = text
                status = "decoded_full_UTF8_not_executed"
            except UnicodeDecodeError:
                status = "unreadable_binary_not_executed"
            members.append(
                {
                    "path": name,
                    "bytes": len(data),
                    "sha256": sha_bytes(data),
                    "text_extraction_status": status,
                    "source_family": family_of(name),
                    "license_status": "not_established_no_redistribution",
                }
            )
            if "/param/" in name and name.endswith(".data") and name in texts:
                topologies[str(Path(name).parent)] = data_topology(text)
        for member in members:
            name = member["path"]
            if name not in texts or "/param/" not in name:
                continue
            folder = str(Path(name).parent)
            context = {}
            for other, text in texts.items():
                if str(Path(other).parent) != folder or not other.endswith(".init"):
                    continue
                for number, line in enumerate(text.splitlines(), 1):
                    words = line.split("#", 1)[0].split()
                    if words and words[0] in {
                        "units",
                        "bond_style",
                        "angle_style",
                        "dihedral_style",
                        "improper_style",
                        "pair_style",
                        "pair_modify",
                        "special_bonds",
                    }:
                        context.setdefault(words[0], []).append(
                            {"file": other, "line": number, "tokens": words[1:]}
                        )
            if name.endswith((".settings", ".setting")):
                member_source = {
                    **source,
                    "relative_path": source["relative_path"] + "::" + name,
                    "sha256": member["sha256"],
                    "container_sha256": source["sha256"],
                    "source_family": family_of(name),
                    "version": revision or "archive_version_unspecified",
                }
                values, counts = parse_commands(
                    texts[name], member_source, topologies.get(folder, {}), context
                )
                records.extend(values)
                member["coefficient_commands_parsed"] = counts
                member["instantiated_parameter_records"] = len(values)
        return records, {
            "members": members,
            "archive_comment_revision": revision,
            "topology_facts": topologies,
            "license_file_present": any(
                re.search(r"(?:^|/)(?:license|copying)(?:\.|$)", n, re.IGNORECASE)
                for n in names
            ),
            "execution": "all archive members read as data only; no extraction to disk or command execution",
        }


def parse_pdf(raw):
    # No PDF library is imported by replay, tests or the ordinary package API.
    import pypdf

    try:
        return pdf_facts(pypdf, raw)
    except pypdf.errors.PyPdfError as error:
        raise ValueError("unreadable PDF: " + str(error)) from error


def pdf_facts(pypdf, raw):
    reader = pypdf.PdfReader(io.BytesIO(raw), strict=False)
    pages, dois = [], set()
    for number, page in enumerate(reader.pages, 1):
        text = page.extract_text() or ""
        pages.append(
            {
                "page": number,
                "extracted_characters": len(text),
                "text_sha256": sha_bytes(text.encode()),
                "status": "text_extracted"
                if text.strip()
                else "no_extractable_text_no_OCR",
            }
        )
        dois.update(re.findall(r"10\.\d{4,9}/[A-Za-z0-9._;()/:-]+", text))
    return {
        "parser": "pypdf " + pypdf.__version__,
        "pages": pages,
        "metadata": {str(k): str(v) for k, v in (reader.metadata or {}).items()},
        "DOI_candidates_not_authority": sorted(dois)[:30],
        "status": "all_pages_text_parsed"
        if pages and all(p["status"] == "text_extracted" for p in pages)
        else "partial_text_extraction_no_OCR",
        "reading_claim": "text parsing only; no new literature adjudication or manual full-reading claim",
    }


def read_parents(repo):
    return {n: read(repo / n) for n in PARENTS}


def capture(repo, root, directory):
    scope()
    require(
        not (root / "source-file-inventory.json").exists(), "do_not_overwrite_capture"
    )
    directory = directory.resolve()
    require(
        not root.resolve().is_relative_to(directory),
        "evidence_output_inside_readonly_source",
    )
    p = read_parents(repo)
    raw_inventory = inventory(directory)
    root.mkdir(parents=True, exist_ok=True)
    # Freeze original file identities before opening any format parser.
    write(
        root / "source-inputs.json",
        signed(
            {
                "schema": "island.j74-source-input-freeze.v1",
                "declared_source_directory": DECLARED,
                "declared_directory_exists": Path(DECLARED).is_dir(),
                "actual_readonly_directory": str(directory),
                "files": raw_inventory,
            }
        ),
    )
    historical = {r["sha256"]: r for r in p[PARENTS[5]]["files"]}
    old_records = p[PARENTS[4]]["records"]
    rows, files, detail = [], [], {}
    seen_content = {}
    for entry in raw_inventory:
        name, sha = entry["relative_path"], entry["sha256"]
        old = historical.get(sha)
        source = {
            **entry,
            "source_family": family_of(name, old),
            "citation": {
                k: old.get(k) for k in ("title", "authors", "journal_year", "doi")
            }
            if old
            else {"status": "local_file_only_not_author_authenticated"},
            "version": "not_established",
            "license_status": "not_established_no_source_assets_redistributed",
            "redistribution_status": "not_granted_in_this_intake; source asset excluded from git",
            "exact_derived_copied_unknown": "copied_origin_unverified",
            "source_role": "contextual_literature"
            if name.endswith(".pdf")
            else "external_implementation"
            if name.endswith(".zip")
            else "native_named_force_field_file",
            "source_authority": "local_copy_not_author_primary",
            "primary_artifact_supplied": False,
            "native_source_identity_complete": False,
            "parser_used": "not_parsed",
            "text_extraction_status": "not_parsed",
            "duplicate_content_of": seen_content.get(sha),
        }
        raw = (directory / name).read_bytes()
        require(sha_bytes(raw) == sha, "source_hash_changed_during_capture")
        new = []
        try:
            if name.endswith(".pdf"):
                meta = parse_pdf(raw)
                source["parser_used"] = meta["parser"]
                source["text_extraction_status"] = meta["status"]
                detail[name] = meta
                if old and sha not in seen_content:
                    for original in old_records:
                        if original["source_file_sha256"] != sha:
                            continue
                        new.append(
                            make_record(
                                source,
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
                source["parameter_extraction_status"] = (
                    "inherited_exact_hash_J71_facts_only"
                    if old
                    else "no_reviewed_numeric_extraction; no_automatic_table_inference"
                )
            elif name.endswith(".frc"):
                text = raw.decode("utf-8")
                source["version"] = version_lines(text)
                source["citation"] = {**source["citation"], **frc_provenance(text)}
                source["parser_used"] = "J74 literal FRC parser v1"
                source["text_extraction_status"] = "decoded_full_UTF8"
                if sha == PCFF_SHA and name == "pcff.frc":
                    citation = next(
                        r["citation"]
                        for r in p[PARENTS[3]]["rows"].values()
                        if r.get("source_file") == "pcff.frc"
                    )
                    source["citation"] = {**citation, **frc_provenance(text)}
                    source["native_source_identity_complete"] = True
                    source["source_authority"] = (
                        "byte_exact_copy_of_retained_LAMMPS_PCFF_release; not_new_author_primary"
                    )
                    source["exact_derived_copied_unknown"] = (
                        "exact_copy_matches_retained_upstream_hash"
                    )
                new, detail[name] = parse_frc(text, source)
            elif name.endswith(".zip"):
                source["parser_used"] = (
                    "Python zipfile + J74 literal LAMMPS/text parser v1"
                )
                source["text_extraction_status"] = "all_members_read_as_data"
                new, detail[name] = parse_archive(raw, source)
                source["version"] = detail[name]["archive_comment_revision"]
            else:
                text = raw.decode("utf-8")
                source["parser_used"] = "Python UTF-8 text reader"
                source["text_extraction_status"] = "decoded_full_UTF8"
                detail[name] = {
                    "lines": len(text.splitlines()),
                    "text_sha256": sha_bytes(text.encode()),
                    "templates_not_executed": True,
                }
        except (
            OSError,
            ValueError,
            ImportError,
            EOFError,
            zipfile.BadZipFile,
            RecursionError,
        ) as error:
            # Unreadable files are explicit retained blockers, never zero rows.
            source["text_extraction_status"] = "parse_failed"
            source["parser_used"] = (
                source["parser_used"]
                if source["parser_used"] != "not_parsed"
                else "attempted_format_parser"
            )
            detail[name] = {
                "error_type": type(error).__name__,
                "error": str(error)[:400],
                "status": "unreadable_no_OCR_or_substitution",
            }
        rows.extend(new)
        files.append(signed({k: v for k, v in source.items() if k != "identity"}))
        seen_content.setdefault(sha, name)
        print(
            json.dumps(
                {
                    "file": name,
                    "status": source["text_extraction_status"],
                    "records": len(new),
                }
            ),
            flush=True,
        )
    exact(inventory(directory), raw_inventory, "source_tree_changed_during_capture")
    write(
        root / "source-file-inventory.json",
        signed(
            {
                "schema": "island.j74-source-inventory.v1",
                "declared_directory": DECLARED,
                "actual_directory": str(directory),
                "files": files,
                "source_assets_committed": False,
            }
        ),
    )
    write(
        root / "source-parsing-facts.json",
        signed(
            {
                "schema": "island.j74-parsing-facts.v1",
                "files": detail,
                "scientific_evaluations": 0,
            }
        ),
    )
    write(
        root / "extracted-parameter-records.json",
        signed(
            {
                "schema": "island.j74-extracted-records.v1",
                "records": rows,
                "normalization_or_parameter_generation_performed": False,
            }
        ),
    )
    tree = subprocess.run(
        ["git", "ls-tree", "-r", "--full-tree", BASE],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    history = {}
    for line in tree.splitlines():
        metadata, name = line.split("\t", 1)
        mode = metadata.split()[0]
        path = repo / name
        history[name] = (
            sha_bytes(str(path.readlink()).encode()) if mode == "120000" else pin(path)
        )
    write(root / "historical-files.json", history)
    frozen = signed(
        {
            "schema": "island.j74-frozen-inputs.v1",
            "base_commit": BASE,
            "previous_J73_commit": J73,
            "fetched_origin_main": MAIN,
            "branch": BRANCH,
            "source_facts": {
                name: pin(root / name)
                for name in (
                    "source-inputs.json",
                    "source-file-inventory.json",
                    "source-parsing-facts.json",
                    "extracted-parameter-records.json",
                    "historical-files.json",
                )
            },
            "parent_files": {name: pin(repo / name) for name in PARENTS},
            "scope": SCOPE,
            "initial_readonly_inspection_preceded_J74_fence_preflight": True,
            "inspection_deviation": "Initial filesystem listing, FRC headers and ZIP member listing were read-only but preceded explicit J74 preflight. Subsequent capture and all terminal verification use C1. No scientific executable or evaluator ran.",
            "bounded_capture": {
                "one_worker": True,
                "max_file_bytes": MAX_FILE,
                "max_archive_member_bytes": MAX_MEMBER,
                "max_expanded_archive_bytes": MAX_FILE,
                "scientific_solver_calls": 0,
            },
        }
    )
    write(root / "frozen-inputs.json", frozen)
    print(
        json.dumps(
            {
                "frozen_inputs_sha256": pin(root / "frozen-inputs.json"),
                "files": len(files),
                "records": len(rows),
                "historical_files": len(history),
            }
        ),
        flush=True,
    )


def inputs(repo, root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "resigned_frozen_source_contract",
    )
    frozen = read(root / "frozen-inputs.json")
    verify_signed(frozen)
    exact(frozen["scope"], SCOPE, "frozen_scope")
    require(
        frozen["base_commit"] == BASE
        and frozen["previous_J73_commit"] == J73
        and frozen["fetched_origin_main"] == MAIN,
        "base_ancestry",
    )
    require(set(frozen["parent_files"]) == set(PARENTS), "explicit_parent_inventory")
    for name, expected in frozen["parent_files"].items():
        require(pin(repo / name) == expected, "historical_source_hash_mutation:" + name)
    for name, expected in frozen["source_facts"].items():
        require(
            pin(root / name) == expected, "fabricated_or_resigned_source_facts:" + name
        )
    return (
        frozen,
        read(root / "source-file-inventory.json"),
        read(root / "extracted-parameter-records.json"),
        read(root / "source-parsing-facts.json"),
        read_parents(repo),
    )


def validate_inventory(document):
    verify_signed(document)
    paths = set()
    for source in document["files"]:
        verify_signed(source)
        require(
            {
                "relative_path",
                "sha256",
                "license_status",
                "redistribution_status",
                "primary_artifact_supplied",
                "native_source_identity_complete",
            }
            <= set(source),
            "source_inventory_schema",
        )
        require(source["relative_path"] not in paths, "duplicate_source_identity")
        paths.add(source["relative_path"])
        require(
            bool(re.fullmatch(r"[0-9a-f]{64}", source["sha256"])),
            "source_hash_required",
        )
        require(
            source["license_status"] and source["redistribution_status"],
            "license_metadata_required",
        )
        require(
            source["primary_artifact_supplied"] is False,
            "local_copy_is_not_author_primary",
        )
        if source["native_source_identity_complete"]:
            require(
                source["relative_path"] == "pcff.frc"
                and source["sha256"] == PCFF_SHA
                and source["citation"]["revision"]
                == "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
                "native_PCFF_authority_identity",
            )


def compatible(record, target):
    """No cross-family authority, type analogy, or inferred convention closure."""
    family = record["source_family"]
    types = record["ordered_atom_types"]
    type_match = types in target["symmetry_equivalent_types"]
    family_match = family == "native_PCFF_CFF91_distribution_candidate"
    contract = target["required_evidence"]["equation_and_units"][
        "historical_unresolved_contract"
    ]
    unit_match = (
        record["units"].get("status") == "complete_source_contract"
        and record["units"].get("symbols") == contract["unit_symbols"]
    )
    form_match = (
        record["analytic_equation"].get("status") == "complete_source_contract"
        and record["analytic_equation"].get("scalar") == contract["scalar_energy"]
    )
    convention_complete = all(
        record[k].get("status") == "complete_source_contract"
        for k in (
            "phase_or_multiplicity",
            "orientation_order_convention",
            "equilibrium_dependencies",
        )
    )
    return {
        "exact_type_match": type_match,
        "exact_graph_match": record["exact_graph_applicability"] is True,
        "function_form_match": record["interaction_family"] == target["family"]
        and family_match
        and form_match,
        "unit_match": unit_match and family_match,
        "convention_match": convention_complete and family_match,
        "equilibrium_dependency_match": convention_complete and family_match,
        "source_family_match": family_match,
        "eligibility": False,
        "reason": "no complete native PCFF contract; external models, missing types and unresolved conventions cannot resolve these requests",
    }


def no_profile(document):
    exact(document["scope"], SCOPE, "profile_or_scientific_scope_activation")


def reconstruct(repo, root):
    _frozen, sources, facts, parsing, p = inputs(repo, root)
    validate_inventory(sources)
    verify_signed(facts)
    rows = facts["records"]
    require(
        len(rows) == len({r["record_id"] for r in rows}),
        "duplicate_extracted_record_identity",
    )
    for row in rows:
        verify_signed(row)
        require(
            row["native_PCFF_eligible"] is False
            and row["used_to_fill_missing_row"] is False,
            "analog_or_external_promotion",
        )
    native_alternatives = defaultdict(list)
    for row in rows:
        if row["source_family"] == "native_PCFF_CFF91_distribution_candidate":
            native_alternatives[
                (row["interaction_family"], tuple(row["ordered_atom_types"]))
            ].append(row)
    distinct_release_alternatives = [
        {
            "interaction_family": family,
            "ordered_types": list(types),
            "records": [
                {
                    "record_id": r["record_id"],
                    "source_sha256": r["source_sha256"],
                    "source_path": r["source_path"],
                    "location": r["location"],
                    "raw_values": r["coefficients_raw"],
                }
                for r in members
            ],
            "status": "distinct_release_or_source_rows_preserved_without_authority_selection",
            "numeric_comparison_performed": False,
        }
        for (family, types), members in sorted(native_alternatives.items())
        if len({r["source_sha256"] for r in members}) > 1
    ]
    requests = p[PARENTS[0]]["requests"]
    require(
        len(requests) == len({r["id"] for r in requests}) == 378,
        "incomplete_378_parent_coverage",
    )
    by_old = {
        r["historical_J71_record_id"]: r
        for r in rows
        if "historical_J71_record_id" in r
    }
    native = [
        r
        for r in rows
        if r["source_family"] == "native_PCFF_CFF91_distribution_candidate"
    ]
    resolutions = []
    for target in requests:
        candidates = {}
        native_queries = []
        for source in sources["files"]:
            if source["source_family"] != "native_PCFF_CFF91_distribution_candidate":
                continue
            source_facts = parsing["files"][source["relative_path"]]
            equivalences = source_facts["ordinary_equivalences"]
            role = "angle" if target["family"] == "angle-angle" else "torsion"
            mapped = [
                [equivalences.get(t, {}).get(role, t) for t in types]
                for types in target["symmetry_equivalent_types"]
            ]
            matches = [
                r
                for r in native
                if r["source_path"] == source["relative_path"]
                and r["interaction_family"] == target["family"]
                and r["ordered_atom_types"] in mapped
                and not r.get("wildcard")
                and not r.get("auto_equivalence")
            ]
            for r in matches:
                candidates[r["record_id"]] = r
            native_queries.append(
                {
                    "source_path": source["relative_path"],
                    "sha256": source["sha256"],
                    "native_source_identity_complete": source[
                        "native_source_identity_complete"
                    ],
                    "mapped_ordinary_equivalence_types": mapped,
                    "equivalence_role": role,
                    "matched_records": [r["record_id"] for r in matches],
                    "auto_equivalence_or_wildcard_used": False,
                }
            )
        for record in native:
            if (
                record["interaction_family"] == target["family"]
                and record["ordered_atom_types"] in target["symmetry_equivalent_types"]
            ):
                candidates[record["record_id"]] = record
        for alternative in target["unresolved_alternatives"]:
            if alternative["record_id"] in by_old:
                record = by_old[alternative["record_id"]]
                candidates[record["record_id"]] = record
        state = STATUS_MAP[target["J72_status"]]
        resolutions.append(
            signed(
                {
                    "id": target["id"],
                    "J70_identity": target["J70_identity"],
                    "J72_identity": target["J72_identity"],
                    "J73_identity": target["identity"],
                    "ordered_atom_types": target["ordered_types"],
                    "site_ids": target["stable_sites"],
                    "local_graph_identity": p[PARENTS[0]]["graph_identity"],
                    "chemical_environment": target["chemical_environment"],
                    "interaction_family": target["family"],
                    "candidate_source_records": [
                        {
                            "record_id": r["record_id"],
                            "source_family_classification": r["source_family"],
                            "evidence_reference": {
                                "path": r["source_path"],
                                "sha256": r["source_sha256"],
                                "location": r["location"],
                            },
                            **compatible(r, target),
                        }
                        for r in sorted(
                            candidates.values(), key=lambda x: x["record_id"]
                        )
                    ],
                    "native_exact_match_count": sum(
                        1
                        for r in candidates.values()
                        if r["source_family"]
                        == "native_PCFF_CFF91_distribution_candidate"
                    ),
                    "native_source_queries": native_queries,
                    "external_implementation_family_catalog_reference": "family-coverage-report.json#/external_implementation_candidate_catalog/"
                    + target["family"],
                    "native_source_lookup": target["existing_native_source_lookup"],
                    "resolution_status": state,
                    "statuses": {s: s == state for s in STATES},
                    "author_confirmed_reported": "28 historical user attestations, no invented per-row mapping",
                    "author_confirmed_audited": False,
                    "blocking_reason": "Exact native source row absent; no external/analog substitution"
                    if state != "convention_blocked"
                    else "Exact native Phi/AAT scalar, units, phase and dependency authority remain unresolved; local copies do not adjudicate this",
                    "required_evidence": target["required_evidence"],
                    "historical_alternative_record_ids": [
                        c["record_id"] for c in target["unresolved_alternatives"]
                    ],
                    "new_primary_artifact": None,
                }
            )
        )
    # No same-contract missing native row was supplied by these frozen assets.
    # Refuse an unexpected discovery rather than silently promote it.
    require(
        all(r["native_exact_match_count"] == 0 for r in resolutions),
        "unexpected_native_match_requires_separate_review",
    )
    counts = dict(Counter(r["resolution_status"] for r in resolutions))
    exact(
        counts,
        {"analog_only": 294, "convention_blocked": 60, "source_missing": 24},
        "unchanged_378_statuses",
    )
    ids = [r["id"] for r in resolutions]
    missing_families = dict(Counter(r["interaction_family"] for r in resolutions))
    blockers = [
        "No complete authoritative same-contract evidence for all 378 native PCFF instances",
        "Native Phi/AAT and cross-term dependencies remain unresolved",
        "Local copies/benchmark implementations do not establish author-primary authority",
        "Redistribution rights are unverified; no source assets included",
    ]
    if not read(root / "source-inputs.json")["declared_directory_exists"]:
        blockers.insert(
            0,
            "Declared /island_ffref/pcff absent; discovered checkout copy inventoried explicitly without aliasing",
        )
    unreadable = [
        s["relative_path"]
        for s in sources["files"]
        if s["text_extraction_status"] == "parse_failed"
    ]
    if unreadable:
        blockers.append(
            "Unreadable source files remain explicit blockers: " + ", ".join(unreadable)
        )
    archive_members = [
        {
            "archive_path": source["relative_path"],
            "archive_sha256": source["sha256"],
            "relative_path": member["path"],
            "sha256": member["sha256"],
            "bytes": member["bytes"],
            "MIME_type": mimetypes.guess_type(member["path"])[0] or "text/plain",
            "parser_used": "Python zipfile and literal text/topology parser v1; no execution",
            "text_extraction_status": member["text_extraction_status"],
            "source_family": member["source_family"],
            "citation": source["citation"],
            "version": parsing["files"][source["relative_path"]][
                "archive_comment_revision"
            ],
            "license_status": member["license_status"],
            "exact_derived_copied_unknown": "local_archive_copy_origin_unverified",
            "source_role": "external_implementation_or_contextual_input",
            "primary_artifact_supplied": False,
        }
        for source in sources["files"]
        if source["relative_path"].endswith(".zip")
        for member in parsing["files"][source["relative_path"]]["members"]
    ]
    ledger = signed(
        {
            "schema": "island.j74-artifact-ledger.v1",
            "status": "blocked",
            "blockers": blockers,
            "source_files": [
                {
                    k: s[k]
                    for k in (
                        "relative_path",
                        "sha256",
                        "source_family",
                        "citation",
                        "version",
                        "source_authority",
                        "primary_artifact_supplied",
                        "native_source_identity_complete",
                        "license_status",
                    )
                }
                for s in sources["files"]
            ],
            "author_confirmed_reported": 28,
            "author_confirmed_audited": 0,
            "new_trusted_author_primary_artifacts": 0,
            "local_copy_is_not_primary_authority": True,
            "archive_member_inventory": archive_members,
            "scope": SCOPE,
        }
    )
    reports = []
    selectors = {
        "native_PCFF": lambda f: f == "native_PCFF_CFF91_distribution_candidate",
        "native_COMPASS": lambda f: f == "native_COMPASS_distribution_candidate",
        "external_Huang_PDMS": lambda f: f == "external_Huang_PDMS" or "Huang" in f,
        "external_Xiang_PDMS_benchmark": lambda f: f.startswith("external_Xiang_"),
        "other_external_or_analog": lambda f: (
            f
            not in {
                "native_PCFF_CFF91_distribution_candidate",
                "native_COMPASS_distribution_candidate",
                "external_Huang_PDMS",
            }
            and not f.startswith("external_Xiang_")
        ),
    }
    for name, selector in selectors.items():
        subset = [r for r in rows if selector(r["source_family"])]
        reports.append(
            {
                "family": name,
                "extracted_facts": len(subset),
                "exact_rows_supplied": 0,
                "exact_rows_unresolved": 378,
                "unresolved_row_ids": ids,
                "count_scope": "exact J73 native PCFF contract, not the completeness of a different external model",
                "convention_fields_unresolved": [
                    "native equation/units",
                    "phase and multiplicity",
                    "ordering",
                    "equilibrium dependencies",
                    "charge origin",
                    "exclusion and LJ contract",
                    "exact graph applicability",
                ],
                "conflicting_rows": "all original alternatives preserved; see source-conflicts.json",
                "unsupported_graph_environments": sorted(
                    {" ".join(r["ordered_atom_types"]) for r in resolutions}
                ),
                "missing_Class_II_cross_terms": missing_families,
                "automatic_assignment": False,
            }
        )
    return {
        "source-artifact-ledger.json": ledger,
        "row-resolution-matrix.json": signed(
            {
                "schema": "island.j74-row-resolution.v1",
                "instances": resolutions,
                "counts": counts,
                "scope": SCOPE,
            }
        ),
        "convention-resolution.json": signed(
            {
                "schema": "island.j74-conventions.v1",
                "resolved_exact_rows": 0,
                "convention_resolved_rows": 0,
                "executable_bundle_qualified": False,
                "instances": [
                    {
                        "id": r["id"],
                        "J73_identity": r["J73_identity"],
                        "convention_resolved": False,
                        "missing_fields": r["required_evidence"][
                            "J72_remaining_missing_conventions"
                        ],
                        "scope": "request only; no cross-family convention inferred",
                    }
                    for r in resolutions
                ],
                "FRC_literal_header_evidence": {
                    n: d.get("section_headers", {})
                    for n, d in parsing["files"].items()
                    if n.endswith(".frc")
                },
                "scope": SCOPE,
            }
        ),
        "source-conflicts.json": signed(
            {
                "schema": "island.j74-conflicts.v1",
                "historical_J72_matrix_sha256": pin(repo / PARENTS[6]),
                "historical_equal_duplicate_equivalences": 16,
                "historical_unresolved_unequal_COMPASS_assignments": 2,
                "new_conflicts_adjudicated": 0,
                "distinct_native_release_alternatives": distinct_release_alternatives,
                "historical_alternatives_preserved": True,
                "resolution_by_numeric_average_or_residual": False,
                "separate_source_families": sorted({r["source_family"] for r in rows}),
                "J73_conflict_references": [
                    {
                        "id": t["id"],
                        "J73_identity": t["identity"],
                        "categories": t["inherited_conflict_categories"],
                        "alternative_record_ids": [
                            c["record_id"] for c in t["unresolved_alternatives"]
                        ],
                    }
                    for t in requests
                ],
                "scope": SCOPE,
            }
        ),
        "family-coverage-report.json": signed(
            {
                "schema": "island.j74-family-coverage.v1",
                "reports": reports,
                "interaction_family_extraction_inventory": {
                    family: {
                        "record_count": sum(
                            r["interaction_family"] == family for r in rows
                        ),
                        "absence_means": "no extracted exact source record, never coefficient zero",
                    }
                    for family in sorted({f[0] for f in FAMILIES.values()})
                },
                "external_implementation_candidate_catalog": {
                    family: [
                        {
                            "record_id": r["record_id"],
                            "source_path": r["source_path"],
                            "ordered_types": r["ordered_atom_types"],
                            "source_family": r["source_family"],
                            "exact_native_PCFF_match": False,
                            "exact_HMDS_graph_match": False,
                            "blocking_reason": "different source/model contract; no numeric atom-type or element analogy establishes a native PCFF match",
                        }
                        for r in rows
                        if r["interaction_family"] == family
                        and "::" in r["source_path"]
                    ]
                    for family in sorted(missing_families)
                },
                "scope": SCOPE,
            }
        ),
        "license-and-redistribution-report.json": signed(
            {
                "schema": "island.j74-license-report.v1",
                "source_assets_committed": False,
                "files": [
                    {
                        "path": s["relative_path"],
                        "sha256": s["sha256"],
                        "license_status": s["license_status"],
                        "redistribution_status": s["redistribution_status"],
                    }
                    for s in sources["files"]
                ],
                "policy": "Only hashes, citation metadata and extracted facts retained. No PDF/FRC/ZIP/source asset redistributed; no license inferred from local presence or upstream project license.",
                "scope": SCOPE,
            }
        ),
    }


def validate_resolution(document, expected):
    verify_signed(document)
    no_profile(document)
    require(
        len(document["instances"])
        == len({r["id"] for r in document["instances"]})
        == 378,
        "incomplete_or_duplicate_378_coverage",
    )
    for row in document["instances"]:
        verify_signed(row)
        require(
            row["resolution_status"] in STATES and set(row["statuses"]) == set(STATES),
            "separate_status_schema",
        )
        require(
            row["statuses"]["executable_bundle_qualified"] is False,
            "premature_profile_qualification",
        )
    exact(document, expected, "resigned_row_or_convention_promotion")


def replay(repo, root, create=False):
    documents = reconstruct(repo, root)
    for name, expected in documents.items():
        if create:
            require(not (root / name).exists(), "do_not_overwrite_resolution_attempt")
            write(root / name, expected)
        actual = read(root / name)
        if name == "row-resolution-matrix.json":
            validate_resolution(actual, expected)
        else:
            exact(actual, expected, "resigned_derived_evidence:" + name)
    return {
        "validation_status": "passed",
        "source_resolution_status": "blocked",
        "instances": 378,
        "resolved_rows": 0,
        "counts": documents["row-resolution-matrix.json"]["counts"],
        "blockers": documents["source-artifact-ledger.json"]["blockers"],
        "scope": SCOPE,
    }


def originals(root, directory):
    scope()
    expected = read(root / "source-inputs.json")["files"]
    exact(
        inventory(directory),
        expected,
        "local_source_hash_or_recursive_inventory_mutation",
    )
    return {
        "source_hash_verification": "passed",
        "files": len(expected),
        "source_files_modified": False,
    }


def main():
    try:
        scope()
    except RuntimeError as error:
        print(json.dumps({"status": "blocked", "reason": str(error), "scope": SCOPE}))
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--source", type=Path, default=Path(DECLARED))
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--capture", action="store_true")
    mode.add_argument("--build", action="store_true")
    mode.add_argument("--verify-originals", action="store_true")
    args = parser.parse_args()
    try:
        result = (
            capture(args.repo.resolve(), args.directory.resolve(), args.source)
            if args.capture
            else originals(args.directory, args.source)
            if args.verify_originals
            else replay(
                args.repo.resolve(), args.directory.resolve(), create=args.build
            )
        )
        if result is not None:
            print(json.dumps(result, sort_keys=True))
    except (IntakeError, OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps({"status": "blocked", "reason": str(error), "scope": SCOPE}))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
