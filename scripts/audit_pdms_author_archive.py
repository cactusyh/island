"""Audit the user-supplied author ZIP offline, without extracting or executing it.

This is a source-copy comparison, never an execution or licensing authority.
The expected tree is the previously audited, independently pinned API response.
"""

import argparse
import hashlib
import io
import json
import shlex
import zipfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath

TREE_SHA256 = "94af420dc0eb7148a3b97f8b459f66213f8167e20ddf1d828300731ea9316b16"
COMMIT = "0053fda8449797e7e2b0790426cfc52dcd023e3e"
TREE = "5b3bcc406226a9f5af9fcaee4abf5a465e9abc26"
ROOT = "PDMS_benchmark-main"


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def commands(raw):
    return [
        {"line": i, "raw_text": line, "tokens": tokens}
        for i, line in enumerate(raw.decode().splitlines(), 1)
        if (tokens := shlex.split(line, comments=True))
    ]


def data_summary(raw):
    """Read only the author's full-style atom records; preserve charge literals."""
    charges, masses, molecules = defaultdict(Counter), {}, set()
    lines = raw.decode().splitlines()
    header, counts, section = lines[0], {}, None
    for line in lines:
        tokens = line.split("#", 1)[0].split()
        if not tokens:
            continue
        if not tokens[0].isdigit():
            section = " ".join(tokens)
        elif section == "Masses":
            masses[tokens[0]] = tokens[1]
        elif section == "Atoms":
            require(len(tokens) in (7, 10), "unexpected_full_style_atom_record")
            charges[tokens[2]][tokens[3]] += 1
            molecules.add(tokens[1])
        elif len(tokens) == 2 and tokens[1] in (
            "atoms",
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
        ):
            counts[tokens[1]] = int(tokens[0])
    atom_count = sum(sum(v.values()) for v in charges.values())
    require(atom_count == counts["atoms"], "atom_inventory_count_mismatch")
    return {
        "header": {"line": 1, "raw_text": header},
        "interaction_counts": counts,
        "atoms_read": atom_count,
        "molecules_read": len(molecules),
        "used_type_masses": {k: masses[k] for k in sorted(charges)},
        "type_charge_literal_histograms": dict(charges),
    }


def audit(archive, expected_sha256, tree_path):
    raw = Path(archive).read_bytes()
    require(sha256(raw) == expected_sha256, "archive_sha256_mismatch")
    tree_raw = Path(tree_path).read_bytes()
    require(sha256(tree_raw) == TREE_SHA256, "historical_tree_sha256_mismatch")
    tree = json.loads(tree_raw)
    require(tree["sha"] == TREE and tree["truncated"] is False, "tree_identity")
    expected = {e["path"]: e for e in tree["tree"] if e["type"] == "blob"}
    inventory, data, inputs, names = [], {}, {}, set()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        require(z.comment.decode() == COMMIT, "unexpected_archive_commit_comment")
        entry_count = len(z.infolist())
        # read() verifies each member CRC. No file is extracted or executed.
        for member in z.infolist():
            if member.is_dir():
                continue
            path = PurePosixPath(member.filename)
            require(
                not path.is_absolute()
                and ".." not in path.parts
                and path.parts[0] == ROOT,
                "unexpected_archive_path",
            )
            name = str(PurePosixPath(*path.parts[1:]))
            require(name not in names, "duplicate_archive_path")
            names.add(name)
            content = z.read(member)
            git_blob = hashlib.sha1(
                b"blob " + str(len(content)).encode() + b"\0" + content
            ).hexdigest()
            match = name in expected and (
                expected[name]["sha"] == git_blob
                and expected[name]["size"] == len(content)
            )
            require(match, "archive_differs_from_pinned_tree:" + name)
            inventory.append(
                {
                    "path": name,
                    "size": len(content),
                    "sha256": sha256(content),
                    "git_blob_sha1": git_blob,
                    "crc32": f"{member.CRC:08x}",
                    "matches_pinned_tree": match,
                }
            )
            if name.endswith(".data"):
                data[name] = data_summary(content)
            elif "/" in name and not name.endswith(".charges"):
                inputs[name] = commands(content)
    require(names == set(expected), "archive_missing_expected_files")
    license_files = [
        name
        for name in sorted(names)
        if PurePosixPath(name)
        .name.upper()
        .startswith(("LICENSE", "LICENCE", "COPYING", "COPYRIGHT", "NOTICE"))
    ]
    model_templates = [name for name in sorted(names) if name.endswith(".lt")]
    build_code = [
        name
        for name in sorted(names)
        if PurePosixPath(name).suffix.lower() in (".py", ".c", ".cpp", ".sh")
    ]
    include_inventory = []
    for name, rows in inputs.items():
        for row in rows:
            if row["tokens"][0] not in ("include", "read_data"):
                continue
            target = str(PurePosixPath(name).parent / row["tokens"][1])
            include_inventory.append(
                {
                    "source_file": name,
                    **row,
                    "target_relative_to_script_directory": target,
                    "target_present_in_that_directory": target in names,
                    "scope": "static path audit; LAMMPS resolves from execution cwd",
                }
            )
    return {
        "schema": "island_user_pdms_author_archive_audit_v1",
        "paper_doi": "10.1021/acs.jpcb.4c08471",
        "repository_commit": COMMIT,
        "repository_tree": TREE,
        "tree_response_sha256": TREE_SHA256,
        "archive": {
            "sha256": expected_sha256,
            "size": len(raw),
            "comment": COMMIT,
            "entries": entry_count,
            "files": len(inventory),
            "uncompressed_bytes": sum(e["size"] for e in inventory),
            "all_member_crcs_verified": True,
        },
        "files": inventory,
        "all_files_identical_to_previously_audited_tree": True,
        "new_or_missing_files": [],
        "license_files": license_files,
        "moltemplate_templates": model_templates,
        "build_code_files": build_code,
        "separate_charge_files": sorted(n for n in names if n.endswith(".charges")),
        "opls_active_init": inputs["param/OPLS/opls_system.in.init"],
        "static_include_inventory": include_inventory,
        "data_files": data,
        "interpretation": {
            "source_copy_verified": True,
            "author_model_version_changed": False,
            "OPLS_Qeq_input_closure_newly_supplied": False,
            "new_repository_license_supplied": False,
            "license_scope": "No new grant; article CC-BY does not by itself establish repository redistribution rights. Local research use is a separate question.",
            "model": "Previously audited DFT/NPA-charge CHARMM-switched-LJ OPLS input branch; not confirmed paper Qeq main-result input closure.",
            "no_native_forcefield_changes": True,
            "execution_profile_created": False,
            "numerical_PDMS_execution": False,
        },
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.archive, args.sha256, args.tree)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(
        json.dumps(
            {
                "archive_verified": True,
                "files_identical": len(result["files"]),
                "execution_profile_created": False,
                "numerical_PDMS_execution": False,
            }
        )
    )


if __name__ == "__main__":
    main()
