"""Explicit developer-only author-source fetch and offline convention audit."""

import argparse
import base64
import json
import urllib.request
from hashlib import sha1, sha256
from pathlib import Path

from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_author import (
    COMMIT,
    PINS,
    audit_author_opls,
    load_author_opls_audit,
)
from island.forcefields.oplsaa.pdms_literature import (
    build_opls_pdms_literature,
)
from island.forcefields.pdms_literature import INPUT_ERRORS


def fetch(directory):
    # An explicit fetch command only. This function is never called by ISLAND preparation.
    root = Path(directory)
    if root.exists():
        raise ValueError("fetch_destination_already_exists")
    headers = {
        "User-Agent": "ISLAND-author-source-audit",
        "Accept": "application/vnd.github+json",
    }
    base = "https://api.github.com/repos/COMO-lab/PDMS_benchmark"
    tree_sha = "5b3bcc406226a9f5af9fcaee4abf5a465e9abc26"
    allowed = (
        base + "/",
        "https://api.github.com/repos/WuGroup-XJTLU/PDMS_benchmark/",
        "https://api.github.com/repositories/902160948/",
    )

    def read(url):
        with urllib.request.urlopen(
            urllib.request.Request(url, headers=headers), timeout=45
        ) as response:
            if not response.url.startswith(allowed):
                raise ValueError("unexpected_source_redirect")
            return json.loads(response.read(25000000)), response.url

    commit, commit_url = read(base + "/git/commits/" + COMMIT)
    if commit["sha"] != COMMIT or commit["tree"]["sha"] != tree_sha:
        raise ValueError("wrong_author_commit_tree")
    tree, tree_url = read(base + "/git/trees/" + tree_sha + "?recursive=1")
    if tree["truncated"] or tree["sha"] != tree_sha:
        raise ValueError("wrong_or_incomplete_author_tree")
    entries = {e["path"]: e for e in tree["tree"]}
    files, provenance = {}, []
    for name, expected in PINS.items():
        entry = entries[name]
        url = base + "/git/blobs/" + entry["sha"]
        blob, final_url = read(url)
        raw = base64.b64decode(blob["content"], validate=False)
        if len(raw) != entry["size"] or sha256(raw).hexdigest() != expected:
            raise ValueError("downloaded_author_source_hash_mismatch")
        if (
            sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
            != entry["sha"]
        ):
            raise ValueError("downloaded_git_blob_mismatch")
        files[name] = raw
        provenance.append(
            {
                "path": name,
                "sha256": expected,
                "size": len(raw),
                "git_blob_sha1": entry["sha"],
                "requested_url": url,
                "final_url": final_url,
            }
        )
    root.mkdir(parents=True)
    for name, raw in files.items():
        dest = root / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
    result = {
        "commit": COMMIT,
        "commit_api_url": commit_url,
        "tree": tree_sha,
        "tree_api_url": tree_url,
        "files": provenance,
    }
    (root / "download-provenance.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch-author-source-to", type=Path)
    parser.add_argument("--sources", type=Path)
    parser.add_argument("--engine-docs", type=Path)
    parser.add_argument("--references", type=Path)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reconstruct", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    if args.fetch_author_source_to:
        if any((args.sources, args.output, args.reconstruct, args.require_executable)):
            parser.error(
                "--fetch-author-source-to is a separate explicit download operation"
            )
    elif not all(
        (args.sources, args.engine_docs, args.references, args.selection, args.catalog)
    ):
        parser.error(
            "offline audit requires --sources, --engine-docs, --references, --selection, --catalog"
        )
    elif bool(args.output) == bool(args.reconstruct) or (
        args.reconstruct and not args.sha256
    ):
        parser.error("use --output or --reconstruct with an independent --sha256")
    try:
        if args.fetch_author_source_to:
            print(json.dumps(fetch(args.fetch_author_source_to)))
            return 0
        supplement = build_opls_pdms_literature(
            DomainSelection(args.selection.read_text(), args.catalog.read_text())
        )
        options = {
            "sources": args.sources,
            "supplement": supplement,
            "references": args.references,
            "engine_docs": args.engine_docs,
        }
        if args.reconstruct:
            result = load_author_opls_audit(args.reconstruct, args.sha256, **options)
        else:
            result = audit_author_opls(**options)
        p = result.payload
        report = {
            k: p[k]
            for k in (
                "identity",
                "data_identity",
                "status",
                "stages",
                "missing_dependencies",
                "profile_identity",
                "execution_profiles",
                "backend_executions",
                "acceptance_gate",
                "limitations",
            )
        }
        report["record_sha256"] = sha256(result.json_text.encode()).hexdigest()
        if args.output:
            args.output.mkdir(parents=True, exist_ok=False)
            (args.output / "record.json").write_text(result.json_text)
            (args.output / "audit.json").write_text(
                json.dumps(p, ensure_ascii=False, indent=2) + "\n"
            )
            (args.output / "report.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
        print(json.dumps(report))
        return 1 if args.require_executable else 0
    except INPUT_ERRORS as error:
        print(
            json.dumps(
                {
                    "status": getattr(
                        error,
                        "status",
                        "asset_unavailable"
                        if isinstance(error, OSError)
                        else "source_incompatible",
                    ),
                    "diagnostics": getattr(
                        error, "diagnostics", [{"reason": str(error)}]
                    ),
                    "execution_profiles": [],
                    "acceptance_gate": 1,
                }
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
