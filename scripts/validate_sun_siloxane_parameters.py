"""Validate Sun parameter evidence, reconstructed graphs and historical bytes."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree

from island.forcefields.external_pdms import LIMITS, digest, require
from scripts.assess_sun_siloxane_parameters import run


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def authenticated(path):
    value = json.loads(Path(path).read_text())
    require(
        value["checksum"]
        == digest({k: v for k, v in value.items() if k != "checksum"}),
        "evidence_checksum_changed",
    )
    return value


def extracted_archive(directory, name, hashes, destination):
    with tarfile.open(directory / name) as archive:
        members = archive.getmembers()
        require(
            len({m.name for m in members}) == len(members), "duplicate_archive_member"
        )
        require({m.name for m in members} == set(hashes), "archive_inventory_changed")
        for member in members:
            require(
                member.isfile() and Path(member.name).name == member.name,
                "unsafe_archive_member",
            )
            data = archive.extractfile(member).read()
            require(
                sha256(data).hexdigest() == hashes[member.name],
                "archive_member_changed",
            )
            (destination / member.name).write_bytes(data)


def validate(directory, *, references, repo=None):
    repo = Path(__file__).resolve().parents[1] if repo is None else Path(repo)
    directory = Path(directory)
    manifest = authenticated(directory / "manifest.json")
    receipt = authenticated(directory / "receipt.json")
    require(manifest["schema"] == "island.sun-siloxane-manifest.v1", "manifest_schema")
    require(receipt["schema"] == "island.sun-siloxane-receipt.v1", "receipt_schema")
    require(
        receipt["manifest_sha256"] == file_hash(directory / "manifest.json"),
        "manifest_changed",
    )
    require(
        receipt["execution_profiles"] == manifest["profile_identities"] == []
        and receipt["backend_executions"] == []
        and receipt["acceptance_gate"] == receipt["native_pcff_acceptance_gate"] == 1
        and receipt["limitations"] == LIMITS,
        "execution_or_readiness_promoted",
    )
    historical = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", manifest["historical_commit"]],
        cwd=repo,
        text=True,
    ).splitlines()
    require(
        set(historical) == set(manifest["historical_files"]),
        "historical_inventory_changed",
    )
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "empty_manifest_group")
        for name, expected in manifest[group].items():
            require(file_hash(repo / name) == expected, "changed:" + name)
    for name, expected in manifest["evidence_files"].items():
        require(file_hash(directory / name) == expected, "changed_evidence:" + name)
    for name, expected in manifest["historical_local_sources"].items():
        require(
            file_hash(Path(references) / name) == expected,
            "changed_historical_source:" + name,
        )
    with tempfile.TemporaryDirectory(prefix="island-sun-evidence-") as temp:
        root = Path(temp)
        extracted_archive(directory, "study.tar.xz", manifest["study_files"], root)
        report = run(
            directory / "parameters.json",
            manifest["catalog_sha256"],
            references,
            root,
            reconstruct=True,
        )
        require(receipt["study_report"] == report, "study_report_changed")
        logs = root / "validation"
        logs.mkdir()
        extracted_archive(
            directory, "validation.tar.xz", manifest["validation_files"], logs
        )
        for name in ("focused", "full"):
            counts = {
                k: sum(
                    int(s.attrib.get(k, 0))
                    for s in ElementTree.parse(logs / (name + ".xml")).getroot()
                )
                for k in ("tests", "failures", "errors", "skipped")
            }
            require(counts == receipt["test_counts"][name], "test_counts_changed")
            require(counts["errors"] == counts["failures"] == 0, "test_failures")
        for name, expected in receipt["verification"].items():
            command = json.loads((logs / (name + "-command.json")).read_text())
            require(
                command == expected and command["exit_code"] == 0,
                "failed_or_changed_verification:" + name,
            )
    require(
        receipt["source_reconstruction"]["parameter_rows_reconstructed"] == 105
        and receipt["source_reconstruction"]["partial_rows_not_zero_filled"] == 8,
        "source_counts_changed",
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files_unchanged": len(historical),
        "historical_local_sources_unchanged": len(manifest["historical_local_sources"]),
        "graph_reconstructions": len(report["fixtures"]),
        "acceptance_gate": 1,
        "execution_profiles": [],
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--references", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.directory, references=args.references), indent=2))


if __name__ == "__main__":
    main()
