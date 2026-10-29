"""Replay relocated PDMS evidence and verify historical bytes and test receipts."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree

from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_source_completion import file_hash
from scripts.audit_pdms_source_completion import audit
from scripts.complete_pdms_source_parameters import run
from scripts.validate_sun_siloxane_parameters import extracted_archive


def authenticated(path):
    data = json.loads(Path(path).read_text())
    assert data["checksum"] == digest(
        {k: v for k, v in data.items() if k != "checksum"}
    )
    return data


def validate(directory, distribution, references):
    directory = Path(directory)
    repo = Path(__file__).resolve().parents[1]
    manifest = authenticated(directory / "manifest.json")
    receipt = authenticated(directory / "receipt.json")
    assert receipt["manifest_sha256"] == file_hash(directory / "manifest.json")
    assert receipt["limitations"] == LIMITS
    assert receipt["execution_profiles"] == manifest["profile_identities"] == []
    assert receipt["acceptance_gate"] == receipt["native_pcff_acceptance_gate"] == 1
    assert receipt["backend_MD_executions"] == []
    tracked = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", manifest["base_commit"]],
        cwd=repo,
        text=True,
    ).splitlines()
    assert set(tracked) == set(manifest["historical_tracked_files"])
    for group in (
        "historical_tracked_files",
        "preserved_prior_supplement",
        "implementation_files",
    ):
        for path, pin in manifest[group].items():
            assert file_hash(repo / path) == pin, path
    for name, pin in manifest["evidence_files"].items():
        assert file_hash(directory / name) == pin, name
    for name, pin in manifest["source_hashes"].items():
        assert file_hash(Path(distribution) / name) == pin
    for name, pin in manifest["historical_local_sources"].items():
        assert file_hash(Path(references) / name) == pin
    with tempfile.TemporaryDirectory(prefix="island-pdms-source-v2-") as temp:
        root = Path(temp)
        extracted_archive(directory, "study.tar.xz", manifest["study_files"], root)
        report = run(
            distribution,
            repo / "docs/evidence/sun_siloxane_parameters/parameters.json",
            references,
            root,
            reconstruct=True,
        )
        assert report == receipt["study_report"]
        # The standalone ASCII audit is ordinary JSON: integer atom-map keys
        # are serialized as strings. Compare that declared representation.
        assert json.loads(json.dumps(audit(root, distribution))) == json.loads(
            (directory / "independent-source-audit.json").read_text()
        )
        logs = root / "validation"
        logs.mkdir()
        extracted_archive(
            directory, "validation.tar.xz", manifest["validation_files"], logs
        )
        for name, counts in receipt["test_counts"].items():
            actual = {
                k: sum(
                    int(s.attrib.get(k, 0))
                    for s in ElementTree.parse(logs / (name + ".xml")).getroot()
                )
                for k in ("tests", "failures", "errors", "skipped")
            }
            assert counts == actual and counts["failures"] == counts["errors"] == 0
        for name, expected in receipt["verification"].items():
            cmd = json.loads((logs / (name + "-command.json")).read_text())
            assert cmd == expected and cmd["exit_code"] == 0
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_tracked_files_unchanged": len(tracked),
        "prior_supplement_files_unchanged": len(manifest["preserved_prior_supplement"]),
        "fresh_graph_reconstructions": len(report["fixtures"]),
        "execution_profiles": [],
        "acceptance_gate": 1,
        **LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("directory", "distribution", "references"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            validate(args.directory, args.distribution, args.references), indent=2
        )
    )


if __name__ == "__main__":
    main()
