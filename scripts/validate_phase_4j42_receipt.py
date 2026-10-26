"""Verify J42 blocked-boundary archives, exact reconstruction and historical bytes."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree

from scripts.phase_4j41_c2 import read
from scripts.phase_4j41_c5 import require
from scripts.phase_4j42 import BASE, BRANCH, reconstruct
from scripts.validate_phase_4j41_c3_receipt import extract
from scripts.validate_phase_4j41_receipt import sha, signed


def validate(directory, *, repo=None):
    directory = Path(directory).resolve()
    repo = Path.cwd() if repo is None else Path(repo)
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(receipt == signed(receipt), "Receipt checksum changed")
    require(manifest == signed(manifest), "Manifest checksum changed")
    require(receipt["schema"] == "island.j42-receipt.v1", "Wrong receipt schema")
    require(
        receipt["base_commit"] == BASE and receipt["branch"] == BRANCH,
        "Wrong base/branch",
    )
    require(
        receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Manifest changed",
    )
    require(
        receipt["acceptance_gate"] == receipt["native_pcff_acceptance_gate"] == 1
        and receipt["status"] == "asset_unavailable"
        and receipt["execution_profiles"] == []
        and receipt["backend_executions"] == []
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(
        set(manifest["archives"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Archive set changed",
    )
    for name, expected in manifest["archives"].items():
        require(sha(directory / name) == expected, "Archive changed: " + name)
    historical = set(
        subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", BASE, "docs"], cwd=repo, text=True
        ).splitlines()
    )
    require(
        historical == set(manifest["historical_files"]),
        "Historical inventory incomplete",
    )
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty preservation inventory")
        for name, expected in manifest[group].items():
            require(sha(repo / name) == expected, "Changed " + name)
    with tempfile.TemporaryDirectory(prefix="island-j42-receipt-") as temp:
        root = Path(temp)
        for archive in manifest["archives"]:
            extract(directory / archive, root)
        for folder, inventory in (
            ("study", "study_files"),
            ("validation", "validation_files"),
        ):
            files = {
                str(p.relative_to(root / folder)): sha(p)
                for p in (root / folder).rglob("*")
                if p.is_file()
            }
            require(files == manifest[inventory], "Archive membership/bytes changed")
        require(
            receipt["reconstruction"] == reconstruct(root / "study"),
            "Reconstruction differs",
        )
        report = read(root / "study/report.json")
        require(report == receipt["study_report"], "Study report changed")
        require(
            len(report["cases"]) == 5
            and all(
                c["status"] == "asset_unavailable"
                and c["required_interactions"] == c["missing_interactions"]
                for c in report["cases"]
            ),
            "Incomplete blocked fixture inventory",
        )
        validation = root / "validation"
        for name, summary in receipt["verification"].items():
            command = json.loads((validation / f"{name}-command.json").read_text())
            require(
                command == summary and command["exit_code"] == 0,
                "Failed check: " + name,
            )
            require(
                command["tested_implementation_commit"]
                == receipt["tested_implementation_commit"],
                "Different tested commit",
            )
        for name in ("focused", "full"):
            counts = {
                k: sum(
                    int(s.attrib.get(k, 0))
                    for s in ElementTree.parse(validation / f"{name}.xml").getroot()
                )
                for k in ("tests", "failures", "errors", "skipped")
            }
            require(
                counts == receipt["test_counts"][name]
                and counts["failures"] == counts["errors"] == 0,
                "Test counts changed",
            )
        require(
            json.loads((validation / "regressions.json").read_text())
            == receipt["regressions"],
            "Regression inventory changed",
        )
        execution = json.loads((validation / "execution-command.json").read_text())
        require(
            execution == receipt["full_source_execution"]
            and execution["exit_code"] == 1,
            "Blocked acceptance changed",
        )
        require(
            "Traceback" not in (validation / "execution.log").read_text(),
            "Unclassified acceptance failure",
        )
    return {
        "receipt_valid": True,
        "historical_files": len(historical),
        "acceptance_gate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    print(json.dumps(validate(parser.parse_args().directory)))


if __name__ == "__main__":
    main()
