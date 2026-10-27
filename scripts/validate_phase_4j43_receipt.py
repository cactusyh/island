"""Validate J43 blocked evidence, exact reconstruction and historical bytes."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree

from island.forcefields.external_pdms import LIMITS, require
from scripts.phase_4j41_c2 import read
from scripts.phase_4j43 import BASE, BRANCH, reconstruct
from scripts.validate_phase_4j41_c3_receipt import extract
from scripts.validate_phase_4j41_receipt import sha, signed


def validate(directory, *, repo=None):
    directory = Path(directory).resolve()
    repo = Path.cwd() if repo is None else Path(repo)
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(receipt == signed(receipt), "Receipt checksum changed")
    require(manifest == signed(manifest), "Manifest checksum changed")
    require(receipt["schema"] == "island.j43-receipt.v1", "Wrong receipt schema")
    require(manifest["schema"] == "island.j43-manifest.v1", "Wrong manifest schema")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "Wrong base or branch",
    )
    require(
        receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Changed manifest",
    )
    require(
        receipt["status"] == "asset_unavailable"
        and receipt["acceptance_gate"] == 1
        and receipt["execution_profiles"] == []
        and receipt["backend_executions"] == []
        and receipt["limitations"] == LIMITS
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    names = set(
        subprocess.check_output(
            [
                "git",
                "ls-tree",
                "-r",
                "--name-only",
                BASE,
                "docs",
                "src",
                "scripts",
                "tests",
            ],
            cwd=repo,
            text=True,
        ).splitlines()
    )
    require(names == set(manifest["historical_files"]), "Historical inventory changed")
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty source inventory")
        for name, expected in manifest[group].items():
            require(sha(repo / name) == expected, "Changed " + name)
    require(
        set(manifest["archives"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Archive set changed",
    )
    with tempfile.TemporaryDirectory(prefix="island-j43-receipt-") as temp:
        root = Path(temp)
        for name, expected in manifest["archives"].items():
            require(sha(directory / name) == expected, "Changed archive " + name)
            extract(directory / name, root)
        for folder, group in (
            ("study", "study_files"),
            ("validation", "validation_files"),
        ):
            files = {
                str(p.relative_to(root / folder)): sha(p)
                for p in (root / folder).rglob("*")
                if p.is_file()
            }
            require(files == manifest[group], "Archive inventory or bytes changed")
        require(
            receipt["reconstruction"] == reconstruct(root / "study"),
            "Changed reconstructed result",
        )
        require(
            receipt["study_report"] == read(root / "study/report.json"), "Changed study"
        )
        for name, command in receipt["verification"].items():
            stored = json.loads(
                (root / "validation" / (name + "-command.json")).read_text()
            )
            require(
                stored == command
                and command["tested_implementation_commit"]
                == receipt["tested_implementation_commit"],
                "Wrong tested command",
            )
            require(
                command["exit_code"] == (1 if name == "execution" else 0),
                "Failed verification " + name,
            )
        for name in ("focused", "full", "regressions"):
            counts = {
                k: sum(
                    int(s.attrib.get(k, 0))
                    for s in ElementTree.parse(
                        root / "validation" / (name + ".xml")
                    ).getroot()
                )
                for k in ("tests", "failures", "errors", "skipped")
            }
            require(
                counts == receipt["test_counts"][name]
                and counts["failures"] == counts["errors"] == 0,
                "Changed test counts",
            )
    return {
        "receipt_valid": True,
        "historical_files": len(names),
        "acceptance_gate": 1,
        "execution_profiles": [],
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    print(json.dumps(validate(parser.parse_args().directory)))


if __name__ == "__main__":
    main()
