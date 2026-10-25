"""Strict C14 report, manifest, source and historical byte validation."""

import argparse
import json
from pathlib import Path

from scripts.phase_4j40_c11_numeric import require, sha, signed
from scripts.phase_4j40_c14 import C13_RECEIPT, build_report


def validate_report(report, source):
    require(
        report == build_report(source),
        "Changed backend report, source bindings, residuals or readiness",
    )
    return True


def validate(directory, *, repo=None):
    repo = Path(repo) if repo is not None else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    require(
        receipt["schema"] == "island.j40-c14-receipt.v1" and receipt == signed(receipt),
        "Invalid C14 receipt",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        manifest["schema"] == "island.j40-c14-manifest.v1"
        and manifest == signed(manifest),
        "Invalid C14 manifest",
    )
    require(
        sha(directory / "manifest.json") == receipt["manifest_sha256"],
        "Changed manifest",
    )
    require(
        set(manifest["files"]) == {"comparison-report.json", "validation.tar.xz"},
        "Incomplete evidence manifest",
    )
    for name, digest in manifest["files"].items():
        require(sha(directory / name) == digest, "Changed C14 evidence: " + name)
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty preservation manifest")
        for name, digest in manifest[group].items():
            require(sha(repo / name) == digest, "Changed " + group + ": " + name)
    c13 = repo / "docs/evidence/phase_4j40_c13"
    require(
        manifest["historical_files"].get("docs/evidence/phase_4j40_c13/receipt.json")
        == C13_RECEIPT,
        "Missing pinned C13 preservation",
    )
    report = json.loads((directory / "comparison-report.json").read_text())
    validate_report(report, c13)
    require(receipt["report_identity"] == report["identity"], "Changed report identity")
    require(
        receipt["source_counts"]
        == {
            k: report[k]
            for k in ("native_evaluations", "comparison_count", "passing_comparisons")
        },
        "Changed C13 counts",
    )
    require(
        receipt["backend_states"]
        == {cell["cell"]: cell["decision"] for cell in report["cells"]},
        "Changed backend decisions",
    )
    require(
        receipt["scientific_acceptance_exit_code"] == 1
        and type(receipt["scientific_acceptance_exit_code"]) is int,
        "Promoted scientific acceptance",
    )
    require(
        receipt["historical_byte_hashes_preserved"]
        == len(manifest["historical_files"]),
        "Changed preservation count",
    )
    return {
        "integrity_exit_code": 0,
        "scientific_acceptance_exit_code": 1,
        "native_evaluations_preserved": 300,
        "comparisons_preserved": 108,
        "passing_comparisons_preserved": 13,
        "historical_hashes_verified": len(manifest["historical_files"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate(args.directory), indent=2))


if __name__ == "__main__":
    main()
