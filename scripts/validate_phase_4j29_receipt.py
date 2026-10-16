"""Validate the J29 evidence receipt against its reproducibility manifest."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from island.exceptions import AnalysisError

SHA1 = re.compile(r"^[0-9a-f]{40}$")
RECEIPT_SCHEMA = "island.phase-4j29-evidence.v1"


def validate_receipt(receipt_path: str | Path, manifest_path: str | Path) -> dict:
    receipt = json.loads(Path(receipt_path).read_text())
    manifest = json.loads(Path(manifest_path).read_text())
    if receipt.get("schema") != RECEIPT_SCHEMA:
        raise AnalysisError("J29 receipt schema mismatch")
    required = ("branch", "base_commit", "validated_commit")
    if any(name not in receipt for name in required):
        raise AnalysisError("J29 receipt validated provenance fields are incomplete")
    if not isinstance(receipt["branch"], str) or not receipt["branch"].startswith("codex/"):
        raise AnalysisError("J29 receipt branch is invalid")
    for name in ("base_commit", "validated_commit"):
        if not SHA1.fullmatch(receipt[name]):
            raise AnalysisError(f"J29 receipt {name} must be a lowercase 40-character SHA-1")
    if receipt.get("historical_preservation", {}).get("readiness") != {"production_validated": False, "simulation_readiness": "not_established"}:
        raise AnalysisError("J29 receipt readiness claims are invalid")
    phase = manifest.get("phase_4j29")
    if not isinstance(phase, dict):
        raise AnalysisError("Manifest phase_4j29 metadata is missing")
    for name in required:
        if phase.get(name) != receipt[name]:
            raise AnalysisError(f"J29 receipt/manifest mismatch: {name}")
    if manifest.get("branch") != receipt["branch"] or manifest.get("base_commit") != receipt["base_commit"]:
        raise AnalysisError("J29 top-level manifest provenance mismatch")
    if manifest.get("production_validated") is not False or manifest.get("simulation_readiness") != "not_established":
        raise AnalysisError("J29 manifest readiness claims are invalid")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    validate_receipt(args.receipt, args.manifest)
    print("J29 receipt validation: PASS")


if __name__ == "__main__":
    main()
