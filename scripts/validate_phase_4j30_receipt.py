"""Validate the J30 multi-chain evidence receipt against the manifest."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from island.exceptions import AnalysisError

SHA256 = re.compile(r"^[0-9a-f]{64}$")
SHA1 = re.compile(r"^[0-9a-f]{40}$")


def _sha(value: object, label: str) -> None:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        raise AnalysisError(f"J30 {label} must be a lowercase 64-character SHA-256")


def validate_receipt(receipt_path: str | Path, manifest_path: str | Path) -> dict:
    receipt = json.loads(Path(receipt_path).read_text())
    manifest = json.loads(Path(manifest_path).read_text())
    if receipt.get("schema") != "island.phase-4j30-evidence.v1":
        raise AnalysisError("J30 receipt schema mismatch")
    for name in ("base_commit", "validated_commit"):
        if not isinstance(receipt.get(name), str) or not SHA1.fullmatch(receipt[name]):
            raise AnalysisError(f"J30 receipt {name} must be a lowercase 40-character SHA-1")
    if not isinstance(receipt.get("branch"), str) or not receipt["branch"].startswith("codex/"):
        raise AnalysisError("J30 receipt branch is invalid")
    fixture = receipt.get("fixture")
    if not isinstance(fixture, dict):
        raise AnalysisError("J30 fixture metadata is missing")
    for name in ("fixture_hash", "final_graph_identity", "packing_identity", "force_field_source_identity", "typed_graph_identity", "charge_identity", "workflow_identity", "checkpoint_identity", "trajectory_identity", "atom_mapping_identity", "selection_identity"):
        _sha(fixture.get(name), f"fixture {name}")
    if fixture.get("chains") != 2 or fixture.get("atoms") != 33 or fixture.get("frames") != 7 or fixture.get("chain_lengths_atoms") != [20, 13] or fixture.get("unequal_chain_lengths") is not True:
        raise AnalysisError("J30 fixture chain, length, atom, or frame metadata is inconsistent")
    aggregation = receipt.get("aggregation")
    if not isinstance(aggregation, dict) or aggregation.get("policies") != ["equal_weight", "mass_weighted"]:
        raise AnalysisError("J30 aggregation policies are invalid")
    definitions = aggregation.get("definitions", {})
    if set(definitions) != {"equal_weight", "mass_weighted"}:
        raise AnalysisError("J30 aggregation definitions are incomplete")
    independent = aggregation.get("independent_first_frame", {})
    if set(independent) != {"equal_weight", "mass_weighted"}:
        raise AnalysisError("J30 independent aggregate controls are incomplete")
    if receipt.get("readiness") != {"production_validated": False, "simulation_readiness": "not_established"}:
        raise AnalysisError("J30 receipt readiness claims are invalid")
    phase = manifest.get("phase_4j30")
    if not isinstance(phase, dict):
        raise AnalysisError("Manifest phase_4j30 metadata is missing")
    for name in ("branch", "base_commit", "validated_commit"):
        if phase.get(name) != receipt.get(name):
            raise AnalysisError(f"J30 receipt/manifest mismatch: {name}")
    for name in ("fixture_hash", "chains", "chain_lengths_atoms", "atoms", "frames", "atom_mapping_identity", "final_graph_identity", "packing_identity", "workflow_identity", "checkpoint_identity", "trajectory_identity"):
        if phase.get(name) != fixture.get(name):
            raise AnalysisError(f"J30 receipt/manifest mismatch: {name}")
    if phase.get("aggregation_policies") != aggregation.get("policies") or phase.get("aggregation_definitions") != definitions or phase.get("independent_first_frame") != independent:
        raise AnalysisError("J30 receipt/manifest aggregation mismatch")
    counts = receipt.get("test_counts")
    if not isinstance(counts, dict) or phase.get("test_counts") != counts:
        raise AnalysisError("J30 receipt/manifest test-count mismatch")
    focused = receipt.get("focused_tests")
    if not isinstance(focused, dict) or focused.get("result") != counts.get("focused"):
        raise AnalysisError("J30 focused test count is stale")
    if phase.get("focused_tests") != counts.get("focused") or phase.get("compatibility_analysis_tests") != counts.get("analysis_compatibility") or phase.get("j22_j23_tests") != counts.get("j22_j23") or phase.get("full_supported_suite") != counts.get("full_supported"):
        raise AnalysisError("J30 manifest test counts are stale")
    integration = receipt.get("integration")
    if not isinstance(integration, dict) or phase.get("integration") != integration:
        raise AnalysisError("J30 receipt/manifest integration mismatch")
    if manifest.get("production_validated") is not False or manifest.get("simulation_readiness") != "not_established":
        raise AnalysisError("J30 manifest readiness claims are invalid")
    for name in ("focused_tests", "compatibility_analysis_tests", "j22_j23_tests", "full_supported_suite"):
        if name in phase and not isinstance(phase[name], str):
            raise AnalysisError(f"J30 manifest test count is malformed: {name}")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    validate_receipt(args.receipt, args.manifest)
    print("J30 receipt validation: PASS")


if __name__ == "__main__":
    main()
