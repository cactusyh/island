import json
from pathlib import Path

import pytest

from island.exceptions import AnalysisError
from scripts.validate_phase_4j30_receipt import validate_receipt

RECEIPT = Path("docs/evidence/phase_4j30.json")
MANIFEST = Path("docs/paper/reproducibility_manifest.json")


def test_j30_receipt_validator_accepts_real_files():
    receipt = validate_receipt(RECEIPT, MANIFEST)
    assert receipt["fixture"]["chains"] == 2
    assert receipt["fixture"]["unequal_chain_lengths"] is True
    assert receipt["test_counts"]["full_supported"] == "2107 passed, 10 skipped, 0 failed, 0 errors"


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (lambda receipt, manifest: {**receipt, "fixture": {**receipt["fixture"], "chains": 3}}, "chain"),
        (lambda receipt, manifest: {**receipt, "aggregation": {**receipt["aggregation"], "policies": ["equal_weight"]}}, "aggregation"),
        (lambda receipt, manifest: {**receipt, "test_counts": {**receipt["test_counts"], "focused": "25 passed"}}, "test-count"),
        (lambda receipt, manifest: {**receipt, "readiness": {"production_validated": True, "simulation_readiness": "established"}}, "readiness"),
    ),
)
def test_j30_receipt_validator_rejects_mutated_receipt_files(tmp_path, mutation, message):
    receipt = json.loads(RECEIPT.read_text())
    manifest = json.loads(MANIFEST.read_text())
    mutated = tmp_path / "receipt.json"
    mutated.write_text(json.dumps(mutation(receipt, manifest)))
    with pytest.raises(AnalysisError, match=message):
        validate_receipt(mutated, MANIFEST)


def test_j30_receipt_validator_rejects_mutated_manifest_files(tmp_path):
    manifest = json.loads(MANIFEST.read_text())
    mutated = tmp_path / "manifest.json"
    changed = dict(manifest, phase_4j30=dict(manifest["phase_4j30"], fixture_hash="0" * 64))
    mutated.write_text(json.dumps(changed))
    with pytest.raises(AnalysisError, match="fixture_hash"):
        validate_receipt(RECEIPT, mutated)
