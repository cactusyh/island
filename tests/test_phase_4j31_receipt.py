"""Receipt mutation controls using copies of retained J31 evidence."""

import pytest

from island.analysis.network import FLAGS, _signed
from island.exceptions import AnalysisError
from scripts.validate_phase_4j31_receipt import (
    validate_structure,
)


def example():
    # Positive acceptance uses the real retained six-success/two-rejection record.
    import json
    from pathlib import Path

    return (
        json.loads(Path("docs/evidence/phase_4j31.json").read_text()),
        json.loads(Path("docs/paper/reproducibility_manifest.json").read_text()),
    )


def test_receipt_structure_accepts_retained_j31_evidence():
    validate_structure(*example())


@pytest.mark.parametrize(
    "defect",
    (
        "chains",
        "graph_identity",
        "cycle_rank",
        "crosslink_count",
        "test_counts",
        "readiness",
        "blocked_execution",
        "missing_check",
        "missing_family",
        "schema",
    ),
)
def test_resigned_semantic_mutations_rejected(defect):
    receipt, manifest = example()
    if defect == "chains":
        receipt["fixtures"]["crosslinked"]["chains"] = 3
    elif defect == "graph_identity":
        receipt["fixtures"]["crosslinked"]["final_graph_identity"] = "invalid"
    elif defect in {"cycle_rank", "crosslink_count"}:
        receipt["fixtures"]["crosslinked"]["independent_reference"][defect] = 999
    elif defect == "test_counts":
        receipt["validation"]["suites"]["focused"]["counts"]["passed"] = 999
    elif defect == "readiness":
        receipt["readiness"] = {**FLAGS, "production_validated": True}
    elif defect == "blocked_execution":
        receipt["force_field_source_completeness"]["matrix"]["pcff/crosslinked"][
            "dynamics_executed"
        ] = True
    elif defect == "missing_check":
        receipt["validation"]["checks"].pop("pip_check")
    elif defect == "missing_family":
        receipt["force_field_source_completeness"]["matrix"].pop("gaff/crosslinked")
    else:
        receipt["schema"] = "unknown"
    receipt = _signed({k: v for k, v in receipt.items() if k != "identity"})
    with pytest.raises(AnalysisError):
        validate_structure(receipt, manifest)


def test_manifest_mismatch_and_corrupted_receipt_rejected():
    receipt, manifest = example()
    manifest["phase_4j31"]["base_commit"] = "f" * 40
    with pytest.raises(AnalysisError, match="manifest"):
        validate_structure(receipt, manifest)
    receipt, manifest = example()
    receipt["branch"] = "another-branch"
    with pytest.raises(AnalysisError, match="checksum"):
        validate_structure(receipt, manifest)
