"""Structural receipt mutation controls with explicitly synthetic test records.

The real final receipt additionally validates source commits, JUnit files,
acceptance artifacts and reconstructed fixtures through the command-line tool.
"""

from copy import deepcopy

import pytest

from island.analysis.network import FLAGS, _signed
from island.exceptions import AnalysisError
from scripts.validate_phase_4j31_receipt import (
    BASE,
    BRANCH,
    CHECKS,
    SCHEMA,
    SUITES,
    manifest_section,
    validate_structure,
)


def example():
    fixtures = {}
    for kind in ("uncrosslinked", "branched", "crosslinked"):
        fixtures[kind] = {
            "chains": 2,
            **{
                name: "a" * 64
                for name in (
                    "final_graph_identity",
                    "crosslink_plan_identity",
                    "packing_identity",
                    "selection_identity",
                    "atom_mapping_identity",
                    "frame_content_identity",
                )
            },
            "independent_reference": {
                "component_count": 2 if kind == "uncrosslinked" else 1,
                "crosslink_count": {
                    "uncrosslinked": 0,
                    "branched": 1,
                    "crosslinked": 2,
                }[kind],
                "cycle_rank": int(kind == "crosslinked"),
            },
        }
    # These are unit-test data, not claims about a real suite or force field.
    counts = {"tests": 1, "passed": 1, "failures": 0, "errors": 0, "skipped": 0}
    receipt = _signed(
        {
            "schema": SCHEMA,
            "branch": BRANCH,
            "base_commit": BASE,
            "validated_commit": BASE,
            "readiness": FLAGS,
            "fixtures": fixtures,
            "validation": {
                "suites": {name: {"counts": deepcopy(counts)} for name in SUITES},
                "checks": {name: {"exit_code": 0} for name in CHECKS},
            },
            "force_field_source_completeness": {
                "matrix": {
                    f"{family}/{kind}": {
                        "status": "blocked",
                        "reason": "Synthetic receipt unit test; no source examined",
                        "dynamics_executed": False,
                        **FLAGS,
                    }
                    for family in ("pcff", "oplsaa", "gaff", "gaff2")
                    for kind in ("branched", "crosslinked")
                }
            },
        }
    )
    manifest = {"phase_4j31": manifest_section(receipt), **FLAGS}
    return receipt, manifest


def test_receipt_structure_accepts_explicit_unit_fixture():
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
