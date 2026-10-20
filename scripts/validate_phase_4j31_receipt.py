"""Validate J31 receipt, real evidence files, fixture reconstruction and manifest.

Use a trusted expected receipt identity to detect wholesale consistently re-signed
replacement. Historical receipts are also compared with the exact reviewed base.
"""

import argparse
import json
import re
import subprocess
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

from island.analysis.network import FLAGS, _body, _digest, _require
from island.exceptions import AnalysisError
from scripts.phase_4j31_contract import expected_typing, require_matrix, typing_evidence

BASE = "aa3014648ddda4b032d6b07cdea381792330f99d"
BRANCH = "codex/phase-4j31-crosslinked-periodic-network-analysis"
SCHEMA = "island.phase-4j31-evidence.v1"
SUITES = {"focused", "j26_j30", "j22_j23", "full_supported"}
CHECKS = {"ruff", "compileall", "pip_check", "foyer_pip_check", "git_diff_check"}


def digest_file(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def junit_counts(path):
    root = ET.parse(path).getroot()
    rows = [root] if root.tag == "testsuite" else root.findall("testsuite")
    result = {
        key: sum(int(row.get(key, 0)) for row in rows)
        for key in ("tests", "failures", "errors", "skipped")
    }
    result["passed"] = (
        result["tests"] - result["failures"] - result["errors"] - result["skipped"]
    )
    return result


def manifest_section(receipt):
    return {
        "branch": receipt["branch"],
        "base_commit": receipt["base_commit"],
        "validated_commit": receipt["validated_commit"],
        "receipt_identity": receipt["identity"],
        "fixture_identities": {
            name: row["selection_identity"] for name, row in receipt["fixtures"].items()
        },
        "test_counts": {
            name: row["counts"] for name, row in receipt["validation"]["suites"].items()
        },
        "source_completeness": {
            name: row["status"]
            for name, row in receipt["force_field_source_completeness"][
                "matrix"
            ].items()
        },
        **FLAGS,
    }


def validate_structure(receipt, manifest):
    _require(
        receipt.get("schema") == SCHEMA
        and receipt.get("identity") == _digest(_body(receipt)),
        "J31 receipt checksum/schema mismatch",
    )
    _require(
        receipt.get("base_commit") == BASE and receipt.get("branch") == BRANCH,
        "J31 branch/base mismatch",
    )
    _require(
        type(receipt.get("validated_commit")) is str
        and re.fullmatch(r"[0-9a-f]{40}", receipt["validated_commit"]),
        "J31 validated commit must be full SHA",
    )
    _require(
        receipt.get("readiness") == FLAGS
        and type(receipt["readiness"].get("production_validated")) is bool,
        "J31 readiness mismatch",
    )
    fixtures = receipt.get("fixtures", {})
    _require(
        set(fixtures) == {"uncrosslinked", "branched", "crosslinked"},
        "J31 fixture coverage mismatch",
    )
    for name, row in fixtures.items():
        _require(
            row.get("chains") == 2 and type(row["chains"]) is int,
            "J31 fixture chain count mismatch",
        )
        for key in (
            "final_graph_identity",
            "crosslink_plan_identity",
            "packing_identity",
            "selection_identity",
            "atom_mapping_identity",
            "frame_content_identity",
        ):
            _require(
                type(row.get(key)) is str and re.fullmatch(r"[0-9a-f]{64}", row[key]),
                f"J31 fixture identity missing: {key}",
            )
        reference = row.get("independent_reference", {})
        _require(
            reference.get("component_count") == (2 if name == "uncrosslinked" else 1),
            "J31 component count mismatch",
        )
        _require(
            reference.get("crosslink_count")
            == {"uncrosslinked": 0, "branched": 1, "crosslinked": 2}[name],
            "J31 crosslink count mismatch",
        )
        _require(
            reference.get("cycle_rank") == (1 if name == "crosslinked" else 0),
            "J31 cycle rank mismatch",
        )
    validation = receipt.get("validation", {})
    suites = validation.get("suites", {})
    _require(set(suites) == SUITES, "J31 suite coverage mismatch")
    for name, row in suites.items():
        counts = row.get("counts", {})
        _require(
            set(counts) == {"tests", "passed", "failures", "errors", "skipped"}
            and all(type(v) is int and v >= 0 for v in counts.values()),
            "J31 malformed test counts",
        )
        _require(
            counts["tests"]
            == sum(counts[k] for k in ("passed", "failures", "errors", "skipped"))
            and counts["passed"] > 0
            and counts["failures"] == counts["errors"] == 0,
            f"J31 test counts inconsistent/failed: {name}",
        )
    checks = validation.get("checks", {})
    _require(
        set(checks) == CHECKS
        and all(row.get("exit_code") == 0 for row in checks.values()),
        "J31 required checks incomplete",
    )
    matrix = receipt.get("force_field_source_completeness", {}).get("matrix", {})
    require_matrix(matrix, legacy_pcff=True)
    _require(
        manifest.get("phase_4j31") == manifest_section(receipt),
        "J31 receipt/manifest mismatch",
    )
    _require(
        manifest.get("production_validated") is False
        and manifest.get("simulation_readiness") == "not_established",
        "J31 manifest readiness mismatch",
    )


def fixture_summaries():
    from scripts.phase_4j31_network_fixture import (
        build_network_fixture,
        graph_control_trajectory,
        independent_topology,
    )

    result = {}
    for name in ("uncrosslinked", "branched", "crosslinked"):
        fixture = build_network_fixture(name)
        trajectory, selection = graph_control_trajectory(fixture)
        result[name] = {
            "chains": len(set(trajectory.molecule_membership.values())),
            "atoms": len(trajectory.metadata.atom_ids),
            "frames": trajectory.frame_count,
            "final_graph_identity": fixture.packed.graph.identity,
            "crosslink_plan_identity": fixture.crosslinked.plan.identity,
            "packing_identity": fixture.packed.plan.identity,
            "selection_identity": selection.identity,
            "atom_mapping_identity": trajectory.mapping_identity,
            "frame_content_identity": selection.payload["frame_content_identity"],
            "independent_reference": independent_topology(fixture.packed.graph),
        }
    return result


def validate_receipt(
    receipt_path, manifest_path, *, root=Path("."), expected_identity=None
):
    root = Path(root)
    receipt = json.loads(Path(receipt_path).read_text())
    manifest = json.loads(Path(manifest_path).read_text())
    validate_structure(receipt, manifest)
    if expected_identity is not None:
        _require(
            receipt["identity"] == expected_identity, "Unexpected J31 receipt identity"
        )
    _require(
        receipt["validated_commit"] == "ba0e2abd951278958a852964da790126f49bbf65",
        "J31 historical implementation commit mismatch",
    )
    frozen = json.loads(
        subprocess.check_output(
            [
                "git",
                "show",
                "6606895d643806c731b52cf0fafa00e338ce3adb:docs/evidence/phase_4j31.json",
            ],
            cwd=root,
        )
    )
    for field in (
        "implementation_files",
        "artifact_files",
        "historical_receipts",
        "source_matrix_path",
        "acceptance_directory",
    ):
        _require(
            receipt.get(field) == frozen[field],
            f"J31 frozen evidence inventory mismatch: {field}",
        )
    for path, checksum in receipt["implementation_files"].items():
        # Historical source is always checked at its actual validated commit.
        # Four acceptance-only files intentionally differ in C1; current bytes
        # are covered by the additive C1 receipt, without rewriting J31 hashes.
        corrected = {
            "scripts/validate_phase_4j31_network.py",
            "scripts/validate_phase_4j31_receipt.py",
            "tests/test_phase_4j31_execution.py",
            "tests/test_phase_4j31_receipt.py",
        }
        if path not in corrected:
            _require(
                digest_file(root / path) == checksum,
                f"J31 changed implementation file: {path}",
            )
        committed = subprocess.check_output(
            ["git", "show", f"{receipt['validated_commit']}:{path}"], cwd=root
        )
        _require(
            sha256(committed).hexdigest() == checksum,
            f"J31 validated commit source mismatch: {path}",
        )
    _require(bool(receipt["implementation_files"]), "J31 implementation hashes missing")
    for path, checksum in receipt["historical_receipts"].items():
        original = subprocess.check_output(["git", "show", f"{BASE}:{path}"], cwd=root)
        _require(
            sha256(original).hexdigest() == digest_file(root / path) == checksum,
            f"J31 historical receipt modified: {path}",
        )
    expected_history = {
        str(path.relative_to(root))
        for path in (root / "docs/evidence").glob("phase_4j*.json")
        if path.name != "phase_4j31.json"
    }
    _require(
        set(receipt["historical_receipts"]) == expected_history,
        "J31 historical receipt inventory mismatch",
    )
    for path, checksum in receipt["artifact_files"].items():
        _require(
            digest_file(root / path) == checksum,
            f"J31 acceptance artifact hash mismatch: {path}",
        )
    for name, row in receipt["validation"]["suites"].items():
        _require(
            junit_counts(root / row["junit"]) == row["counts"],
            f"J31 JUnit/test-count mismatch: {name}",
        )
    matrix = json.loads((root / receipt["source_matrix_path"]).read_text())
    _require(
        matrix == receipt["force_field_source_completeness"],
        "J31 source matrix artifact mismatch",
    )
    from island.analysis import load_network_analysis_bundle, load_verified_trajectory
    from island.analysis.network import _frame_data
    from island.periodic_md import (
        compare_periodic_md_workflows,
        load_periodic_md_workflow,
    )

    acceptance = root / receipt["acceptance_directory"]
    for name, row in receipt["force_field_source_completeness"]["matrix"].items():
        directory = acceptance / name.replace("/", "-")
        _require(
            json.loads((directory / "result.json").read_text()) == row,
            "J31 retained case result/matrix mismatch",
        )
        if row["family"] == "pcff":
            from scripts.phase_4j31_network_fixture import FRC, build_network_fixture

            evidence = typing_evidence(build_network_fixture(row["fixture"]), FRC)
            _require(
                expected_typing(row["fixture"], evidence),
                "J31 legacy PCFF source/domain replay mismatch",
            )
            continue
        workflow = load_periodic_md_workflow(
            directory / "full", expected_backend_identity=row["backend_identity"]
        )
        _require(
            workflow.payload["manifest"]["backend"]["config"]["source"]
            == row["source"],
            "J31 archived source binding mismatch",
        )
        comparison = compare_periodic_md_workflows(
            directory / "full", directory / "relocated"
        )
        _require(
            comparison == row["execution"]["comparison"],
            "J31 archived split/continuous workflow mismatch",
        )
        _, selection, record = load_network_analysis_bundle(
            directory / "analysis-relocated",
            expected_identity=row["execution"]["analysis_record_identity"],
        )
        _require(
            all(
                selection.payload[key] == row[key]
                for key in (
                    "final_graph_identity",
                    "crosslink_plan_identity",
                    "packing_identity",
                    "typed_graph_identity",
                    "charge_identity",
                    "assignment_identity",
                )
            ),
            "J31 archived graph/source/assignment binding mismatch",
        )
        _require(
            selection.identity == row["execution"]["selection_identity"],
            "J31 archived selection mismatch",
        )
        restored = load_verified_trajectory(directory / "full")
        for key in (
            "workflow_identity",
            "checkpoint_identity",
            "final_graph_identity",
            "packing_identity",
            "force_field_source_identity",
            "typed_graph_identity",
            "charge_identity",
            "assignment_identity",
        ):
            _require(
                getattr(restored.metadata, key) == selection.payload[key],
                f"J31 archived workflow binding mismatch: {key}",
            )
        _require(
            restored.metadata.force_field_family
            == {"gaff": "GAFF", "gaff2": "GAFF2", "oplsaa": "OPLS-AA"}[row["family"]],
            "J31 archived force-field family mismatch",
        )
        if row["family"] in {"gaff", "gaff2"}:
            _require(
                digest_file(directory / "amber-source-artifacts/result.prmtop")
                == row["source"]["sha256"],
                "J31 archived Amber source mismatch",
            )
        _require(
            _digest(_frame_data(restored))
            == selection.payload["frame_content_identity"],
            "J31 archived trajectory/frame mismatch",
        )
        _require(
            record.payload["observables"] == row["execution"]["network_observables"],
            "J31 archived observable mismatch",
        )
    _require(
        json.loads(json.dumps(fixture_summaries())) == receipt["fixtures"],
        "J31 reconstructed fixture/reference mismatch",
    )
    # Every pre-existing manifest entry is historical; add only phase_4j31.
    original_manifest = json.loads(
        subprocess.check_output(
            ["git", "show", f"{BASE}:docs/paper/reproducibility_manifest.json"],
            cwd=root,
        )
    )
    _require(
        {k: v for k, v in manifest.items() if k != "phase_4j31"} == original_manifest,
        "J31 modified historical manifest entries",
    )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--expected-identity")
    args = parser.parse_args()
    validate_receipt(
        args.receipt, args.manifest, expected_identity=args.expected_identity
    )
    print("J31 receipt validation: PASS")


if __name__ == "__main__":
    try:
        main()
    except AnalysisError as error:
        raise SystemExit(str(error)) from error
