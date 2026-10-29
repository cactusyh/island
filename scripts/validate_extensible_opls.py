"""Validate the independently pinned external OPLS software-integration receipt."""

import argparse
import json
import lzma
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.pdms_literature import pinned_bytes
from island.workflows.storage import child


def validate(directory, expected, repo):
    directory, repo = Path(directory), Path(repo)
    manifest = json.loads(pinned_bytes(directory / "manifest.json", expected))
    require(
        manifest["schema"] == "island_external_opls_integration_manifest_v1",
        "manifest_schema",
    )
    require(
        manifest["identity"]
        == digest({k: v for k, v in manifest.items() if k != "identity"}),
        "manifest_identity",
    )
    for name, pin in manifest["artifacts"].items():
        prefix = "docs/evidence/opls_extensible_v1/"
        # Verify the files actually read below, including a relocated evidence
        # directory; the original repository copy cannot vouch for a mutation.
        path = (
            child(directory, name[len(prefix) :])
            if name.startswith(prefix)
            else child(repo, name)
        )
        pinned_bytes(path, pin)
    baseline = json.loads(
        lzma.decompress((directory / "protected-file-hashes.json.xz").read_bytes())
    )
    require(
        len(baseline) == manifest["historical_files_protected"],
        "historical_inventory_changed",
    )
    for name, pin in baseline.items():
        pinned_bytes(child(repo, name), pin)
    receipt = json.loads((directory / "receipt.json").read_text())
    require(
        receipt["identity"]
        == digest({k: v for k, v in receipt.items() if k != "identity"}),
        "receipt_identity",
    )
    require(
        receipt["limitations"] == LIMITS
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "promoted_readiness",
    )
    require(
        receipt["software_integration_passed"] is True
        and receipt["physical_PDMS_profile_identity"] is None
        and receipt["PDMS_executable"] is False
        and receipt["native_pcff_acceptance_gate"]
        == receipt["PDMS_acceptance_gate"]
        == 1
        and receipt["new_physical_parameter_profiles"] == [],
        "forged_execution_authority",
    )
    runtime = json.loads((directory / "runtime.json").read_text())
    require(
        runtime["limitations"] == LIMITS
        and runtime["PDMS_executable"] is False
        and runtime["software_integration_passed"] is True,
        "runtime_readiness_changed",
    )
    require(
        runtime["xml_sha256"]
        == manifest["control_xml_sha256"]
        == "c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4",
        "wrong_control_source",
    )
    require(
        len(runtime["controls"]) == 2 and len(runtime["PDMS_cases"]) == 10,
        "changed_fixture_inventory",
    )
    for r in runtime["controls"]:
        require(
            r["status"] == "prepared"
            and r["source_row_reconstruction"][
                "native_atom_type_order_charge_mass_mapping_equal"
            ]
            is True
            and r["source_row_reconstruction"][
                "native_interaction_counts_and_selected_source_rows_equal"
            ]
            is True,
            "control_mapping_failed",
        )
        require(
            r["native_energy_absolute_difference"] < 1e-8
            and r["native_force_max_absolute_difference"] < 1e-8
            and max(r["independent_component_errors"].values()) < 1e-8
            and max(r["independent_finite_difference_force_max_errors"].values())
            < 1e-4,
            "numerical_verification_failed",
        )
        require(
            r["minimization"]["converged"] is True
            and r["minimization"]["final_verified"] is True
            and r["split_continuation_state_exact"] is True
            and r["split_continuation_RNG_exact"] is True
            and r["checkpoint_roundtrip_bytes_exact"] is True,
            "workflow_verification_failed",
        )
    for r in runtime["PDMS_cases"]:
        require(
            r["status"] in {"asset_unavailable", "typing_blocked"}
            and r["execution_profile"] is None
            and r["MD_run"] is False
            and r["numerical_engine_instantiated"] is False
            and r["graph_valid"] is True,
            "PDMS_boundary_bypassed",
        )
    reports = {}
    for name in ("focused-tests.xml", "full-suite.xml"):
        root = ET.parse(directory / name).getroot()
        suites = [root] if root.tag == "testsuite" else list(root)
        stats = {
            k: sum(int(s.attrib.get(k, 0)) for s in suites)
            for k in ("tests", "failures", "errors", "skipped")
        }
        require(
            stats["tests"] > 0 and stats["errors"] == stats["failures"] == 0,
            "failed_test_report:" + name,
        )
        if name == "focused-tests.xml":
            require(
                stats["tests"] >= 58 and stats["skipped"] == 0,
                "focused_coverage_missing",
            )
        reports[name] = stats
    quality = json.loads((directory / "quality-checks.json").read_text())
    require(
        all(c["exit_code"] == 0 for c in quality["checks"]), "failed_quality_checks"
    )
    require(
        len([c for c in quality["checks"] if c["name"].startswith("pip_")]) == 5,
        "five_pip_checks_required",
    )
    return {
        "manifest_valid": True,
        "receipt_valid": True,
        "historical_files_unchanged": len(baseline),
        "software_integration_passed": True,
        "PDMS_executable": False,
        "PDMS_acceptance_gate": 1,
        "test_reports": reports,
        "limitations": LIMITS,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--sha256", required=True)
    p.add_argument("--repo", type=Path, default=Path("."))
    a = p.parse_args()
    print(json.dumps(validate(a.directory, a.sha256, a.repo)))


if __name__ == "__main__":
    main()
