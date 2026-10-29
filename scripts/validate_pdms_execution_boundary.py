"""Validate the independently pinned PDMS boundary evidence, including relocation."""

import argparse
import json
import lzma
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.pdms_literature import pinned_bytes
from island.workflows.storage import child

PREFIX = "docs/evidence/opls_pdms_execution_boundary_v4/"


def validate(directory, expected, repo):
    directory, repo = Path(directory), Path(repo)
    manifest = json.loads(pinned_bytes(directory / "manifest.json", expected))
    require(
        manifest["schema"] == "island_PDMS_execution_boundary_manifest_v4",
        "manifest_schema",
    )
    require(
        manifest["identity"]
        == digest({k: v for k, v in manifest.items() if k != "identity"}),
        "manifest_identity",
    )
    for name, pin in manifest["artifacts"].items():
        path = (
            child(directory, name[len(PREFIX) :])
            if name.startswith(PREFIX)
            else child(repo, name)
        )
        pinned_bytes(path, pin)
    baseline = json.loads(
        lzma.decompress((directory / "protected-file-hashes.json.xz").read_bytes())
    )
    require(len(baseline) == manifest["historical_files_protected"], "historical_count")
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
        receipt["PDMS_executable"] is False
        and receipt["PDMS_MD_run"] is False
        and receipt["execution_profiles"] == []
        and receipt["acceptance_gate"] == 1,
        "forged_execution_authority",
    )
    rights = json.loads((directory / "publication-rights.json").read_text())
    require(
        rights["grant"][0]["raw_text"] == "This article is licensed under CC-BY 4.0",
        "paper_grant",
    )
    require(
        rights["sha256"]
        == "04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc",
        "paper_pin",
    )
    graph = json.loads((directory / "mapping/report.json").read_text())
    require(
        graph["source_rights"]["author_repository"]["use_permission"]
        == "not_established",
        "forged_asset_permission",
    )
    require(
        graph["PDMS_executable"] is False and graph["PDMS_MD_run"] is False,
        "mapping_execution_changed",
    )
    for name in ("hmds", "tms_capped_dp2", "tms_capped_dp4"):
        row = graph["cases"][name]
        require(row["status"] == "asset_unavailable", "blocked_status_changed")
        report = json.loads(
            (directory / ("mapping/" + name + "-independent.json")).read_text()
        )
        require(
            report["independent_row_mapping_equal"] is True
            and report["numerical_PDMS_execution"] is False,
            "source_reconstruction",
        )
    for name in ("tms_capped_dp2_unoriented", "tms_capped_dp4_unoriented"):
        require(
            graph["cases"][name]["status"] == "typing_ambiguous",
            "ambiguous_graph_promoted",
        )
    require(
        graph["cases"]["unsupported_sih3"]["status"] == "typing_blocked",
        "unsupported_Si_promoted",
    )
    backend = json.loads((directory / "longrange/report.json").read_text())
    require(
        backend["software_equations_verified"] is True
        and backend["author_model_qualified"] is False
        and backend["physical_PDMS_executable"] is False
        and backend["periodic_scientific_closure"] is False,
        "backend_scope_promoted",
    )
    require(
        backend["controls"][0]["status"] == "scientifically_unqualified"
        and backend["author_default_accuracy_control_verified"] is False,
        "author_failure_erased",
    )
    tight = backend["controls"][-1]
    require(
        tight["status"] == "verified"
        and tight["max_force_difference"] < 3e-4
        and max(tight["component_absolute_errors"].values()) < 3e-4,
        "backend_comparison",
    )
    require(
        backend["independent_Ewald_alpha_and_shell_convergence"] < 1e-8,
        "Ewald_reference_convergence",
    )
    for row in tight["finite_difference_errors"].values():
        require(max(row.values()) < 1e-4, "independent_force_failure")
    stats = {}
    for name in ("focused-tests.xml", "regressions.xml", "backend-final-tests.xml"):
        root = ET.parse(directory / name).getroot()
        suites = [root] if root.tag == "testsuite" else list(root)
        count = {
            k: sum(int(s.attrib.get(k, 0)) for s in suites)
            for k in ("tests", "failures", "errors", "skipped")
        }
        require(
            count["failures"] == count["errors"] == count["skipped"] == 0,
            "test_failure_or_skip:" + name,
        )
        stats[name] = count
    return {
        "valid": True,
        "historical_files_unchanged": len(baseline),
        "tests": stats,
        "PDMS_executable": False,
        "PDMS_acceptance_gate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.directory, args.sha256, args.repo), indent=2))


if __name__ == "__main__":
    main()
