"""Validate J38 artifacts; --diagnostic-only never means benchmark acceptance."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

from island.analysis.npt import _digest, frame_measurements
from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.graph.final import _system_from_payload
from scripts.phase_4j38_contract import (
    BASE,
    require,
    stable,
    validate_identity,
    validate_matrix,
)
from scripts.phase_4j38_fixture import BENCHMARKS, build
from scripts.phase_4j38_physics import check


def read(path):
    return json.loads(Path(path).read_text())


def checksum(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def validate_structure(receipt, manifest, matrix):
    for value in (receipt, manifest):
        validate_identity(value)
    validate_matrix(matrix)
    require(
        receipt["schema"] == "island.j38-receipt.v1"
        and manifest["schema"] == "island.j38-reproducibility.v1",
        "Wrong receipt/manifest schema",
    )
    require(receipt["base_commit"] == manifest["base_commit"] == BASE, "Changed base")
    require(
        receipt["manifest_identity"] == manifest["identity"]
        and receipt["matrix_identity"]
        == matrix["identity"]
        == manifest["matrix_identity"],
        "Receipt/manifest/matrix mismatch",
    )
    require(
        manifest["benchmarks"] == [asdict(b) for b in BENCHMARKS],
        "Changed benchmark specification",
    )
    require(
        receipt["cases"] == matrix["cases"],
        "Receipt claims differ from executed matrix",
    )
    require(
        receipt["acceptance_met"] is False
        and receipt["software_acceptance_met"] is matrix["software_acceptance_met"],
        "Incorrect acceptance claim",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Wrong readiness",
    )
    require(
        matrix["pcff_control"]["acceptance_met"] is True
        and matrix["pcff_control"]["dynamics_executed"] is False
        and matrix["pcff_control"]["diagnostic"]["classification"]
        == "expected_typing_rejection",
        "Missing PCFF negative control",
    )
    require(
        matrix["native_am1bcc"]["status"] == "blocked"
        and matrix["native_am1bcc"]["acceptance_met"] is False,
        "Native AM1-BCC must remain blocked",
    )
    require(
        set(receipt["tests"]) == {"focused", "j22-j37", "full_supported"},
        "Test inventory mismatch",
    )
    for counts in receipt["tests"].values():
        require(
            counts["passed"] > 0 and counts["failed"] == counts["errors"] == 0,
            "Failed or missing required tests",
        )
    return True


def validate_case_bundle(root, row):
    """Recompute scientific data bindings and execution claims, not just JSON hashes."""
    key = row["benchmark"]
    fixture, evidence = build(key)
    require(
        read(root / "fixture.json") == json.loads(json.dumps(evidence)),
        "Changed PSMILES/packing/mapping/fixture evidence",
    )
    from scripts.validate_phase_4j38 import input_provenance

    require(
        row["input_provenance"]
        == json.loads(json.dumps(input_provenance(fixture, evidence))),
        "Changed benchmark input provenance",
    )
    if not row["execution_met"]:
        require(row.get("diagnostic"), "Missing blocked diagnostic")
        return
    workflow, restart = load_continuation(root / "uninterrupted")
    backend = read(root / "uninterrupted/backend.json")
    comparison = compare_continuations(
        root / "uninterrupted", root / "resumed", root / "relocated"
    )
    require(
        comparison["passed"] and comparison == read(root / "comparison-j38.json"),
        "Continuation comparison changed",
    )
    bindings = {
        "final_graph_identity": fixture.packed.graph.identity,
        "packing_identity": fixture.packed.plan.identity,
        "source_identity": workflow["manifest"]["source_identity"],
        "assignment_identity": backend["assignment_identity"],
        "charge_identity": backend["charge_identity"],
        "typed_graph_identity": backend["typed_graph_identity"],
        "workflow_identity": workflow["manifest"]["identity"],
        "trajectory_identity": _digest(workflow["frames"]),
        "checkpoint_identity": restart["checkpoint_identity"],
        "comparison_identity": comparison["identity"],
        "frame_count": len(workflow["frames"]),
        "density_blocks": workflow["analysis"]["density"]["blocks"],
    }
    require(
        all(row[k] == v for k, v in bindings.items()),
        "Workflow/trajectory/source/assignment binding changed",
    )
    require(
        backend["final_graph_identity"] == fixture.packed.graph.identity
        and backend["config"]["family"]
        == ("OPLS-AA" if row["family"] == "oplsaa" else "GAFF2"),
        "Wrong backend graph/family",
    )
    require(row["setup"] == evidence["setup"], "Copied setup metadata")
    require(
        row["gates"]["finite_time_stability"]["status"]
        == ("passed" if stable(workflow["frames"], workflow["analysis"]) else "failed"),
        "Changed stability claim",
    )
    require(
        len({f["volume_nm3"] for f in workflow["frames"] if f["stage"] == "npt"}) > 1,
        "NPT box did not change",
    )
    frames = frame_measurements(
        workflow["frames"], workflow["manifest"]["total_mass_dalton"]
    )
    observable = read(root / "observables.json")
    validate_identity(observable)
    require(
        observable["identity"] == row["analysis_identity"]
        and observable["frame_content_identity"] == bindings["trajectory_identity"],
        "Changed analysis/frame content",
    )
    require(
        [f["density_g_cm3"] for f in observable["frames"]]
        == [f["density_g_cm3"] for f in frames if f["stage"] == "npt"],
        "Copied density metadata",
    )
    from scripts.validate_phase_4j38 import ANALYSIS

    child = subprocess.run(
        [
            str(ANALYSIS),
            "-c",
            "import json,sys; from scripts.phase_4j38_observables import analyze; print(json.dumps(analyze(sys.argv[1])))",
            str(root),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    require(
        child.returncode == 0 and json.loads(child.stdout) == observable,
        "Fresh-process observable reconstruction failed",
    )
    import numpy as np

    with np.load(root / "trajectory.npz", allow_pickle=False) as export:
        require(
            np.array_equal(
                export["positions_nm"],
                np.array([f["positions_nm"] for f in workflow["frames"]]),
            )
            and np.array_equal(
                export["box_vectors_nm"],
                np.array([f["box_vectors_nm"] for f in workflow["frames"]]),
            )
            and export["atom_ids"].tolist() == workflow["manifest"]["atom_ids"],
            "Trajectory export changed",
        )
    if row["family"] == "oplsaa":
        require(
            check(
                SimpleNamespace(
                    system=_system_from_payload(backend["system"]), payload=backend
                ),
                workflow["frames"],
            )
            == read(root / "physics.json"),
            "Independent reference changed",
        )
    else:
        from scripts.phase_4j37_physics import check as amber_check

        require(
            amber_check(root) == read(root / "physics.json"),
            "Independent Amber reference changed",
        )


def validate(path):
    import xml.etree.ElementTree as ET

    path = Path(path)
    directory = path.parent
    receipt = read(path)
    manifest = read(directory / "manifest.json")
    matrix = read(directory / "matrix.json")
    validate_structure(receipt, manifest, matrix)
    for relative, expected in manifest["files"].items():
        require(
            not Path(relative).is_absolute() and ".." not in Path(relative).parts,
            "Unsafe manifest path",
        )
        require(
            checksum(directory / relative) == expected, "Changed artifact: " + relative
        )
    sha = receipt["implementation_commit"]
    require(
        len(sha) == 40
        and subprocess.check_output(["git", "rev-parse", sha], text=True).strip()
        == sha,
        "Missing implementation commit",
    )
    for name, expected in manifest["implementation_files"].items():
        require(
            sha256(
                subprocess.check_output(["git", "show", f"{sha}:{name}"])
            ).hexdigest()
            == expected
            and checksum(name) == expected,
            "Changed tested implementation: " + name,
        )
    for group, counts in receipt["tests"].items():
        suites = (
            ET.parse(directory / "validation" / (group + ".xml"))
            .getroot()
            .findall("testsuite")
        )
        actual = {
            k: sum(int(s.get(attr, 0)) for s in suites)
            for k, attr in [
                ("total", "tests"),
                ("failed", "failures"),
                ("errors", "errors"),
                ("skipped", "skipped"),
            ]
        }
        actual["passed"] = (
            actual["total"] - actual["failed"] - actual["errors"] - actual["skipped"]
        )
        require(actual == counts, "Test counts changed")
        command = read(directory / "validation" / (group + "-command.json"))
        require(
            command["implementation_commit"] == sha and command["exit_code"] == 0,
            "Wrong tested commit or failed test gate",
        )
    for name, expected_code in receipt["validation_gates"].items():
        command = read(directory / "validation" / (name + "-command.json"))
        require(
            command["implementation_commit"] == sha
            and command["exit_code"] == expected_code,
            "Unmet validation gate: " + name,
        )
    for name, expected in manifest["historical_files"].items():
        require(
            checksum(name)
            == expected
            == sha256(
                subprocess.check_output(["git", "show", BASE + ":" + name])
            ).hexdigest(),
            "Historical bytes changed: " + name,
        )
    with tempfile.TemporaryDirectory(prefix="j38-offline-") as tmp:
        root = Path(tmp)
        with tarfile.open(directory / "artifacts/workflows.tar.gz") as archive:
            archive.extractall(root, filter="data")
        runs = root / "runs"
        require(read(runs / "matrix.json") == matrix, "Archived matrix differs")
        from scripts.validate_phase_4j37 import validate_bundle

        for kind in ("linear", "network"):
            validate_bundle(runs / "verified-resp", kind)
        from scripts.phase_4j31_contract import typing_evidence
        from scripts.phase_4j32_assignment import prepare
        from scripts.phase_4j32_real_fixture import FRC
        from scripts.validate_phase_4j32 import pcff_contract

        fixture, _ = build("network")
        control = matrix["pcff_control"]
        d = control["diagnostic"]
        require(
            typing_evidence(fixture, FRC) == d["typing_evidence"],
            "Changed PCFF typing evidence",
        )
        try:
            prepare(fixture, "pcff", root)
        except Exception as error:  # noqa: BLE001 -- require exact native typing traceback
            require(
                pcff_contract(fixture, d["typing_evidence"], error)
                and str(error) == d["message"]
                and type(error).__name__ == d["exception_type"],
                "Changed expected PCFF rejection",
            )
        else:
            raise ValueError("Unexpected PCFF preparation success")
        for key, row in matrix["cases"].items():
            require(
                read(runs / key.replace("/", "-") / "result.json") == row,
                "Archived case differs",
            )
            validate_case_bundle(runs / key.replace("/", "-"), row)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("receipt", type=Path)
    p.add_argument("--diagnostic-only", action="store_true")
    args = p.parse_args()
    try:
        result = validate(args.receipt)
    except Exception as error:  # noqa: BLE001 -- diagnostic failure is nonzero
        print("INVALID: " + str(error))
        return 1
    print("VALID DIAGNOSTIC; acceptance=" + str(result["acceptance_met"]))
    return 0 if args.diagnostic_only else int(not result["acceptance_met"])


if __name__ == "__main__":
    raise SystemExit(main())
