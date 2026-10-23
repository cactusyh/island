"""Validate J39 model provenance and actual histories; native OPLS stays rejected."""

import argparse
import json
import subprocess
import tarfile
import tempfile
import xml.etree.ElementTree as ET
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.analysis.npt import FLAGS, _digest
from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.forcefields.oplsaa.resp import CHARGES, OPLSRESPBackend
from scripts.phase_4j38_c1_audit import require, sha
from scripts.phase_4j38_c2 import aromatic_decision
from scripts.phase_4j38_c3_reference import audit
from scripts.phase_4j38_contract import stable, validate_identity
from scripts.validate_phase_4j39 import BASE, matrix

SOURCES = (
    "src/island/forcefields/oplsaa/resp.py",
    "src/island/periodic_md.py",
    "scripts/phase_4j39_physics.py",
    "scripts/validate_phase_4j39.py",
    "scripts/validate_phase_4j39_receipt.py",
    "tests/test_phase_4j39.py",
)


def read(path):
    return json.loads(Path(path).read_text())


def validate_bundle(root):
    root = Path(root)
    backend = OPLSRESPBackend((root / "backend.json").read_text())
    backend.validate_integrity()
    p = backend.payload
    evidence = read(root / "fixture.json")
    require(
        evidence["final_graph"]["graph_identity"] == p["final_graph_identity"],
        "Wrong fixture graph",
    )
    native_audit = read(root / "native-source-audit.json")
    # Raw source inventory remains pinned independently of the alternate model.
    from xml.etree import ElementTree

    xml = p["parameter_assignment"]["native_payload"]["record"]["source_xml"]
    require(
        native_audit["source_sha256"] == sha256(xml.encode()).hexdigest(),
        "Changed native source hash",
    )
    require(
        native_audit["pin_sha256"] == sha("src/island/forcefields/oplsaa/pin.json"),
        "Changed source pin",
    )
    table = {
        a.attrib["type"]: a.attrib
        for a in ElementTree.fromstring(xml).find("NonbondedForce")
    }
    for row, site in zip(
        native_audit["rows"], evidence["final_graph"]["sites"], strict=True
    ):
        t = "opls_145" if site["element"] == "C" else "opls_146"
        require(
            row["atom_id"] == site["id"]
            and row["source_type"] == t
            and row["source_charge_decimal"] == table[t]["charge"],
            "Changed original source charge row",
        )
    decision = aromatic_decision(native_audit, evidence)
    require(
        decision == read(root / "native-source-decision.json"),
        "Native incompatibility changed",
    )
    raw = root / "raw"
    anchors = read("scripts/phase_4j38_c3_anchors.json")
    require(anchors == read(raw / "anchors.json"), "Changed trusted raw anchors")
    review = audit(
        raw,
        "aromatic",
        read(raw / "aromatic-charges.json"),
        read("docs/evidence/phase_4j36/environment.json"),
        anchors,
    )
    require(
        review == read(root / "raw-review.json") == read(raw / "aromatic-review.json"),
        "Changed raw independent review",
    )
    require(
        read(raw / "aromatic-charges.json")
        == p["parameter_assignment"]["native_payload"]["record"]["resp_charges"],
        "Wrong RESP source",
    )
    workflow, _ = load_continuation(root / "uninterrupted")
    require(read(root / "uninterrupted/backend.json") == p, "Workflow model differs")
    require(
        read(root / "resume.command.json")["exit_code"] == 0
        and read(root / "observables.command.json")["exit_code"] == 0,
        "Child process failed",
    )
    comparison = compare_continuations(
        root / "uninterrupted", root / "resumed", root / "relocated"
    )
    require(
        comparison == read(root / "comparison.json") and comparison["passed"],
        "Changed comparison",
    )
    require(
        comparison["checkpoint_hashes"] == comparison["resumed_checkpoint_hashes"]
        and comparison["rng_state_hashes"] == comparison["resumed_rng_state_hashes"],
        "Checkpoint/RNG mismatch",
    )
    frames = workflow["frames"]
    require(
        len(frames) == 321
        and stable(frames, workflow["analysis"])
        and workflow["analysis"]["density"]["blocks"] == 4,
        "Insufficient or unstable trajectory",
    )
    require(
        len({f["volume_nm3"] for f in frames if f["stage"] == "npt"}) > 1,
        "Fixed NPT volume",
    )
    require(
        all(f["measured_pressure_bar"] is None for f in frames), "Invented pressure"
    )
    with np.load(root / "trajectory.npz") as trajectory:
        for key in ("positions_nm", "velocities_nm_ps", "box_vectors_nm"):
            require(
                np.array_equal(trajectory[key], np.array([f[key] for f in frames])),
                "Trajectory content changed: " + key,
            )
        require(
            trajectory["atom_ids"].tolist() == workflow["manifest"]["atom_ids"],
            "Changed trajectory mapping",
        )
    physics = read(root / "physics.json")
    require(
        physics["passed"] is True
        and physics["charge_origin"] == "RESP/Psi4"
        and physics["native_opls_charges"] is False,
        "Wrong reference model",
    )
    require(
        all(
            s["energy_error_kj_mol"] <= 1e-5
            and s["force_error_kj_mol_angstrom"] <= 1e-5
            and max(s["component_errors_kj_mol"].values()) <= 1e-5
            for s in physics["snapshots"]
        ),
        "Reference mismatch",
    )
    observables = read(root / "observables.json")
    validate_identity(observables)
    require(
        observables["equilibrium_status"] == "insufficient_sampling",
        "False equilibrium claim",
    )
    result = read(root / "result.json")
    validate_identity(result)
    require(
        result["label"] == "OPLS-compatible parameters with RESP/Psi4 charges"
        and result["charge_origin"] == "RESP/Psi4"
        and result["native_opls_charges"] is result["native_am1bcc_success"] is False,
        "Misrepresented charge model",
    )
    expected = dict(
        charge_vector_identity=CHARGES,
        raw_review_identity=review["identity"],
        final_graph_identity=p["final_graph_identity"],
        packing_identity=evidence["packing_plan"]["identity"],
        backend_identity=p["identity"],
        workflow_identity=workflow["manifest"]["identity"],
        trajectory_identity=_digest(frames),
        analysis_identity=observables["identity"],
        comparison_identity=comparison["identity"],
        frame_count=321,
        split_boundary=128,
        comparison_tolerances="exact",
        checkpoint_byte_equal=True,
        rng_equal=True,
        status="passed",
        software_acceptance_met=True,
        dynamics_executed=True,
        equilibration_status="insufficient_sampling",
        measured_pressure_status="unavailable",
        **FLAGS,
    )
    require(
        all(
            type(result.get(k)) is type(v) and result[k] == v
            for k, v in expected.items()
        ),
        "Contradictory result binding",
    )
    require(
        read(root / "matrix.json") == matrix(result, decision),
        "Changed matrix or native success claim",
    )
    return result


def validate_receipt(path):
    path = Path(path)
    receipt = read(path)
    validate_identity(receipt)
    require(
        receipt["base_commit"] == BASE
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Changed base/readiness",
    )
    root = path.parent
    manifest = read(root / "manifest.json")
    require(
        receipt["manifest_sha256"] == sha(root / "manifest.json"), "Changed manifest"
    )
    for name, value in manifest.items():
        require(sha(root / name) == value, "Changed artifact: " + name)
    for name in SOURCES:
        raw = subprocess.check_output(
            ["git", "show", receipt["implementation_commit"] + ":" + name]
        )
        require(
            raw == Path(name).read_bytes(), "Tested implementation differs: " + name
        )
    historical = read(root / "historical-hashes.json")
    from scripts.validate_phase_4j38_c2_receipt import historical_hashes

    require(historical == historical_hashes(BASE), "Historical evidence changed")
    for name in ("focused", "j22-j38-c3", "full_supported"):
        command = read(root / "validation" / (name + "-command.json"))
        require(
            command["exit_code"] == 0
            and command["implementation_commit"] == receipt["implementation_commit"],
            "Wrong tested commit/result",
        )
        cases = list(ET.parse(root / "validation" / (name + ".xml")).iter("testcase"))
        actual = {
            k: sum(1 for case in cases if case.find(k) is not None)
            for k in ("failure", "error", "skipped")
        }
        actual["tests"] = len(cases)
        require(
            actual == receipt["tests"][name]
            and actual["failure"] == actual["error"] == 0,
            "JUnit counts changed",
        )
    for name in (
        "ruff",
        "compileall",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
        "diff_check",
        "historical_hashes",
        "historical_receipts",
        "environment_inventory",
        "independent_raw_audit",
        "relocated_bundle",
    ):
        require(
            read(root / "validation" / (name + "-command.json"))["exit_code"] == 0,
            "Unmet validation command: " + name,
        )
    with tempfile.TemporaryDirectory(prefix="j39-receipt-") as temporary:
        with tarfile.open(root / "artifacts.tar.gz") as archive:
            archive.extractall(temporary, filter="data")
        result = validate_bundle(Path(temporary) / "run")
    require(
        receipt["result_identity"] == result["identity"]
        and receipt["software_acceptance_met"] is True
        and receipt["native_opls_acceptance_met"] is False,
        "Incorrect acceptance claim",
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "receipt",
        nargs="?",
        type=Path,
        default=Path("docs/evidence/phase_4j39/receipt.json"),
    )
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--software-only", action="store_true")
    args = parser.parse_args()
    try:
        validate_bundle(args.bundle) if args.bundle else validate_receipt(args.receipt)
        print(
            "J39 software evidence valid; native OPLS source_incompatible; equilibrium not established"
        )
        return 0 if args.software_only else 1
    except Exception as error:  # noqa: BLE001 -- invalid evidence must return nonzero with diagnostic
        print(type(error).__name__ + ": " + str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
