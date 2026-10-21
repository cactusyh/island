"""Re-signed software controls against genuine retained J32 workflow bundles."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.analysis.condensed import files
from island.analysis.network import _signed
from scripts.phase_4j32_receipt import extract_archive
from scripts.validate_phase_4j32 import ANALYSIS_PYTHON


@pytest.fixture(scope="module")
def retained(tmp_path_factory):
    candidate = os.environ.get("ISLAND_J32_TEST_RUN")
    if candidate:
        return Path(candidate)
    root = tmp_path_factory.mktemp("j32-bundles")
    extract_archive(Path("docs/evidence/phase_4j32/artifacts/acceptance.tar.gz"), root)
    return root


def test_readonly_real_bundle_reconstruction(retained):
    p = subprocess.run(
        [
            str(ANALYSIS_PYTHON),
            "-m",
            "scripts.phase_4j32_receipt",
            "--runs",
            str(retained),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode == 0, p.stdout + p.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "source",
        "source_missing",
        "frame",
        "checkpoint",
        "npt",
        "thermodynamic",
        "tg",
        "mechanical",
        "charge",
        "relocation",
        "fresh_process",
        "frames_short",
    ],
)
def test_bundle_mutations_fail_with_diagnostics(retained, tmp_path, mutation):
    target = tmp_path / "bundle"
    shutil.copytree(retained / "oplsaa-linear/bundle", target)
    path = target / "analysis.json"
    record = json.loads(path.read_text())
    if mutation == "source":
        with (target / "source-assets/oplsaa.xml").open("a") as f:
            f.write("\n<!-- software negative control -->\n")
    elif mutation == "source_missing":
        (target / "source-assets/oplsaa.xml").unlink()
    elif mutation in {"npt", "thermodynamic"}:
        execution_path = target / "execution.json"
        e = json.loads(execution_path.read_text())
        if mutation == "npt":
            e["npt"]["status"] = "complete"
        else:
            e["thermodynamic_origin"] = "inherited from J31"
        execution_path.write_text(json.dumps(e))
    elif mutation in {"tg", "mechanical"}:
        key = "glass_transition" if mutation == "tg" else "mechanical_properties"
        record["observables"][key].update(status="ok", mean=300)
    elif mutation == "charge":
        assignment_path = target / "assignment.json"
        assignment = json.loads(assignment_path.read_text())
        assignment["atom_types_and_charges"][0]["charge"] += 0.01
        assignment_path.write_text(json.dumps(assignment))
    elif mutation in {"relocation", "fresh_process"}:
        shutil.rmtree(target / "continuation")
    else:
        # Mutate real archived content, keeping the top-level input manifest
        # consistently re-signed. J23's own verification must still reject it.
        actual = next((target / "workflow").rglob("record.json"))
        envelope = json.loads(actual.read_text())
        if mutation == "frames_short":
            envelope["payload"]["frames"] = envelope["payload"]["frames"][:4]
        elif mutation == "checkpoint":
            envelope["payload"]["frames"][-1]["checkpoint_sha256"] = "0" * 64
        else:
            envelope["payload"]["frames"][-1]["positions_nm"][0][0] += 0.01
        actual.write_text(json.dumps(envelope))
    record["input_files"] = files(target)
    path.write_text(json.dumps(_signed(record)))
    code = "from island.analysis.condensed import load_bundle;import sys;load_bundle(sys.argv[1])"
    p = subprocess.run(
        [str(ANALYSIS_PYTHON), "-c", code, str(target)],
        capture_output=True,
        text=True,
        check=False,
    )
    (tmp_path / "failure.log").write_text(p.stdout + p.stderr)
    assert p.returncode != 0
    assert "Traceback" in p.stderr
    assert (tmp_path / "failure.log").stat().st_size > 0


def test_archive_traversal_rejected(tmp_path):
    import io
    import tarfile

    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        info = tarfile.TarInfo("../escape")
        info.size = 1
        stream.addfile(info, io.BytesIO(b"x"))
    with pytest.raises(ValueError, match="Unsafe"):
        extract_archive(archive, tmp_path / "destination")


def test_inventory_cli_does_not_silently_accept_missing_cases(tmp_path):
    (tmp_path / "matrix.json").write_text(
        json.dumps(_signed({"cases": {}, "acceptance_met": True}))
    )
    p = subprocess.run(
        [sys.executable, "-m", "scripts.phase_4j32_receipt", "--runs", str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0
    assert "inventory" in p.stderr


def test_record_publication_failure_is_transactional(tmp_path, monkeypatch):
    import island.analysis.condensed as module

    monkeypatch.setattr(module, "analyze", lambda _: {"label": "software control only"})

    def failure(*args):
        raise OSError("injected publication failure")

    monkeypatch.setattr(os, "link", failure)
    with pytest.raises(OSError, match="publication"):
        module.save_bundle_record(tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("field", ["force_field_family", "packing_identity"])
def test_case_bundle_bindings_cannot_be_resigned_away(retained, field):
    from scripts.phase_4j32_receipt import validate_bundle_binding

    matrix = json.loads((retained / "matrix.json").read_text())
    row = matrix["cases"]["oplsaa/linear"]["result"]
    record = json.loads((retained / "oplsaa-linear/bundle/analysis.json").read_text())
    validate_bundle_binding("oplsaa/linear", row, record)
    record["bindings"][field] = "GAFF2" if field == "force_field_family" else "0" * 64
    with pytest.raises(ValueError, match="Case/bundle"):
        validate_bundle_binding("oplsaa/linear", row, _signed(record))


def test_generation_cli_requires_full_validation_before_success(tmp_path):
    code = """
import sys
import scripts.phase_4j32_receipt as receipt
receipt.generate = lambda *args: None
def rejected(*args):
    raise ValueError('injected receipt-validation failure')
receipt.validate = rejected
sys.argv = ['receipt', str(sys.argv[1]), '--generate', 'software-control']
receipt.main()
"""
    p = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path / "receipt.json")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0
    assert "receipt-validation failure" in p.stderr
    assert "PASS" not in p.stdout
