"""Exact-graph, charge-origin and periodic history controls for J39."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt import _signed
from island.exceptions import IslandError
from island.forcefields.oplsaa.resp import CHARGES, OPLSRESPBackend, validate_backend
from island.periodic_md import _backend
from scripts.phase_4j38_c2 import SourceIncompatibleError, require_native_source
from scripts.phase_4j38_c3_reference import audit
from scripts.validate_phase_4j39 import matrix
from scripts.validate_phase_4j39_receipt import validate_bundle

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    if os.environ.get("ISLAND_J39_RUN"):
        return Path(os.environ["ISLAND_J39_RUN"])
    root = tmp_path_factory.mktemp("j39-retained")
    with tarfile.open(ROOT / "docs/evidence/phase_4j39/artifacts.tar.gz") as archive:
        archive.extractall(root, filter="data")
    return root / "run"


@pytest.fixture(scope="module")
def backend(bundle):
    return read(bundle / "backend.json")


def test_accepted_model_is_explicit_resp_not_native(bundle, backend):
    OPLSRESPBackend(json.dumps(backend)).validate_integrity()
    result = validate_bundle(bundle)
    assert result["charge_origin"] == "RESP/Psi4"
    assert result["charge_vector_identity"] == CHARGES
    assert result["native_opls_charges"] is result["native_am1bcc_success"] is False
    assert result["checkpoint_byte_equal"] and result["rng_equal"]
    assert result["frame_count"] == 321
    assert result["production_validated"] is False


def test_native_residual_and_original_values_unchanged(bundle):
    decision = read(bundle / "native-source-decision.json")
    with pytest.raises(SourceIncompatibleError):
        require_native_source(decision)
    a = read(bundle / "native-source-audit.json")
    assert [c["residual_decimal_e"] for c in a["components"]] == ["-0.230", "-0.230"]
    assert a["total_decimal_e"] == "-0.460"
    assert {r["charge_e"] for r in a["rows"] if r["element"] == "C"} == {-0.115}
    assert {r["charge_e"] for r in a["rows"] if r["element"] == "H"} == {0.115}
    assert a["dynamics_executed"] is False


@pytest.mark.parametrize(
    "mutation",
    [
        "origin",
        "label",
        "native",
        "schema",
        "source",
        "vector",
        "unknown_id",
        "duplicate_id",
        "order",
        "nonfinite",
        "boolean",
        "unit",
        "component",
        "graph",
        "executable",
        "basis",
        "policy",
        "rounding",
        "different_molecule",
        "normalization",
    ],
)
def test_invalid_exact_charge_provenance_rejected(backend, mutation):
    p = deepcopy(backend)
    r = p["parameter_assignment"]["native_payload"]["record"]
    q = r["resp_charges"]
    if mutation == "origin":
        r["charge_origin"] = "native_nonbonded"
    elif mutation == "label":
        r["label"] = "native OPLS-AA"
    elif mutation == "native":
        r["native_opls_charges"] = True
    elif mutation == "schema":
        r["schema"] = "native"
    elif mutation == "source":
        r["source_xml"] += " "
    elif mutation == "vector":
        q["rows"][0]["charge_e"] += 0.001
    elif mutation == "unknown_id":
        q["rows"][0]["atom_id"] = 9999
    elif mutation == "duplicate_id":
        q["rows"][1]["atom_id"] = q["rows"][0]["atom_id"]
    elif mutation == "order":
        q["rows"].reverse()
    elif mutation == "nonfinite":
        q["rows"][0]["charge_e"] = float("nan")
    elif mutation == "boolean":
        q["rows"][0]["charge_e"] = True
    elif mutation == "unit":
        q["units"] = "coulomb"
    elif mutation == "component":
        q["source_audit"]["components"][0]["formal_total_e"] = 1
    elif mutation == "graph":
        q["final_graph_identity"] = "0" * 64
    elif mutation == "executable":
        q["runtime"]["files"]["bin/python"] = "0" * 64
    elif mutation == "basis":
        q["source_audit"]["basis"] = "STO-3G"
    elif mutation == "policy":
        q["policy"]["system"]["allowed_residual_e"] = 1.0
    elif mutation == "rounding":
        q["rows"][0]["charge_e"] = round(q["rows"][0]["charge_e"], 3)
    elif mutation == "different_molecule":
        q["charge_origin"] = "external_oplsaa"
    elif mutation == "normalization":
        q["normalization"] = True
    # Re-sign all public envelopes to prove semantic/source validation is required.
    from island.forcefields.ambertools.provided import reference_digest

    with pytest.raises((ValueError, RuntimeError, IslandError)):
        q["identity"] = reference_digest(
            {k: v for k, v in q.items() if k != "identity"}
        )
        p = _signed({k: v for k, v in p.items() if k != "identity"})
        validate_backend(p)


@pytest.mark.parametrize("key", ["charge", "sigma", "epsilon", "atom_type"])
def test_changed_numerical_source_assignment_rejected(backend, key):
    p = deepcopy(backend)
    p["numerical_data"]["sites"][0][key] = "opls_147" if key == "atom_type" else 0.123
    with pytest.raises((ValueError, IslandError), match="contradicts"):
        validate_backend(p)


def test_serialized_loader_enforces_resp_marker(backend):
    p = deepcopy(backend)
    p["parameter_assignment"]["native_payload"]["record"]["schema"] = "native"
    # Nested re-signing cannot disable detection by the retained typing profile.
    with pytest.raises((ValueError, IslandError)):
        _backend(p)


def test_no_charge_modification_or_serialization_loss(backend):
    q = backend["parameter_assignment"]["native_payload"]["record"]["resp_charges"]
    expected = {r["atom_id"]: r["ieee754_hex"] for r in q["rows"]}
    assert all(
        r["charge"].hex() == expected[r["id"]]
        for r in backend["numerical_data"]["sites"]
    )
    assert all(c["accepted"] for c in q["policy"]["components"])
    assert q["policy"]["system"]["accepted"]
    assert q["normalization"] is q["redistribution"] is False


def test_fresh_offline_child_reconstructs_relocated_model(bundle, tmp_path):
    shutil.copy(bundle / "backend.json", tmp_path / "backend.json")
    command = [
        sys.executable,
        "-c",
        'from pathlib import Path; from island.forcefields.oplsaa.resp import OPLSRESPBackend; OPLSRESPBackend(Path(__import__("sys").argv[1]).read_text()).validate_integrity()',
        str(tmp_path / "backend.json"),
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": ".:src"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_raw_artifacts_audited_without_production_importer(bundle):
    root = bundle / "raw"
    result = audit(
        root,
        "aromatic",
        read(root / "aromatic-charges.json"),
        read("docs/evidence/phase_4j36/environment.json"),
        read("scripts/phase_4j38_c3_anchors.json"),
    )
    assert result == read(bundle / "raw-review.json")
    assert result["import_bit_exact"] and result["raw_provenance_verified"]


@pytest.mark.parametrize(
    "field",
    [
        "positions_nm",
        "velocities_nm_ps",
        "box_vectors_nm",
        "potential_energy_kj_mol",
        "measured_temperature_kelvin",
        "stage",
        "checkpoint_sha256",
    ],
)
def test_trajectory_mutation_rejected(bundle, tmp_path, field):
    # Only the modified history is copied; other artifacts are read-only symlinks.
    for path in bundle.iterdir():
        (tmp_path / path.name).symlink_to(
            path.resolve(), target_is_directory=path.is_dir()
        )
    (tmp_path / "resumed").unlink()
    shutil.copytree(bundle / "resumed", tmp_path / "resumed")
    path = tmp_path / "resumed/frames.json"
    value = read(path)
    f = value["frames"][-1]
    if field in ("positions_nm", "velocities_nm_ps", "box_vectors_nm"):
        f[field][0][0] += 0.01
    elif field == "stage":
        f[field] = "nvt"
    elif field == "checkpoint_sha256":
        f[field] = "0" * 64
    else:
        f[field] += 1.0
    path.write_text(
        json.dumps(_signed({k: v for k, v in value.items() if k != "identity"}))
    )
    with pytest.raises((ValueError, IslandError)):
        validate_bundle(tmp_path)


def test_matrix_cannot_relabel_native_block_as_success(bundle):
    r = read(bundle / "result.json")
    m = matrix(r, read(bundle / "native-source-decision.json"))
    assert m == read(bundle / "matrix.json")
    assert m["paths"]["native_oplsaa"]["status"] == "source_incompatible"
    assert m["paths"]["native_oplsaa"]["dynamics_executed"] is False
    assert m["paths"]["opls_compatible_provided"]["status"] == "asset_unavailable"
    assert m["paths"]["opls_compatible_resp_psi4"]["status"] == "passed"
    assert m["acceptance_met"] is m["native_opls_acceptance_met"] is False
