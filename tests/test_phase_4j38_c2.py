"""C2 source chemistry and provenance controls; no synthetic acceptance evidence."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt import _signed
from scripts.phase_4j38_c2 import (
    SourceIncompatibleError,
    aromatic_decision,
    require_native_source,
    validate_matrix,
)

ROOT = Path(__file__).resolve().parents[1]
PYTHON = ROOT.parent / "island-validation/phase4g1-env/bin/python"
ENV = {
    **os.environ,
    "PYTHONPATH": ".:src",
    "OPENMM_CPU_THREADS": "1",
    "OMP_NUM_THREADS": "1",
}
C1 = ROOT / "docs/evidence/phase_4j38_c1/run"


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value))


def cli(module, *args):
    return subprocess.run(
        [str(PYTHON), "-m", module, *map(str, args)],
        cwd=ROOT,
        env=ENV,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    root = tmp_path_factory.mktemp("j38-c2") / "output"
    result = cli("scripts.validate_phase_4j38_c2", root)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "acceptance unmet" in result.stdout
    return root


def test_native_source_chemically_incompatible_not_ionic():
    value = aromatic_decision(
        read(C1 / "aromatic-audit.json"), read(C1 / "aromatic-fixture.json")
    )
    assert value["native_system_total_e"] == "-0.460"
    assert value["formal_system_total_e"] == 0
    assert [c["residual_e"] for c in value["components"]] == ["-0.230"] * 2
    assert all(
        c["chemically_intended_ionic_repeat"] is False for c in value["components"]
    )
    assert (
        value["native_path_allowed"]
        is value["normalization"]
        is value["redistribution"]
        is False
    )
    with pytest.raises(SourceIncompatibleError, match="neutral C12H10"):
        require_native_source(value)


@pytest.mark.parametrize(
    "mutation",
    [
        "zero_charge",
        "rounded",
        "boolean",
        "nonfinite",
        "order",
        "component",
        "molecule",
        "formal",
        "units",
        "graph",
        "normalize",
    ],
)
def test_native_source_changes_rejected(mutation):
    audit, evidence = (
        read(C1 / "aromatic-audit.json"),
        read(C1 / "aromatic-fixture.json"),
    )
    if mutation == "zero_charge":
        audit["rows"][3]["charge_e"] = 0.0
        audit["rows"][3]["source_charge_decimal"] = "0.0"
    elif mutation == "rounded":
        audit["rows"][0]["charge_e"] = -0.12
    elif mutation == "boolean":
        audit["rows"][0]["charge_e"] = False
    elif mutation == "nonfinite":
        audit["rows"][0]["charge_e"] = float("nan")
    elif mutation == "order":
        audit["rows"].reverse()
    elif mutation == "component":
        audit["rows"][0]["component"] = 1
    elif mutation == "molecule":
        audit["rows"][0]["molecule"] = "other"
    elif mutation == "formal":
        evidence["final_graph"]["sites"][0]["formal_charge"] = -1
    elif mutation == "units":
        audit["units"] = "coulomb"
    elif mutation == "graph":
        audit["final_graph_identity"] = "0" * 64
    else:
        audit["normalized"] = True
    with pytest.raises(ValueError):
        aromatic_decision(audit, evidence)


@pytest.mark.parametrize(
    "origin", ["external_provided", "external_oplsaa", "RESP/Psi4", "AM1-BCC"]
)
def test_no_foreign_origin_can_be_native_opls(origin):
    value = aromatic_decision(
        read(C1 / "aromatic-audit.json"), read(C1 / "aromatic-fixture.json")
    )
    value["charge_origin"] = origin
    value = _signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(ValueError, match="External charges"):
        require_native_source(value)


def test_precise_matrix_origins_and_blocked_gates(bundle):
    matrix = read(bundle / "matrix.json")
    assert len(matrix["cases"]) == 8
    assert matrix["new_md_executions"] == 0
    for key, row in matrix["cases"].items():
        assert row["new_dynamics_executed"] is False
        assert row["charge_paths"]["externally_provided"]["enabled"] is False
        assert row["charge_paths"]["externally_provided"]["status"] == "not_selected"
        assert row["stages"]["scientific_validation"]["status"] == "blocked"
        if key == "oplsaa/aromatic":
            assert row["selected_charge_origin"] == "native_oplsaa"
            assert row["charge_source_status"] == "source_incompatible"
        elif key in {"gaff2/pe", "gaff2/aromatic"}:
            assert row["charge_source_status"] == "asset_unavailable"
            assert (
                row["charge_paths"]["resp_psi4_exact_graph"]["exact_graph_verified"]
                is False
            )
            assert (
                row["stages"]["parameter_completeness"]["status"]
                == row["stages"]["executable_backend"]["status"]
                == "not_executed"
            )
        elif key.startswith("gaff2/"):
            assert row["charge_paths"]["resp_psi4_exact_graph"]["review_identity"]
    assert not list(bundle.rglob("checkpoint*"))


@pytest.mark.parametrize(
    "mutation",
    [
        "origin",
        "external_as_native",
        "charge_status",
        "false_resp_match",
        "false_review",
        "wrong_graph",
        "dynamics",
        "parameter",
        "equilibrium",
        "acceptance",
        "inventory",
    ],
)
def test_resigned_matrix_mutations_rejected(bundle, mutation):
    matrix = deepcopy(read(bundle / "matrix.json"))
    row = matrix["cases"]["gaff2/pe"]
    if mutation == "origin":
        row["selected_charge_origin"] = "AM1-BCC"
    elif mutation == "external_as_native":
        matrix["cases"]["oplsaa/aromatic"]["charge_paths"]["externally_provided"][
            "enabled"
        ] = True
    elif mutation == "charge_status":
        row["charge_source_status"] = "retained_verified"
    elif mutation == "false_resp_match":
        row["charge_paths"]["resp_psi4_exact_graph"]["exact_graph_verified"] = True
    elif mutation == "false_review":
        row["charge_paths"]["resp_psi4_exact_graph"]["review_identity"] = "1" * 64
    elif mutation == "wrong_graph":
        row["final_graph_identity"] = matrix["cases"]["gaff2/peo"][
            "final_graph_identity"
        ]
    elif mutation == "dynamics":
        row["new_dynamics_executed"] = True
    elif mutation == "parameter":
        row["stages"]["parameter_completeness"]["status"] = "passed"
    elif mutation == "equilibrium":
        row["stages"]["equilibration"]["status"] = "passed"
    elif mutation == "acceptance":
        matrix["acceptance_met"] = True
    else:
        del matrix["cases"]["oplsaa/aromatic"]
    matrix = _signed({k: v for k, v in matrix.items() if k != "identity"})
    with pytest.raises(ValueError):
        validate_matrix(
            matrix,
            read(bundle / "c1-recheck/matrix.json"),
            read(bundle / "native-source-decision.json"),
            read(bundle / "c1-recheck/resp-inventory.json"),
        )


@pytest.mark.parametrize(
    "mutation", ["source_hash", "residual", "source_compatible", "origin", "graph"]
)
def test_signed_source_decision_mutations_rejected_fresh_process(
    bundle, tmp_path, mutation
):
    root = tmp_path / "changed"
    shutil.copytree(bundle, root)
    value = read(root / "native-source-decision.json")
    if mutation == "source_hash":
        value["source_sha256"] = "0" * 64
    elif mutation == "residual":
        value["components"][0]["residual_e"] = "0"
    elif mutation == "source_compatible":
        value["source_compatible"] = True
    elif mutation == "origin":
        value["charge_origin"] = "external_provided"
    else:
        value["final_graph_identity"] = "0" * 64
    value = _signed({k: v for k, v in value.items() if k != "identity"})
    write(root / "native-source-decision.json", value)
    result = cli(
        "scripts.validate_phase_4j38_c2_receipt", root, "--bundle", "--diagnostic-only"
    )
    assert result.returncode != 0
    assert "Changed native source compatibility" in result.stderr
    assert (root / "native-rejection.json").exists()


def test_relocation_and_child_reconstruction_are_diagnostic_only(bundle, tmp_path):
    root = tmp_path / "relocated"
    shutil.copytree(bundle, root)
    result = cli(
        "scripts.validate_phase_4j38_c2_receipt", root, "--bundle", "--diagnostic-only"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    result = cli("scripts.validate_phase_4j38_c2_receipt", root, "--bundle")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "acceptance=False" in result.stdout
