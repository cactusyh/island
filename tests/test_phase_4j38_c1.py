"""Real raw-source negative controls; blocked fixtures are never positive MD evidence."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt import _signed
from scripts.phase_4j38_c1_audit import components
from scripts.phase_4j38_c1_contract import FAILED, STAGES, validate

ROOT = Path(__file__).resolve().parents[1]
PYTHON = ROOT.parent / "island-validation/phase4g1-env/bin/python"
ENV = {
    **os.environ,
    "PYTHONPATH": ".:src",
    "OPENMM_CPU_THREADS": "1",
    "OMP_NUM_THREADS": "1",
}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value))


def cli(module, *args):
    return subprocess.run(
        [str(PYTHON), "-m", module, *map(str, args)],
        cwd=ROOT,
        env=ENV,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    path = tmp_path_factory.mktemp("j38-c1") / "diagnostic"
    result = cli("scripts.validate_phase_4j38_c1", path)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "acceptance unmet" in result.stdout
    return path


def test_independent_raw_source_reproduction(bundle):
    audit = read(bundle / "aromatic-audit.json")
    assert [c["source_total_decimal_e"] for c in audit["components"]] == [
        "-0.230",
        "-0.230",
    ]
    assert [c["formal_total_e"] for c in audit["components"]] == [0, 0]
    assert audit["total_decimal_e"] == "-0.460"
    assert audit["normalized"] is audit["redistributed"] is False
    assert len(audit["rows"]) == 44
    assert {r["source_type"] for r in audit["rows"] if r["element"] == "C"} == {
        "opls_145"
    }
    assert [c["carbon_without_hydrogen"] for c in audit["components"]] == [
        [4, 7],
        [26, 29],
    ]
    native = read(bundle / "native-rejection.json")
    assert native["diagnostic"] == "Native component neutrality failed: -0.23 e"
    assert native["dynamics_executed"] is False


def test_exact_resp_assets_unavailable(bundle):
    report = read(bundle / "resp-inventory.json")
    assert set(report["requested"]) == {"gaff2/pe", "gaff2/aromatic"}
    retained = {r["final_graph_identity"] for r in report["inventory"]}
    for cell in report["requested"].values():
        assert cell["status"] == "asset_unavailable"
        assert cell["final_graph_identity"] not in retained
        assert cell["exact_graph_match"] is cell["dynamics_executed"] is False
    assert not list(bundle.rglob("checkpoint*"))
    assert read(bundle / "run-state.json")["acceptance_met"] is False


def test_nine_stages_and_no_inherited_new_execution(bundle):
    matrix = read(bundle / "matrix.json")
    assert len(matrix["cases"]) == 8
    for key, row in matrix["cases"].items():
        assert set(row["stages"]) == set(STAGES)
        assert row["new_dynamics_executed"] is False
        assert row["execution_origin"] == (
            "not_executed" if key in FAILED else "retained_j38"
        )
        if key in FAILED:
            assert row["stages"]["executable_backend"]["status"] == "not_executed"
            assert row["stages"]["parameter_completeness"]["status"] == "not_executed"
    assert (
        matrix["cases"]["oplsaa/aromatic"]["stages"]["charge_completeness"]["status"]
        == "source_incomplete"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_cell",
        "key_identity",
        "inherited_md",
        "charge_complete",
        "equilibrium",
        "parameter_complete",
        "wrong_graph",
        "wrong_packing",
        "readiness",
        "all_success",
    ],
)
def test_resigned_matrix_cannot_hide_failed_stages(bundle, mutation):
    value = read(bundle / "matrix.json")
    row = value["cases"]["oplsaa/aromatic"]
    if mutation == "missing_cell":
        del value["cases"]["gaff2/pe"]
    elif mutation == "key_identity":
        row["case"] = "gaff2/aromatic"
    elif mutation == "inherited_md":
        row["new_dynamics_executed"] = True
    elif mutation in {"charge_complete", "parameter_complete", "equilibrium"}:
        stage = {
            "charge_complete": "charge_completeness",
            "parameter_complete": "parameter_completeness",
            "equilibrium": "equilibration",
        }[mutation]
        row["stages"][stage]["status"] = "passed"
    elif mutation == "wrong_graph":
        row["final_graph_identity"] = "0" * 64
    elif mutation == "wrong_packing":
        row["packing_identity"] = "0" * 64
    elif mutation == "readiness":
        value["production_validated"] = True
    else:
        for r in value["cases"].values():
            for stage in r["stages"].values():
                stage["status"] = "passed"
    value = _signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(ValueError):
        validate(
            value,
            read(ROOT / "docs/evidence/phase_4j38/matrix.json"),
            read(bundle / "aromatic-audit.json"),
            read(bundle / "resp-inventory.json"),
        )


@pytest.mark.parametrize(
    "mutation",
    [
        "charge",
        "rounding",
        "order",
        "component",
        "source",
        "units",
        "formal",
        "provenance",
        "executable",
        "wrong_resp_graph",
        "wrong_resp_match",
        "resp_dynamics",
    ],
)
def test_raw_audit_and_resp_mutations_rejected_in_fresh_process(
    bundle, tmp_path, mutation
):
    path = tmp_path / "changed"
    shutil.copytree(bundle, path)
    if mutation.startswith("wrong_resp") or mutation == "resp_dynamics":
        target = path / "resp-inventory.json"
        report = read(target)
        cell = report["requested"]["gaff2/pe"]
        if mutation == "wrong_resp_graph":
            cell["final_graph_identity"] = report["inventory"][0][
                "final_graph_identity"
            ]
        elif mutation == "wrong_resp_match":
            cell["exact_graph_match"] = True
        else:
            cell["dynamics_executed"] = True
    else:
        target = path / "aromatic-audit.json"
        report = read(target)
        if mutation == "charge":
            report["rows"][0]["charge_e"] += 0.23
        elif mutation == "rounding":
            report["rows"][0]["source_charge_decimal"] = "-0.12"
        elif mutation == "order":
            report["rows"].reverse()
        elif mutation == "component":
            report["rows"][0]["component"] = 1
        elif mutation == "source":
            report["source_sha256"] = "0" * 64
        elif mutation == "units":
            report["units"] = "coulomb"
        elif mutation == "formal":
            report["components"][0]["formal_total_e"] = -1
        elif mutation == "provenance":
            report["rows"][0]["metadata"]["repeat_unit_index"] = 1
        else:
            report["runtime"]["executable_sha256"] = "0" * 64
    write(target, report)
    result = cli(
        "scripts.validate_phase_4j38_c1_receipt", path, "--bundle", "--diagnostic-only"
    )
    assert result.returncode != 0
    assert "Changed" in result.stderr
    assert target.exists()  # Diagnostic retained on rejection.


def test_relocation_fresh_process_and_nonzero_acceptance(bundle, tmp_path):
    moved = tmp_path / "new-path"
    shutil.copytree(bundle, moved)
    result = cli(
        "scripts.validate_phase_4j38_c1_receipt", moved, "--bundle", "--diagnostic-only"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    result = cli("scripts.validate_phase_4j38_c1_receipt", moved, "--bundle")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "acceptance remains unmet" in result.stdout


@pytest.mark.parametrize(
    "mutation", ["duplicate", "reorder", "component", "membership", "bond"]
)
def test_independent_connectivity_controls(bundle, mutation):
    graph = deepcopy(read(bundle / "aromatic-fixture.json")["final_graph"])
    if mutation == "duplicate":
        graph["sites"].append(graph["sites"][0])
    elif mutation == "reorder":
        graph["sites"].reverse()
    elif mutation == "component":
        graph["components"][0].pop()
    elif mutation == "membership":
        graph["component_membership"]["1"] = 1
    else:
        graph["bonds"].append(graph["bonds"][0])
    with pytest.raises(ValueError):
        components(graph)


@pytest.mark.parametrize(
    "mutation",
    [
        "atom_id",
        "molecule",
        "end_group",
        "coordinates",
        "formal",
        "packing",
        "repeat",
        "omitted_atom",
    ],
)
def test_exact_final_graph_controls(bundle, tmp_path, mutation):
    path = tmp_path / "changed-graph"
    shutil.copytree(bundle, path)
    target = path / "aromatic-fixture.json"
    evidence = read(target)
    graph = evidence["final_graph"]
    if mutation == "atom_id":
        graph["sites"][0]["id"] = 99
    elif mutation == "molecule":
        graph["molecule_membership"]["1"] = "other-chain"
    elif mutation == "end_group":
        evidence["endpoints"][0][0] = 4
    elif mutation == "coordinates":
        evidence["packed_system"]["coordinates"]["1"][0] += 0.01
    elif mutation == "formal":
        graph["sites"][0]["formal_charge"] = 1
    elif mutation == "packing":
        evidence["packing_plan"]["site_id_remapping"][0]["1"] = 2
    elif mutation == "repeat":
        graph["sites"][0]["metadata"]["repeat_unit_index"] = 1
    else:
        graph["sites"].pop()
    write(target, evidence)
    result = cli(
        "scripts.validate_phase_4j38_c1_receipt", path, "--bundle", "--diagnostic-only"
    )
    assert result.returncode != 0
    assert "Changed exact graph" in result.stderr
