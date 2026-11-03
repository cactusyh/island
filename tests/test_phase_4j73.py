"""Evidence-only tests: every executable/import probe uses marker-only fake code."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import phase_4j72_c1_fence as fence
from scripts.phase_4j73_source_intake import (
    COUNTS,
    INPUTS,
    ROOT,
    ROW_STATES,
    SCOPE,
    IntakeError,
    artifact_contract,
    artifact_preflight,
    exact,
    parents,
    pin,
    read,
    reconstruct,
    render,
    replay,
    signed,
    validate_requests,
    write,
)
from scripts.validate_phase_4j73_receipt import semantic

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def expected():
    return reconstruct(REPO, REPO / ROOT)


def resign(value):
    return signed({k: v for k, v in value.items() if k != "identity"})


def fake(path, marker):
    path.write_text(
        "#!/usr/bin/env python3\n# ISLAND_C1_SYNTHETIC_ONLY\n"
        + f"open({str(marker)!r},'w').write('fake-only')\n"
    )
    path.chmod(0o755)
    return path


def child(code, *arguments, **kwargs):
    return subprocess.run(
        [sys.executable, "-c", code, *map(str, arguments)],
        capture_output=True,
        text=True,
        check=False,
        **kwargs,
    )


def test_exact_378_coverage_and_unchanged_statuses(expected):
    p = expected["source-intake.json"]
    validate_requests(p, p)
    assert p["unresolved_counts"] == COUNTS
    assert len({r["id"] for r in p["requests"]}) == 378
    assert p["family_counts"] == {
        "angle-angle": 78,
        **{
            f: 60
            for f in (
                "angle-angle-torsion_1",
                "angle-torsion_3",
                "bond-bond_1_3",
                "end_bond-torsion_3",
                "middle_bond-torsion_3",
            )
        },
    }
    for r in p["requests"]:
        assert r["statuses"] == ROW_STATES and not r["source_records"]
        assert [s["id"] for s in r["chemical_environment"]] == r["stable_sites"]
        assert [s["type"] for s in r["chemical_environment"]] == r["ordered_types"]
        assert r["existing_native_source_lookup"]["status"] == "missing"
        assert len(r["existing_native_source_lookup"]["source_sha256"]) == 64
        assert set(r["interaction_context"]) == {
            "bond",
            "angle",
            "torsion",
            "improper",
            "Class_II",
        }
    groups = expected["source-artifact-requests.json"]["groups"]
    ids = [i for g in groups for i in g["instance_ids"]]
    assert len(ids) == len(set(ids)) == 378
    assert set(ids) == {r["id"] for r in p["requests"]}


@pytest.mark.parametrize(
    "operation",
    [
        "missing",
        "duplicate_id",
        "duplicate_identity",
        "hide_conflict",
        "analog_as_native",
        "unresolved_convention",
        "equilibrium",
        "unit",
        "source_hash",
        "site",
        "profile",
    ],
)
def test_resigned_row_mutations_rejected(expected, operation):
    p = copy.deepcopy(expected["source-intake.json"])
    row = p["requests"][0]
    if operation == "missing":
        p["requests"].pop()
    elif operation == "duplicate_id":
        p["requests"][1]["id"] = row["id"]
        p["requests"][1] = resign(p["requests"][1])
    elif operation == "duplicate_identity":
        p["requests"][1] = copy.deepcopy(row)
    elif operation == "hide_conflict":
        row["unresolved_alternatives"].pop()
    elif operation == "analog_as_native":
        row["source_records"] = [
            {"native_PCFF": True, "origin": "analog", "sha256": "1" * 64}
        ]
    elif operation == "unresolved_convention":
        row["statuses"]["convention_resolved"] = True
        row["required_evidence"]["equation_and_units"][
            "historical_unresolved_contract"
        ]["phase"] = None
    elif operation == "equilibrium":
        row["required_evidence"]["equilibria"]["historical_dependency_context_only"][0][
            "equilibrium"
        ] += 0.01
    elif operation == "unit":
        row["required_evidence"]["equation_and_units"][
            "historical_unresolved_contract"
        ]["source_units"] = ["unknown"]
    elif operation == "source_hash":
        row["existing_native_source_lookup"]["source_sha256"] = "0" * 64
    elif operation == "site":
        row["ordered_types"][0] = "analog_Si"
    else:
        p["no_execution_gate"]["profile"] = "premature_pcff"
    p["requests"][0] = resign(row)
    p = resign(p)
    with pytest.raises(IntakeError):
        validate_requests(p, expected["source-intake.json"])


@pytest.mark.parametrize("state", list(ROW_STATES))
def test_attestation_is_not_primary_row_authority(expected, state):
    p = copy.deepcopy(expected["source-intake.json"])
    r = p["requests"][0]
    r["statuses"][state] = True
    p["requests"][0] = resign(r)
    with pytest.raises(IntakeError):
        validate_requests(resign(p), expected["source-intake.json"])
    statuses = expected["status-separation.json"]
    assert statuses["author_confirmed_reported"]["count"] == 28
    assert statuses["author_confirmed_audited"]["count"] == 0
    assert statuses["primary_artifact_supplied"]["count"] == 0


def submission(tmp_path):
    path = tmp_path / "synthetic-author-file.txt"
    path.write_text(
        "ISLAND_C1_SYNTHETIC_ONLY: fabricated fixture, no scientific values\n"
    )
    c = artifact_contract()
    return {
        "original_filename": path.name,
        "sha256": pin(path),
        "citation": {k: "synthetic-only" for k in c["required_citation_fields"]},
        "license_status": "synthetic",
        "redistribution_status": "synthetic",
        "origin": "original_author_file",
        "provenance_chain": "self-asserted synthetic fixture",
        "source_version": "synthetic-v1",
        "locations": [{k: "synthetic location" for k in c["required_location_fields"]}],
        "target_ids": ["synthetic-control-only"],
    }


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_hash",
        "wrong_hash",
        "renamed_file",
        "missing_location",
        "analog_origin",
        "complete_fabricated_file",
    ],
)
def test_fabricated_primary_artifact_never_bootstraps_authority(tmp_path, mutation):
    r = submission(tmp_path)
    if mutation == "missing_hash":
        r["sha256"] = None
    elif mutation == "wrong_hash":
        r["sha256"] = "0" * 64
    elif mutation == "renamed_file":
        r["original_filename"] = "changed.txt"
    elif mutation == "missing_location":
        r["locations"][0]["table"] = None
    elif mutation == "analog_origin":
        r["origin"] = "COMPASS_analog"
    with pytest.raises(IntakeError):
        artifact_preflight(r, tmp_path)


def receipt_pair():
    # A synthetic envelope is enough to isolate semantic checks. Full real
    # receipts additionally reconstruct all inputs and execution ledgers.
    r = signed(
        {
            "schema": "island.j73-receipt.v1",
            "base_commit": "73eb64e02ab5d970e02f444e1b41c8183c49034d",
            "branch": "codex/phase-4j73-authoritative-pcff-source-intake",
            "scope": SCOPE.copy(),
            "manifest_sha256": "1" * 64,
            "historical_files_preserved": 0,
            "verification": {},
            "offline_reconstruction": {
                "offline": True,
                "instances": 378,
                "unresolved_counts": COUNTS,
                "scope": SCOPE.copy(),
            },
        }
    )
    m = signed(
        {
            "schema": "island.j73-manifest.v1",
            "base_commit": r["base_commit"],
            "scope": SCOPE.copy(),
            "implementation_files": {},
            "artifacts": {},
            "frozen_inputs_sha256": "2" * 64,
        }
    )
    return r, m


@pytest.mark.parametrize("key", list(SCOPE))
@pytest.mark.parametrize("where", ["receipt", "manifest", "nested"])
def test_resigned_scope_semantics_never_erase_violation_or_enable_execution(
    tmp_path, key, where
):
    r, m = receipt_pair()
    target = (
        m["scope"]
        if where == "manifest"
        else r["offline_reconstruction"]["scope"]
        if where == "nested"
        else r["scope"]
    )
    value = target[key]
    target[key] = (
        not value
        if type(value) is bool
        else value + 1
        if type(value) is int
        else "promoted"
    )
    m = resign(m)
    write(tmp_path / "manifest.json", m)
    r["manifest_sha256"] = pin(tmp_path / "manifest.json")
    r = resign(r)
    with pytest.raises(IntakeError):
        semantic(r, m)


def test_integer_counts_are_not_booleans():
    with pytest.raises(IntakeError):
        exact({"count": False}, {"count": 0}, "typed_counts")


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "antechamber",
        "sander",
        "pmemd",
        "lammps",
        "msi2lmp",
        "fit_hmds",
        "parameter_generation",
        "reserved_validation",
    ],
)
def test_fake_scientific_processes_never_launch(tmp_path, name):
    marker = tmp_path / "ran"
    exe = fake(tmp_path / name, marker)
    argv = [str(exe), "original-argument"]
    with pytest.raises(fence.ScopeDenied):
        subprocess.run(argv, check=False)
    assert argv == [str(exe), "original-argument"] and not marker.exists()
    result = child(
        "import subprocess,sys\ntry: subprocess.run([sys.argv[1]])\nexcept RuntimeError: pass\nelse: raise AssertionError('launched')",
        exe,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


@pytest.mark.parametrize(
    "name", ["fit_hmds", "parameter_generation", "reserved_validation"]
)
def test_python_scientific_driver_rejected_before_child(tmp_path, name):
    marker = tmp_path / "ran"
    driver = fake(tmp_path / (name + ".py"), marker)
    with pytest.raises(fence.ScopeDenied, match="J73_scientific_driver_forbidden"):
        subprocess.run([sys.executable, str(driver)], check=False)
    assert not marker.exists()


@pytest.mark.parametrize("name", ["openmm", "lammps", "psi4"])
def test_fake_module_discovery_allowed_execution_blocked(tmp_path, name):
    marker = tmp_path / "ran"
    fake(tmp_path / (name + ".py"), marker)
    result = child(
        "import sys,importlib.util\nfrom scripts.phase_4j73_scope import install\ninstall()\nsys.path.insert(0,sys.argv[1])\nassert importlib.util.find_spec(sys.argv[2])\ntry: __import__(sys.argv[2])\nexcept RuntimeError: pass\nelse: raise AssertionError('executed')",
        tmp_path,
        name,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


def test_no_fence_stops_before_source_access(monkeypatch):
    monkeypatch.setattr(
        fence,
        "assert_installed",
        lambda: (_ for _ in ()).throw(RuntimeError("C1 fence absent")),
    )
    with pytest.raises(RuntimeError, match="C1 fence absent"):
        parents(Path("nonexistent-input"), Path("nonexistent-evidence"))


def test_historical_source_and_freeze_resigning_rejected(tmp_path):
    freeze = read(REPO / ROOT / "frozen-inputs.json")
    freeze["inputs"][INPUTS[0]] = "0" * 64
    write(tmp_path / "frozen-inputs.json", resign(freeze))
    with pytest.raises(IntakeError, match="frozen_parent_contract_mutation"):
        parents(REPO, tmp_path)
    shutil.copyfile(REPO / ROOT / "frozen-inputs.json", tmp_path / "frozen-inputs.json")
    for name in INPUTS:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    source = read(REPO / INPUTS[0])
    source["source_complete"] = True
    write(tmp_path / INPUTS[0], resign(source))
    with pytest.raises(IntakeError, match="historical_source_hash_mutation"):
        parents(tmp_path, tmp_path)


def test_duplicate_json_keys_fail_closed(tmp_path):
    path = tmp_path / "duplicates.json"
    path.write_text('{"source_id":"first","source_id":"hidden-conflict"}')
    with pytest.raises(IntakeError, match="duplicate_JSON_key"):
        read(path)


def test_relocation_and_offline_replay_with_no_scientific_modules(tmp_path):
    relocated = tmp_path / "repository"
    (relocated / ROOT).mkdir(parents=True)
    shutil.copytree(
        REPO / "scripts",
        relocated / "scripts",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    for name in INPUTS:
        target = relocated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    for name in (
        "frozen-inputs.json",
        "historical-files.json",
        "source-intake.json",
        "artifact-contract.json",
        "source-artifact-requests.json",
        "status-separation.json",
        "no-execution-gate.json",
        "evidence-request-report.md",
    ):
        shutil.copyfile(REPO / ROOT / name, relocated / ROOT / name)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j73_receipt",
            str(ROOT),
            "--replay",
        ],
        cwd=relocated,
        env={**os.environ, "PYTHONPATH": str(relocated)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert json.loads(result.stdout)["instances"] == 378
    assert not any(n in sys.modules for n in ("openmm", "lammps", "psi4"))


def test_deterministic_report_and_replay(expected):
    assert render(expected) == (REPO / ROOT / "evidence-request-report.md").read_text()
    assert replay(REPO, REPO / ROOT)["unresolved_counts"] == COUNTS


@pytest.mark.parametrize(
    "control",
    [
        "QM",
        "MD",
        "OpenMM",
        "LAMMPS",
        "AmberTools",
        "msi2lmp_converter",
        "fitting",
        "parameter_generation",
        "reserved_validation",
    ],
)
@pytest.mark.skip(
    reason="J73 evidence-only allocation: real scientific execution explicitly forbidden; synthetic rejection controls run instead"
)
def test_live_scientific_control_scope_skip(control):
    raise AssertionError("must never reach a live scientific fixture")
