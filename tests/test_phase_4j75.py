"""Offline evidence tests. Process probes use marker-only synthetic files."""

import copy
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import phase_4j72_c1_fence as fence
from scripts import phase_4j75_source_path as j75
from scripts.validate_phase_4j75_receipt import FILES, artifacts

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / j75.ROOT


@pytest.fixture(scope="module", autouse=True)
def evidence_scope():
    j75.scope()


@pytest.fixture(scope="module")
def expected():
    return j75.derive(REPO, ROOT)


def resign(document):
    return j75.signed({k: v for k, v in document.items() if k != "identity"})


def test_exact_declared_absolute_path_and_frozen_existence():
    j75.check_declared(j75.AUTHORITATIVE)
    path = j75.read(ROOT / "source-path-report.json")
    assert path["requested_absolute_path"] == j75.AUTHORITATIVE
    assert path["path_exists"] and path["directory_status"]
    assert path["resolved_realpath"] == j75.AUTHORITATIVE
    assert path["inside_git_checkout"] and not path["is_symlink"]
    assert path["fence_installed_before_source_inspection"]
    assert not path["fallback_used"]


@pytest.mark.parametrize(
    "path",
    [
        "/island_ffref/pcff",
        "island_ffref/pcff",
        j75.AUTHORITATIVE + "/",
        j75.AUTHORITATIVE + "/../pcff",
        "/tmp/pcff",
    ],
)
def test_incorrect_path_rejected_before_inspection(monkeypatch, path):
    def forbidden(*args, **kwargs):
        pytest.fail("an incorrect source path was inspected")

    monkeypatch.setattr(Path, "exists", forbidden)
    with pytest.raises(j75.IntakeError, match="incorrect_authoritative"):
        j75.path_report(path, REPO)


def test_absent_exact_path_creates_no_source_facts_or_fallback(tmp_path, monkeypatch):
    absent = str(tmp_path / "synthetic-unsupplied")
    monkeypatch.setattr(j75, "AUTHORITATIVE", absent)
    output = tmp_path / "evidence"
    j75.capture(REPO, output, absent)
    report = j75.read(output / "source-path-report.json")
    assert not report["path_exists"] and not report["fallback_used"]
    assert (
        j75.read(output / "frozen-inputs.json")["decision"] == "source_path_unavailable"
    )
    assert j75.read(output / "source-root-manifest.json")["files"] == []
    assert j75.read(output / "source-file-inventory.json")["files"] == []
    assert j75.read(output / "extracted-parameter-records.json")["records"] == []


def test_symlink_identity_is_reported_and_never_parsed(tmp_path, monkeypatch):
    real = tmp_path / "synthetic-source"
    real.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(real, target_is_directory=True)
    monkeypatch.setattr(j75, "AUTHORITATIVE", str(alias))
    report = j75.path_report(str(alias), REPO)
    assert report["is_symlink"] and not report["root_identity_usable"]
    assert report["resolved_realpath"] == str(real)
    with pytest.raises(j75.IntakeError, match="identity_blocked"):
        j75.capture(REPO, tmp_path / "evidence", str(alias))
    assert not (tmp_path / "evidence/source-root-manifest.json").exists()


def test_deterministic_tree_permissions_and_hash_mutation(tmp_path):
    (tmp_path / "nested").mkdir()
    (tmp_path / "nested/a.txt").write_text("synthetic a")
    (tmp_path / "z.txt").write_text("synthetic z")
    first = j75.tree_inventory(tmp_path)
    assert first == j75.tree_inventory(tmp_path)
    assert [r["relative_path"] for r in first] == ["nested/a.txt", "z.txt"]
    assert all(r["readable"] and r["permissions"] for r in first)
    (tmp_path / "nested/a.txt").write_text("changed synthetic")
    assert j75.tree_inventory(tmp_path)[0]["sha256"] != first[0]["sha256"]
    (tmp_path / "alias").symlink_to("nested", target_is_directory=True)
    with pytest.raises(j75.IntakeError, match="symlink"):
        j75.tree_inventory(tmp_path)


def test_root_manifest_complete_and_frozen_before_parsing():
    manifest = j75.read(ROOT / "source-root-manifest.json")
    assert manifest["recursive_file_count"] == len(manifest["files"]) == 27
    assert manifest["total_bytes"] == sum(f["bytes"] for f in manifest["files"])
    assert manifest["frozen_before_format_parsing"]
    assert manifest["no_symlinks_followed"]


@pytest.mark.parametrize("field", ["sha256", "bytes"])
def test_checkout_external_mismatch_not_hidden(field):
    current = [{"relative_path": "synthetic", "sha256": "1" * 64, "bytes": 2}]
    old = copy.deepcopy(current)
    old[0][field] = "2" * 64 if field == "sha256" else 3
    report = j75.read(ROOT / "source-path-report.json")
    result = j75.compare_copies(current, old, report, j75.AUTHORITATIVE)
    assert not result["files"][0]["exact_duplicate"]
    assert not result["files"][0]["source_authority_inferred"]


def test_actual_copy_comparison_is_same_root_not_independent_evidence(expected):
    result = expected["checkout-comparison.json"]
    assert all(r["exact_duplicate"] for r in result["files"])
    assert result["source_root_identity"]["J74_recorded_root_same_path"]
    assert not result["source_root_identity"]["independent_external_root"]
    assert not result["hash_equality_proves_authorship"]


@pytest.mark.parametrize(
    "name,historical,family",
    [
        ("pcff.frc", None, 0),
        ("cff91.frc", None, 0),
        ("compass_published.frc", None, 1),
        ("PDMS_benchmark-main/param/COMPASS/file", None, 2),
        ("PDMS_benchmark-main/param/OPLS/file", None, 2),
        ("PDMS_benchmark-main/param/HUANG/file", None, 3),
        ("paper.pdf", {"force_field": "Smith2004_PDMS"}, 4),
        (
            "paper.pdf",
            {
                "authors": "Sun and Rigby",
                "force_field": "PCFF",
                "title": "Polysiloxanes: force field",
            },
            5,
        ),
        (
            "paper.pdf",
            {
                "authors": "Sun and Rigby",
                "force_field": "PCFF_related",
                "title": "Poly(ethylene oxide)",
            },
            6,
        ),
        (
            "paper.pdf",
            {"authors": "Sun", "force_field": "CFF93", "title": "Polycarbonates"},
            6,
        ),
        ("paper.pdf", {"force_field": "PCFF_COMPASS_comparison"}, 6),
        (
            "paper.pdf",
            {
                "force_field": "multiple_benchmark",
                "discussed_models": ["Huang", "COMPASS"],
            },
            6,
        ),
        ("analog.pdf", None, 6),
    ],
)
def test_source_families_never_merged(name, historical, family):
    assert j75.source_family(name, historical) == j75.FAMILIES[family]


@pytest.mark.parametrize(
    "mutation", ["duplicate", "hash", "license", "author", "COMPASS"]
)
def test_inventory_identity_license_authority_rejection(mutation):
    inventory = j75.read(ROOT / "source-file-inventory.json")
    s = next(
        s for s in inventory["files"] if s["relative_path"] == "compass_published.frc"
    )
    if mutation == "duplicate":
        inventory["files"].append(copy.deepcopy(s))
    elif mutation == "hash":
        s["sha256"] = ""
    elif mutation == "license":
        del s["license_status"]
    elif mutation == "author":
        s["primary_artifact_supplied"] = True
    else:
        s["source_family"] = j75.FAMILIES[0]
    inventory["files"] = [resign(s) for s in inventory["files"]]
    with pytest.raises(j75.IntakeError):
        j75.validate_inventory(resign(inventory))


@pytest.mark.parametrize(
    "mutation",
    ["sha", "family", "equation", "unit", "convention", "zero", "record_id", "root"],
)
def test_resigned_fabricated_source_record_rejected(tmp_path, mutation):
    freeze = j75.read(ROOT / "frozen-inputs.json")
    shutil.copyfile(ROOT / "frozen-inputs.json", tmp_path / "frozen-inputs.json")
    for name in freeze["capture_files"]:
        shutil.copyfile(ROOT / name, tmp_path / name)
    facts = j75.read(tmp_path / "extracted-parameter-records.json")
    r = facts["records"][0]
    key, value = {
        "sha": ("source_sha256", "0" * 64),
        "family": ("source_family", j75.FAMILIES[0]),
        "equation": (
            "analytic_equation",
            {"status": "complete_source_contract", "scalar": "fabricated"},
        ),
        "unit": (
            "units",
            {"status": "complete_source_contract", "symbols": "invented"},
        ),
        "convention": ("convention_resolved", True),
        "zero": ("coefficients_raw", ["0"]),
        "record_id": ("record_id", "self-signed-source"),
        "root": ("source_root", "/island_ffref/pcff"),
    }[mutation]
    r[key] = value
    facts["records"][0] = resign(r)
    j75.write(tmp_path / "extracted-parameter-records.json", resign(facts))
    # Even recomputing the outer freeze cannot confer source authority.
    freeze["capture_files"]["extracted-parameter-records.json"] = j75.pin(
        tmp_path / "extracted-parameter-records.json"
    )
    j75.write(tmp_path / "frozen-inputs.json", resign(freeze))
    with pytest.raises(j75.IntakeError, match="resigned_capture_freeze"):
        j75.inputs(REPO, tmp_path)


def test_missing_conventions_and_cross_terms_cannot_resolve_rows(expected):
    resolutions = expected["row-resolution-matrix.json"]
    assert len({r["id"] for r in resolutions["instances"]}) == 378
    assert resolutions["counts"] == {
        "analog_only": 294,
        "convention_blocked": 60,
        "source_missing": 24,
    }
    for r in resolutions["instances"]:
        assert not any(q["matches"] for q in r["native_source_queries"])
        assert not r["statuses"]["exact_row_resolved"]
        assert not r["statuses"]["convention_resolved"]
        assert not r["statuses"]["executable_bundle_qualified"]
        for c in r["candidate_source_records"]:
            assert not c["eligibility"]
            assert not c["exclusion_match"] and not c["charge_convention_match"]


@pytest.mark.parametrize(
    "mutation",
    ["omit", "duplicate", "promote", "profile", "phase", "candidate_duplicate"],
)
def test_resigned_coverage_and_profile_mutations_rejected(tmp_path, expected, mutation):
    # Replay reconstructs all scientific declarations from the pinned freeze.
    for name in (*j75.CAPTURE_FILES, "frozen-inputs.json"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    for name, document in expected.items():
        j75.write(tmp_path / name, document)
    rows = j75.read(tmp_path / "row-resolution-matrix.json")
    if mutation == "omit":
        rows["instances"].pop()
    elif mutation == "duplicate":
        rows["instances"].append(copy.deepcopy(rows["instances"][0]))
    elif mutation == "profile":
        rows["scope"]["automatic_assignment"] = True
    else:
        r = rows["instances"][0]
        if mutation == "promote":
            r["resolution_status"] = "exact_row_resolved"
        elif mutation == "phase":
            r["candidate_source_records"][0]["phase_multiplicity_match"] = True
        else:
            r["candidate_source_records"].append(
                copy.deepcopy(r["candidate_source_records"][0])
            )
        rows["instances"][0] = resign(r)
    j75.write(tmp_path / "row-resolution-matrix.json", resign(rows))
    with pytest.raises(j75.IntakeError, match="resigned_derived_evidence"):
        j75.replay(REPO, tmp_path)


def test_duplicate_candidate_rejected():
    with pytest.raises(j75.IntakeError, match="duplicate_candidate"):
        j75.candidate_set([{"record_id": "synthetic"}, {"record_id": "synthetic"}])


@pytest.mark.parametrize(
    "mutation", ["unit", "equation", "COMPASS", "Huang", "analog", "claimed_convention"]
)
def test_explicit_source_contract_mismatches_remain_ineligible(mutation):
    target = j75.read(REPO / j75.old.PARENTS[0])["requests"][0]
    record = copy.deepcopy(
        j75.read(ROOT / "extracted-parameter-records.json")["records"][0]
    )
    if mutation == "unit":
        record["units"] = {
            "status": "complete_source_contract",
            "symbols": {"E": "invented kJ"},
        }
        key = "unit_match"
    elif mutation == "equation":
        record["analytic_equation"] = {
            "status": "complete_source_contract",
            "scalar": "invented prefactor",
        }
        key = "function_form_match"
    elif mutation in {"COMPASS", "Huang", "analog"}:
        record["source_family"] = mutation
        record["source_model_identity"] = mutation
        key = "source_family_match"
    else:
        record["source_model_identity"] = "native_PCFF_CFF91_distribution_candidate"
        record["source_family"] = j75.FAMILIES[0]
        for field in (
            "phase_or_multiplicity",
            "orientation_order_convention",
            "equilibrium_dependencies",
            "applicable_exclusions",
            "charge_origin",
            "LJ_equation_and_mixing",
        ):
            record[field] = {"status": "complete_source_contract"}
        key = "convention_match"
    result = j75.match(record, target)
    assert not result[key] and not result["eligibility"]
    assert not result["target_convention_authority_resolved"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("evaluator_executions", 1),
        ("automatic_assignment", True),
        ("J72_zero_execution_scope_met", True),
        ("accepted_coefficients", 1),
    ],
)
def test_resigned_receipt_cannot_hide_scope_or_promote_profile(field, value):
    from scripts.validate_phase_4j75_receipt import semantic

    scope = copy.deepcopy(j75.SCOPE)
    scope[field] = value
    receipt = j75.signed(
        {
            "schema": "island.j75-receipt.v1",
            "base_commit": j75.BASE,
            "fetched_origin_main": j75.MAIN,
            "scope": scope,
        }
    )
    manifest = j75.signed(
        {
            "schema": "island.j75-manifest.v1",
            "base_commit": j75.BASE,
            "fetched_origin_main": j75.MAIN,
            "scope": j75.SCOPE,
        }
    )
    with pytest.raises(j75.IntakeError, match="resigned_scientific_scope"):
        semantic(receipt, manifest)


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "psi4",
        "gaussian",
        "antechamber",
        "resp",
        "lammps",
        "msi2lmp",
        "converter",
        "fit_parameters",
        "evaluate_energy",
        "reserved_validation",
    ],
)
def test_synthetic_native_name_cannot_execute(tmp_path, name):
    marker = tmp_path / "executed"
    fake = tmp_path / name
    fake.write_text(
        "#!/usr/bin/env python\nfrom pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic marker')\n"
    )
    fake.chmod(0o755)
    with pytest.raises((fence.ScopeDenied, PermissionError)):
        subprocess.run([str(fake)], check=True)
    assert not marker.exists()


@pytest.mark.parametrize(
    "driver",
    [
        "run_qm",
        "run_md",
        "run_converter",
        "evaluate_energy",
        "evaluate_forces",
        "evaluate_pressure",
        "generate_parameters",
        "optimize_coefficients",
        "reserved_validation",
    ],
)
def test_synthetic_python_driver_cannot_execute(tmp_path, driver):
    marker = tmp_path / "executed"
    fake = tmp_path / (driver + ".py")
    fake.write_text(
        "from pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic')\n"
    )
    with pytest.raises(fence.ScopeDenied):
        subprocess.run([sys.executable, str(fake)], check=True)
    assert not marker.exists()


def test_child_inherits_fence_before_synthetic_forbidden_launch(tmp_path):
    marker = tmp_path / "executed"
    fake = tmp_path / "sqm"
    fake.write_text("#!/usr/bin/env python\nraise RuntimeError('must never launch')\n")
    fake.chmod(0o755)
    code = (
        "from scripts.phase_4j75_source_path import scope; scope(); import subprocess; from scripts.phase_4j72_c1_fence import ScopeDenied,assert_installed; assert_installed();\ntry: subprocess.run(["
        + repr(str(fake))
        + "])\nexcept ScopeDenied: print('synthetic rejection before execution')\nelse: raise RuntimeError('fence bypass')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert "synthetic rejection before execution" in result.stdout
    assert not marker.exists()


def test_no_fence_stops_before_source_access(monkeypatch):
    def absent():
        raise RuntimeError("absent fence")

    monkeypatch.setattr(fence, "assert_installed", absent)
    with pytest.raises(RuntimeError, match="absent fence"):
        j75.path_report(j75.AUTHORITATIVE, REPO)


@pytest.mark.parametrize("module", ["psi4", "openmm", "lammps"])
def test_synthetic_module_discovery_allowed_execution_denied(
    tmp_path, monkeypatch, module
):
    marker = tmp_path / "executed"
    (tmp_path / (module + ".py")).write_text(
        "from pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic')\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    import importlib.util

    spec = importlib.util.find_spec(module)
    assert spec is not None
    with pytest.raises(fence.ScopeDenied):
        importlib.import_module(module)
    assert not marker.exists()


def test_fresh_process_replay_has_no_optional_scientific_imports():
    code = "from scripts.phase_4j75_source_path import replay,ROOT; from pathlib import Path; import sys; r=replay(Path.cwd(),ROOT); assert r['resolved_native_rows']==0; assert not any(n.split('.')[0] in {'openmm','lammps','psi4','pypdf'} for n in sys.modules); print('offline replay passed')"
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert "offline replay passed" in result.stdout


def test_relocated_fresh_process_replay_needs_no_source_assets(tmp_path):
    relocated = tmp_path / "relocated"
    for name in {*j75.PARENTS, *j75.EXTRA_PARENTS, *j75.REPLAY_FILES}:
        target = relocated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    for name in (*j75.CAPTURE_FILES, "frozen-inputs.json"):
        target = relocated / j75.ROOT / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    documents = j75.derive(REPO, ROOT)
    for name, document in documents.items():
        j75.write(relocated / j75.ROOT / name, document)
    import os

    result = subprocess.run(
        [sys.executable, "-m", "scripts.phase_4j75_source_path", "replay"],
        cwd=relocated,
        env={**os.environ, "PYTHONPATH": str(relocated)},
        text=True,
        capture_output=True,
        check=True,
    )
    assert json.loads(result.stdout)["decision"] == "parameter_incomplete"
    assert not (relocated / "island_ffref").exists()


def test_explicit_implementation_inventory_and_readiness(expected):
    assert len(FILES) == 4
    assert expected["capability-report.json"]["decision"] == "parameter_incomplete"
    assert j75.SCOPE["J75_new_QM_attempts"] == j75.SCOPE["J75_MD_steps"] == 0
    assert j75.SCOPE["J72_unintended_external_QM_attempts"] == 4
    assert j75.SCOPE["J72_zero_execution_scope_met"] is False
    assert j75.SCOPE["automatic_assignment"] is False


def test_source_assets_not_redistributed(tmp_path):
    (tmp_path / "synthetic.frc").write_text("synthetic")
    with pytest.raises(j75.IntakeError, match="redistribution"):
        artifacts(tmp_path)


@pytest.mark.parametrize(
    "control",
    [
        "QM",
        "MD",
        "converter",
        "fitting",
        "energy_force_evaluator",
        "parameter_generation",
        "reserved_validation",
    ],
)
def test_live_controls_outside_j75_scope(control):
    pytest.skip(control + " forbidden in J75; no live executable or fixture launched")
