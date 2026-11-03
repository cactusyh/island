"""Offline source facts only; all process probes are marker-only fake files."""

import copy
import io
import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from scripts import phase_4j72_c1_fence as fence
from scripts.phase_4j74_source_intake import (
    PARENTS,
    ROOT,
    SCOPE,
    STATES,
    IntakeError,
    compatible,
    family_of,
    inputs,
    inventory,
    originals,
    parse_archive,
    parse_frc,
    read,
    reconstruct,
    replay,
    scope,
    signed,
    validate_inventory,
    validate_resolution,
    write,
)
from scripts.validate_phase_4j74_receipt import (
    FILES,
    PARSER_PIP_DIAGNOSTIC,
    artifacts,
    semantic,
)

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", autouse=True)
def evidence_only_scope():
    scope()


@pytest.fixture(scope="module")
def expected():
    return reconstruct(REPO, REPO / ROOT)


def resign(value):
    return signed({k: v for k, v in value.items() if k != "identity"})


def test_deterministic_recursive_inventory_and_byte_mutation(tmp_path):
    (tmp_path / "nested").mkdir()
    (tmp_path / "z.txt").write_text("synthetic-only z\n")
    (tmp_path / "nested/a.txt").write_text("synthetic-only a\n")
    first = inventory(tmp_path)
    assert first == inventory(tmp_path)
    assert [s["relative_path"] for s in first] == ["nested/a.txt", "z.txt"]
    (tmp_path / "nested/a.txt").write_text("synthetic-only changed\n")
    second = inventory(tmp_path)
    assert second[0]["sha256"] != first[0]["sha256"]
    assert second[1] == first[1]
    assert inventory(tmp_path / "not-supplied") == []


def test_symlink_source_not_silently_followed(tmp_path):
    (tmp_path / "copy.txt").write_text("synthetic")
    (tmp_path / "alias.txt").symlink_to("copy.txt")
    with pytest.raises(IntakeError, match="symlink"):
        inventory(tmp_path)


@pytest.mark.parametrize(
    "mutation",
    ["duplicate", "no_license", "no_redistribution", "author_upgrade", "native_hash"],
)
def test_source_identity_license_and_authority_rejected(mutation):
    p = read(REPO / ROOT / "source-file-inventory.json")
    if mutation == "duplicate":
        p["files"].append(copy.deepcopy(p["files"][0]))
    else:
        r = next(r for r in p["files"] if r["relative_path"] == "pcff.frc")
        if mutation == "no_license":
            del r["license_status"]
        elif mutation == "no_redistribution":
            r["redistribution_status"] = ""
        elif mutation == "author_upgrade":
            r["primary_artifact_supplied"] = True
        else:
            r["sha256"] = "0" * 64
        p["files"] = [resign(x) for x in p["files"]]
    with pytest.raises(IntakeError):
        validate_inventory(resign(p))


@pytest.mark.parametrize(
    "mutation",
    ["sha", "unit", "equation", "family", "row_id", "location", "zero_fill", "profile"],
)
def test_resigned_fabricated_source_facts_rejected(tmp_path, mutation):
    frozen = read(REPO / ROOT / "frozen-inputs.json")
    shutil.copyfile(REPO / ROOT / "frozen-inputs.json", tmp_path / "frozen-inputs.json")
    for name in frozen["source_facts"]:
        shutil.copyfile(REPO / ROOT / name, tmp_path / name)
    facts = read(tmp_path / "extracted-parameter-records.json")
    row = facts["records"][0]
    if mutation == "sha":
        row["source_sha256"] = "0" * 64
    elif mutation == "unit":
        row["units"] = {"status": "complete", "raw": "invented kJ"}
    elif mutation == "equation":
        row["analytic_equation"] = {
            "status": "complete",
            "raw": "invented half-prefactor",
        }
    elif mutation == "family":
        row["source_family"] = "native_PCFF_CFF91_distribution_candidate"
    elif mutation == "row_id":
        row["record_id"] = "self-signed-authority"
    elif mutation == "location":
        row["location"] = {"page": "invented"}
    elif mutation == "zero_fill":
        row["coefficients_raw"] = ["0"]
        row["used_to_fill_missing_row"] = True
    else:
        row["native_PCFF_eligible"] = True
    facts["records"][0] = resign(row)
    write(tmp_path / "extracted-parameter-records.json", resign(facts))
    with pytest.raises(IntakeError, match="fabricated_or_resigned_source_facts"):
        inputs(REPO, tmp_path)


@pytest.mark.parametrize(
    "name,expected_family",
    [
        ("pcff.frc", "native_PCFF_CFF91_distribution_candidate"),
        ("cff91.frc", "native_PCFF_CFF91_distribution_candidate"),
        ("compass_published.frc", "native_COMPASS_distribution_candidate"),
        ("PDMS_benchmark-main/param/HUANG/H_system.in.setting", "external_Huang_PDMS"),
        (
            "PDMS_benchmark-main/param/COMPASS/system.in.settings",
            "external_Xiang_COMPASS_implementation",
        ),
        (
            "PDMS_benchmark-main/param/UA/UA_system.in.setting",
            "external_Xiang_united_atom_implementation",
        ),
        (
            "PDMS_benchmark-main/param/OPLS/opls_system.in.settings",
            "external_Xiang_OPLS_implementation",
        ),
    ],
)
def test_families_never_merged_by_local_presence(name, expected_family):
    assert family_of(name) == expected_family


def synthetic_source():
    return {
        "relative_path": "synthetic.frc",
        "sha256": "1" * 64,
        "citation": {"status": "synthetic"},
        "source_family": "native_PCFF_CFF91_distribution_candidate",
        "source_authority": "synthetic_control_not_real_source",
    }


def test_frc_raw_values_order_missing_cross_terms_and_no_zero_fill():
    text = """#atom_types cff91
1.0 1 sio 28.086 Si 4 synthetic
#equivalence cff91
1.0 1 c3 c c c c c
#quartic_bond cff91
> E = K2*(r-r0)^2 + K3*(r-r0)^3 + K4*(r-r0)^4
> units kcal/mol and angstrom
1.0 1 osi sio 1.6562 306.1232 -517.3424 673.7067
#end_bond-torsion_3 cff91
> coefficients absent in this synthetic control
"""
    records, facts = parse_frc(text, synthetic_source())
    assert len(records) == 1
    assert records[0]["ordered_atom_types"] == ["osi", "sio"]
    assert records[0]["coefficients_raw"] == [
        "1.6562",
        "306.1232",
        "-517.3424",
        "673.7067",
    ]
    assert facts["ordinary_equivalences"]["c3"]["bond"] == "c"
    assert "end_bond-torsion_3:cff91" in facts["section_headers"]
    assert not any(r["interaction_family"] == "end_bond-torsion_3" for r in records)
    assert records[0]["normalized_values"] is None
    assert records[0]["used_to_fill_missing_row"] is False


@pytest.mark.parametrize(
    "mutation", ["unit", "equation", "COMPASS", "Huang", "types", "graph"]
)
def test_explicit_contract_comparison_rejects_incompatible_metadata(mutation):
    target = read(REPO / PARENTS[0])["requests"][0]
    contract = target["required_evidence"]["equation_and_units"][
        "historical_unresolved_contract"
    ]
    record = {
        "source_family": "native_PCFF_CFF91_distribution_candidate",
        "ordered_atom_types": target["ordered_types"],
        "interaction_family": target["family"],
        "exact_graph_applicability": True,
        "units": {
            "status": "complete_source_contract",
            "symbols": contract["unit_symbols"],
        },
        "analytic_equation": {
            "status": "complete_source_contract",
            "scalar": contract["scalar_energy"],
        },
        **{
            k: {"status": "complete_source_contract"}
            for k in (
                "phase_or_multiplicity",
                "orientation_order_convention",
                "equilibrium_dependencies",
            )
        },
    }
    assert compatible(record, target)["unit_match"]
    assert compatible(record, target)["function_form_match"]
    if mutation == "unit":
        record["units"]["symbols"] = {"E": "kJ/mol"}
        key = "unit_match"
    elif mutation == "equation":
        record["analytic_equation"]["scalar"] = "different prefactor/function"
        key = "function_form_match"
    elif mutation in ("COMPASS", "Huang"):
        record["source_family"] = mutation
        key = "source_family_match"
    elif mutation == "types":
        record["ordered_atom_types"] = ["Si", "O", "Si", "C"]
        key = "exact_type_match"
    else:
        record["exact_graph_applicability"] = False
        key = "exact_graph_match"
    result = compatible(record, target)
    assert not result[key] and not result["eligibility"]


def test_archive_duplicate_identity_and_traversal_rejected():
    for name in ("../outside.txt", "/absolute.txt"):
        data = io.BytesIO()
        with zipfile.ZipFile(data, "w") as z:
            z.writestr(name, "synthetic only")
        with pytest.raises(IntakeError, match="archive_member_path"):
            parse_archive(data.getvalue(), synthetic_source())
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as z:
        z.writestr("synthetic.txt", "first")
        with pytest.warns(UserWarning, match="Duplicate name"):
            z.writestr("synthetic.txt", "hidden alternative")
    with pytest.raises(IntakeError, match="duplicate_archive_source_identity"):
        parse_archive(data.getvalue(), synthetic_source())


def test_real_coverage_family_catalog_and_unchanged_requests(expected):
    matrix = expected["row-resolution-matrix.json"]
    assert matrix["counts"] == {
        "analog_only": 294,
        "convention_blocked": 60,
        "source_missing": 24,
    }
    assert len(matrix["instances"]) == 378
    for row in matrix["instances"]:
        assert set(row["statuses"]) == set(STATES)
        assert not row["statuses"]["primary_artifact_supplied"]
        assert not row["statuses"]["exact_row_resolved"]
        assert not row["statuses"]["convention_resolved"]
        assert not row["statuses"]["executable_bundle_qualified"]
        assert row["native_source_queries"]
        assert all(
            not q["matched_records"] and not q["auto_equivalence_or_wildcard_used"]
            for q in row["native_source_queries"]
        )
    coverage = expected["family-coverage-report.json"]
    assert len(coverage["reports"]) == 5
    assert all(
        r["exact_rows_supplied"] == 0 and not r["automatic_assignment"]
        for r in coverage["reports"]
    )
    assert "torsion-torsion" in coverage["interaction_family_extraction_inventory"]


@pytest.mark.parametrize(
    "mutation",
    [
        "missing",
        "duplicate",
        "identity",
        "primary",
        "resolved",
        "convention",
        "native_analog",
        "profile",
    ],
)
def test_resigned_resolution_coverage_and_promotion_rejected(expected, mutation):
    matrix = copy.deepcopy(expected["row-resolution-matrix.json"])
    row = matrix["instances"][0]
    if mutation == "missing":
        matrix["instances"].pop()
    elif mutation == "duplicate":
        matrix["instances"][1] = copy.deepcopy(row)
    elif mutation == "identity":
        row["J73_identity"] = "0" * 64
    elif mutation in ("primary", "resolved", "convention"):
        key = {
            "primary": "primary_artifact_supplied",
            "resolved": "exact_row_resolved",
            "convention": "convention_resolved",
        }[mutation]
        row["statuses"][key] = True
        row["resolution_status"] = key
    elif mutation == "native_analog":
        row["native_exact_match_count"] = 1
    else:
        matrix["scope"]["profile"] = "forbidden_PCFF_profile"
    matrix["instances"][0] = resign(row)
    with pytest.raises(IntakeError):
        validate_resolution(resign(matrix), expected["row-resolution-matrix.json"])


def envelope():
    from scripts.phase_4j74_source_intake import BASE, BRANCH, J73, MAIN

    common = {
        "base_commit": BASE,
        "previous_J73_commit": J73,
        "fetched_origin_main": MAIN,
        "scope": SCOPE.copy(),
    }
    r = signed(
        {
            **common,
            "schema": "island.j74-receipt.v1",
            "branch": BRANCH,
            "status": "blocked",
            "manifest_sha256": "0" * 64,
            "historical_files_preserved": 0,
            "source_resolution": {
                "source_resolution_status": "blocked",
                "resolved_rows": 0,
                "counts": {
                    "analog_only": 294,
                    "convention_blocked": 60,
                    "source_missing": 24,
                },
                "scope": SCOPE.copy(),
            },
            "verification": {
                "initial_readonly_inspection_before_J74_preflight": True,
                "initial_sequence_compliant": False,
                "parser_environment_pip_check": {
                    "returncode": 1,
                    "status": "existing_torch_dependency_conflicts_retained; no_dependency_changes",
                    "diagnostic": PARSER_PIP_DIAGNOSTIC,
                },
            },
        }
    )
    m = signed(
        {
            **common,
            "schema": "island.j74-manifest.v1",
            "implementation_files": dict.fromkeys(FILES, "0" * 64),
            "artifacts": {},
            "frozen_inputs_sha256": "0" * 64,
        }
    )
    return r, m


@pytest.mark.parametrize("key", list(SCOPE))
@pytest.mark.parametrize("where", ["receipt", "manifest", "nested"])
def test_resigned_receipts_cannot_erase_history_or_enable_execution(key, where):
    r, m = envelope()
    target = (
        m["scope"]
        if where == "manifest"
        else r["source_resolution"]["scope"]
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
    with pytest.raises(IntakeError):
        semantic(resign(r), resign(m))


def test_receipt_cannot_erase_preflight_sequence_deviation():
    r, m = envelope()
    r["verification"]["initial_sequence_compliant"] = True
    with pytest.raises(IntakeError, match="deviation"):
        semantic(resign(r), m)


def test_receipt_cannot_hide_parser_environment_failure():
    r, m = envelope()
    r["verification"]["parser_environment_pip_check"]["returncode"] = 0
    with pytest.raises(IntakeError, match="failure_must_not_be_erased"):
        semantic(resign(r), m)


@pytest.mark.parametrize("suffix", [".pdf", ".frc", ".zip", ".dat", ".data"])
def test_unlicensed_source_assets_cannot_enter_manifest(tmp_path, suffix):
    (tmp_path / ("synthetic-placeholder" + suffix)).write_text(
        "synthetic-only, no source data"
    )
    with pytest.raises(IntakeError, match="redistribution_forbidden"):
        artifacts(tmp_path)


def fake(path, marker):
    path.write_text(
        "#!/usr/bin/env python3\n# ISLAND_C1_SYNTHETIC_ONLY\n"
        + f"open({str(marker)!r},'w').write('fake')\n"
    )
    path.chmod(0o755)
    return path


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "antechamber",
        "psi4",
        "lammps",
        "msi2lmp",
        "fit_hmds",
        "reserved_validation",
        "evaluate_energy",
        "evaluate_forces",
        "evaluate_pressure",
    ],
)
def test_fake_native_science_cannot_launch_in_parent_or_child(tmp_path, name):
    marker = tmp_path / "ran"
    executable = fake(tmp_path / name, marker)
    with pytest.raises(fence.ScopeDenied):
        subprocess.run([str(executable)], check=False)
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import subprocess,sys\nfrom scripts.phase_4j74_source_intake import scope\nscope()\ntry: subprocess.run([sys.argv[1]])\nexcept RuntimeError: pass\nelse: raise AssertionError('launched')",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


@pytest.mark.parametrize(
    "name",
    [
        "fit_hmds",
        "parameter_generation",
        "reserved_validation",
        "evaluate_energy",
        "evaluate_forces",
        "evaluate_pressure",
        "run_converter",
    ],
)
def test_fake_python_scientific_driver_cannot_launch(tmp_path, name):
    marker = tmp_path / "ran"
    executable = fake(tmp_path / (name + ".py"), marker)
    with pytest.raises(fence.ScopeDenied):
        subprocess.run([sys.executable, str(executable)], check=False)
    assert not marker.exists()


def test_missing_fence_blocks_before_source_reads(monkeypatch):
    monkeypatch.setattr(
        fence,
        "assert_installed",
        lambda: (_ for _ in ()).throw(RuntimeError("absent fence")),
    )
    with pytest.raises(RuntimeError, match="absent fence"):
        inventory(Path("must-not-access"))


def test_cli_missing_fence_returns_blocked_receipt_before_intake(monkeypatch, capsys):
    from scripts.validate_phase_4j74_receipt import main

    monkeypatch.setattr(
        fence,
        "assert_installed",
        lambda: (_ for _ in ()).throw(RuntimeError("absent fence")),
    )
    assert main() == 2
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "blocked" and result["scope"] == SCOPE


def test_original_source_byte_mutation_rejected(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    file = source / "synthetic.txt"
    file.write_text("original synthetic-only bytes")
    root = tmp_path / "evidence"
    write(root / "source-inputs.json", {"files": inventory(source)})
    file.write_text("changed synthetic-only bytes")
    with pytest.raises(IntakeError, match="inventory_mutation"):
        originals(root, source)


def test_relocation_and_offline_replay_without_sources(tmp_path):
    relocated = tmp_path / "repository"
    shutil.copytree(
        REPO / "scripts",
        relocated / "scripts",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    for name in PARENTS:
        target = relocated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    shutil.copytree(
        REPO / ROOT,
        relocated / ROOT,
        ignore=shutil.ignore_patterns("validation", "manifest.json", "receipt.json"),
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j74_receipt",
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
    assert not (relocated / "island_ffref").exists()


def test_replay_and_original_hash_mutation_check(expected, tmp_path):
    assert replay(REPO, REPO / ROOT)["source_resolution_status"] == "blocked"
    with pytest.raises(IntakeError, match="inventory_mutation"):
        originals(REPO / ROOT, tmp_path)


@pytest.mark.parametrize(
    "control",
    [
        "QM",
        "MD",
        "OpenMM",
        "LAMMPS",
        "AmberTools",
        "converter",
        "fitting",
        "parameter_generation",
        "reserved_validation",
        "energy_force_pressure_evaluator",
    ],
)
@pytest.mark.skip(
    reason="J74 read-only intake: live science and scientific evaluator controls forbidden; only synthetic rejection probes run"
)
def test_live_scientific_control_scope_skip(control):
    raise AssertionError("live scientific fixture must not execute")
