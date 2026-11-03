"""J72 source controls; synthetic transformations never count as HMDS evidence."""

import copy
import json
import math
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import pytest

from island.pcff_completion import CompletionError, signed
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.source_adjudication import (
    CONTRACT_KEYS,
    compare,
    native_authority,
    normalized_cell,
    profile_gate,
)
from scripts.adjudicate_phase_4j72 import (
    BASE,
    BRANCH,
    DECLARATIONS,
    PARENT,
    REQUIRED,
    ROOT,
    inputs,
    replay,
    reports,
)
from scripts.check_phase_4j72_sources import numeric_row
from scripts.validate_phase_4j72_receipt import offline_declarations, semantics

REPO = Path(__file__).resolve().parents[1]
E = REPO / ROOT


def resign(value):
    return signed({k: v for k, v in value.items() if k != "identity"})


@pytest.fixture(scope="module")
def retained():
    return reports(REPO, E)


def control(value="1", unit="kcal/mol"):
    return signed(
        {
            "record_id": "synthetic-control",
            "fields": ["k"],
            "contract": {k: "same_complete_control_contract" for k in CONTRACT_KEYS},
            "comparison_values": [normalized_cell(value, unit)],
        }
    )


def test_synthetic_same_row_numeric_conflict():
    result = compare(control("1"), control("2"))
    assert result["numerically_compared"]
    assert result["categories"] == ["same_row_numeric_conflict"]
    assert result["normalized_difference"] == ["-1"]
    assert result["authority_selected"] is None and not result["values_averaged"]


def test_independent_PDF_numeric_row_parser_rejects_ambiguity_and_missing_cells():
    assert numeric_row("A 1.2 −3.4\nB 5 6", "A", 2) == ["1.2", "−3.4"]
    for text in ("A 1.2", "A unclear −3.4", "A 1 2\nA 1 2"):
        with pytest.raises(ValueError, match="ambiguous_or_missing"):
            numeric_row(text, "A", 2)


def test_synthetic_unit_equivalent_transformation_has_no_native_authority():
    result = compare(control("4.184", "kJ/mol"), control("1", "kcal/mol"))
    assert result["numerically_compared"] and not result["categories"]
    assert result["equal_with_declared_arithmetic_allowance"]
    assert not result["native_PCFF_resolution"]
    assert Decimal(normalized_cell("1", "nm")["value"]) == 10


def test_synthetic_documented_degree_radian_potential_and_force_transformation():
    # Same scalar harmonic potential, including coefficient angular dimensions.
    theta_deg, theta0_deg, k_deg = 120.0, 109.5, 0.07
    theta = float(normalized_cell(str(theta_deg), "degree")["value"])
    theta0 = float(normalized_cell(str(theta0_deg), "degree")["value"])
    k_rad = float(normalized_cell(str(k_deg), "kcal/(mol*degree^2)")["value"])
    assert k_deg * (theta_deg - theta0_deg) ** 2 == pytest.approx(
        k_rad * (theta - theta0) ** 2
    )
    # Independent chain rule and finite difference in the radian coordinate.
    force = -2 * k_deg * (theta_deg - theta0_deg) * 180 / math.pi
    step = 1e-6
    energy = lambda t: k_rad * (t - theta0) ** 2
    assert force == pytest.approx(
        -(energy(theta + step) - energy(theta - step)) / (2 * step), rel=1e-9
    )


@pytest.mark.parametrize("key,category", list(CONTRACT_KEYS.items()))
def test_synthetic_incompatible_metadata_prevents_numerical_comparison(key, category):
    a, b = control(), control("999999")
    b["contract"][key] = "different"
    result = compare(a, resign(b))
    assert category in result["categories"]
    assert not result["numerically_compared"]
    assert result["normalized_difference"] is None


@pytest.mark.parametrize("key", list(CONTRACT_KEYS))
def test_synthetic_missing_metadata_not_resolved_by_equal_values(key):
    a, b = control(), control()
    b["contract"][key] = None
    result = compare(a, resign(b))
    assert result["categories"] == ["insufficient_evidence"]
    assert not result["numerically_compared"]


def test_synthetic_unknown_unit_and_missing_cell_are_unresolved():
    for unit, value in (("unidentified", "1"), ("kcal/mol", None)):
        result = compare(control(value, unit), control())
        assert "insufficient_evidence" in result["categories"]
        assert not result["numerically_compared"]
    assert "unit_conflict" in compare(control("1", "angstrom"), control())["categories"]


def authority_control():
    return {
        "selection_basis": "source_evidence",
        "authority_level": 1,
        "source_family": "PCFF",
        "resolution": "all_atom",
        "exact_environment": True,
        "same_source_contract": True,
        "complete_equation_conventions": True,
        "unresolved_conflicts": [],
        "zero_fill": False,
        "wildcard": False,
        "source_sha256": "a" * 64,
        "source_location": "manufactured source:1",
        "source_row_id": "manufactured-row",
        "equation_metadata": {
            k: "manufactured complete metadata"
            for k in (
                "energy_equation",
                "force_equation",
                "units",
                "phase",
                "prefactor",
                "multiplicity",
                "orientation",
                "equilibrium_dependencies",
            )
        },
    }


@pytest.mark.parametrize("level", [1, 2, 3])
def test_synthetic_source_hierarchy_necessary_complete_contract(level):
    row = authority_control()
    row["authority_level"] = level
    assert native_authority(row)


@pytest.mark.parametrize(
    "field,value",
    [
        ("authority_level", 4),
        ("authority_level", 5),
        ("authority_level", 6),
        ("authority_level", True),
        ("source_family", "COMPASS"),
        ("source_family", "CFF93"),
        ("source_family", "CFF91"),
        ("source_family", "OPLS"),
        ("source_family", "OPLS-AA"),
        ("source_family", "GAFF"),
        ("source_family", "GAFF2"),
        ("source_family", "PDMS_ClassII"),
        ("source_family", "QM_research_correction"),
        ("resolution", "united_atom"),
        ("resolution", "coarse_grained"),
        ("exact_environment", False),
        ("same_source_contract", False),
        ("complete_equation_conventions", False),
        ("equation_metadata", {}),
        ("unresolved_conflicts", ["source_authority_conflict"]),
        ("zero_fill", True),
        ("wildcard", True),
        ("source_location", None),
        ("source_row_id", None),
        ("selection_basis", "best_HMDS_residual"),
        ("selection_basis", "best_MD_performance"),
    ],
)
def test_synthetic_native_authority_rejects_overlay_analog_or_residual_selection(
    field, value
):
    row = authority_control()
    row[field] = value
    with pytest.raises(CompletionError):
        native_authority(row)


def test_real_J71_reconstruction_and_raw_identity_preservation(retained):
    original = read(REPO / PARENT / "parameter-records.json")["records"]
    normalized = retained["comparison-normalized-records.json"]["records"]
    assert len(original) == len(normalized) == 467
    assert sum(v is not None for r in original for v in r["raw_value"]) == 1470
    for old, new in zip(original, normalized, strict=True):
        assert old["record_id"] == new["record_id"]
        assert old["identity"] == new["J71_identity"]
        for key in (
            "raw_value",
            "raw_unit",
            "location",
            "citation",
            "source_file_sha256",
        ):
            assert new[key] == old[key]
        assert new["authority_level_for_J70_target"] == 5
        assert not new["native_pcff_authority"]


def test_real_duplicate_equivalences_are_not_numeric_repairs_or_PCFF_rows(retained):
    result = retained["mathematical-resolutions.json"]
    assert (
        len(result["resolutions"]) == 16 and result["unequal_conflicts_resolved"] == []
    )
    ids = []
    for row in result["resolutions"]:
        ids.extend(row["alternatives"])
        assert row["comparison"]["equal_with_declared_arithmetic_allowance"]
        assert row["native_PCFF_rows_resolved"] == row["target_instances_resolved"] == 0
        assert len(row["citations"]) == 2
    assert len(set(ids)) == 32


def test_real_two_unequal_angle_assignments_remain_unselected(retained):
    conflicts = retained["source-conflict-matrix.json"][
        "documented_source_disagreements"
    ]
    assert len(conflicts) == 2
    for row in conflicts:
        assert row["status"] == "source_conflict_unresolved"
        assert row["selected_record"] is None
        assert row["comparison"]["numerically_compared"]
        assert "same_row_numeric_conflict" in row["categories"]
        assert "source_authority_conflict" in row["categories"]


def test_real_378_instances_keep_exact_mapping_dependencies_and_all_candidates(
    retained,
):
    old = read(REPO / "docs/evidence/phase_4j70/completeness-matrix.json")["instances"]
    j71 = read(REPO / PARENT / "mapping-report.json")["instances"]
    rows = retained["row-adjudication.json"]["instances"]
    assert len(rows) == len({r["id"] for r in rows}) == 378
    for source, screen, row in zip(old, j71, rows, strict=True):
        assert row["id"] == source["id"]
        assert row["J70_identity"] == source["identity"]
        assert row["equilibrium_dependencies"] == source["dependencies"]
        assert row["stable_sites"] == source["stable_sites"]
        assert row["ordered_types"] == source["ordered_types"]
        assert row["symmetry_equivalent_types"] == source["symmetry_equivalent_types"]
        assert {c["record_id"] for c in row["candidate_records"]} == {
            m["record_id"] for m in screen["candidate_mappings"]
        }
        assert row["selected_source"] is None and not row["resolved"]
    assert retained["decision.json"]["summary"]["row_status_counts"] == {
        "analog_unqualified": 294,
        "convention_unresolved": 60,
        "source_missing": 24,
    }


def test_real_family_catalogs_remain_separate_partial_and_nonexecutable(retained):
    catalogs = retained["candidate-source-bundles.json"]["catalogs"]
    ids = [i for c in catalogs for i in c["record_ids"]]
    assert len(ids) == len(set(ids)) == 467
    for catalog in catalogs:
        assert catalog["parameter_completeness"] == "incomplete"
        assert not catalog["coherent_complete_contract"]
        assert catalog["execution_profile"] is None
        assert catalog["native_rows"] == catalog["external_complete_rows"] == []


def test_synthetic_profile_gate_exact_coverage_and_converter_requirements():
    targets = [
        {
            "id": str(i),
            "identity": str(i),
            "stable_sites": [i],
            "ordered_types": ["synthetic"],
        }
        for i in range(378)
    ]
    rows = [
        {
            **authority_control(),
            "target_id": t["id"],
            "target_identity": t["identity"],
            "stable_sites": t["stable_sites"],
            "ordered_types": t["ordered_types"],
            "exact": True,
            "complete": True,
            "contract_id": "synthetic-only",
        }
        for t in targets
    ]
    assert profile_gate(rows, targets, strict_conversion=True, static_closure=True)
    for conversion, closure in ((False, False), (True, False), (False, True)):
        with pytest.raises(CompletionError):
            profile_gate(
                rows, targets, strict_conversion=conversion, static_closure=closure
            )
    for field, value in (
        ("target_id", "missing"),
        ("target_identity", "bad"),
        ("stable_sites", [0]),
        ("ordered_types", ["bad"]),
        ("complete", False),
        ("exact", False),
        ("source_family", "COMPASS"),
        ("zero_fill", True),
        ("wildcard", True),
        ("contract_id", "other"),
        ("authority_level", True),
        ("same_source_contract", False),
        ("exact_environment", False),
        ("complete_equation_conventions", False),
        ("resolution", "united_atom"),
        ("equation_metadata", {}),
    ):
        mutated = copy.deepcopy(rows)
        mutated[-1][field] = value
        with pytest.raises(CompletionError):
            profile_gate(mutated, targets, strict_conversion=True, static_closure=True)
    with pytest.raises(CompletionError):
        profile_gate(rows[:-1], targets, strict_conversion=True, static_closure=True)
    rows[-1] = rows[0]
    with pytest.raises(CompletionError):
        profile_gate(rows, targets, strict_conversion=True, static_closure=True)


def test_real_strict_conversion_static_closure_and_validation_stay_inactive(
    retained, monkeypatch
):
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **kw: pytest.fail(
            "replay must not launch QM, fitting, MD or converter"
        ),
    )
    assert reports(REPO, E) == retained
    for name in ("strict-converter-report.json", "static-closure-report.json"):
        r = retained[name]
        assert not r["eligible"] and not r["executed"] and not r["complete"]
    d = retained["decision.json"]
    assert d["new_QM"] == 4 and d["HMDS_new_QM"] == 0
    assert d["QM_fits"] == d["accepted_coefficients"] == 0
    assert d["MD_steps"] is None
    assert d["MD_steps_verified_minimum"] == 46
    assert d["MD_steps_possible_maximum"] == 48
    assert not d["zero_execution_scope_met"]
    assert (
        not d["reserved_validation_activated"] and not d["reserved_validation_consumed"]
    )


def metadata_control():
    manifest = signed(
        {
            "schema": "island.j72-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {},
            "historical_files": {},
            "artifacts": {},
            "scientific_scope": DECLARATIONS,
        }
    )
    receipt = signed(
        {
            "schema": "island.j72-receipt.v1",
            **DECLARATIONS,
            "base_commit": BASE,
            "required_J71": REQUIRED,
            "branch": BRANCH,
            "manifest_sha256": "synthetic",
            "software_counts": {"tests": 1, "failure": 0, "error": 0, "skipped": 0},
            "terminal_unique_tests": 1,
            "historical_files_preserved": 1,
            "offline_reconstruction": offline_declarations(),
        }
    )
    return receipt, manifest


@pytest.mark.parametrize("field", list(DECLARATIONS))
@pytest.mark.parametrize("where", ["receipt", "manifest", "nested"])
def test_resigned_scientific_declarations_rejected(field, where):
    receipt, manifest = copy.deepcopy(metadata_control())
    semantics(receipt, manifest)
    target = (
        receipt
        if where == "receipt"
        else manifest["scientific_scope"]
        if where == "manifest"
        else receipt["offline_reconstruction"]
    )
    old = target[field]
    target[field] = (
        not old if type(old) is bool else old + 1 if type(old) is int else "forged"
    )
    with pytest.raises(CompletionError):
        semantics(resign(receipt), resign(manifest))


@pytest.mark.parametrize(
    "field", ["automatic_assignment", "new_QM", "accepted_coefficients"]
)
def test_resigned_bool_integer_semantics_are_strict(field):
    receipt, manifest = copy.deepcopy(metadata_control())
    receipt[field] = int(receipt[field]) if type(receipt[field]) is bool else False
    with pytest.raises(CompletionError):
        semantics(resign(receipt), manifest)


@pytest.mark.parametrize(
    "name,field",
    [
        ("comparison-normalized-records.json", "raw_unit"),
        ("comparison-normalized-records.json", "source_file_sha256"),
        ("comparison-normalized-records.json", "raw_value"),
        ("comparison-normalized-records.json", "location"),
        ("comparison-normalized-records.json", "contract"),
        ("comparison-normalized-records.json", "J71_identity"),
        ("source-conflict-matrix.json", "summary"),
        ("row-adjudication.json", "summary"),
        *[
            ("row-adjudication.json", k)
            for k in (
                "stable_sites",
                "ordered_types",
                "source_lookup",
                "equilibrium_dependencies",
                "final_J72_status",
                "selected_source",
            )
        ],
        ("candidate-source-bundles.json", "parameter_completeness"),
        ("strict-converter-report.json", "eligible"),
        ("decision.json", "automatic_assignment"),
    ],
)
def test_offline_replay_rejects_resigned_derived_semantic_mutation(
    tmp_path, name, field
):
    # Keep nested and outer identities valid: rejection must use reconstruction.
    shutil.copytree(E / "scope-deviation", tmp_path / "scope-deviation")
    for p in E.glob("*"):
        if p.is_file():
            shutil.copyfile(p, tmp_path / p.name)
    value = read(tmp_path / name)
    target = (
        value["records"][0]
        if "records" in value
        else value["instances"][0]
        if name == "row-adjudication.json" and field != "summary"
        else value
    )
    target[field] = "forged"
    if target is not value:
        value["records" if "records" in value else "instances"][0] = resign(target)
    (tmp_path / name).write_text(json.dumps(resign(value)))
    with pytest.raises(CompletionError, match="semantic_reconstruction"):
        replay(REPO, tmp_path)


def test_frozen_contract_resigning_cannot_change_parent_evidence(tmp_path):
    value = read(E / "contract.json")
    value["scientific_scope"]["new_QM"] = 1
    (tmp_path / "contract.json").write_text(json.dumps(resign(value)))
    with pytest.raises(CompletionError, match="frozen_contract_mutation"):
        inputs(REPO, tmp_path)


def test_fresh_process_relocation_without_git_pdfs_or_solver(tmp_path):
    for name in ("src", "scripts"):
        shutil.copytree(
            REPO / name, tmp_path / name, ignore=shutil.ignore_patterns("__pycache__")
        )
    shutil.copytree(E, tmp_path / ROOT, ignore=shutil.ignore_patterns("validation"))
    shutil.copytree(
        REPO / PARENT, tmp_path / PARENT, ignore=shutil.ignore_patterns("validation")
    )
    # Explicit parent pins, not a historical broad implementation glob.
    names = set(read(E / "contract.json")["parent_pins"]) | set(
        read(REPO / PARENT / "contract.json")["parent_pins"]
    )
    for name in names:
        source, destination = REPO / name, tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    assert not (tmp_path / "island_ffref").exists() and not (tmp_path / ".git").exists()
    env = {**os.environ, "PYTHONPATH": str(tmp_path) + ":" + str(tmp_path / "src")}
    p = subprocess.run(
        [
            sys.executable,
            "scripts/adjudicate_phase_4j72.py",
            "--repo",
            str(tmp_path),
            "--replay",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert p.returncode == 0, p.stdout + p.stderr
    assert "source_conflict_unresolved" in p.stdout and "'new_QM': 4" in p.stdout


def test_canonical_api_remains_fail_closed_with_exact_378_missing_instances():
    from island.pcff_pdms import PCFFPDMSAssets, prepare_polymer

    assets = (
        REPO.parent
        / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
    )
    if not assets.is_dir():
        pytest.skip(
            "declared native source assets unavailable; offline capability still checked"
        )
    system, prepared = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="pcff",
        dp=1,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        assets=PCFFPDMSAssets(assets),
    )
    result = prepared.validate_integrity(system)
    assert len(result["missing_terms"]) == 378
    assert not result["parameter_complete"] and result["selected_profile"] is None
    with pytest.raises(Exception, match="parameter_complete"):
        prepared.require()


@pytest.mark.parametrize(
    "field,value",
    [
        ("failure", False),
        ("error", False),
        ("tests", True),
        ("skipped", -1),
        ("unexpected", 0),
    ],
)
def test_resigned_software_count_schema_mutations(field, value):
    receipt, manifest = copy.deepcopy(metadata_control())
    receipt["software_counts"][field] = value
    with pytest.raises(CompletionError, match="software_count_schema"):
        semantics(resign(receipt), manifest)


@pytest.mark.parametrize(
    "artifact", ["source-conflict-matrix.json", "source-conflict-report.md"]
)
def test_strict_validation_rejects_semantics_after_resigning_entire_outer_chain(
    tmp_path, artifact
):
    from island.pcff_completion import pin, write
    from scripts.validate_phase_4j64_receipt import artifacts
    from scripts.validate_phase_4j72_receipt import FILES, validate

    shutil.copytree(E / "scope-deviation", tmp_path / "scope-deviation")
    for p in E.iterdir():
        if p.is_file() and p.name not in ("manifest.json", "receipt.json"):
            shutil.copyfile(p, tmp_path / p.name)
    if artifact.endswith(".json"):
        changed = read(tmp_path / artifact)
        finding = changed["documented_source_disagreements"][0]
        finding["selected_record"] = finding["alternatives"][0]
        changed["documented_source_disagreements"][0] = resign(finding)
        write(tmp_path / artifact, resign(changed))
    else:
        with (tmp_path / artifact).open("a") as stream:
            stream.write("\nForged claim: parameters are now accepted.\n")
    # Synthetic verification header; real retained source evidence is the target
    # of the mutation. Every outer identity/inventory is recomputed correctly.
    receipt, manifest = metadata_control()
    manifest.update(
        {
            "implementation_files": {n: pin(REPO / n) for n in FILES},
            "artifacts": artifacts(tmp_path),
        }
    )
    write(tmp_path / "manifest.json", resign(manifest))
    receipt["manifest_sha256"] = pin(tmp_path / "manifest.json")
    write(tmp_path / "receipt.json", resign(receipt))
    assert read(tmp_path / "manifest.json")["artifacts"] == artifacts(tmp_path)
    assert read(tmp_path / "receipt.json")["manifest_sha256"] == pin(
        tmp_path / "manifest.json"
    )
    with pytest.raises(
        CompletionError, match="semantic_reconstruction|human_report_mutation"
    ):
        validate(REPO, tmp_path)


def test_real_convention_findings_separate_proved_mismatch_from_suspected_error(
    retained,
):
    rows = retained["convention-adjudications.json"]["findings"]
    assert len(rows) == 5
    ids = {
        r["record_id"]
        for r in retained["comparison-normalized-records.json"]["records"]
    }
    for row in rows:
        assert row["candidate_record_ids"] and set(row["candidate_record_ids"]) <= ids
        assert row["status"] == "unresolved" and row["selected_record"] is None
        assert row["native_PCFF_rows_resolved"] == 0
        assert not row["numeric_comparison_performed"]
    by_id = {r["id"]: r for r in rows}
    assert by_id["PolCA_OH_H_printed_length"]["categories"] == ["insufficient_evidence"]
    assert by_id["PolCA_OH_H_printed_length"]["candidate_record_ids"] == ["cell-0133"]
    assert "equation_conflict" in by_id["native_Phi_versus_cosine_AAT"]["categories"]
    assert by_id["Si_C_H_reference_order"]["candidate_record_ids"] == [
        "cell-0241",
        "cell-0242",
    ]


@pytest.mark.parametrize(
    "argv",
    [
        ["python", "-m", "scripts.validate_phase_4j31_network"],
        ["/native/sqm"],
        ["antechamber"],
        ["python", "-S", "file.py"],
        ["python", "-c", "resume_periodic_md_workflow()"],
    ],
)
def test_child_execution_entry_points_are_blocked(argv):
    from scripts.phase_4j72_execution_fence import denied

    assert denied(argv)
    assert not denied(["python", "-m", "scripts.phase_4j37_reference", "--help"])


def test_child_fence_blocks_integrator_and_native_command_before_execution(tmp_path):
    if "ISLAND_J72_FENCE_DIRECTORY" not in os.environ:
        pytest.skip("requires the J72 child-fenced verification environment")
    probe = tmp_path / "probe.py"
    (tmp_path / "psi4.py").write_text(
        "raise AssertionError('module code must never execute')\n"
    )
    probe.write_text(
        "import openmm,subprocess,importlib.util\nblocked=0\n"
        "assert importlib.util.find_spec('psi4') is not None\n"
        "try:\n import psi4\nexcept RuntimeError as e:\n assert 'J72 child fence' in str(e);blocked+=1\n"
        "try:\n openmm.VerletIntegrator(0.001).step(1)\nexcept RuntimeError as e:\n assert 'J72 child fence' in str(e);blocked+=1\n"
        "try:\n subprocess.run(['sqm'])\nexcept RuntimeError as e:\n assert 'J72 child fence' in str(e);blocked+=1\n"
        "assert blocked==3\nprint('all blocked before execution')\n"
    )
    result = subprocess.run(
        [sys.executable, str(probe)], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "all blocked before execution" in result.stdout
    (tmp_path / "psi4.py").unlink()
    probe.write_text(
        "import importlib.util; assert importlib.util.find_spec('psi4') is None\n"
    )
    result = subprocess.run(
        [sys.executable, str(probe)], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_actual_scope_violation_counts_reconstruct_without_scientific_execution():
    from island.pcff_qm_pilot.source_adjudication import execution_deviation

    result = execution_deviation(E)
    assert result["new_QM"] == 4 and result["new_HMDS_QM"] == 0
    assert result["MD_steps"] is None
    assert result["MD_steps_verified_minimum"] == 46
    assert result["MD_steps_possible_maximum"] == 48
    assert (
        not result["zero_execution_scope_met"]
        and not result["use_for_source_adjudication"]
    )


@pytest.mark.parametrize(
    "artifact", ["execution-deviation.json", "quarantined-j31-regressions.tar.xz"]
)
def test_quarantined_execution_evidence_mutation_rejected(tmp_path, artifact):
    from island.pcff_qm_pilot.source_adjudication import execution_deviation

    shutil.copytree(E / "scope-deviation", tmp_path / "scope-deviation")
    path = tmp_path / "scope-deviation" / artifact
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(CompletionError, match="deviation_.*_mutation"):
        execution_deviation(tmp_path)


def test_terminal_correction_retains_raw_failures_without_double_counting(tmp_path):
    from xml.etree import ElementTree as ET

    from island.pcff_completion import write
    from scripts.validate_phase_4j72_receipt import (
        terminal_sources,
        terminal_xml_closure,
    )

    raw = tmp_path / "validation/guarded-final"
    raw.mkdir(parents=True)
    nodes = [
        "tests/control.py::old",
        "tests/test_phase_4j35_c1.py::affected",
        "tests/test_phase_4j72.py::new",
    ]
    for i in range(8):
        suite = ET.Element("testsuite")
        if i == 0:
            ET.SubElement(suite, "testcase", classname="tests.control", name="old")
            c = ET.SubElement(
                suite, "testcase", classname="tests.test_phase_4j35_c1", name="affected"
            )
            ET.SubElement(c, "failure", message="retained old fence failure")
        ET.ElementTree(suite).write(raw / f"full-{i}.xml")
    suite = ET.Element("testsuite")
    ET.SubElement(
        suite, "testcase", classname="tests.test_phase_4j35_c1", name="affected"
    )
    ET.SubElement(suite, "testcase", classname="tests.test_phase_4j72", name="new")
    ET.ElementTree(suite).write(tmp_path / "validation/guard-corrections-terminal.xml")
    selected = terminal_sources(tmp_path, nodes)
    assert (
        len(selected) == 3
        and len(list(ET.parse(raw / "full-0.xml").iter("failure"))) == 1
    )
    merged = ET.Element("testsuite")
    for _, case in selected.values():
        merged.append(copy.deepcopy(case))
    ET.ElementTree(merged).write(tmp_path / "validation/terminal-results.xml")
    write(tmp_path / "validation/junit-files.json", ["validation/terminal-results.xml"])
    write(
        tmp_path / "validation/terminal-selection.json",
        {n: p for n, (p, _) in selected.items()},
    )
    terminal_xml_closure(tmp_path, nodes)
    merged[0].set("time", "forged")
    ET.ElementTree(merged).write(tmp_path / "validation/terminal-results.xml")
    with pytest.raises(CompletionError, match="terminal_result_mutation"):
        terminal_xml_closure(tmp_path, nodes)
