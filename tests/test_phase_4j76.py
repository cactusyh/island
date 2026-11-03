"""Exact-row source controls, kept separate from real HMDS source evidence."""

import copy
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import phase_4j72_c1_fence as fence
from scripts import phase_4j76_source_resolution as j76

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / j76.ROOT


@pytest.fixture(scope="module", autouse=True)
def inherited_scope():
    j76.scope()


@pytest.fixture(scope="module")
def expected():
    return j76.reconstruct(REPO, ROOT)


def resign(data):
    return j76.signed({k: v for k, v in data.items() if k != "identity"})


def synthetic():
    contract = {f: "synthetic explicit " + f for f in j76.FIELDS}
    candidate = {
        "id": "synthetic-native-row",
        "source_family": "native_PCFF",
        "exact_types": True,
        "exact_environment": True,
        "authority_verified": True,
        "contract": copy.deepcopy(contract),
        "raw_values": ["1.25"],
        "value_origin": "source_literal",
        "closure": {
            f: {
                "status": "explicit_source_supported",
                "value": contract[f],
                "source_references": [
                    {"file": "synthetic-only", "sha256": "1" * 64, "line": 1}
                ],
            }
            for f in j76.FIELDS
        },
    }
    return contract, candidate


def test_manufactured_complete_exact_native_row_is_resolved():
    contract, c = synthetic()
    assert (
        j76.adjudicate([c], "source_missing", contract) == "exact_native_row_resolved"
    )


@pytest.mark.parametrize("field", j76.FIELDS)
def test_each_missing_or_mismatched_convention_blocks_exact_row(field):
    contract, c = synthetic()
    c["closure"].pop(field)
    assert (
        j76.adjudicate([c], "source_missing", contract)
        == "exact_native_row_convention_blocked"
    )
    _, c = synthetic()
    c["contract"][field] = "different units/equation/sign/dependency"
    assert (
        j76.adjudicate([c], "source_missing", contract)
        == "exact_native_row_convention_blocked"
    )


def test_unclosed_target_or_unproven_evidence_cannot_be_self_declared_complete():
    contract, c = synthetic()
    assert (
        j76.adjudicate([c], "source_missing") == "exact_native_row_convention_blocked"
    )
    c["closure"]["Class_II_dependencies"]["source_references"] = []
    assert (
        j76.adjudicate([c], "source_missing", contract)
        == "exact_native_row_convention_blocked"
    )


@pytest.mark.parametrize(
    "family,status",
    [
        ("COMPASS", "compass_only"),
        ("external_PDMS", "external_only"),
        ("analog", "source_missing"),
        ("CFF91_unqualified", "source_missing"),
    ],
)
def test_complete_non_native_rows_never_become_native(family, status):
    contract, c = synthetic()
    c["source_family"] = family
    assert j76.adjudicate([c], "source_missing", contract) == status


@pytest.mark.parametrize(
    "field", ["exact_types", "exact_environment", "authority_verified"]
)
def test_type_analogy_or_untrusted_artifact_is_not_exact(field):
    contract, c = synthetic()
    c[field] = False
    assert j76.adjudicate([c], "analog_only", contract) == "analog_only"


@pytest.mark.parametrize("difference", ["raw_values", "equation", "units"])
def test_unequal_rows_or_incompatible_contracts_remain_conflicting(difference):
    contract, c = synthetic()
    other = copy.deepcopy(c)
    other["id"] = "synthetic-alternative"
    if difference == "raw_values":
        other["raw_values"] = ["7.5"]
    else:
        field = (
            "functional_form" if difference == "equation" else "units_and_conversion"
        )
        other["contract"][field] = "incompatible; numeric comparison prohibited"
    before = copy.deepcopy([c, other])
    assert j76.adjudicate([c, other], "analog_only", contract) == "conflicting_source"
    assert [c, other] == before


@pytest.mark.parametrize(
    "basis", ["residual_improvement", "average", "fitted_value", "generic_default"]
)
def test_no_selection_by_fit_average_or_default(basis):
    contract, c = synthetic()
    c["selection_basis"] = basis
    with pytest.raises(j76.IntakeError, match="unqualified_parameter_selection"):
        j76.adjudicate([c], "source_missing", contract)


@pytest.mark.parametrize(
    "origin",
    ["zero_fill", "wildcard_fallback", "guessed_mixing_rule", "parameter_generation"],
)
def test_no_synthetic_missing_value_substitution(origin):
    contract, c = synthetic()
    c["value_origin"] = origin
    with pytest.raises(j76.IntakeError, match="missing_value_substitution"):
        j76.adjudicate([c], "source_missing", contract)


def test_duplicate_candidate_identity_rejected():
    contract, c = synthetic()
    with pytest.raises(j76.IntakeError, match="duplicate_candidate"):
        j76.adjudicate([c, copy.deepcopy(c)], "source_missing", contract)


def fake_source():
    return {
        "relative_path": "synthetic.frc",
        "sha256": "1" * 64,
        "source_family": "synthetic",
        "native_source_identity_complete": False,
    }


def fake_target():
    return {
        "family": "angle-angle",
        "symmetry_equivalent_types": [["a", "b", "c", "d"], ["d", "b", "c", "a"]],
    }


def test_independent_parser_covers_unretained_rows_without_zero_fill():
    text = """#angle-angle cff91
> E = K * (Theta - Theta0) * (Theta' - Theta0')
1.0 1 a b c d 1.25
1.0 1 a b c d 2.5
1.0 1 * b c d 0.0
#end_bond-torsion_3 cff91
"""
    parsed = j76.parse_frc(text, fake_source(), [])
    assert len(parsed["target_family_row_index"]) == 3
    assert parsed["target_family_row_index"][0]["retained_J75_fact_ids"] == []
    assert (
        len(j76.query_source(parsed, fake_target())["exact_numeric_row_matches"]) == 2
    )
    assert not any(
        r["section"] == "end_bond-torsion_3" for r in parsed["target_family_row_index"]
    )


def test_ordinary_equivalence_is_explicit_and_AA_reversal_is_restricted():
    text = """#equivalence cff91
1.0 1 a na na na na na
#equivalence cff91_auto
1.0 1 a wrong wrong wrong wrong wrong
#angle-angle cff91
1.0 1 d c b a 1.25
1.0 1 d b c a 2.5
"""
    parsed = j76.parse_frc(text, fake_source(), [])
    assert parsed["ordinary_equivalences"]["a"]["angle"] == "na"
    found = j76.query_source(parsed, fake_target())["exact_numeric_row_matches"]
    assert len(found) == 1 and found[0]["ordered_types"] == ["d", "b", "c", "a"]


@pytest.mark.parametrize(
    "section,row,alternative",
    [
        ("atom_types", "1.0 1 a 12.0 C 4 synthetic", "1.0 1 a 12.0 C 4 synthetic"),
        ("equivalence", "1.0 1 a na na na na na", "1.0 1 a na different na na na"),
    ],
)
def test_duplicate_typing_and_equivalence_do_not_silently_override(
    section, row, alternative
):
    with pytest.raises(j76.IntakeError, match="duplicate_or_conflicting"):
        j76.parse_frc(
            "#" + section + " cff91\n" + row + "\n" + alternative + "\n",
            fake_source(),
            [],
        )


def test_identical_equivalence_alternatives_retain_both_locations_and_versions():
    parsed = j76.parse_frc(
        "#equivalence cff91\n2.1 11 a na na na na na\n1.0 1 a na na na na na\n",
        fake_source(),
        [],
    )
    declarations = parsed["ordinary_equivalences"]["a"]["source_declarations"]
    assert [d["line"] for d in declarations] == [2, 3]
    assert [d["version"] for d in declarations] == ["2.1", "1.0"]
    assert declarations[0]["mapping"] == declarations[1]["mapping"]


def test_real_all_378_have_independent_native_queries_and_no_false_exact_row(expected):
    rows = expected["row-resolution-matrix.json"]
    assert len(rows["instances"]) == 378
    assert rows["counts"] == {
        **dict.fromkeys(j76.STATES, 0),
        "analog_only": 294,
        "unresolved": 60,
        "source_missing": 24,
    }
    assert rows["family_counts"] == {
        "angle-angle": 78,
        **dict.fromkeys(sorted(j76.TARGET_FAMILIES - {"angle-angle"}), 60),
    }
    for r in rows["instances"]:
        assert not r["exact_native_row_present"] and r["closed_convention_fields"] == 0
        assert len(r["source_queries"]) == 3
        assert not any(q["exact_numeric_row_matches"] for q in r["source_queries"])
        assert all(c["eligibility"] is False for c in r["candidate_records"])
        if r["J75_status"] == "convention_blocked":
            assert r["status"] == "unresolved"


def test_source_inventory_same_root_and_source_family_separation(expected):
    capture = j76.inputs(REPO, ROOT)
    assert capture["exact_source_root"] == j76.j75.AUTHORITATIVE
    assert len(capture["stable_graph"]) == 27
    assert capture["J75_retained_source_facts"] == 1012
    assert expected["native-pcff-coverage.json"]["retained_facts"] == 358
    for name in ("native-pcff", "compass", "external-pdms"):
        assert expected[name + "-coverage.json"]["exact_native_requests_resolved"] == 0
        assert (
            expected[name + "-coverage.json"][
                "COMPASS_or_external_values_promoted_to_PCFF"
            ]
            is False
        )


def test_closure_distinguishes_family_forms_from_missing_rows(expected):
    g = expected["convention-closure-matrix.json"]["global_conventions"]
    assert g["bond_function"]["evidence"] and g["Lennard_Jones_form"]["evidence"]
    assert all(s["element_and_valence_match"] for s in g["atom_typing"]["sites"])
    for k in (
        "charges",
        "one_two_rules",
        "one_three_rules",
        "one_four_rules",
        "long_range_electrostatics",
        "cutoff",
        "tail",
    ):
        assert "unresolved" in g[k]["status"]


def test_conflicts_and_all_evidence_requests_are_preserved(expected):
    c = expected["source-conflicts.json"]
    assert c["equal_duplicate_equivalences"] == 16
    assert c["unequal_unresolved_COMPASS_assignments"] == 2
    assert c["new_resolutions"] == 0
    requests = expected["remaining-evidence-requests.json"]["requests"]
    assert len({r["id"] for r in requests}) == 378
    assert requests == sorted(requests, key=lambda r: r["id"])


@pytest.mark.parametrize(
    "path",
    [
        "/island_ffref/pcff",
        "island_ffref/pcff",
        j76.j75.AUTHORITATIVE + "/",
        "/tmp/pcff",
    ],
)
def test_path_substitution_rejected_before_source_access(tmp_path, path, monkeypatch):
    monkeypatch.setattr(
        j76.j75, "originals", lambda *_: pytest.fail("wrong source root accessed")
    )
    with pytest.raises(j76.IntakeError, match="incorrect_authoritative"):
        j76.capture(REPO, tmp_path / "out", path)


@pytest.mark.parametrize(
    "mutation",
    [
        "source_hash",
        "row_values",
        "row_location",
        "family",
        "source_path",
        "profile",
        "units",
        "dependency",
    ],
)
def test_resigned_capture_or_authority_mutation_rejected(tmp_path, mutation):
    for n in ("source-capture.json", "historical-files.json", "frozen-inputs.json"):
        shutil.copyfile(ROOT / n, tmp_path / n)
    c = j76.read(tmp_path / "source-capture.json")
    s = c["FRC_sources"][0]
    if mutation == "source_hash":
        s["source"]["sha256"] = "0" * 64
    elif mutation == "source_path":
        c["exact_source_root"] = "/island_ffref/pcff"
    elif mutation == "profile":
        c["scope"]["executable_profile_created"] = True
    else:
        s["target_family_row_index"][0][mutation] = "fabricated"
    j76.write(tmp_path / "source-capture.json", resign(c))
    freeze = j76.read(tmp_path / "frozen-inputs.json")
    freeze["captured_files"]["source-capture.json"] = j76.pin(
        tmp_path / "source-capture.json"
    )
    j76.write(tmp_path / "frozen-inputs.json", resign(freeze))
    with pytest.raises(j76.IntakeError, match="resigned_input_freeze"):
        j76.inputs(REPO, tmp_path)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_row",
        "duplicate_row",
        "native_promotion",
        "convention_promotion",
        "profile",
    ],
)
def test_resigned_derived_evidence_mutation_rejected(tmp_path, expected, mutation):
    for n in ("source-capture.json", "historical-files.json", "frozen-inputs.json"):
        shutil.copyfile(ROOT / n, tmp_path / n)
    for name, document in expected.items():
        j76.write(tmp_path / name, document)
    data = j76.read(tmp_path / "row-resolution-matrix.json")
    if mutation == "missing_row":
        data["instances"].pop()
    elif mutation == "duplicate_row":
        data["instances"].append(copy.deepcopy(data["instances"][0]))
    elif mutation == "profile":
        data["scope"]["automatic_assignment_enabled"] = True
    else:
        r = data["instances"][0]
        if mutation == "native_promotion":
            r["status"] = "exact_native_row_resolved"
        else:
            r["convention_closure"]["units_and_conversion"]["status"] = (
                "explicit_source_supported"
            )
        data["instances"][0] = resign(r)
    j76.write(tmp_path / "row-resolution-matrix.json", resign(data))
    with pytest.raises(j76.IntakeError, match="resigned_scientific_evidence"):
        j76.replay(REPO, tmp_path)


def test_fresh_relocated_process_replay_has_no_source_assets(tmp_path):
    relocated = tmp_path / "offline"
    for name in {*j76.PARENTS, *j76.REPLAY_FILES}:
        target = relocated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    for name in (
        "source-capture.json",
        "historical-files.json",
        "frozen-inputs.json",
        *j76.reconstruct(REPO, ROOT),
    ):
        target = relocated / j76.ROOT / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    p = subprocess.run(
        [sys.executable, "-m", "scripts.phase_4j76_source_resolution", "replay"],
        cwd=relocated,
        env={**os.environ, "PYTHONPATH": str(relocated)},
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(p.stdout)["status_counts"]["exact_native_row_resolved"] == 0
    assert not (relocated / "island_ffref").exists()


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "psi4",
        "lammps",
        "msi2lmp",
        "antechamber",
        "converter",
        "evaluate_energy",
        "generate_parameters",
        "reserved_validation",
    ],
)
def test_marker_only_fake_native_executables_cannot_launch(tmp_path, name):
    marker = tmp_path / "ran"
    fake = tmp_path / name
    fake.write_text(
        "#!/usr/bin/env python\nfrom pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic')\n"
    )
    fake.chmod(0o755)
    with pytest.raises(fence.ScopeDenied):
        subprocess.run([str(fake)], check=True)
    assert not marker.exists()


@pytest.mark.parametrize(
    "driver",
    [
        "run_qm",
        "run_md",
        "run_converter",
        "evaluate_energy",
        "generate_parameters",
        "optimize_coefficients",
        "reserved_validation",
    ],
)
def test_marker_only_fake_python_drivers_cannot_launch(tmp_path, driver):
    fake = tmp_path / (driver + ".py")
    fake.write_text("raise RuntimeError('synthetic; must never execute')\n")
    with pytest.raises(fence.ScopeDenied):
        subprocess.run([sys.executable, str(fake)], check=True)


def test_child_inherits_fence_and_rejects_before_fake_execution(tmp_path):
    fake = tmp_path / "sqm"
    fake.write_text(
        "#!/usr/bin/env python\nraise RuntimeError('synthetic must never run')\n"
    )
    fake.chmod(0o755)
    code = (
        "from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; from scripts.phase_4j72_c1_fence import ScopeDenied,assert_installed; assert_installed();\ntry: subprocess.run(["
        + repr(str(fake))
        + "])\nexcept ScopeDenied: print('rejected synthetic child')\nelse: raise RuntimeError('bypass')"
    )
    p = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert "rejected synthetic child" in p.stdout


def test_missing_fence_prevents_source_access(monkeypatch, tmp_path):
    def missing():
        raise RuntimeError("missing fence")

    monkeypatch.setattr(fence, "assert_installed", missing)
    with pytest.raises(RuntimeError, match="missing fence"):
        j76.capture(REPO, tmp_path / "out", j76.j75.AUTHORITATIVE)


@pytest.mark.parametrize("module", ["openmm", "psi4", "lammps"])
def test_fake_module_discovery_allowed_but_import_cannot_execute(
    module, tmp_path, monkeypatch
):
    marker = tmp_path / "ran"
    (tmp_path / (module + ".py")).write_text(
        "from pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic')\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    assert importlib.util.find_spec(module) is not None
    with pytest.raises(fence.ScopeDenied):
        __import__(module)
    assert not marker.exists()


@pytest.mark.parametrize(
    "control",
    [
        "QM",
        "MD",
        "LAMMPS_OpenMM",
        "converter",
        "fitting_generation",
        "evaluator",
        "reserved_validation",
    ],
)
def test_live_controls_are_explicit_scope_skips(control):
    pytest.skip(control + " prohibited in J76; no scientific fixture launched")
