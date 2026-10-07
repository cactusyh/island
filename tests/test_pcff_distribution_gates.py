"""Software gate controls; no converter or physical calculations in these tests."""

import importlib

import pytest

from island.workflows import storage


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend("scripts")
    return (
        importlib.import_module("check_pcff_distribution"),
        importlib.import_module("validate_pcff_distribution"),
        importlib.import_module("pcff_distribution_reference"),
    )


@pytest.mark.parametrize(
    "names", [[], ["a"], ["a", "a"], ["a", "b", "b"], ["a", "other"]]
)
def test_required_case_sets_fail_closed(modules, names):
    check, _, _ = modules
    with pytest.raises(ValueError, match="mandatory"):
        check.checked_set(
            [{"case": n} for n in names], {"cases": [{"name": "a"}, {"name": "b"}]}
        )


def test_order_independence_and_invalid_declaration(modules):
    check, _, _ = modules
    check.checked_set(
        [{"case": "b"}, {"case": "a"}], {"cases": [{"name": "a"}, {"name": "b"}]}
    )
    with pytest.raises(ValueError):
        check.checked_set(
            [{"case": "a"}, {"case": "a"}], {"cases": [{"name": "a"}, {"name": "a"}]}
        )


def test_artifact_binding_and_traversal(modules, tmp_path):
    check, _, _ = modules
    p = tmp_path / "log"
    p.write_text("original")
    hashes = {"log": storage.checksum(p.read_bytes())}
    check.artifact_check(tmp_path, hashes)
    p.write_text("changed")
    with pytest.raises(ValueError, match="stale"):
        check.artifact_check(tmp_path, hashes)
    for bad in ({}, {"../outside": "x"}, {"missing": "x"}):
        with pytest.raises(ValueError):
            check.artifact_check(tmp_path, bad)


def test_evidence_stages_are_distinct_not_success_booleans(modules):
    check, _, _ = modules
    cases = [
        {
            "case": "polymer",
            "errors": [],
            "polymer_target": True,
            "native_model_complete": False,
            "acceptance_passed": True,
            "independent_numerical_artifacts_validated": True,
        }
    ]
    args = {
        "audit_validated": True,
        "total_source_labels": 133,
        "verified_source_labels": 0,
    }
    gates = check.derive_gates(cases, {"polymer"}, **args)
    assert gates["converter_boundary_gate"]
    assert not gates["executable_polymer_gate"] and not gates["full_source_complete"]
    cases[0]["native_model_complete"] = True
    assert not check.derive_gates(cases, set(), **args)["executable_polymer_gate"]
    assert check.derive_gates(cases, {"polymer"}, **args)["executable_polymer_gate"]
    cases[0]["errors"] = ["wrong reason"]
    assert not check.derive_gates(cases, {"polymer"}, **args)["converter_boundary_gate"]


@pytest.mark.parametrize(
    "code,text",
    [
        (1, "RuntimeError unrelated"),
        (0, "Unable to find bond data for a b"),
        (13, "Unable to find bond data for a b"),
        (12, "Unable to find invented data for a b"),
    ],
)
def test_wrong_converter_failure_reason_or_stage(modules, code, text):
    _, _, reference = modules
    with pytest.raises(ValueError):
        reference.missing_request({}, code, text)


def test_raw_classifications_preserve_missing_and_auto_boundary(modules, monkeypatch):
    _, _, reference = modules
    monkeypatch.setattr(
        reference,
        "resolve_raw",
        lambda *args, **kwargs: {"status": "missing", "family": "quartic_bond"},
    )
    assert (
        reference.missing_request({}, 12, "Unable to find bond data for a b")[
            "classification"
        ]
        == "source_parameter_missing"
    )
    monkeypatch.setattr(
        reference,
        "resolve_raw",
        lambda *args, **kwargs: {"status": "assigned", "family": "quadratic_bond"},
    )
    assert (
        reference.missing_request({}, 12, "Unable to find bond data for a b")[
            "classification"
        ]
        == "automatic_fallback_not_implemented_by_converter"
    )


def test_unexpected_success_runtime_or_wrong_reason(modules):
    _, runner, _ = modules
    expected = {
        "x": {
            "role": "negative_unchanged_v5",
            "expected": {
                "typing_complete": False,
                "public_exception": {
                    "type": "PCFFError",
                    "message": "Incomplete source graph typing",
                },
            },
        }
    }
    c = {"name": "x", "role": "negative_unchanged_v5"}
    good = {
        "typing_complete": False,
        "public_exception": expected["x"]["expected"]["public_exception"],
        "converter": {
            "state": "blocked_by_native_typing_or_charge",
            "not_executed": True,
        },
    }
    assert runner.check_case(good, c, expected) == []
    for change in (
        {"public_preparation_complete": True},
        {"execution_error": {"type": "RuntimeError"}},
        {"public_exception": {"type": "PCFFError", "message": "wrong source"}},
        {"converter": {"returncode": 0}},
    ):
        assert runner.check_case({**good, **change}, c, expected)


def test_converter_decimal_serialization_is_exact(modules):
    check, _, _ = modules
    original = [[2.3241877043067753, -0.14516730121201168, -0.23796651175723083]]
    atoms = [["1", "1", "1", "0", "2.324187704", "-0.145167301", "-0.237966512"]]
    result = check.check_coordinates(atoms, original)
    assert result["exact_serialized_frame"]
    assert 0 < result["maximum_difference_from_original_angstrom"] < 5e-10
    atoms[0][4] = "2.324187705"
    with pytest.raises(AssertionError):
        check.check_coordinates(atoms, original)
