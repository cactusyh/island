"""Synthetic software tests; no synthetic fixture is real-source acceptance."""

from copy import deepcopy

import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make, row  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.exceptions import PCFFError
from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    special_pair_policy,
)
from island.forcefields.pcff.fallbacks import POLICY, lookup


def test_trace_does_not_change_resolution():
    catalog = {"a": row("bond-angle", ["a", "b", "c"], [2.0, 7.0], namespace="cff91")}
    for labels in (["c", "b", "a"], ["a", "d", "c"]):
        old = lookup("bond-angle", labels, "cff91", catalog, [])
        traced = lookup("bond-angle", labels, "cff91", catalog, [], trace=True)
        assert traced["searches"][0]["position_roles"] == ["angle"] * 3
        stripped = {
            k: v
            for k, v in traced.items()
            if k not in ("searches", "missing_equivalence_types")
        }
        assert stripped == old
        assert traced["searches"][0]["records_examined"] == 1
    catalog["b"] = row(
        "bond-angle", ["a", "b", "c"], [4.0, 9.0], rid="b", namespace="cff91"
    )
    assert (
        lookup("bond-angle", ["a", "b", "c"], "cff91", catalog, [], trace=True)[
            "status"
        ]
        == "ambiguous"
    )


def test_owned_inspection_preserves_model_and_dependencies(fallback_source):  # noqa: F811
    system, typing, charges, assignment, model = make(fallback_source)
    original = deepcopy(system)
    kwargs = {
        "resolution_policy": POLICY,
        "special_pairs": model.payload["special_pairs"],
    }
    result = inspect_pcff_operational_support(system, typing, **kwargs)
    assert result["diagnostic_only"]
    assert result["charge_identity"] == charges.identity
    assert result["assignment_identity"] == assignment.identity
    assert result["model_identity"] == model.identity
    assert result["model_complete"]
    coupling = next(
        q for q in result["interaction_queries"] if q["family"] == "bond-angle"
    )
    assert len(coupling["equilibrium_dependencies"]) == 3
    assert all(d["status"] == "assigned" for d in coupling["equilibrium_dependencies"])
    assert any(
        "quadratic_bond:cff91_auto:" in r
        for d in coupling["equilibrium_dependencies"]
        for r in d["source_rows"]
    )
    result["native_charge_record"].clear()
    result["source_rows"].clear()
    typing.validate_integrity(original)
    from island.workflows.bundle import system_data
    from island.workflows.storage import json_bytes

    assert json_bytes(system_data(system)) == json_bytes(system_data(original))
    assert (
        inspect_pcff_operational_support(system, typing, **kwargs)["charge_identity"]
        == charges.identity
    )
    changed = deepcopy(system)
    changed.topology.sites[11].mass += 1
    with pytest.raises(PCFFError):
        inspect_pcff_operational_support(changed, typing, **kwargs)
    kwargs["special_pairs"]["invented"] = True
    with pytest.raises(PCFFError):
        inspect_pcff_operational_support(system, typing, **kwargs)


def test_failed_native_charge_never_enters_parameterization(
    fallback_source,  # noqa: F811
    monkeypatch,
):
    from island.forcefields.pcff import operational

    system, typing, charge, _, _ = make(fallback_source)

    class FailedNativeControl:
        """Injected software outcome only; not a publishable native record."""

        complete = False
        payload = charge.payload
        identity = charge.identity

    monkeypatch.setattr(
        operational,
        "assign_automatic_pcff_charges",
        lambda *a, **kw: FailedNativeControl(),
    )

    def forbidden(*args, **kwargs):
        pytest.fail("Failed charge reached scientific parameterization")

    monkeypatch.setattr(operational, "assign_pcff_parameters", forbidden)
    monkeypatch.setattr(operational, "define_pcff_model", forbidden)
    result = inspect_pcff_operational_support(
        system,
        typing,
        resolution_policy=POLICY,
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    assert not result["charges_complete"] and not result["model_complete"]
    assert result["model_identity"] is None and result["interaction_queries"]
    assert "native_charge_incomplete" in result["model_diagnostics"][0]["reason"]


def test_operational_gate_is_not_a_diagnostic_success(monkeypatch):
    from pathlib import Path

    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "scripts"))
    from validate_pcff_operational_support import operational_gate

    keys = (
        "typed",
        "charges_complete",
        "model_complete",
        "source_verified",
        "independent_selection_verified",
        "whole_system_verified",
        "relocated_bundle_verified",
    )
    passed = {"case": "required", **dict.fromkeys(keys, True)}
    assert operational_gate([passed], ["required"])
    assert not operational_gate([], ["required"])
    assert not operational_gate([passed, passed], ["required"])
    for key in keys:
        assert not operational_gate([dict(passed, **{key: False})], ["required"])
    assert not operational_gate(
        [dict(passed, error="publication failed")], ["required"]
    )


def test_inspection_needs_no_optional_science(fallback_source, monkeypatch):  # noqa: F811
    import builtins

    system, typing, _, _, model = make(fallback_source)
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"openmm", "scipy", "rdkit", "parmed", "foyer"}:
            raise AssertionError("Forbidden scientific import: " + name)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    inspection = inspect_pcff_operational_support(
        system,
        typing,
        resolution_policy=POLICY,
        special_pairs=model.payload["special_pairs"],
    )
    assert inspection["model_identity"] == model.identity


def test_independent_raw_oracle_roles_and_source_missing(monkeypatch):
    from pathlib import Path

    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "scripts"))
    from pcff_j5_reference import match, read_source, resolve_raw

    raw = b"""#bond-angle cff91
1 1 a b c 2 7
#angle-angle cff91
1 1 a b c d 3
#quadratic_angle cff91_auto
1 1 *2 center x 100 2
1 1 y center *4 120 8
"""
    rows = read_source(raw)
    assert match(rows, "bond-angle", "cff91", ["c", "b", "a"])["values"] == [
        7 * 4.184,
        2 * 4.184,
    ]
    assert (
        match(rows, "angle-angle", "cff91", ["d", "b", "c", "a"])["status"]
        == "assigned"
    )
    assert (
        match(rows, "angle-angle", "cff91", ["d", "c", "b", "a"])["status"] == "missing"
    )
    assert (
        match(rows, "quadratic_angle", "cff91_auto", ["y", "center", "x"])["status"]
        == "ambiguous"
    )
    missing = resolve_raw(rows, "bond-angle", ["x", "y", "z"])
    assert missing["status"] == "missing"
    assert all(s["namespace"] != "cff91_auto" for s in missing["searches"])
