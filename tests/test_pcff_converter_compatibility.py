"""Synthetic software controls; real chemistry is in separate J13 receipts."""

from copy import deepcopy

import pytest
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make, row  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack
from island.exceptions import PCFFError
from island.forcefields import PCFFOptions, adopt_forcefield
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
)
from island.forcefields.pcff.fallbacks import (
    COMPATIBILITY_POLICY,
    CONVERTER_BB13,
    MSI_POLICY,
    lookup,
    policy_evidence,
)
from island.forcefields.pcff.model import bb13_policy_applies


def test_policy_is_explicit_and_missing_only():
    assert (
        PCFFOptions(
            "local.frc",
            (0, 0, 1),
            (0, 0, 1),
            resolution_policy=COMPATIBILITY_POLICY,
            typing_profile="island_pcff_source_graph_v4",
        ).resolution_policy
        == COMPATIBILITY_POLICY
    )
    deps = [{"status": "assigned"}] * 2

    def applies(
        f="bond-bond_1_3",
        labels=("hc", "c2", "oc", "c2"),
        dependencies=deps,
        status="missing",
        policy=COMPATIBILITY_POLICY,
    ):
        return bb13_policy_applies(
            f, labels, dependencies, resolution_policy=policy, status=status
        )

    assert applies()
    assert not applies(policy=MSI_POLICY)
    for status in ("assigned", "ambiguous", "not_applicable"):
        assert not applies(status=status)
    for k in range(4):
        labels = ["c2"] * 4
        labels[k] = "cp"
        assert not applies(labels=labels)
    assert not applies(dependencies=[{"status": "missing"}] * 2)
    assert not applies(dependencies=[])
    assert not applies(f="bond-angle")
    assert policy_evidence(MSI_POLICY)["cross_terms"].endswith(
        "no policy-derived BB13 zeros"
    )
    assert policy_evidence(COMPATIBILITY_POLICY)["converter_bb13"] == CONVERTER_BB13


def test_guarded_matcher_conflicts_stay_conflicts():
    rows = [
        row("quadratic_angle", ["*", "b", "*"], [1.0, 2.0], rid="a"),
        row("quadratic_angle", ["a", "b", "*1"], [1.0, 3.0], rid="b"),
    ]
    catalog = {r["id"]: dict(r, line=i + 1) for i, r in enumerate(rows)}
    for policy in (MSI_POLICY, COMPATIBILITY_POLICY):
        assert (
            lookup(
                "quadratic_angle",
                ["a", "b", "c"],
                "cff91_auto",
                catalog,
                [],
                policy=policy,
            )["status"]
            == "ambiguous"
        )


def test_native_model_distinguishes_converter_zeros_and_rejects_resigning(
    fallback_source,  # noqa: F811
):
    system, typing, _, _, old = make(fallback_source)
    before = system.to_dict()
    charge = assign_automatic_pcff_charges(
        system, typing, resolution_policy=COMPATIBILITY_POLICY
    )
    assignment = assign_pcff_parameters(
        system, typing, charge, resolution_policy=COMPATIBILITY_POLICY
    )
    model = define_pcff_model(assignment, special_pairs=old.payload["special_pairs"])
    assert model.payload["model_definition_complete"]
    assert not model.payload["raw_parameter_coverage_complete"]
    zeros = [
        t for t in model.payload["terms"] if t["origin"] == "converter_derived_zero"
    ]
    assert zeros and all(
        t["source_rows"] == [] and t["converter_evidence"] == CONVERTER_BB13
        for t in zeros
    )
    prepared = adopt_forcefield(system, "pcff", model)
    assert prepared.native_result.identity == model.identity
    for field, value in [
        ("origin", "source_row"),
        ("coefficients", [1.0]),
        ("source_rows", [{"invented": True}]),
    ]:
        payload = deepcopy(model.payload)
        next(t for t in payload["terms"] if t["origin"] == "converter_derived_zero")[
            field
        ] = value
        with pytest.raises(PCFFError):
            type(model)(pack(payload), assignment).validate_integrity(system)
    assert system.to_dict() == before


def test_assigned_zero_ambiguous_and_unrelated_missing_never_collapsed(fallback_source):  # noqa: F811
    from island.forcefields.pcff.model import definition

    system, typing, _, _, old = make(fallback_source)
    charge = assign_automatic_pcff_charges(
        system, typing, resolution_policy=COMPATIBILITY_POLICY
    )
    assignment = assign_pcff_parameters(
        system, typing, charge, resolution_policy=COMPATIBILITY_POLICY
    )
    data = assignment.payload
    bb13 = next(a for a in data["assignments"] if a["family"] == "bond-bond_1_3")
    # Direct internal numerical fixtures test assembly decisions, not signed source data.
    bb13.update(
        status="assigned",
        normalized_values=[0.0],
        selected=[{"record_id": "synthetic-explicit-zero"}],
    )
    output = definition(data, old.payload["special_pairs"])
    t = next(t for t in output["terms"] if t["assignment_id"] == bb13["id"])
    assert t["origin"] == "source_row" and t["source_rows"] == bb13["selected"]
    bb13.update(status="ambiguous", reason="synthetic conflict")
    output = definition(data, old.payload["special_pairs"])
    assert not output["model_definition_complete"]
    assert any(d["id"] == bb13["id"] for d in output["diagnostics"])
    bb13["status"] = "missing"
    other = next(a for a in data["assignments"] if a["family"] == "bond-angle")
    other.update(status="missing", reason="synthetic absent coupling")
    assert not definition(data, old.payload["special_pairs"])[
        "model_definition_complete"
    ]
