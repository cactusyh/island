"""Synthetic matcher/integrity controls; not independent chemical acceptance."""

from itertools import product

import pytest
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make  # noqa: F401
from test_pcff_fallbacks import row as synthetic_row
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields import ForceFieldRequest, PCFFOptions
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_parameters,
    save_pcff_parameters,
)
from island.forcefields.pcff.fallbacks import MSI_POLICY, POLICY, lookup
from island.forcefields.pcff.model import bb13_policy_applies
from island.forcefields.pcff.resolution import classify
from island.forcefields.pcff.source import PCFFSource, digest


def row(*args, **kwargs):
    return synthetic_row(*args, **kwargs, namespace="cff91")


def query(family, types, rows, **kwargs):
    catalog = {r["id"]: dict(r, line=i + 1) for i, r in enumerate(rows)}
    return lookup(family, types, "cff91", catalog, [], policy=MSI_POLICY, **kwargs)


def test_exact_precedes_wildcard_but_wildcard_specificity_is_not_authority():
    f = "quartic_angle"
    rows = [
        row(f, ["*1", "b", "*2"], [1, 2]),
        row(f, ["a", "b", "*3"], [1, 3], rid="specific"),
    ]
    for table in (rows, rows[::-1]):
        result = query(f, ["a", "b", "c"], table)
        assert result["status"] == "ambiguous"
        assert {c["record_id"] for c in result["decision_trace"]["candidates"]} == {
            "x",
            "specific",
        }
        assert query(f, ["c", "b", "a"], table)["status"] == "ambiguous"
        exact = row(f, ["a", "b", "c"], [1, 4], rid="exact")
        r = query(f, ["a", "b", "c"], [*table, exact])
        assert r["status"] == "assigned" and r["normalized_values"] == [1, 4]
        assert r["decision_trace"]["tier"] == "exact"
        assert {c["decision"] for c in r["decision_trace"]["candidates"]} == {
            "lower_match_tier",
            "coefficient_agreement_required",
        }


def test_versions_zero_and_agreeing_wildcards():
    f = "bond-angle"
    old = row(f, ["a", "b", "*"], [3, 9], version="1", rid="old")
    new = row(f, ["a", "b", "*"], [0, 0], version="2", rid="new")
    broad = row(f, ["*1", "b", "*2"], [0, 0], rid="broad")
    r = query(f, ["a", "b", "c"], [old, new, broad])
    assert r["normalized_values"] == [0, 0]
    assert {c["record_id"] for c in r["selected"]} == {"new", "broad"}
    assert (
        r["decision_trace"]["candidates"][0]["decision"] == "superseded_source_version"
    )
    assert query(f, ["a", "z", "c"], [old, new, broad])["status"] == "missing"


@pytest.mark.parametrize(
    "family,values,expected",
    [
        ("bond-angle", [2, 7], [7, 2]),
        ("end_bond-torsion_3", [1, 2, 3, 7, 8, 9], [7, 8, 9, 1, 2, 3]),
        ("angle-torsion_3", [1, 2, 3, 7, 8, 9], [7, 8, 9, 1, 2, 3]),
    ],
)
def test_directional_roles(family, values, expected):
    labels = ["a", "b", "c"] if len(values) == 2 else ["a", "b", "c", "d"]
    r = query(family, labels[::-1], [row(family, labels, values)])
    assert r["normalized_values"] == expected
    assert r["selected"][0]["permutation"] == list(reversed(range(len(labels))))


def test_improper_center_and_shared_arm_cannot_move():
    rows = [row("angle-angle", ["a", "b", "c", "d"], [7])]
    assert query("angle-angle", ["d", "b", "c", "a"], rows)["status"] == "assigned"
    assert query("angle-angle", ["d", "c", "b", "a"], rows)["status"] == "missing"
    rows = [row("wilson_out_of_plane", ["a", "b", "c", "d"], [7, 0])]
    assert (
        query("wilson_out_of_plane", ["c", "b", "d", "a"], rows)["status"] == "assigned"
    )
    assert (
        query("wilson_out_of_plane", ["b", "a", "c", "d"], rows)["status"] == "missing"
    )


def test_repeated_labels_cannot_hide_directional_conflict():
    assert (
        query(
            "bond-angle", ["a", "b", "a"], [row("bond-angle", ["a", "b", "a"], [3, 4])]
        )["status"]
        == "ambiguous"
    )


def test_source_required_policy_does_not_reinterpret_historical_zeros():
    deps = [{"status": "assigned"}]
    assert bb13_policy_applies("bond-bond_1_3", ["hc", "c", "c", "hc"], deps)
    assert not bb13_policy_applies(
        "bond-bond_1_3", ["hc", "c", "c", "hc"], deps, resolution_policy=MSI_POLICY
    )
    r = query("bond-bond_1_3", ["hc", "c", "c", "hc"], [])
    classified = classify(
        {
            "resolution": r,
            "family": "bond-bond_1_3",
            "supplied_types": ["hc", "c", "c", "hc"],
            "equilibrium_dependencies": deps,
            "dependencies_complete": True,
            "request_id": "q",
            "sites": [1, 2, 3, 4],
        },
        {},
    )
    assert classified["structural_model_blocker"]
    assert (
        classified["classification"]
        == "source_parameter_missing_under_declared_searches"
    )


def test_complete_pipeline_and_trace_tampering(fallback_source, monkeypatch, tmp_path):  # noqa: F811
    from island.forcefields.pcff import expanded, model, source

    raw = (
        fallback_source.raw.replace(b"\n#end\n", b"\n")
        + b"\n#bond-bond_1_3 cff91\n"
        + "\n".join(
            "1 1 " + " ".join(labels) + " 0" for labels in product(("c", "h"), repeat=4)
        ).encode()
        + b"\n#end\n"
    )
    for mapping, key in [
        (source.PIN, "sha256"),
        (expanded.PROFILE, "source_sha256"),
        (model.PROFILE, "frc_sha256"),
    ]:
        monkeypatch.setitem(mapping, key, digest(raw))
    src = PCFFSource(raw, digest(raw))
    m, t, _, _, _ = make(src)
    before = m.to_dict()
    q = assign_automatic_pcff_charges(m, t, resolution_policy=MSI_POLICY)
    a = assign_pcff_parameters(m, t, q, resolution_policy=MSI_POLICY)
    p = define_pcff_model(
        a, special_pairs=model.special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert p.payload["model_definition_complete"]
    assert p.payload["term_origins"] == {"source_row": len(p.payload["terms"])}
    assert any(x["family"] == "quadratic_bond" for x in a.payload["assignments"])
    save_pcff_parameters(a, tmp_path / "assignment.json")
    assert (
        load_pcff_parameters(tmp_path / "assignment.json", src, system=m).identity
        == a.identity
    )
    changed = unpack(a.json_text)
    next(x for x in changed["assignments"] if "decision_trace" in x)["decision_trace"][
        "tier"
    ] = "invented"
    with pytest.raises(PCFFError):
        type(a)(pack(changed), src).validate_integrity(m)
    assert m.to_dict() == before
    assert a.payload["resolution_policy"] == MSI_POLICY
    assert ForceFieldRequest(
        "pcff",
        PCFFOptions(
            "external.frc",
            (0, 0, 1),
            (0, 0, 1),
            typing_profile="island_pcff_source_graph_v4",
            resolution_policy=MSI_POLICY,
        ),
    )


def test_historical_specificity_behavior_stays_unchanged():
    f = "quartic_angle"
    rows = [
        row(f, ["*1", "b", "*2"], [1, 2]),
        row(f, ["a", "b", "*3"], [1, 3], rid="specific"),
    ]
    r = lookup(
        f, ["a", "b", "c"], "cff91", {r["id"]: r for r in rows}, [], policy=POLICY
    )
    assert r["status"] == "assigned" and r["normalized_values"] == [1, 3]
    assert "decision_trace" not in r and "searches" not in r


def test_direct_wildcard_before_ordinary_equivalence_and_positional_roles():
    eq = [
        {"id": "e" + t, "version": "1", "data": {"type": t, "families": {"angle": v}}}
        for t, v in zip(("A", "B", "C"), ("a", "b", "c"), strict=True)
    ]
    rows = [
        dict(row("quartic_angle", ["*", "B", "*"], [1, 2]), line=1),
        dict(row("quartic_angle", ["a", "b", "c"], [1, 9], rid="equiv"), line=2),
    ]
    r = lookup(
        "quartic_angle",
        ["A", "B", "C"],
        "cff91",
        {x["id"]: x for x in rows},
        eq,
        policy=MSI_POLICY,
    )
    assert r["path"] == "direct" and r["normalized_values"] == [1, 2]
    r = lookup(
        "quartic_angle",
        ["C", "B", "A"],
        "cff91",
        {"equiv": rows[1]},
        eq,
        policy=MSI_POLICY,
    )
    assert r["path"] == "ordinary_family_equivalence"
    assert r["resolved_types"] == ["c", "b", "a"]
    assert r["position_roles"] == ["angle"] * 3
    assert len(r["equivalence_evidence"]) == 3
    autoeq = [
        {
            "id": "ae" + t,
            "version": "1",
            "data": {
                "type": t,
                "families": {"angle_end": "end_" + t, "angle_apex": "center_" + t},
            },
        }
        for t in ("A", "B", "C")
    ]
    auto = dict(
        row("quadratic_angle", ["end_A", "center_B", "end_C"], [1, 3]),
        namespace="cff91_auto",
        line=4,
    )
    r = lookup(
        "quadratic_angle",
        ["C", "B", "A"],
        "cff91_auto",
        {"auto": auto},
        autoeq,
        policy=MSI_POLICY,
    )
    assert r["status"] == "assigned" and r["path"] == "automatic_position_equivalence"
    assert r["position_roles"] == ["angle_end", "angle_apex", "angle_end"]
    assert r["selected"][0]["permutation"] == [2, 1, 0]


def test_missing_bb13_blocks_public_native_adoption(fallback_source):  # noqa: F811
    from island.forcefields import PreparedForceFieldError, adopt_forcefield

    m, t, _, _, old = make(fallback_source)
    assert old.payload["model_definition_complete"]
    q = assign_automatic_pcff_charges(m, t, resolution_policy=MSI_POLICY)
    a = assign_pcff_parameters(m, t, q, resolution_policy=MSI_POLICY)
    p = define_pcff_model(a, special_pairs=old.payload["special_pairs"])
    assert not p.payload["model_definition_complete"]
    assert p.payload["term_origins"].get("policy_derived_zero", 0) == 0
    assert any(d["id"].startswith("bond-bond_1_3") for d in p.payload["diagnostics"])
    with pytest.raises(PreparedForceFieldError):
        adopt_forcefield(m, "pcff", p)


def test_unused_equivalence_conflict_cannot_preempt_direct_match():
    eq = [
        {"id": str(i), "version": "1", "data": {"type": "a", "families": {"angle": x}}}
        for i, x in enumerate(("x", "y"))
    ]
    rows = {"r": dict(row("quartic_angle", ["a", "b", "c"], [1, 2]), line=1)}
    r = lookup("quartic_angle", ["a", "b", "c"], "cff91", rows, eq, policy=MSI_POLICY)
    assert r["status"] == "assigned" and r["path"] == "direct"
    with pytest.raises(PCFFError, match="Conflicting source candidates"):
        lookup("quartic_angle", ["a", "b", "c"], "cff91", {}, eq, policy=MSI_POLICY)
