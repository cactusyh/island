"""Synthetic software contracts for source/model adjudication; no new physics."""

import builtins
from copy import deepcopy

import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFResolutionAssessment,
    assess_pcff_resolution,
    load_pcff_resolution_assessment,
    save_pcff_resolution_assessment,
)
from island.forcefields.pcff.fallbacks import POLICY
from island.forcefields.pcff.model import bb13_policy_applies
from island.forcefields.pcff.resolution import classify


def query(family="bond-bond_1_3", labels=("hc", "c3", "c3", "hc"), status="missing"):
    return {
        "request_id": family + ":1,2,3,4",
        "family": family,
        "sites": [1, 2, 3, 4],
        "supplied_types": list(labels),
        "resolution": {"status": status, "normalized_values": [0, 0]},
        "equilibrium_dependencies": [{"status": "assigned"}, {"status": "assigned"}],
        "dependencies_complete": True,
    }


def test_existing_non_cp_policy_all_roles_and_dependencies():
    q = query()
    result = classify(q, {})
    assert result["classification"] == "existing_policy_derived_zero"
    assert not result["structural_model_blocker"]
    assert result["raw_source_status"] == "missing"
    assert "non_cp" in result["policy"] and not result["selected_source_rows"]
    for role in range(4):
        changed = deepcopy(q)
        changed["supplied_types"][role] = "cp"
        assert classify(changed, {})["structural_model_blocker"]
    for dep in ("missing", "ambiguous", "unresolved_dependencies"):
        changed = deepcopy(q)
        changed["equilibrium_dependencies"][0]["status"] = dep
        assert classify(changed, {})["structural_model_blocker"]
    assert not bb13_policy_applies(
        "bond-angle", q["supplied_types"], q["equilibrium_dependencies"]
    )


def test_no_automatic_coupling_or_wildcard_resolution():
    for family in ("bond-angle", "bond-bond", "angle-angle", "middle_bond-torsion_3"):
        assert (
            classify(query(family), {})["classification"]
            == "source_parameter_missing_under_declared_searches"
        )
    ambiguous = query("quartic_angle", status="ambiguous")
    assert classify(ambiguous, {})["classification"] == "interpretation_unresolved"
    explicit_zero = query("bond-bond", status="assigned")
    assert classify(explicit_zero, {})["classification"] == "source_assigned"
    explicit_zero["dependencies_complete"] = False
    assert classify(explicit_zero, {})["structural_model_blocker"]
    wilson = query("wilson_out_of_plane", status="assigned")
    wilson["resolution"]["normalized_values"] = [1.0, 0.2]
    assert (
        classify(wilson, {})["classification"] == "nonzero_wilson_semantics_unresolved"
    )


def test_assessment_owned_persistence_and_rechecksummed_contradictions(
    fallback_source,  # noqa: F811
    tmp_path,
    monkeypatch,
):
    system, typing, _, _, model = make(fallback_source)
    result = assess_pcff_resolution(
        system,
        typing,
        resolution_policy=POLICY,
        special_pairs=model.payload["special_pairs"],
    )
    from island.forcefields import PreparedForceFieldError, adopt_forcefield

    with pytest.raises(PreparedForceFieldError):
        adopt_forcefield(system, "pcff", result, source=fallback_source)
    p = result.payload
    assert p["diagnostic_only"] and p["native_model_complete"]
    assert p["structural_blocker_count"] == 0
    assert p["counts"]["existing_policy_derived_zero"] == 9
    assert p["native_identities"]["model_identity"] == model.identity
    path = tmp_path / "assessment.json"
    save_pcff_resolution_assessment(result, path)
    original = path.read_bytes()
    importer = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"openmm", "scipy", "rdkit", "foyer", "parmed"}:
            raise AssertionError(name)
        return importer(name, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(builtins, "__import__", blocked)
        restored = load_pcff_resolution_assessment(path, fallback_source, system=system)
        assert restored.json_text == result.json_text
    assert path.read_bytes() == original
    p["counts"].clear()
    assert result.payload["counts"]
    for mutate in ("policy", "classification", "charge"):
        data = unpack(result.json_text)
        if mutate == "policy":
            data["entries"][0]["policy"] = "invented"
        elif mutate == "classification":
            data["entries"][0]["structural_model_blocker"] = True
        else:
            data["charges_complete"] = False
        bad = PCFFResolutionAssessment(pack(data), fallback_source)
        with pytest.raises(PCFFError, match="Contradictory"):
            save_pcff_resolution_assessment(bad, tmp_path / f"bad-{mutate}.json")
        assert not (tmp_path / f"bad-{mutate}.json").exists()
    changed = deepcopy(system)
    changed.topology.sites[11].mass += 1
    with pytest.raises(PCFFError):
        result.validate_integrity(changed)
    with pytest.raises(PCFFError):
        save_pcff_resolution_assessment(result, path)
    assert path.read_bytes() == original
