"""Content-cache software tests; no claimed chemical/reference coverage."""

from dataclasses import replace

import pytest
from test_pcff_operational_profile import ring, synthetic_profile_source  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield
from island.forcefields.pcff import PCFFModelSpecification
from island.forcefields.pcff import validation_cache as cache


@pytest.fixture
def cached_model(synthetic_profile_source):  # noqa: F811
    path, selection = synthetic_profile_source
    cache.clear_pcff_validation_cache()
    system = ring()
    options = PCFFOptions(
        path,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        source_profile=selection,
    )
    prepared = prepare_forcefield(system, ForceFieldRequest("pcff", options))
    return system, prepared.native_result


def test_hot_cache_keeps_live_graph_binding_and_coordinates(cached_model):
    system, model = cached_model
    model.validate_integrity(system)
    info = cache.pcff_validation_cache_info()
    original = model.identity
    for _ in range(5):
        model.validate_integrity(system)
    assert cache.pcff_validation_cache_info()["hits"] > info["hits"]
    for changed_field in ("element", "mass", "stereo", "bond", "id"):
        changed = system.copy()
        atom = changed.topology.sites[11]
        if changed_field == "element":
            atom.element, atom.atomic_number = "N", 7
        elif changed_field == "mass":
            atom.mass += 1
        elif changed_field == "stereo":
            atom.metadata["chiral_tag"] = "CHI_TETRAHEDRAL_CW"
        elif changed_field == "bond":
            changed.topology.remove_bond(11, 23)
        else:
            changed.topology.remove_site(11)
        with pytest.raises(PCFFError):
            model.validate_integrity(changed)
    moved = system.copy()
    moved.coordinates.translate([2, 3, 4])
    model.validate_integrity(moved)
    assert model.identity == original


def test_rechecksummed_content_and_nested_source_changes_miss_cache(cached_model):
    system, model = cached_model
    model.validate_integrity(system)
    data = unpack(model.json_text)
    data["terms"][0]["coefficients"][0] += 1
    with pytest.raises(PCFFError, match="Contradictory"):
        replace(model, json_text=pack(data)).validate_integrity(system)
    source = replace(model.assignment.source, expected_sha256="0" * 64)
    altered = replace(model, assignment=replace(model.assignment, source=source))
    with pytest.raises(PCFFError):
        altered.validate_integrity(system)
    nested = unpack(model.assignment.json_text)
    nested["assignments"][0]["normalized_values"][0] += 1
    with pytest.raises(PCFFError):
        replace(
            model, assignment=replace(model.assignment, json_text=pack(nested))
        ).validate_integrity(system)


def test_trusted_rule_function_changes_invalidate_success(cached_model, monkeypatch):
    system, model = cached_model
    from island.forcefields.pcff import model as module

    model.validate_integrity(system)
    old = module.definition

    def changed(*args, **kwargs):
        result = old(*args, **kwargs)
        result["model_definition_complete"] = False
        return result

    monkeypatch.setattr(module, "definition", changed)
    with pytest.raises(PCFFError, match="Contradictory"):
        model.validate_integrity(system)


def test_summary_ownership_and_byte_budget(cached_model, monkeypatch):
    system, model = cached_model
    first = cache.profile_validation_summary(model)
    first["model"]["terms"].clear()
    assert cache.profile_validation_summary(model)["model"]["terms"]
    monkeypatch.setattr(cache, "MAX_BYTES", 1)
    cache.clear_pcff_validation_cache()
    model.validate_integrity(system)
    assert cache.pcff_validation_cache_info()["entries"] == 0
    assert model.identity


def test_historical_identity_and_source_inventory_ownership(cached_model):
    system, model = cached_model
    from island.forcefields.pcff.charges import identity

    expected = identity(unpack(model.json_text))
    model.validate_integrity(system)
    assert model.identity == expected
    source = model.assignment.source
    inv = source.inventory
    inv["sections"].clear()
    assert source.inventory["sections"]
    assert (
        PCFFModelSpecification(model.json_text, model.assignment).identity == expected
    )
