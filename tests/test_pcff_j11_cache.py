"""Synthetic cache controls; no new chemical-domain or numerical evidence."""

import pytest
from test_pcff_validation_cache import (  # noqa: F401
    cached_model,
    synthetic_profile_source,
)

from island.exceptions import PCFFError
from island.forcefields.pcff import model as model_module
from island.forcefields.pcff import source as source_module


def test_in_place_rule_code_change_invalidates_success(cached_model, monkeypatch):  # noqa: F811
    system, model = cached_model
    model.validate_integrity(system)

    def changed(assignment, policy):
        raise ValueError("J11 changed rule executed")

    monkeypatch.setattr(model_module.definition, "__code__", changed.__code__)
    with pytest.raises(PCFFError, match="J11 changed rule executed"):
        model.validate_integrity(system)


@pytest.mark.parametrize("kind", ["defaults", "kwdefaults", "closure"])
def test_mutable_rule_state_invalidates_success(cached_model, monkeypatch, kind):  # noqa: F811
    system, model = cached_model
    original = model_module.definition
    state = {"reject": False}

    def default_rule(assignment, policy, state=state):
        if state["reject"]:
            raise ValueError("J11 mutable rule changed")
        return original(assignment, policy)

    def keyword_rule(assignment, policy, *, state=state):
        if state["reject"]:
            raise ValueError("J11 mutable rule changed")
        return original(assignment, policy)

    def closure_rule(assignment, policy):
        if state["reject"]:
            raise ValueError("J11 mutable rule changed")
        return original(assignment, policy)

    rule = {
        "defaults": default_rule,
        "kwdefaults": keyword_rule,
        "closure": closure_rule,
    }[kind]
    monkeypatch.setattr(model_module, "definition", rule)
    model.validate_integrity(system)
    state["reject"] = True
    with pytest.raises(PCFFError, match="J11 mutable rule changed"):
        model.validate_integrity(system)


def test_in_place_parser_code_change_invalidates_source(cached_model, monkeypatch):  # noqa: F811
    _, model = cached_model
    source = model.assignment.source
    source.validate_integrity()

    def changed(token):
        raise ValueError("J11 changed parser executed")

    monkeypatch.setattr(source_module.number, "__code__", changed.__code__)
    with pytest.raises(PCFFError, match="J11 changed parser executed"):
        source.validate_integrity()


def test_nested_policy_change_invalidates_success(cached_model, monkeypatch):  # noqa: F811
    system, model = cached_model
    model.validate_integrity(system)
    monkeypatch.setitem(model_module.PROFILE["equations"], "quartic_bond", "invalid")
    with pytest.raises(PCFFError, match="Contradictory"):
        model.validate_integrity(system)


def test_rehashed_source_bytes_cannot_reuse_pinned_validation(cached_model):  # noqa: F811
    from dataclasses import replace
    from hashlib import sha256

    system, model = cached_model
    model.validate_integrity(system)
    original = model.assignment.source
    raw = original.raw.replace(b"cp h -0.1 0.1", b"cp h -0.2 0.2")
    assert raw != original.raw
    changed = replace(original, raw=raw, expected_sha256=sha256(raw).hexdigest())
    altered = replace(model, assignment=replace(model.assignment, source=changed))
    with pytest.raises(PCFFError):
        altered.validate_integrity(system)


def test_callable_state_snapshots_cycles_without_user_methods():
    from island.forcefields.pcff._validation_state import callable_state

    class Untrusted:
        def __repr__(self):
            raise AssertionError("Must not call user repr during cache lookup")

    state = {"external": Untrusted()}
    state["cycle"] = state

    def rule():
        return state

    original = callable_state(rule)
    assert callable_state(rule) == original
    state["policy"] = "new"
    assert callable_state(rule) != original
