"""Software-only driver lifecycle controls; do not replace native validation."""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from types import SimpleNamespace

import pytest

from island.charge_references.records import pack
from island.exceptions import PCFFError
from island.forcefields.pcff.charges import identity


@pytest.fixture
def driver(monkeypatch):
    monkeypatch.syspath_prepend("scripts")
    spec = spec_from_file_location(
        "variant_driver", "scripts/validate_pcff_variants.py"
    )
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_driver_uses_validated_loaders_once_without_identity_replay(
    driver, monkeypatch
):
    payload = {"explicitly_synthetic_software_fixture": True}

    class LoadedImmutableFixture:
        json_text = pack(payload)

        @property
        def identity(self):
            raise AssertionError("Redundant semantic validation after validated loader")

    calls = []

    def comparison(path, *, variants):
        calls.append(("comparison", path, variants))
        return LoadedImmutableFixture()

    def audit(path, *, assessment, variants):
        calls.append(("audit", path, assessment, variants))
        return LoadedImmutableFixture()

    monkeypatch.setattr(driver, "load_pcff_source_comparison", comparison)
    monkeypatch.setattr(driver, "load_pcff_source_variant_audit", audit)
    sources = ("explicit checked sources",)
    original = SimpleNamespace(label="synthetic owned assessment")
    assert driver._load_comparison_identity("comparison.json", sources) == identity(
        payload
    )
    assert driver._load_audit_identity("audit.json", original, sources) == identity(
        payload
    )
    assert calls == [
        ("comparison", "comparison.json", sources),
        ("audit", "audit.json", original, sources),
    ]


def test_driver_propagates_semantic_loader_failure(driver, monkeypatch):
    def rejected(*a, **kw):
        raise PCFFError("semantic rejection")

    monkeypatch.setattr(driver, "load_pcff_source_comparison", rejected)
    monkeypatch.setattr(driver, "load_pcff_source_variant_audit", rejected)
    with pytest.raises(PCFFError, match="semantic rejection"):
        driver._load_comparison_identity("bad.json", ())
    with pytest.raises(PCFFError, match="semantic rejection"):
        driver._load_audit_identity("bad.json", None, ())
