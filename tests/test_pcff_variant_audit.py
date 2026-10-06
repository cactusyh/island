"""Synthetic cross-source audit contracts; no operational alternate-source claim."""

import builtins

import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFSourceVariant,
    PCFFSourceVariantAudit,
    assess_pcff_resolution,
    audit_pcff_source_variants,
    load_pcff_source_variant_audit,
    save_pcff_source_variant_audit,
)
from island.forcefields.pcff.fallbacks import POLICY
from island.forcefields.pcff.source import digest


def provenance():
    return {
        "repository": "synthetic test",
        "revision": "fixture",
        "source_path": "synthetic.frc",
        "family": "PCFF",
        "license_status": "test-only",
        "citation": "software fixture, no scientific evidence",
    }


def test_audit_preserves_existing_model_and_rejects_alternate_semantics(
    fallback_source,  # noqa: F811
):
    system, typing, _, _, model = make(fallback_source)
    assessment = assess_pcff_resolution(
        system,
        typing,
        resolution_policy=POLICY,
        special_pairs=model.payload["special_pairs"],
    )
    original = PCFFSourceVariant(
        fallback_source.raw, fallback_source.expected_sha256, pack(provenance())
    )
    raw = fallback_source.raw + b"\n"
    alternate = PCFFSourceVariant(raw, digest(raw), pack(provenance()))
    result = audit_pcff_source_variants(assessment, [original, alternate])
    p = result.payload
    assert (
        p["original_assessment"]["native_identities"]["model_identity"]
        == model.identity
    )
    pinned = next(
        r for r in p["variant_results"] if r["candidate_model_semantics_authorized"]
    )
    candidate = next(
        r for r in p["variant_results"] if not r["candidate_model_semantics_authorized"]
    )
    assert pinned["native_model_complete"]
    assert pinned["classification_counts"]["policy_derived_zero"] == 9
    assert candidate["classification_counts"].get("policy_derived_zero", 0) == 0
    assert not candidate["native_model_complete"]
    assert not p["operational_gate"] and not p["full_source_complete"]
    assert candidate["charge_audit"]["source_charge_calculation_complete"]
    assert not candidate["charge_audit"]["native_charge_validation_authorized"]
    assert (
        assessment.payload["native_model_complete"]
        and model.identity
        == p["original_assessment"]["native_identities"]["model_identity"]
    )


def test_audit_persistence_tampering_and_lazy_loading(
    fallback_source,  # noqa: F811
    tmp_path,
    monkeypatch,
):
    system, typing, _, _, model = make(fallback_source)
    assessment = assess_pcff_resolution(
        system,
        typing,
        resolution_policy=POLICY,
        special_pairs=model.payload["special_pairs"],
    )
    source = PCFFSourceVariant(
        fallback_source.raw, fallback_source.expected_sha256, pack(provenance())
    )
    result = audit_pcff_source_variants(assessment, [source])
    path = tmp_path / "audit.json"
    save_pcff_source_variant_audit(result, path)
    before = path.read_bytes()
    importer = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"rdkit", "openmm", "scipy", "foyer", "parmed"}:
            raise AssertionError(name)
        return importer(name, *args, **kwargs)

    with monkeypatch.context() as m:
        m.setattr(builtins, "__import__", blocked)
        restored = load_pcff_source_variant_audit(
            path, assessment=assessment, variants=[source]
        )
        assert restored.identity == result.identity
    for mutate in ["charge", "dependency", "policy"]:
        p = unpack(result.json_text)
        if mutate == "charge":
            p["variant_results"][0]["charge_audit"]["charges_e"][11] = "1"
        elif mutate == "dependency":
            p["variant_results"][0]["entries"][-1]["dependencies"] = []
        else:
            p["variant_results"][0]["candidate_model_semantics_authorized"] = False
        bad = PCFFSourceVariantAudit(pack(p), assessment, (source,))
        with pytest.raises(PCFFError, match="Contradictory"):
            save_pcff_source_variant_audit(bad, tmp_path / (mutate + ".json"))
    with pytest.raises(PCFFError):
        save_pcff_source_variant_audit(result, path)
    assert path.read_bytes() == before
