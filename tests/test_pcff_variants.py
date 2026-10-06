"""Synthetic source-variant software fixtures, never scientific parameter data."""

import builtins
from copy import deepcopy
from dataclasses import replace

import pytest

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFSourceComparison,
    PCFFVariantPolicy,
    PCFFVariantResolution,
    PCFFVariantSelection,
    compare_pcff_sources,
    load_pcff_source_comparison,
    load_pcff_source_variant,
    load_pcff_variant_resolution,
    resolve_pcff_variant_record,
    save_pcff_source_comparison,
    save_pcff_variant_resolution,
    select_native_pcff_source,
)
from island.forcefields.pcff import source as historical
from island.forcefields.pcff.source import PCFFSource, digest
from island.forcefields.pcff.variants import CANDIDATE_PROFILE

RAW = b"""!BIOSYM forcefield 1
#version synthetic.frc 1.0 synthetic
#atom_types cff91
1 1 c 12 C 4 synthetic carbon
1 1 h 1 H 1 synthetic hydrogen
1 1 o 16 O 2 synthetic oxygen
#equivalence cff91
1 1 c c c c c c
1 1 h h h h h h
1 1 o o o o o o
#auto_equivalence cff91_auto
1 1 c c c c_ c_ c_ c_ c_ c_ c_
1 1 h h h h_ h_ h_ h_ h_ h_ h_
1 1 o o o o_ o_ o_ o_ o_ o_ o_
#bond_increments cff91_auto
1 1 c h -0.1 0.1
#quartic_bond cff91
1 1 c h 1.1 2 -3 4
#quadratic_angle cff91_auto
1 1 h_ c_ *2 120 37.5
1 1 *4 c_ o_ 120 40
#bond-angle cff91
1 1 h c o 2 3
#angle-angle cff91
1 1 h c h o 1
#wilson_out_of_plane cff91
1 1 h c h o 1 0
#torsion-torsion_1 cff91
#unknown_term alternate
1 1 c h 0.5
#end
"""


def provenance(label):
    return {
        "repository": "synthetic software fixture",
        "revision": label,
        "source_path": label + ".frc",
        "family": "unreviewed PCFF candidate",
        "license_status": "test fixture",
        "citation": "synthetic, no scientific evidence",
    }


def variant(tmp_path, raw=RAW, label="a"):
    path = tmp_path / (label + ".frc")
    path.write_bytes(raw)
    return load_pcff_source_variant(
        path, expected_sha256=digest(raw), provenance=provenance(label)
    )


def selection(v):
    return PCFFVariantSelection(**v.selection)


def test_inventory_differences_namespace_versions_and_ownership(tmp_path):
    a = variant(tmp_path)
    b = variant(
        tmp_path,
        RAW.replace(b"1 1 c h -0.1 0.1", b"1 1 c h -0.2 0.2\n2 1 c o 0 0"),
        "b",
    )
    c = compare_pcff_sources([a, b])
    assert c.identity == compare_pcff_sources([b, a]).identity
    p = c.payload
    assert any(
        s["name"] == "torsion-torsion_1" and s["record_count"] == 0
        for s in a.payload["sections"]
    )
    assert any(
        s["namespace"] == "alternate" and s["record_count"] == 1
        for s in a.payload["sections"]
    )
    d = p["differences"][0]
    assert sorted([len(d["left_only"]), len(d["right_only"])]) == [1, 2]
    p["variants"][0]["counts"].clear()
    assert c.payload["variants"][0]["counts"]
    original = provenance("owned")
    path = tmp_path / "owned.frc"
    path.write_bytes(RAW)
    owned = load_pcff_source_variant(
        path, expected_sha256=digest(RAW), provenance=original
    )
    original["revision"] = "changed"
    assert owned.payload["provenance"]["revision"] == "owned"
    with pytest.raises(PCFFError, match="hash mismatch"):
        load_pcff_source_variant(path, expected_sha256="0" * 64, provenance=original)
    with pytest.raises(PCFFError):
        compare_pcff_sources([a, a])
    # Repeated ordered candidates remain recorded; they are not flattened into one row.
    repeated = variant(
        tmp_path,
        RAW.replace(b"1 1 c h 1.1 2 -3 4", b"1 1 c h 1.1 2 -3 4\n2 1 c h 1.2 3 -4 5"),
        "repeated",
    )
    assert repeated.payload["counts"]["quartic_bond"] == 2


def test_explicit_selection_no_implicit_mixing_and_fallback(tmp_path):
    a = variant(tmp_path)
    b = variant(tmp_path, RAW.replace(b"1 1 c h -0.1 0.1", b"1 1 c o 0.2 -0.2"), "b")
    primary = PCFFVariantPolicy(selection(a))
    missing = resolve_pcff_variant_record(
        [a, b],
        policy=primary,
        family="bond_increments",
        namespace="cff91_auto",
        types=["c", "o"],
        sites=[11, 89],
    )
    assert missing.payload["classification"] == "source_row_absent"
    assert len(missing.payload["attempts"]) == 1
    explicit = PCFFVariantPolicy(selection(a), (selection(b),))
    resolved = resolve_pcff_variant_record(
        [a, b],
        policy=explicit,
        family="bond_increments",
        namespace="cff91_auto",
        types=["c", "o"],
        sites=[11, 89],
    )
    p = resolved.payload
    assert p["classification"] == "source_variant_only"
    assert p["selected_source"] == b.selection
    assert not p["mixed_source_model"] and not p["model_assembly_authorized"]
    reverse = resolve_pcff_variant_record(
        [b, a],
        policy=explicit,
        family="bond_increments",
        namespace="cff91_auto",
        types=["o", "c"],
        sites=[89, 11],
    )
    assert reverse.payload["attempts"][-1]["forward"]["normalized_values"] == [
        -0.2,
        0.2,
    ]
    assert reverse.payload["attempts"][-1]["physical_reversal_invariant"]
    for policy in [
        PCFFVariantPolicy(PCFFVariantSelection(a.expected_sha256, "wrong")),
        PCFFVariantPolicy(selection(a), (selection(a),)),
        PCFFVariantPolicy(
            selection(a), (PCFFVariantSelection("0" * 64, CANDIDATE_PROFILE),)
        ),
    ]:
        with pytest.raises(PCFFError):
            resolve_pcff_variant_record(
                [a, b], policy=policy, family="quartic_bond", types=["c", "h"]
            )
    with pytest.raises(PCFFError, match="semantics not audited"):
        select_native_pcff_source([a, b], selection(b))
    with pytest.raises(PCFFError):
        resolve_pcff_variant_record(
            [a, b], policy=primary, family="quartic_bond", types="ch"
        )
    with pytest.raises(PCFFError):
        resolve_pcff_variant_record(
            [a, b],
            policy=primary,
            family="quartic_bond",
            types=["c", "h"],
            sites=[True, 9],
        )


def test_wildcard_reversal_gate_center_and_asymmetric_coefficients(tmp_path):
    a = variant(tmp_path)
    p = PCFFVariantPolicy(selection(a))
    ambiguous = resolve_pcff_variant_record(
        [a],
        policy=p,
        family="quadratic_angle",
        namespace="cff91_auto",
        types=["h", "c", "o"],
    )
    assert ambiguous.payload["classification"] == "wildcard_ambiguous"
    assert not ambiguous.payload["attempts"][0]["physical_reversal_invariant"]
    with pytest.raises(PCFFError):
        resolve_pcff_variant_record(
            [a],
            policy=replace(p, resolution_policy="numeric_wildcard_priority"),
            family="quadratic_angle",
            namespace="cff91_auto",
            types=["h", "c", "o"],
        )
    directional = resolve_pcff_variant_record(
        [a], policy=p, family="bond-angle", types=["h", "c", "o"]
    )
    assert directional.payload["classification"] == "source_row_present"
    q = directional.payload["attempts"][0]
    assert q["physical_reversal_invariant"]
    assert q["forward"]["normalized_values"] == [2 * 4.184, 3 * 4.184]
    assert q["reverse"]["normalized_values"] == [3 * 4.184, 2 * 4.184]
    moved = resolve_pcff_variant_record(
        [a], policy=p, family="angle-angle", types=["c", "h", "h", "o"]
    )
    assert moved.payload["classification"] == "source_row_absent"
    identical = variant(tmp_path, RAW.replace(b"120 40", b"120 37.5"), "identical")
    same = resolve_pcff_variant_record(
        [identical],
        policy=PCFFVariantPolicy(selection(identical)),
        family="quadratic_angle",
        namespace="cff91_auto",
        types=["h", "c", "o"],
    )
    assert same.payload["classification"] == "source_row_present"
    assert same.payload["attempts"][0]["physical_reversal_invariant"]
    missing = resolve_pcff_variant_record(
        [a], policy=p, family="bond-bond", types=["h", "c", "o"]
    )
    assert missing.payload["classification"] == "source_row_absent"


def test_persistence_offline_and_rechecksummed_semantic_tampering(
    tmp_path, monkeypatch
):
    a = variant(tmp_path)
    comparison = compare_pcff_sources([a])
    query = resolve_pcff_variant_record(
        [a],
        policy=PCFFVariantPolicy(selection(a)),
        family="quartic_bond",
        types=["c", "h"],
    )
    cp = tmp_path / "comparison.json"
    rp = tmp_path / "resolution.json"
    save_pcff_source_comparison(comparison, cp)
    save_pcff_variant_resolution(query, rp)
    old = (cp.read_bytes(), rp.read_bytes())
    importer = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"openmm", "scipy", "rdkit", "foyer", "parmed"}:
            raise AssertionError(name)
        return importer(name, *args, **kwargs)

    with monkeypatch.context() as m:
        m.setattr(builtins, "__import__", blocked)
        assert (
            load_pcff_source_comparison(cp, variants=[a]).identity
            == comparison.identity
        )
        assert load_pcff_variant_resolution(rp, variants=[a]).identity == query.identity
    assert old == (cp.read_bytes(), rp.read_bytes())
    for kind in ["comparison", "query"]:
        data = unpack(comparison.json_text if kind == "comparison" else query.json_text)
        if kind == "comparison":
            data["variants"][0]["rows"][0]["types"] = ["o"]
        else:
            data["attempts"][0]["forward"]["normalized_values"][1] = 999
        bad = (PCFFSourceComparison if kind == "comparison" else PCFFVariantResolution)(
            pack(data), (a,)
        )
        save = (
            save_pcff_source_comparison
            if kind == "comparison"
            else save_pcff_variant_resolution
        )
        with pytest.raises(PCFFError, match="Contradictory"):
            save(bad, tmp_path / (kind + "-bad.json"))
        assert not (tmp_path / (kind + "-bad.json")).exists()
    with pytest.raises(PCFFError):
        save_pcff_source_comparison(comparison, cp)
    assert old[0] == cp.read_bytes()
    for bad in ['{"x":1,"x":2}', '{"x":NaN}', "{"]:
        (tmp_path / "malformed.json").write_text(bad)
        with pytest.raises(PCFFError):
            load_pcff_source_comparison(tmp_path / "malformed.json", variants=[a])
    with pytest.raises(PCFFError):
        load_pcff_source_comparison(
            cp, variants=[variant(tmp_path, RAW + b"\n", "different-bytes")]
        )


def test_known_native_source_identity_unchanged(tmp_path, monkeypatch):
    monkeypatch.setitem(historical.PIN, "sha256", digest(RAW))
    path = tmp_path / "native.frc"
    path.write_bytes(RAW)
    prov = provenance("native")
    prov["family"] = "PCFF"
    v = load_pcff_source_variant(path, expected_sha256=digest(RAW), provenance=prov)
    native = select_native_pcff_source([v], selection(v))
    assert native.identity == PCFFSource(RAW, digest(RAW)).identity
    assert v.payload["native_assignment_authorized"]
    assert deepcopy(v.payload) == v.payload


def test_nonzero_wilson_and_conflicting_charge_cannot_be_promoted(tmp_path):
    wilson = variant(
        tmp_path, RAW.replace(b"1 1 h c h o 1 0", b"1 1 h c h o 1 5"), "wilson"
    )
    q = resolve_pcff_variant_record(
        [wilson],
        policy=PCFFVariantPolicy(selection(wilson)),
        family="wilson_out_of_plane",
        types=["h", "c", "h", "o"],
    )
    assert q.payload["classification"] == "model_incomplete"
    assert q.payload["attempts"][0]["row_present"]
    assert not q.payload["attempts"][0]["physical_reversal_invariant"]
    charge = variant(
        tmp_path,
        RAW.replace(b"1 1 c h -0.1 0.1", b"1 1 c h -0.1 0.1\n1 1 h c 0.2 -0.2"),
        "charge",
    )
    q = resolve_pcff_variant_record(
        [charge],
        policy=PCFFVariantPolicy(selection(charge)),
        family="bond_increments",
        namespace="cff91_auto",
        types=["c", "h"],
    )
    assert q.payload["classification"] == "native_charge_incomplete"
    assert q.payload["selected_source"] is None
    assert q.payload["attempts"][0]["forward"]["candidate_ids"]


def test_cli_inspection_is_not_full_source_acceptance(tmp_path):
    import json
    import os
    import subprocess
    import sys

    a = variant(tmp_path)
    declaration = tmp_path / "sources.json"
    declaration.write_text(
        json.dumps(
            {
                "sources": [
                    {
                        "path": str(tmp_path / "a.frc"),
                        "sha256": a.expected_sha256,
                        "profile": CANDIDATE_PROFILE,
                        "provenance": provenance("a"),
                    }
                ]
            }
        )
    )
    env = dict(os.environ, PYTHONPATH="src")
    args = [
        sys.executable,
        "-m",
        "island.forcefields.pcff.variants_cli",
        "--declaration",
        str(declaration),
        "--output",
        str(tmp_path / "cli"),
    ]
    result = subprocess.run(
        args + ["--require-full-source"],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    assert result.returncode == 1, result.stderr
    report = unpack((tmp_path / "cli/report.json").read_text())
    assert report["inspection_complete"] and not report["full_source_complete"]
    before = (tmp_path / "cli/report.json").read_bytes()
    result = subprocess.run(args, capture_output=True, text=True, env=env, check=False)
    assert result.returncode != 0
    assert (tmp_path / "cli/report.json").read_bytes() == before
