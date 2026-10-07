"""Synthetic auxiliary grammar/integrity tests; local real audit is separate."""

import builtins
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from test_pcff_automatic import synthetic  # noqa: F401

from island.charge_references.records import pack
from island.exceptions import PCFFError
from island.forcefields.pcff.distribution import (
    PCFFDistributionAudit,
    inspect_pcff_distribution,
    inspect_pcff_templates,
    load_pcff_distribution_audit,
    save_pcff_distribution_audit,
)
from island.forcefields.pcff.source import digest

TEMPLATE = b"""! synthetic envelope, not a chemical library
 type: hn
 template: (>H (-N))
 atom_test: 1
 aromaticity: NON_AROMATIC
 end_test
 end_type
 type hn
 template (>H (-N))
 end_type
 precedence:
 (? (hn(hn)))
 end_precedence
"""


@pytest.fixture
def paths(tmp_path, synthetic):  # noqa: F811
    values = {"frc": synthetic.raw, "rlb": b"VERSION\nelib\n", "templates": TEMPLATE}
    paths = {}
    for name, raw in values.items():
        p = tmp_path / name
        p.write_bytes(raw)
        paths[name + "_path"] = p
    return paths, {k: digest(v) for k, v in values.items()}


def test_lossless_structural_inventory_no_promotion():
    r = inspect_pcff_templates(TEMPLATE)
    assert r["counts"] == {"records": 2, "labels": 1, "patterns": 2, "tests": 1}
    assert r["duplicate_labels"] == {"hn": 2}
    assert r["records"][0]["id"] != r["records"][1]["id"]
    assert r["precedence"]["tree"] == [["?", ["hn", ["hn"]]]]
    assert not r["operational_authority"]
    assert r["diagnostics"][0]["reason"] == "non_colon_type_syntax_retained"
    assert not r["records"][1]["patterns"][0]["colon"]


@pytest.mark.parametrize(
    "raw",
    [
        TEMPLATE[:-17],
        b"\x00",
        b"\xff",
        b"precedence:\n()\nend_precedence",
        TEMPLATE.replace(b"atom_test: 1", b"atom_test: 0"),
        TEMPLATE.replace(b"(? (hn(hn)))", b"(? (hn(hn))"),
    ],
)
def test_malformed_auxiliary_boundary(raw):
    with pytest.raises(PCFFError):
        inspect_pcff_templates(raw)


def test_uninterpreted_conditions_are_data_only(tmp_path):
    raw = TEMPLATE.replace(
        b"aromaticity: NON_AROMATIC", b"arbitrary: __import__('os').system('false')"
    )
    r = inspect_pcff_templates(raw)
    assert r["records"][0]["tests"][0]["conditions"][0]["key"] == "arbitrary"
    assert not r["operational_authority"]


def test_owned_persistence_relocation_and_semantic_tampering(paths, tmp_path):
    paths, hashes = paths
    result = inspect_pcff_distribution(**paths, expected_hashes=hashes)
    payload = result.payload
    payload["runtime_authorization"] = True
    assert not result.payload["runtime_authorization"]
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFDistributionAudit(pack(payload)).validate_integrity(**paths)
    payload = result.payload
    payload["templates"]["records"][0]["label"] = "invented"
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFDistributionAudit(pack(payload)).validate_integrity(**paths)
    saved = tmp_path / "audit.json"
    save_pcff_distribution_audit(saved, result, **paths)
    before = saved.read_bytes()
    with pytest.raises(PCFFError):
        save_pcff_distribution_audit(saved, result, **paths)
    assert saved.read_bytes() == before
    moved = tmp_path / "relocated"
    moved.mkdir()
    for name, p in list(paths.items()):
        target = moved / p.name
        p.rename(target)
        paths[name] = target
    assert load_pcff_distribution_audit(saved, **paths).identity == result.identity
    paths["rlb_path"].write_bytes(b"VERSION\nchanged\n")
    with pytest.raises(PCFFError, match="hash mismatch"):
        load_pcff_distribution_audit(saved, **paths)


def test_explicit_dependencies_and_lazy_loading(paths, monkeypatch):
    paths, hashes = paths
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        assert name.split(".")[0] not in {"rdkit", "openmm", "scipy", "foyer", "parmed"}
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    r = inspect_pcff_distribution(**paths, expected_hashes=hashes)
    r.validate_integrity(**paths)
    hashes["rlb"] = "0" * 64
    with pytest.raises(PCFFError, match="hash mismatch"):
        inspect_pcff_distribution(**paths, expected_hashes=hashes)


def test_child_reconstruction_without_scientific_imports(paths, tmp_path, monkeypatch):
    paths, hashes = paths
    from island.forcefields.pcff.source import PIN

    monkeypatch.setitem(
        PIN,
        "sha256",
        "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c",
    )
    record = inspect_pcff_distribution(**paths, expected_hashes=hashes)
    saved = tmp_path / "audit.json"
    save_pcff_distribution_audit(saved, record, **paths)
    script = """import builtins,json,sys
orig=builtins.__import__
def blocked(name,globals=None,locals=None,fromlist=(),level=0):
    assert level or name.split('.')[0] not in {'openmm','rdkit','scipy','foyer','parmed'}
    return orig(name,globals,locals,fromlist,level)
builtins.__import__=blocked
from island.forcefields.pcff.distribution import load_pcff_distribution_audit
r=load_pcff_distribution_audit(sys.argv[1],**json.loads(sys.argv[2]))
print(r.identity)
"""
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            script,
            str(saved),
            json.dumps({k: str(v) for k, v in paths.items()}),
        ],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONPATH": str(Path("src").resolve())},
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == record.identity
