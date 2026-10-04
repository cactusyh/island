"""Synthetic source contract tests, not real PCFF parameter/typing evidence."""

import os
import subprocess
import sys
from copy import deepcopy

import pytest

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFChargeResult,
    assign_pcff_charges,
    assign_pcff_types,
    load_pcff_record,
    load_pcff_source,
    save_pcff_record,
)
from island.forcefields.pcff import source as module
from island.forcefields.pcff.source import PCFFSource, digest, parse_frc, select

# No authentic source data is redistributed by this deliberately synthetic fixture.
SYNTHETIC = b"""!BIOSYM forcefield 1
#version synthetic.frc 1.0 synthetic
#atom_types cff91
1.0 1 c 12 C 4 synthetic
1.0 1 c3 12 C 4 synthetic
1.0 1 c2 12 C 4 synthetic
1.0 1 h 1 H 1 synthetic
1.0 1 hc 1 H 1 synthetic
1.0 1 ho 1 H 1 synthetic
1.0 1 o 16 O 2 synthetic
1.0 1 oh 16 O 2 synthetic
1.0 1 oc 16 O 2 synthetic
#equivalence cff91
1.0 1 c3 c c c c c
1.0 1 c2 c c c c c
1.0 1 hc h h h h h
1.0 1 ho h h h h h
1.0 1 oh o o o o o
1.0 1 oc o o o o o
#auto_equivalence cff91_auto
1.0 1 c3 c c c_ c_ c_ c_ c_ c_ c_
#bond_increments cff91_auto
1.0 1 c c 0 0
1.0 1 c h -0.1 0.1
1.0 1 c o 0.2 -0.2
1.0 1 h o 0.3 -0.3
#future_cross_term cff91
1.0 1 c h h 4.0 5.0
#end
"""


@pytest.fixture
def fixture(monkeypatch):
    monkeypatch.setitem(module.PIN, "sha256", digest(SYNTHETIC))
    return PCFFSource(SYNTHETIC, digest(SYNTHETIC))


def molecule(case="ethanol", reverse=False):
    # Example is an explicit graph builder; no scientific dependency or coordinate generation.
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("pcff_example", "examples/pcff_charges.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.molecule(case, reverse=reverse)


def run(source, case="ethanol", reverse=False):
    system, types = molecule(case, reverse)
    typing = assign_pcff_types(
        system, source, types, provenance="synthetic software fixture"
    )
    return system, typing, assign_pcff_charges(system, typing)


def changed_source(monkeypatch, raw):
    monkeypatch.setitem(module.PIN, "sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


def test_analytical_charge_ownership_orientation(fixture):
    system, typing, result = run(fixture)
    for sid in system.topology.sites:
        system.coordinates.set(sid, [float(sid), 0, 0])
    before = deepcopy(system.to_dict())
    # C(H3)-C(H2)-O-H: (-.3, 0, -.5, five *.1, .3).
    assert result.complete
    assert list(result.charges.values()) == pytest.approx(
        [-0.3, 0, -0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3]
    )
    assert abs(result.payload["total_charge"]) < 1e-12
    _, _, reverse = run(fixture, reverse=True)
    assert reverse.identity == result.identity
    result.charges.clear()
    result.payload["typing"]["types"].clear()
    assert result.complete and len(typing.payload["types"]) == 9
    assert system.to_dict() == before
    result.validate_integrity(system)
    system.topology.sites[10].mass += 1
    with pytest.raises(PCFFError, match="graph mismatch"):
        result.validate_integrity(system)


def test_inventory_and_hash(fixture, tmp_path):
    inv = fixture.inventory
    assert any(
        s["name"] == "future_cross_term"
        and s["interpretation"] == "raw_only"
        and "4.0 5.0" in s["lines"][0]["raw"]
        for s in inv["sections"]
    )
    path = tmp_path / "source.frc"
    path.write_bytes(SYNTHETIC)
    assert load_pcff_source(path).identity == fixture.identity
    with pytest.raises(PCFFError, match="hash mismatch"):
        load_pcff_source(path, expected_sha256="0" * 64)


@pytest.mark.parametrize(
    "change",
    [
        lambda b: b.replace(b"c h -0.1 0.1", b"c h -0.1"),
        lambda b: b.replace(b"c h -0.1 0.1", b"c h nan 0.1"),
        lambda b: b.replace(b"#end", b""),
        lambda b: b.replace(b"c3 c c c c c", b"c3 c c"),
        lambda b: b.replace(b"12 C 4 synthetic", b"12 C"),
    ],
)
def test_malformed_supported_rows(change):
    with pytest.raises(PCFFError):
        parse_frc(change(SYNTHETIC))


def test_versions_conflicts_and_zero(fixture, monkeypatch):
    raw = SYNTHETIC.replace(
        b"1.0 1 c h -0.1 0.1", b"0.5 1 c h -0.7 0.7\n1.0 1 c h -0.1 0.1"
    )
    source = changed_source(monkeypatch, raw)
    _, _, result = run(source, "ethane")
    assert result.complete and result.charges[10] == pytest.approx(-0.3)
    assert result.payload["contributions"][0]["increments"] == [0, 0]
    raw = raw.replace(b"1.0 1 c h -0.1 0.1", b"1.0 1 c h -0.1 0.1\n1.0 1 c h -0.2 0.2")
    _, _, failed = run(changed_source(monkeypatch, raw), "ethane")
    assert not failed.complete
    assert any(
        d["reason"] == "ambiguous_increment" for d in failed.payload["diagnostics"]
    )
    with pytest.raises(PCFFError, match="Incomplete"):
        _ = failed.charges


def test_missing_not_zero_or_auto_fallback(fixture, monkeypatch):
    raw = SYNTHETIC.replace(b"1.0 1 c c 0 0\n", b"")
    _, _, result = run(changed_source(monkeypatch, raw), "ethane")
    assert not result.complete
    assert result.payload["diagnostics"][0]["reason"] == "missing_increment"


def test_component_mismatch(fixture, monkeypatch):
    raw = SYNTHETIC.replace(b"c h -0.1 0.1", b"c h -0.1 0.11")
    _, _, result = run(changed_source(monkeypatch, raw), "ethane")
    assert not result.complete
    assert any(
        d["reason"] == "component_charge_mismatch"
        for d in result.payload["diagnostics"]
    )


@pytest.mark.parametrize(
    "kind", ["missing", "unknown", "element", "parent", "bool_id", "missing_h"]
)
def test_invalid_typing(fixture, kind):
    system, labels = molecule()
    if kind == "missing":
        labels.pop(10)
    elif kind == "unknown":
        labels[10] = "NO_TYPE"
    elif kind == "element":
        labels[10] = "hc"
    elif kind == "parent":
        labels[max(labels)] = "hc"
    elif kind == "bool_id":
        labels[True] = "c3"
    else:
        system.topology.remove_site(max(labels))
        labels.pop(max(labels))
    with pytest.raises(PCFFError):
        assign_pcff_types(system, fixture, labels, provenance="synthetic")


def test_data_only_persistence_tampering(fixture, tmp_path):
    system, typing, result = run(fixture)
    for record in (typing, result):
        path = tmp_path / record.payload["schema"]
        save_pcff_record(record, path)
        assert (
            load_pcff_record(path, fixture, system=system).identity == record.identity
        )
        before = path.read_bytes()
        with pytest.raises(PCFFError):
            save_pcff_record(record, path)
        assert path.read_bytes() == before
    for key, value in [
        ("complete", False),
        ("total_charge", 9),
        ("policy", "repair"),
        ("production_validated", True),
        ("typing_identity", "0" * 64),
    ]:
        p = result.payload
        p[key] = value
        bad = PCFFChargeResult(pack(p), fixture)
        with pytest.raises(PCFFError):
            save_pcff_record(bad, tmp_path / "bad")
        assert not (tmp_path / "bad").exists()
    p = result.payload
    p["partial_charges"][10] += 0.01
    with pytest.raises(PCFFError):
        PCFFChargeResult(pack(p), fixture).validate_integrity()
    path = tmp_path / "invalid"
    path.write_text('{"payload":{},"payload":{}}')
    with pytest.raises(PCFFError):
        load_pcff_record(path, fixture)


def test_reverse_source_orientation(fixture, monkeypatch):
    _, _, a = run(fixture)
    raw = SYNTHETIC.replace(b"c h -0.1 0.1", b"h c 0.1 -0.1")
    _, _, b = run(changed_source(monkeypatch, raw))
    # Source changes, numerical assignment does not.
    assert a.json_text != b.json_text
    assert unpack(a.json_text)["partial_charges"] == b.charges


def test_equal_version_conflict_does_not_use_file_order():
    rows = parse_frc(SYNTHETIC)["sections"][1]["records"]
    a = deepcopy(rows[0])
    b = deepcopy(a)
    b["data"]["mass"] += 1
    for candidates in ([a, b], [b, a]):
        with pytest.raises(PCFFError, match="Conflicting"):
            select(candidates)


def test_no_optional_imports(tmp_path):
    src = tmp_path / "synthetic.frc"
    src.write_bytes(SYNTHETIC)
    code = """
import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self, fullname, *a):
  if fullname.split('.')[0] in {'rdkit','openmm','foyer','parmed','scipy'}:
   raise RuntimeError('forbidden import '+fullname)
sys.meta_path.insert(0,Block())
from island.forcefields.pcff import load_pcff_source
from island.forcefields.pcff.source import digest
from pathlib import Path
p=Path(sys.argv[1]);s=load_pcff_source(p,expected_sha256=digest(p.read_bytes()))
assert s.identity['profile'] is None
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(src)],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    assert result.returncode == 0, result.stderr


def test_cli(fixture, tmp_path):
    from island.forcefields.pcff.__main__ import main

    src = tmp_path / "synthetic.frc"
    src.write_bytes(SYNTHETIC)
    assert main([str(src), "--full"]) == 0
    assert main([str(src), "--sha256", "0" * 64]) == 1


def test_component_cancellation_not_success(monkeypatch):
    # Two components have opposite defects; molecular neutrality cannot hide them.
    raw = SYNTHETIC.replace(b"c h -0.1 0.1", b"c h -0.1 0.11").replace(
        b"#end", b"\n#bond_increments cff91_auto\n1.0 1 c3 hc -0.11 0.1\n#end"
    )
    source = changed_source(monkeypatch, raw)
    left, labels = molecule("ethane")
    right, right_labels = molecule("ethane")
    for site in right.topology.sites.values():
        copied = deepcopy(site)
        copied.id += 1000
        left.topology.add_site(copied)
        labels[copied.id] = "c" if right_labels[site.id] == "c3" else "h"
    for bond in right.topology.bonds.values():
        left.topology.add_bond(bond.site1 + 1000, bond.site2 + 1000, order=1)
    result = assign_pcff_charges(
        left, assign_pcff_types(left, source, labels, provenance="synthetic")
    )
    p = result.payload
    assert abs(p["total_charge"]) < 1e-12
    assert not p["complete"]
    assert (
        sum(d["reason"] == "component_charge_mismatch" for d in p["diagnostics"]) == 2
    )


def test_auto_equivalence_only_blocked(monkeypatch):
    raw = SYNTHETIC.replace(b"1.0 1 c c 0 0", b"1.0 1 X X 0 0").replace(
        b"c3 c c c_ c_ c_ c_ c_ c_ c_", b"c3 c X c_ c_ c_ c_ c_ c_ c_"
    )
    _, _, result = run(changed_source(monkeypatch, raw), "ethane")
    assert not result.complete
    assert (
        result.payload["diagnostics"][0]["auto_equivalence_evidence"][0]["record"][
            "data"
        ]["families"]["bond_increment"]
        == "X"
    )


def test_direct_precedence_and_forward_conflict(monkeypatch):
    raw = SYNTHETIC.replace(
        b"#end", b"#bond_increments cff91_auto\n1.0 1 c3 hc -0.2 0.2\n#end"
    )
    _, _, result = run(changed_source(monkeypatch, raw), "ethane")
    assert result.charges[10] == pytest.approx(-0.6)
    raw = raw.replace(
        b"#end", b"#bond_increments cff91_auto\n1.0 1 hc c3 0.3 -0.3\n#end"
    )
    _, _, result = run(changed_source(monkeypatch, raw), "ethane")
    assert not result.complete
    assert any(
        "forward/reverse" in d.get("message", "") for d in result.payload["diagnostics"]
    )


def test_source_mismatch_and_rechecksummed_contributions(fixture, tmp_path):
    system, labels = molecule()
    other_raw = SYNTHETIC.replace(b"synthetic.frc", b"other.frc")
    other = PCFFSource(other_raw, digest(other_raw))
    with pytest.raises(PCFFError, match="not audited"):
        assign_pcff_types(system, other, labels, provenance="synthetic")
    _, _, result = run(fixture)
    p = result.payload
    p["contributions"][0]["source_matches"][0]["selection"]["record"]["line"] += 1
    path = tmp_path / "tampered"
    path.write_text(pack(p))
    with pytest.raises(PCFFError, match="Contradictory"):
        load_pcff_record(path, fixture)
    with pytest.raises(PCFFError, match="mapping"):
        assign_pcff_types(system, fixture, list(labels.items()), provenance="synthetic")


def test_optional_isolation_assignment_and_load(fixture, tmp_path):
    _, _, result = run(fixture)
    record_path = tmp_path / "record.json"
    save_pcff_record(result, record_path)
    source_path = tmp_path / "source.frc"
    source_path.write_bytes(SYNTHETIC)
    code = """
import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self, fullname, *a):
  if fullname.split('.')[0] in {'rdkit','openmm','foyer','parmed','scipy'}:
   raise RuntimeError('forbidden import '+fullname)
sys.meta_path.insert(0,Block())
from pathlib import Path
from island.forcefields.pcff import load_pcff_source, load_pcff_record
from island.forcefields.pcff.source import PIN, digest
p=Path(sys.argv[1]);PIN['sha256']=digest(p.read_bytes()) # explicitly synthetic test pin
s=load_pcff_source(p)
r=load_pcff_record(sys.argv[2],s)
assert r.complete and len(r.charges)==9
"""
    proc = subprocess.run(
        [sys.executable, "-c", code, str(source_path), str(record_path)],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    assert proc.returncode == 0, proc.stderr


def test_acceptance_cli_failure_preserves_report(tmp_path):
    proc = subprocess.run(
        [
            sys.executable,
            "scripts/validate_pcff_charges.py",
            "--source",
            str(tmp_path / "missing"),
            "--output",
            str(tmp_path / "experiment"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 1, proc.stderr
    import json

    p = json.loads((tmp_path / "experiment/report.json").read_text())
    assert not p["passed"] and len(p["outcomes"]) == 3
    assert all(row["status"] == "failed" for row in p["outcomes"])
