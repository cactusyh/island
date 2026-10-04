"""Synthetic software contracts; real source acceptance is a separate script."""

import json
import os
import subprocess
import sys
from copy import deepcopy

import pytest
from test_pcff import SYNTHETIC, molecule

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFAutomaticChargeResult,
    PCFFAutomaticTypingResult,
    assign_automatic_pcff_charges,
    automatic,
    load_pcff_automatic_record,
    save_pcff_automatic_record,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.source import PCFFSource, digest


@pytest.fixture
def synthetic(monkeypatch):
    raw = SYNTHETIC.replace(
        b"1.0 1 c2 12 C 4 synthetic",
        b"1.0 1 c2 12 C 4 synthetic\n1.0 1 c1 12 C 4 synthetic",
    ).replace(
        b"#equivalence cff91",
        b"#equivalence cff91\n1.0 1 c c c c c c\n1.0 1 c1 c c c c c",
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


def explicit(carbon_bonds, hcounts):
    topology = Topology()
    for i in range(len(hcounts)):
        topology.add_site(
            AtomSite(11 + 5 * i, "C", 12.011, element="C", atomic_number=6)
        )
    for a, b in carbon_bonds:
        topology.add_bond(11 + 5 * a, 11 + 5 * b, order=1)
    next_id = 100
    for i, count in enumerate(hcounts):
        for _ in range(count):
            topology.add_site(
                AtomSite(next_id, "H", 1.008, element="H", atomic_number=1)
            )
            topology.add_bond(11 + 5 * i, next_id, order=1)
            next_id += 3
    return MolecularSystem(topology, Coordinates())


def test_carbon_rules_precedence(synthetic):
    for system, expected in [
        (explicit([(0, 1), (1, 2)], [3, 2, 3]), ["c3", "c2", "c3"]),
        (explicit([(0, 1), (1, 2), (1, 3)], [3, 1, 3, 3]), ["c3", "c1", "c3", "c3"]),
        (
            explicit([(0, 1), (0, 2), (0, 3), (0, 4)], [0, 3, 3, 3, 3]),
            ["c", "c3", "c3", "c3", "c3"],
        ),
        (explicit([], [4]), ["c"]),
    ]:
        result = type_pcff_atoms(system, synthetic)
        assert result.complete
        assert [
            result.assignments[11 + 5 * i] for i in range(len(expected))
        ] == expected
        for sid, label in result.assignments.items():
            entry = result.payload["entries"][sid]
            if label in ("c1", "c2", "c3"):
                assert entry["matched_rules"] == [f"C.H{label[-1]}", "C.sp3"]
                assert entry["overridden_rules"] == ["C.sp3"]
        assert assign_automatic_pcff_charges(system, result).complete


def test_ethanol_bridge_metadata_and_ownership(synthetic):
    system, manual = molecule()
    before = deepcopy(system.topology.sites)
    result = type_pcff_atoms(system, synthetic)
    assert result.assignments == manual
    charges = assign_automatic_pcff_charges(system, result)
    assert charges.complete
    assert charges.charges[10] == pytest.approx(-0.3)
    assert (
        "ISLAND automatic typing"
        in charges.payload["native_charge_record"]["typing"]["provenance"]
    )
    assert charges.payload["automatic_typing_identity"] == result.identity
    result.payload["entries"].clear()
    result.assignments.clear()
    charges.charges.clear()
    assert system.topology.sites == before
    for sid, atom in system.topology.sites.items():
        atom.metadata.update(
            repeat_unit_index=900, coordinate_seed=99, arbitrary="ignored"
        )
        system.coordinates.set(sid, [sid, 2, 3])
    system.metadata["coordinate_history"] = ["irrelevant"]
    assert type_pcff_atoms(system, synthetic).identity == result.identity
    charges.validate_integrity(system)
    system.topology.sites[10].metadata["chiral_tag"] = "CHI_TETRAHEDRAL_CW"
    with pytest.raises(PCFFError, match="graph mismatch"):
        result.validate_integrity(system)


def test_id_remapping_orientation(synthetic):
    system, _ = molecule()
    original = type_pcff_atoms(system, synthetic)
    mapping = {sid: 3000 - 13 * sid for sid in system.topology.sites}
    topology = Topology()
    for sid in reversed(list(system.topology.sites)):
        atom = deepcopy(system.topology.sites[sid])
        atom.id = mapping[sid]
        topology.add_site(atom)
    for bond in reversed(list(system.topology.bonds.values())):
        topology.add_bond(mapping[bond.site2], mapping[bond.site1], order=1)
    other = MolecularSystem(topology, Coordinates())
    changed = type_pcff_atoms(other, synthetic)
    assert changed.assignments == {
        mapping[s]: t for s, t in original.assignments.items()
    }
    q = assign_automatic_pcff_charges(system, original).charges
    assert assign_automatic_pcff_charges(other, changed).charges == {
        mapping[s]: v for s, v in q.items()
    }


@pytest.mark.parametrize(
    "smiles,reason",
    [
        ("C1CCCCC1", "cyclic_component"),
        ("CC=O", "unsaturated_or_aromatic_bond"),
        ("c1ccccc1", "aromaticity"),
        ("COOC", "peroxide"),
        ("O", "water"),
        ("COC(OC)C", "multiple_oxygen_carbon"),
        ("CC(O)OC", "multiple_oxygen_carbon"),
        ("CN", "unsupported_element"),
        ("C[O-]", "charged_site"),
        ("[13CH3]C", "nonstandard_mass_or_isotope"),
        ("[1H]C", "nonstandard_mass_or_isotope"),
        ("[CH3]", "incomplete_hydrogens_or_valence"),
    ],
)
def test_unsupported_before_generic(synthetic, smiles, reason):
    from island.chemistry import from_smiles

    system = from_smiles(smiles, random_seed=2026)
    result = type_pcff_atoms(system, synthetic)
    assert not result.complete
    assert not result.payload["assignments"]  # whole offending component withheld
    assert reason in {d["reason"] for d in result.payload["diagnostics"]}
    with pytest.raises(PCFFError, match="Incomplete"):
        assign_automatic_pcff_charges(system, result)


def test_missing_hydrogens_and_malformed(synthetic):
    system, _ = molecule()
    system.topology.remove_site(max(system.topology.sites))
    result = type_pcff_atoms(system, synthetic)
    assert not result.complete
    system.topology.sites[10].mass = float("nan")
    with pytest.raises(PCFFError, match="mass"):
        type_pcff_atoms(system, synthetic)
    system.topology.sites[10].mass = 12.011
    system.topology.sites[10].atomic_number = 8
    with pytest.raises(PCFFError, match="Inconsistent"):
        type_pcff_atoms(system, synthetic)


def test_typing_complete_charge_incomplete(synthetic, monkeypatch):
    raw = synthetic.raw.replace(b"1.0 1 c o 0.2 -0.2\n", b"")
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", digest(raw))
    src = PCFFSource(raw, digest(raw))
    system, _ = molecule()
    result = type_pcff_atoms(system, src)
    assert result.complete
    charges = assign_automatic_pcff_charges(system, result)
    assert not charges.complete
    assert "missing_increment" in {
        d["reason"] for d in charges.payload["native_charge_record"]["diagnostics"]
    }
    with pytest.raises(PCFFError):
        _ = charges.charges


def test_ambiguity_not_iteration_precedence(synthetic, monkeypatch):
    extra = deepcopy(automatic.PROFILE["rules"][1])
    extra.update(id="conflict", type="c")
    monkeypatch.setitem(
        automatic.PROFILE, "rules", automatic.PROFILE["rules"] + [extra]
    )
    system = explicit([(0, 1), (1, 2), (1, 3)], [3, 1, 3, 3])
    result = type_pcff_atoms(system, synthetic)
    assert not result.complete
    assert result.payload["entries"][16]["status"] == "ambiguous"


def test_integrity_persistence_source_profile(synthetic, tmp_path):
    system, _ = molecule()
    auto = type_pcff_atoms(system, synthetic)
    charge = assign_automatic_pcff_charges(system, auto)
    for result in [auto, charge]:
        path = tmp_path / result.payload["schema"]
        save_pcff_automatic_record(result, path)
        before = path.read_bytes()
        assert (
            load_pcff_automatic_record(path, synthetic, system=system).identity
            == result.identity
        )
        with pytest.raises(PCFFError):
            save_pcff_automatic_record(result, path)
        assert path.read_bytes() == before
    for key in ["assignments", "entries", "profile", "coverage"]:
        p = auto.payload
        if key == "assignments":
            p[key][10] = "c2"
        if key == "entries":
            p[key][10]["environment"]["hydrogens"] = 2
        if key == "profile":
            p[key]["hydroxyl_convention"] = "silently_use_ho2"
        if key == "coverage":
            p[key]["typed"] = 99
        with pytest.raises(PCFFError, match="Contradictory"):
            PCFFAutomaticTypingResult(pack(p), synthetic).validate_integrity()
    p = charge.payload
    p["automatic_typing_identity"] = "0" * 64
    with pytest.raises(PCFFError):
        save_pcff_automatic_record(
            PCFFAutomaticChargeResult(pack(p), synthetic), tmp_path / "bad"
        )
    assert not (tmp_path / "bad").exists()
    with pytest.raises(PCFFError, match="profile"):
        type_pcff_atoms(system, synthetic, profile="generic_pcff")
    bad = PCFFSource(synthetic.raw + b"\n", digest(synthetic.raw + b"\n"))
    with pytest.raises(PCFFError):
        type_pcff_atoms(system, bad)


def test_offline_optional_isolation(synthetic, tmp_path):
    system, _ = molecule()
    record = assign_automatic_pcff_charges(system, type_pcff_atoms(system, synthetic))
    p = tmp_path / "record"
    save_pcff_automatic_record(record, p)
    src = tmp_path / "source"
    src.write_bytes(synthetic.raw)
    profile = tmp_path / "profile"
    profile.write_text(json.dumps(automatic.PROFILE))
    code = """
import sys,importlib.abc,json
from pathlib import Path
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,*args):
  if name.split('.')[0] in {'rdkit','scipy','foyer','parmed','openmm'}:raise RuntimeError(name)
sys.meta_path.insert(0,Block())
from island.forcefields.pcff import source, automatic, load_pcff_source, load_pcff_automatic_record
raw=Path(sys.argv[1]).read_bytes();source.PIN['sha256']=source.digest(raw)
automatic.PROFILE=json.loads(Path(sys.argv[3]).read_text()) # synthetic test profile
r=load_pcff_automatic_record(sys.argv[2],load_pcff_source(sys.argv[1]));assert r.complete
assert len(r.charges)==9
"""
    proc = subprocess.run(
        [sys.executable, "-c", code, str(src), str(p), str(profile)],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    assert proc.returncode == 0, proc.stderr


def test_incomplete_record_roundtrip(synthetic, tmp_path):
    system = explicit([(0, 1), (1, 2), (2, 0)], [2, 2, 2])
    result = type_pcff_atoms(system, synthetic)
    save_pcff_automatic_record(result, tmp_path / "diagnostic")
    restored = load_pcff_automatic_record(tmp_path / "diagnostic", synthetic)
    assert not restored.complete
    assert unpack(result.json_text) == restored.payload


@pytest.mark.parametrize(
    "repeat,dp,counts",
    [
        ("[*:1]CC[*:2]", 1, {"c3": 2, "hc": 6}),
        ("[*:1]CC[*:2]", 3, {"c3": 2, "c2": 4, "hc": 14}),
        ("[*:1]CCO[*:2]", 1, {"c3": 1, "c2": 1, "oh": 1, "hc": 5, "ho": 1}),
        ("[*:1]CCO[*:2]", 3, {"c3": 1, "c2": 5, "oh": 1, "oc": 2, "hc": 13, "ho": 1}),
    ],
)
def test_built_chain_ends_without_repeat_input(synthetic, repeat, dp, counts):
    from collections import Counter

    from island.builders import build_linear_polymer

    system = build_linear_polymer(repeat, dp=dp, generate_3d=False)
    result = type_pcff_atoms(system, synthetic)
    assert dict(Counter(result.assignments.values())) == counts
    for atom in system.topology.sites.values():
        for key in (
            "repeat_unit_index",
            "repeat_unit_type",
            "source_repeat_atom_index",
            "chain_id",
        ):
            atom.metadata.pop(key, None)
    system.metadata.clear()
    assert type_pcff_atoms(system, synthetic).identity == result.identity
    assert assign_automatic_pcff_charges(system, result).complete


def test_cli_report_and_exclusive_output(synthetic, tmp_path, capsys):
    from island.forcefields.pcff.typing_cli import main

    src = tmp_path / "source"
    src.write_bytes(synthetic.raw)
    output = tmp_path / "out"
    assert (
        main(
            [
                "--source",
                str(src),
                "--psmiles",
                "[*:1]CCO[*:2]",
                "--dp",
                "1",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    report = json.loads((output / "report.json").read_text())
    assert report["typing_coverage"]["complete"]
    assert report["charge_coverage"]["complete"]
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    assert main(["--source", str(src), "--smiles", "CCO", "--output", str(output)]) == 1
    assert {p.name: p.read_bytes() for p in output.iterdir()} == before
    assert main(["--source", str(src), "--smiles", "O"]) == 1
    # Last report has valid unsupported typing, no charge coverage.
    output = capsys.readouterr().out
    assert "water" in output and '"charge_coverage": null' in output


def test_no_decoding_of_executable_or_bad_json(synthetic, tmp_path):
    for text in ("not json", "{}", '{"payload":{},"payload":{}}', '{"payload":NaN}'):
        path = tmp_path / "invalid"
        path.write_text(text)
        with pytest.raises(PCFFError):
            load_pcff_automatic_record(path, synthetic)


def test_acceptance_missing_source_nonzero(tmp_path):
    proc = subprocess.run(
        [
            sys.executable,
            "scripts/validate_pcff_automatic.py",
            "--source",
            str(tmp_path / "missing"),
            "--output",
            str(tmp_path / "run"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 1
    report = json.loads((tmp_path / "run/report.json").read_text())
    assert not report["passed"] and report["outcomes"][0]["status"] == "failed"
