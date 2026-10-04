"""Explicitly synthetic Class II software contracts, not PCFF scientific data."""

from copy import deepcopy

import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFClass2Result,
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    automatic,
    inspect_pcff_class2,
    load_pcff_parameters,
    save_pcff_parameters,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.class2 import FAMILIES, conversion, resolve
from island.forcefields.pcff.source import PCFFSource, digest


@pytest.fixture
def class2_source(synthetic, monkeypatch):  # noqa: F811
    # Deliberately sparse: missing cross terms stay diagnostic, never zero-filled.
    extra = b"""#quartic_bond cff91
1.0 1 c c 1.5 100 -200 300
1.0 1 c h 1.1 50 -60 70
#quartic_angle cff91
1.0 1 c c c 110 40 -5 6
1.0 1 c c h 109 30 -4 5
1.0 1 h c h 108 20 -3 4
#torsion_3 cff91
1.0 1 h c c h 1 0 -2 180 3 0
1.0 1 c c c h 1 0 -2 180 3 0
#nonbond(9-6) cff91
1.0 1 c 4 0.1
1.0 1 h 2 0.01
#bond-bond cff91
1.0 1 h c h 0
#bond-angle cff91
1.0 1 h c h 2
#end_bond-torsion_3 cff91
1.0 1 h c c h 1 -2 0
#angle-angle cff91
1.0 1 h c h h -1
1.0 1 c c h h -2
1.0 1 h c c h -3
"""
    raw = synthetic.raw.replace(b"#end", extra + b"#end")
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


def assigned(system, src):
    t = type_pcff_atoms(system, src)
    c = assign_automatic_pcff_charges(system, t)
    return assign_pcff_parameters(system, t, c)


def test_units_and_native_signed_coefficients(class2_source):
    inspected = inspect_pcff_class2(class2_source)
    inspected["conventions"]["mixing"].clear()
    assert inspect_pcff_class2(class2_source)["conventions"]["mixing"]
    catalog = inspect_pcff_class2(class2_source)["records"]
    bond = next(r for r in catalog.values() if r["section"] == "quartic_bond")
    assert bond["normalized_values"] == pytest.approx(
        [1.5, 418.4, -836.8, 1255.2], abs=1e-12
    )
    assert bond["source_units"] == ["angstrom", "E/L^2", "E/L^3", "E/L^4"]
    assert conversion("E/A^3") == 4.184  # NOT multiplied by degrees/radian cubed
    assert conversion("degree") == pytest.approx(0.017453292519943295)
    short = next(r for r in catalog.values() if r["section"] == "end_bond-torsion_3")
    assert short["expanded_values"] == [1, -2, 0, 1, -2, 0]
    assert short["original_values"] == [1, -2, 0]


def test_authoritative_inventory_dependencies_and_ownership(class2_source):
    m = explicit([(0, 1)], [3, 3])
    before = deepcopy(m)
    p = assigned(m, class2_source)
    data = p.payload
    assert len(data["inventories"]["bonds"]) == 7
    assert len(data["inventories"]["angles"]) == 12
    assert len(data["inventories"]["proper_torsions"]) == 9
    assert data["coverage"]["angle-angle"] == {"assigned": 24}
    assert data["coverage"]["wilson_out_of_plane"] == {"not_applicable": 8}
    assert data["coverage"]["bond-bond_1_3"] == {"missing": 9}
    assert (
        not data["parameter_coverage_complete"] and not data["physical_model_complete"]
    )
    zero = next(
        a
        for a in data["assignments"]
        if a["family"] == "bond-bond" and a["status"] == "assigned"
    )
    assert zero["normalized_values"] == [0]
    assert [d["equilibrium_value"] for d in zero["dependencies"]] == [1.1, 1.1]
    assert all(d["source_rows"] for d in zero["dependencies"])
    assert automatic.chemical_graph(m) == automatic.chemical_graph(before)
    assert m.metadata == before.metadata
    data["assignments"].clear()
    assert p.payload["assignments"]
    # Caches and coordinates do not enter assignment.
    m.topology.angles.clear()
    m.topology.dihedrals.clear()
    from island.core.topology import Angle, Dihedral, Improper

    # Deliberately impossible cached interactions must not become authoritative.
    m.topology.angles[(900, 901, 902)] = Angle(900, 901, 902)
    m.topology.dihedrals[(900, 901, 902, 903)] = Dihedral(900, 901, 902, 903)
    m.topology.impropers.append(Improper(900, 901, 902, 903))
    for sid in m.topology.sites:
        m.coordinates.set(sid, (sid, 2, 3))
    m.topology.sites = dict(reversed(list(m.topology.sites.items())))
    m.topology.bonds = dict(reversed(list(m.topology.bonds.items())))
    assert assigned(m, class2_source).identity == p.identity


def row(family, labels, values, rid="row", version="1"):
    return {
        "id": rid,
        "section": family,
        "types": labels,
        "normalized_values": values,
        "version": version,
    }


@pytest.mark.parametrize(
    "family,labels,values,expected",
    [
        ("bond-angle", ["a", "b", "c"], [1, 2], [2, 1]),
        (
            "end_bond-torsion_3",
            ["a", "b", "c", "d"],
            [1, 2, 3, 4, 5, 6],
            [4, 5, 6, 1, 2, 3],
        ),
        (
            "angle-torsion_3",
            ["a", "b", "c", "d"],
            [1, 2, 3, 4, 5, 6],
            [4, 5, 6, 1, 2, 3],
        ),
        ("middle_bond-torsion_3", ["a", "b", "c", "d"], [1, 2, 3], [1, 2, 3]),
    ],
)
def test_asymmetric_reversal(family, labels, values, expected):
    selected = resolve(family, labels[::-1], {"row": row(family, labels, values)}, [])
    assert selected["normalized_values"] == expected
    assert selected["selected"][0]["orientation"] == "reverse"


def test_angle_angle_permutation_keeps_center_and_shared_arm():
    catalog = {"row": row("angle-angle", ["a", "center", "shared", "d"], [-3])}
    assert resolve("angle-angle", ["d", "center", "shared", "a"], catalog, [])[
        "normalized_values"
    ] == [-3]
    assert (
        resolve("angle-angle", ["d", "shared", "center", "a"], catalog, [])["status"]
        == "missing"
    )


def test_repeated_labels_and_conflicts():
    r = row("bond-angle", ["h", "c", "h"], [1, 2])
    assert resolve("bond-angle", ["h", "c", "h"], {"r": r}, [])["status"] == "ambiguous"
    r["normalized_values"] = [1, 1]
    other = row("bond-angle", ["h", "c", "h"], [3, 3], "other")
    assert (
        resolve("bond-angle", ["h", "c", "h"], {"r": r, "o": other}, [])["status"]
        == "ambiguous"
    )
    other["version"] = "2"
    assert resolve("bond-angle", ["h", "c", "h"], {"r": r, "o": other}, [])[
        "normalized_values"
    ] == [3, 3]


def test_strict_persistence_and_rechecksummed_tampering(class2_source, tmp_path):
    m = explicit([(0, 1)], [3, 3])
    p = assigned(m, class2_source)
    path = tmp_path / "parameters.json"
    save_pcff_parameters(p, path)
    assert load_pcff_parameters(path, class2_source, system=m).identity == p.identity
    with pytest.raises(PCFFError):
        save_pcff_parameters(p, path)
    for key in ("coverage", "conventions", "inventories"):
        data = unpack(p.json_text)
        data[key] = {}
        with pytest.raises(PCFFError, match="Contradictory"):
            PCFFClass2Result(pack(data), class2_source).validate_integrity()
    data = unpack(p.json_text)
    data["assignments"][0]["normalized_values"][0] += 1
    with pytest.raises(PCFFError):
        PCFFClass2Result(pack(data), class2_source).validate_integrity()
    m.topology.sites[11].mass += 1
    with pytest.raises(PCFFError):
        p.validate_integrity(m)


def test_truncated_source_record(class2_source):
    raw = class2_source.raw.replace(b"1.5 100 -200 300", b"1.5 100")
    src = PCFFSource(raw, digest(raw))
    with pytest.raises(PCFFError, match="quartic_bond line"):
        inspect_pcff_class2(src)


def test_family_units_explicit():
    assert len(FAMILIES) == 13
    assert all(len(v[2]) > 0 for v in FAMILIES.values())


def test_missing_equilibrium_dependency(class2_source, monkeypatch):
    raw = class2_source.raw.replace(b"1.0 1 c h 1.1 50 -60 70\n", b"")
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", digest(raw))
    src = PCFFSource(raw, digest(raw))
    p = assigned(explicit([(0, 1)], [3, 3]), src).payload
    affected = [
        a
        for a in p["assignments"]
        if a["family"] == "bond-bond" and a.get("normalized_values") == [0]
    ]
    assert affected and all(a["status"] == "missing" for a in affected)
    assert all(a["reason"] == "unresolved equilibrium dependency" for a in affected)


def test_noncontiguous_remapping(class2_source):
    from island import AtomSite, Coordinates, MolecularSystem, Topology

    m = explicit([(0, 1)], [3, 3])
    first = assigned(m, class2_source).payload
    ids = {sid: 1000 - 7 * sid for sid in m.topology.sites}
    t = Topology()
    for sid, atom in reversed(list(m.topology.sites.items())):
        t.add_site(
            AtomSite(
                ids[sid],
                atom.name,
                atom.mass,
                element=atom.element,
                atomic_number=atom.atomic_number,
            )
        )
    for a, b in reversed(list(m.topology.bonds)):
        t.add_bond(ids[b], ids[a], order=1)
    second = assigned(MolecularSystem(t, Coordinates()), class2_source).payload
    assert first["coverage"] == second["coverage"]
    # Physical coefficient multiset is invariant under an arbitrary ID remapping.
    values = lambda p: sorted(
        (a["family"], tuple(a.get("normalized_values", []))) for a in p["assignments"]
    )
    assert values(first) == values(second)


def test_source_and_charge_mismatch(class2_source):
    m = explicit([(0, 1)], [3, 3])
    other = explicit([(0, 1), (1, 2)], [3, 2, 3])
    t = type_pcff_atoms(m, class2_source)
    c = assign_automatic_pcff_charges(other, type_pcff_atoms(other, class2_source))
    with pytest.raises(PCFFError):
        assign_pcff_parameters(m, t, c)
    p = assigned(m, class2_source)
    wrong = PCFFSource(class2_source.raw + b"\n", digest(class2_source.raw + b"\n"))
    with pytest.raises(PCFFError):
        PCFFClass2Result(p.json_text, wrong).validate_integrity()


def test_offline_load_without_scientific_dependencies(class2_source, tmp_path):
    import os
    import subprocess
    import sys

    result = assigned(explicit([(0, 1)], [3, 3]), class2_source)
    path = tmp_path / "parameters.json"
    path.write_text(result.json_text)
    src = tmp_path / "synthetic.frc"
    src.write_bytes(class2_source.raw)
    code = """
import sys
class Block:
    def find_spec(self, name, *args):
        if name.split('.')[0] in {'rdkit','openmm','foyer','parmed','scipy'}:
            raise RuntimeError('unexpected scientific import '+name)
sys.meta_path.insert(0,Block())
from pathlib import Path
from island.forcefields.pcff import source,automatic,load_pcff_parameters
raw=Path(sys.argv[1]).read_bytes()
sha=source.digest(raw)
source.PIN['sha256']=sha
automatic.PROFILE['source_sha256']=sha
p=load_pcff_parameters(sys.argv[2],source.PCFFSource(raw,sha))
assert not p.payload['physical_model_complete']
"""
    done = subprocess.run(
        [sys.executable, "-c", code, str(src), str(path)],
        check=False,
        env=os.environ.copy(),
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, done.stderr


@pytest.mark.parametrize("require_complete,exit_code", [(False, 0), (True, 1)])
def test_cli_preserves_diagnostic_and_reports_coverage_gate(
    class2_source, tmp_path, monkeypatch, require_complete, exit_code
):
    from island.forcefields.pcff import class2_cli, save_pcff_automatic_record

    m = explicit([(0, 1)], [3, 3])
    typing = type_pcff_atoms(m, class2_source)
    charges = assign_automatic_pcff_charges(m, typing)
    path = tmp_path / "charges.json"
    save_pcff_automatic_record(charges, path)
    output = tmp_path / "parameters.json"
    monkeypatch.setattr(class2_cli, "load_pcff_source", lambda path: class2_source)
    args = ["--source", "synthetic", "--charges", str(path), "--output", str(output)]
    if require_complete:
        args.append("--require-complete")
    assert class2_cli.main(args) == exit_code
    assert not load_pcff_parameters(output, class2_source).payload[
        "parameter_coverage_complete"
    ]
    original = output.read_bytes()
    assert class2_cli.main(args) == 1
    assert output.read_bytes() == original
