"""Synthetic contract tests; these do not establish real Foyer coverage."""

import subprocess
import sys
from copy import deepcopy

import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack
from island.exceptions import OPLSAssignmentError, OPLSDependencyError
from island.forcefields.oplsaa import (
    FoyerOPLSSource,
    OPLSChargeResult,
    OPLSTypingResult,
    adapter,
    assign_native_charges,
    load_oplsaa_source,
    source,
    type_atoms,
)
from island.forcefields.oplsaa.models import chemical_graph, typing_data


@pytest.fixture
def synthetic(monkeypatch):
    # Deliberately synthetic rule names/charges. Test-only pin substitution.
    xml = b"""<ForceField name="OPLS-AA" version="0.1.0" combining_rule="geometric">
    <AtomTypes><Type name="synthetic_C" element="C" def="C"/>
    <Type name="synthetic_H" element="H" def="H"/></AtomTypes>
    <NonbondedForce><Atom type="synthetic_C" charge="-0.3"/>
    <Atom type="synthetic_H" charge="0.1"/></NonbondedForce></ForceField>"""
    monkeypatch.setitem(source.PIN, "xml_sha256", source.digest(xml))
    graph = Topology()
    ids = [10, 30, 50, 70, 90, 110, 130, 150]
    for i, sid in enumerate(ids):
        graph.add_site(
            AtomSite(
                sid,
                str(sid),
                12.011 if i < 2 else 1.008,
                element="C" if i < 2 else "H",
                atomic_number=6 if i < 2 else 1,
            )
        )
    for a, b in [
        (10, 30),
        (10, 50),
        (10, 70),
        (10, 90),
        (30, 110),
        (30, 130),
        (30, 150),
    ]:
        graph.add_bond(a, b, order=1)
    system = MolecularSystem(graph, Coordinates({s: (0.0, 0.0, 0.0) for s in ids}))
    src = FoyerOPLSSource(xml)
    p = {
        "schema": "island_foyer_typing_v1",
        "source": src.identity,
        "graph": chemical_graph(system),
        "index_to_site_id": dict(enumerate(ids)),
        "environment": {"foyer": source.PIN["foyer_version"]},
        "matches": {},
        "data": {},
    }
    for sid in ids:
        t = "synthetic_" + graph.sites[sid].element
        p["matches"][sid] = {"whitelist": [t], "blacklist": [], "atomtype": t}
    p["data"] = typing_data(p, system, src)
    return system, src, OPLSTypingResult(pack(p))


def test_native_charge_ownership_and_compatibility(synthetic):
    s, src, typing = synthetic
    before = deepcopy(s)
    result = assign_native_charges(s, typing, src)
    result.validate_integrity(s, src)
    assert abs(result.payload["data"]["total_e"]) < 1e-15
    assert not result.payload["data"]["capabilities"]["energy_evaluation"]
    p = result.payload
    p["charges"].clear()
    assert len(result.payload["charges"]) == 8
    assert deepcopy(result).identity == result.identity
    assert s.to_dict() == before.to_dict()
    moved = deepcopy(s)
    moved.coordinates.set(10, (999.0, 0.0, 0.0))
    moved.topology.sites = dict(reversed(list(moved.topology.sites.items())))
    moved.topology.bonds = dict(reversed(list(moved.topology.bonds.items())))
    result.validate_integrity(moved, src)
    moved.topology.sites[10].metadata["cip_label"] = "R"
    with pytest.raises(OPLSAssignmentError, match="graph"):
        result.validate_integrity(moved, src)


@pytest.mark.parametrize(
    "defect", ["type", "missing", "ambiguous", "blacklist", "index", "source", "flag"]
)
def test_rechecksummed_typing_contradictions(synthetic, defect):
    s, src, t = synthetic
    p = t.payload
    if defect == "type":
        p["matches"][10]["atomtype"] = "synthetic_H"
    if defect == "missing":
        del p["matches"][10]
    if defect == "ambiguous":
        p["matches"][10]["whitelist"].append("synthetic_H")
    if defect == "blacklist":
        p["matches"][10]["blacklist"] = ["synthetic_H"]
    if defect == "index":
        p["index_to_site_id"][0] = 30
    if defect == "source":
        p["source"]["source_revision"] = "0" * 40
    if defect == "flag":
        p["data"]["production_validated"] = True
    with pytest.raises(OPLSAssignmentError):
        OPLSTypingResult(pack(p)).validate_integrity(s, src)


@pytest.mark.parametrize("value", [True, 0.2])
def test_native_charge_tamper(synthetic, value):
    s, src, t = synthetic
    p = assign_native_charges(s, t, src).payload
    p["charges"][50]["charge"] = value
    with pytest.raises(OPLSAssignmentError):
        OPLSChargeResult(pack(p)).validate_integrity(s, src)


def test_missing_record_and_non_neutral_no_repair(synthetic, monkeypatch):
    s, src, t = synthetic
    for replacement in [b"", b'<Atom type="synthetic_H" charge="0.2"/>']:
        xml = src.xml.replace(b'<Atom type="synthetic_H" charge="0.1"/>', replacement)
        monkeypatch.setitem(source.PIN, "xml_sha256", source.digest(xml))
        changed = FoyerOPLSSource(xml)
        p = t.payload
        p["source"] = changed.identity
        with pytest.raises(OPLSAssignmentError, match="Missing native|neutrality"):
            assign_native_charges(s, OPLSTypingResult(pack(p)), changed)


def test_invalid_source_dependency_and_graph(synthetic, tmp_path, monkeypatch):
    s, src, _ = synthetic
    path = tmp_path / "source.xml"
    path.write_bytes(src.xml)
    assert load_oplsaa_source(path).identity == src.identity
    path.write_bytes(src.xml + b" ")
    with pytest.raises(OPLSAssignmentError, match="hash"):
        load_oplsaa_source(path)
    monkeypatch.setattr(adapter.importlib.util, "find_spec", lambda _: None)
    with pytest.raises(OPLSDependencyError):
        type_atoms(s, src)
    s.topology.bonds.pop((10, 50))
    with pytest.raises(OPLSAssignmentError, match="valence"):
        type_atoms(s, src)


def test_mocked_upstream_maps_explicit_ids(synthetic, monkeypatch):
    s, src, _ = synthetic

    class Graph:
        def __init__(self):
            self.atoms = {}
            self.bonds = []

        def add_atom(self, i, name, **kw):
            self.atoms[i] = name

        def add_bond(self, i, j, **kw):
            self.bonds.append((i, j, kw))

    def find(g, provider, **kw):
        assert len(g.bonds) == 7 and all(b[2]["bond_order"] == 1 for b in g.bonds)
        return {
            i: {
                "atomtype": "synthetic_" + e,
                "whitelist": {"synthetic_" + e},
                "blacklist": set(),
            }
            for i, e in g.atoms.items()
        }

    monkeypatch.setattr(
        adapter,
        "check_foyer_installation",
        lambda: (
            lambda *args: args,
            find,
            Graph,
            {"foyer": source.PIN["foyer_version"]},
        ),
    )
    result = type_atoms(s, src)
    assert result.payload["index_to_site_id"][0] == 10
    assert assign_native_charges(s, result, src).payload["charges"][50]["charge"] == 0.1


def test_optional_imports():
    code = """
import sys
for n in ('foyer', 'gmso', 'rdkit', 'openmm', 'parmed', 'scipy'): sys.modules[n] = None
import island
import island.forcefields.oplsaa
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.mark.parametrize("invalid", ["not json", "{}", "NaN", "Infinity"])
def test_malformed_public_results(synthetic, invalid):
    s, src, _ = synthetic
    for cls in (OPLSTypingResult, OPLSChargeResult):
        with pytest.raises(OPLSAssignmentError):
            cls(invalid).validate_integrity(s, src)


def test_exclusive_validated_persistence(synthetic, tmp_path):
    from island.forcefields.oplsaa import load_opls_result, save_opls_result

    s, src, t = synthetic
    result = assign_native_charges(s, t, src)
    path = tmp_path / "charges.json"
    save_opls_result(result, s, src, path)
    assert load_opls_result(path, s, src).identity == result.identity
    original = path.read_bytes()
    with pytest.raises(OPLSAssignmentError):
        save_opls_result(result, s, src, path)
    assert path.read_bytes() == original
    p = result.payload
    p["data"]["production_validated"] = True
    with pytest.raises(OPLSAssignmentError):
        save_opls_result(OPLSChargeResult(pack(p)), s, src, tmp_path / "invalid.json")
    assert not (tmp_path / "invalid.json").exists()


@pytest.mark.parametrize(
    "message", ["Found no types for atom 0", "Found multiple types for atom 0"]
)
def test_structured_upstream_failure(synthetic, monkeypatch, message):
    from island.exceptions import OPLSTypingError

    s, src, _ = synthetic

    class Graph:
        def add_atom(self, *a, **kw):
            pass

        def add_bond(self, *a, **kw):
            pass

    def fail(*a, **kw):
        raise ValueError(message)

    monkeypatch.setattr(
        adapter,
        "check_foyer_installation",
        lambda: (lambda *a: None, fail, Graph, {"foyer": source.PIN["foyer_version"]}),
    )
    with pytest.raises(OPLSTypingError) as caught:
        type_atoms(s, src)
    assert caught.value.diagnostics["index_to_site_id"][0] == 10
    d = caught.value.diagnostics
    d["index_to_site_id"].clear()
    assert caught.value.diagnostics["upstream_error"] == message
    assert caught.value.diagnostics["index_to_site_id"]


def test_known_chemical_graph_changes_rejected(synthetic):
    s, src, t = synthetic
    from dataclasses import replace

    changed = deepcopy(s)
    changed.topology.bonds[(10, 30)] = replace(
        changed.topology.bonds[(10, 30)], order=2
    )
    with pytest.raises(OPLSAssignmentError):
        t.validate_integrity(changed, src)
    changed = deepcopy(s)
    changed.topology.sites[10].mass += 1
    with pytest.raises(OPLSAssignmentError):
        t.validate_integrity(changed, src)
    with pytest.raises(OPLSAssignmentError):
        t.validate_integrity(s, FoyerOPLSSource(src.xml + b" "))
