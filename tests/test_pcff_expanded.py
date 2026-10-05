"""Synthetic software regressions. Real pinned-source acceptance is separate."""

from copy import deepcopy
from itertools import permutations

import numpy as np
import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island import Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_source_types,
    automatic,
    expanded,
    load_pcff_automatic_record,
    save_pcff_automatic_record,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.source import PCFFSource, digest
from island.forcefields.pcff.terms import SourceClass2Term, finite_difference_forces


@pytest.fixture
def expanded_source(synthetic, monkeypatch):  # noqa: F811
    raw = synthetic.raw.replace(
        b"#equivalence cff91",
        b"1.0 1 c3h 12 C 4 synthetic\n1.0 1 c4h 12 C 4 synthetic\n#equivalence cff91\n1.0 1 c3h c c c c c\n1.0 1 c4h c c c c c",
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


def test_ring_expansion_preserves_historical_rejection(expanded_source):
    for size, label in [(3, "c3h"), (4, "c4h"), (6, "c2")]:
        m = explicit([(i, (i + 1) % size) for i in range(size)], [2] * size)
        before = deepcopy(m)
        old = type_pcff_atoms(m, expanded_source)
        assert not old.complete
        t = type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME)
        assert t.complete
        assert all(t.assignments[11 + 5 * i] == label for i in range(size))
        q = assign_automatic_pcff_charges(m, t)
        assert q.complete
        assert sum(q.charges.values()) == pytest.approx(0, abs=1e-12)
        assert m.topology.sites == before.topology.sites
        assert m.topology.bonds == before.topology.bonds


def test_explicit_and_automatic_same_charges_owned_persistence(
    expanded_source, tmp_path
):
    m = explicit([(0, 1), (1, 2), (2, 0)], [2] * 3)
    auto = type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME)
    explicit_ = assign_pcff_source_types(
        m, expanded_source, auto.assignments, provenance="synthetic software labels"
    )
    q = assign_automatic_pcff_charges(m, explicit_)
    assert q.charges == assign_automatic_pcff_charges(m, auto).charges
    assert explicit_.payload["origin"] == "explicit_checked"
    save_pcff_automatic_record(q, tmp_path / "q.json")
    assert (
        load_pcff_automatic_record(
            tmp_path / "q.json", expanded_source, system=m
        ).identity
        == q.identity
    )
    q.charges.clear()
    q.payload["automatic_typing"]["assignments"].clear()
    assert len(q.charges) == 9
    bad = auto.assignments
    bad[11] = "c2"
    with pytest.raises(PCFFError, match="contradict"):
        assign_pcff_source_types(m, expanded_source, bad, provenance="invalid")
    bad = unpack(q.json_text)
    bad["native_charge_record"]["partial_charges"][11] = 100
    with pytest.raises(PCFFError, match="Contradictory"):
        type(q)(pack(bad), expanded_source).validate_integrity()
    bad = unpack(auto.json_text)
    bad["entries"][11]["environment"]["smallest_ring"] = 6
    with pytest.raises(PCFFError, match="Contradictory"):
        type(auto)(pack(bad), expanded_source).validate_integrity()


def test_order_remapping_coordinates_and_metadata(expanded_source):
    m = explicit([(0, 1), (1, 2), (2, 0)], [2] * 3)
    original = type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME)
    topology = Topology()
    mapping = {i: 900 - 7 * i for i in m.topology.sites}
    for i in reversed(list(m.topology.sites)):
        atom = deepcopy(m.topology.sites[i])
        atom.id = mapping[i]
        topology.add_site(atom)
    for bond in reversed(list(m.topology.bonds.values())):
        topology.add_bond(mapping[bond.site2], mapping[bond.site1], order=1)
    moved = MolecularSystem(topology, Coordinates())
    t = type_pcff_atoms(moved, expanded_source, profile=expanded.PROFILE_NAME)
    assert t.assignments == {mapping[i]: v for i, v in original.assignments.items()}
    for i in m.topology.sites:
        m.coordinates.set(i, [i, 1, 2])
        m.topology.sites[i].metadata["repeat_unit_index"] = 777
    assert (
        original.identity
        == type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME).identity
    )
    m.topology.sites[11].mass += 1
    with pytest.raises(PCFFError, match="graph mismatch"):
        original.validate_integrity(m)


def test_unsupported_environment_never_generic_fallback(expanded_source):
    m = explicit([(0, 1)], [2, 3])  # Missing H/radical-like incomplete graph
    t = type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME)
    assert not t.complete and not t.payload["assignments"]
    with pytest.raises(PCFFError, match="Incomplete"):
        assign_automatic_pcff_charges(m, t)


def test_wilson_compiled_synthetic_fd_and_permutations():
    mm = pytest.importorskip("openmm")
    from test_pcff_evaluation import compiled

    xyz = np.array(
        [[1.1, 0.2, 0.18], [0, 0, 0], [-0.5, 1.2, -0.07], [-0.8, -0.9, 0.11]]
    )
    term = SourceClass2Term("wilson_out_of_plane", (0, 1, 2, 3), (8.7, 0.0))

    def data(sites):
        return {
            "schema": "island_pcff_source_model_v1",
            "terms": [
                {
                    "family": term.family,
                    "sites": sites,
                    "coefficients": [8.7, 0.0],
                    "equilibria": [],
                }
            ],
            "nonbonded": [
                {"site": i, "rmin": 1, "epsilon": 0, "charge": 0} for i in range(4)
            ],
            "special_pair_inventory": [],
        }

    energy, forces = compiled(data([0, 1, 2, 3]), xyz)
    assert energy == pytest.approx(term.energy(dict(enumerate(xyz))), abs=1e-12)
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        fd = finite_difference_forces(term, dict(enumerate(xyz)), displacement=h)
        errors.append(np.max(np.abs(forces - np.array(list(fd.values())))))
    assert errors[2] < errors[1] < errors[0] and errors[-1] < 1e-7
    for p in permutations((0, 2, 3)):
        e, f = compiled(data([p[0], 1, p[1], p[2]]), xyz)
        assert e == pytest.approx(energy, abs=1e-12)
        np.testing.assert_allclose(f, forces, atol=1e-11)
    assert mm.__version__


def test_full_catalog_namespaces_and_positional_equivalence(
    expanded_source, monkeypatch
):
    from island.forcefields.pcff import (
        inspect_pcff_full_source,
        resolve_pcff_source_record,
    )

    raw = expanded_source.raw.replace(
        b"#end",
        b"""#quadratic_bond cff91_auto
1.0 1 c_ h_ 1.13 300
2.0 1 c_ h_ 1.14 310
#quadratic_angle cff91_auto
1.0 1 h_ c_ h_ 109.5 40
#torsion_1 cff91_auto
1.0 1 * c_ c_ * 0.4 3 180
#wilson_out_of_plane cff91_auto
1.0 1 h_ c_ h_ h_ 3 0
#torsion-torsion_1 cff91
#end""",
    )
    # The synthetic source's auto roles are inspected independently below.
    inv = expanded_source.inventory
    for sec in inv["sections"]:
        if sec["name"] == "auto_equivalence":
            for line in sec["lines"]:
                if line["raw"].strip() and line["raw"].strip()[0].isdigit():
                    label = line["raw"].split()[2]
                    replacement = (
                        "1.0 1 "
                        + label
                        + " "
                        + " ".join(["c_" if label.startswith("c") else "h_"] * 9)
                    )
                    raw = raw.replace(line["raw"].encode(), replacement.encode())
    raw = raw.replace(
        b"#bond_increments",
        b"1.0 1 c2 c_ c_ c_ c_ c_ c_ c_ c_ c_\n1.0 1 hc h_ h_ h_ h_ h_ h_ h_ h_ h_\n#bond_increments",
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    s = PCFFSource(raw, digest(raw))
    catalog = inspect_pcff_full_source(s)
    assert sum(r["section"] == "quadratic_bond" for r in catalog["records"]) == 2
    selected = resolve_pcff_source_record(
        s, "quadratic_bond", ["c2", "hc"], namespace="cff91_auto"
    )
    assert selected["status"] == "assigned"
    assert selected["normalized_values"] == [1.14, 310 * 4.184]
    assert selected["position_roles"] == ["bond", "bond"]
    angle = resolve_pcff_source_record(
        s, "quadratic_angle", ["hc", "c2", "hc"], namespace="cff91_auto"
    )
    assert angle["status"] == "assigned"
    assert angle["position_roles"] == ["angle_end", "angle_apex", "angle_end"]
    assert angle["normalized_values"][0] == pytest.approx(np.deg2rad(109.5))
    assert angle["normalized_values"][1] == 40 * 4.184
    torsion = resolve_pcff_source_record(
        s, "torsion_1", ["hc", "c2", "c2", "hc"], namespace="cff91_auto"
    )
    assert torsion["status"] == "assigned" and torsion["normalized_values"][1] == 3
    bad = raw.replace(b"#end", b"#quadratic_bond cff91_auto\n2.0 1 c_ h_ 1.5 310\n#end")
    monkeypatch.setitem(source.PIN, "sha256", digest(bad))
    assert (
        resolve_pcff_source_record(
            PCFFSource(bad, digest(bad)),
            "quadratic_bond",
            ["c2", "hc"],
            namespace="cff91_auto",
        )["status"]
        == "ambiguous"
    )


def test_wilson_center_fixed_and_repeated_labels():
    from island.forcefields.pcff.class2 import resolve

    row = {
        "id": "synthetic",
        "section": "wilson_out_of_plane",
        "types": ["h", "cp", "cp", "cp"],
        "version": "1.0",
        "normalized_values": [7.0, 0.0],
    }
    for arms in permutations(["h", "cp", "cp"]):
        query = [arms[0], "cp", arms[1], arms[2]]
        chosen = resolve("wilson_out_of_plane", query, {"synthetic": row}, [])
        assert chosen["status"] == "assigned"
    assert (
        resolve("wilson_out_of_plane", ["cp", "h", "cp", "cp"], {"synthetic": row}, [])[
            "status"
        ]
        == "missing"
    )


def test_offline_optional_dependency_isolation(expanded_source, tmp_path):
    import os
    import subprocess
    import sys

    m = explicit([(0, 1), (1, 2), (2, 0)], [2] * 3)
    t = type_pcff_atoms(m, expanded_source, profile=expanded.PROFILE_NAME)
    (tmp_path / "source.frc").write_bytes(expanded_source.raw)
    save_pcff_automatic_record(t, tmp_path / "typing.json")
    code = """
import sys
class Block:
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in ('rdkit','openmm','scipy','parmed','foyer'): raise AssertionError(fullname)
sys.meta_path.insert(0,Block())
from pathlib import Path
from island.forcefields.pcff import load_pcff_automatic_record, expanded, source
raw=Path(sys.argv[1]).read_bytes(); sha=source.digest(raw)
source.PIN['sha256']=sha; expanded.PROFILE['source_sha256']=sha
r=load_pcff_automatic_record(sys.argv[2],source.PCFFSource(raw,sha))
assert r.complete
"""
    subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(tmp_path / "source.frc"),
            str(tmp_path / "typing.json"),
        ],
        check=True,
        env=os.environ.copy(),
    )


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("c1ccccc1", ["cp"] * 6),
        ("c1ccoc1", ["c5", "c5", "c5", "op", "c5"]),
        ("c1ccsc1", ["c5", "c5", "cs", "sp", "cs"]),
        ("C=C", ["c=", "c="]),
        ("CC=C", ["c3", "c=1", "c="]),
        ("CC=CC", ["c3", "c=2", "c=2", "c3"]),
        ("CC#N", ["c3", "ct", "nt"]),
        ("CC(=O)[O-]", ["c3", "c-", "o-", "o-"]),
        ("C[N+](C)(C)C", ["c3", "n4", "c3", "c3", "c3"]),
        ("CCOC(C)OCC", ["c3", "c2", "oc", "coh", "c3", "oc", "c2", "c3"]),
        ("CS", ["c3", "sh"]),
        ("CSC", ["c3", "sc", "c3"]),
        ("CCl", ["c3", "cl"]),
        ("CC(=O)C", ["c3", "c_0", "o_1", "c3"]),
    ],
)
def test_independent_graph_predicates(smiles, expected):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles
    from island.forcefields.pcff.automatic import chemical_graph

    m = from_smiles(smiles, random_seed=2026)
    labels, _, issues = expanded.recognize(chemical_graph(m))
    assert not issues
    # This is a fixed SMILES fixture, not an atom-order production typing rule.
    heavy = [i for i, a in m.topology.sites.items() if a.element != "H"]
    assert [labels[i] for i in heavy] == expected


def test_ambiguous_groups_and_malformed_aromaticity_are_not_absorbed():
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles
    from island.forcefields.pcff.automatic import chemical_graph

    for smiles in ("CC(=O)N", "Nc1ccccc1", "OO", "[CH3]", "[2H]C"):
        m = from_smiles(smiles, random_seed=2026)
        labels, _, issues = expanded.recognize(chemical_graph(m))
        assert issues and not labels
    g = chemical_graph(from_smiles("c1ccccc1", random_seed=2026))
    next(b for b in g["bonds"] if b["aromatic"])["aromatic"] = False
    labels, _, issues = expanded.recognize(g)
    assert not labels and any(
        i["reason"] == "inconsistent_aromatic_bond_semantics" for i in issues
    )


def test_expanded_model_integrity_and_facade_bundle(parameters, tmp_path, monkeypatch):  # noqa: F811
    from island.forcefields import (
        PreparedForceFieldSources,
        adopt_forcefield,
        load_prepared_forcefield,
        save_prepared_forcefield,
    )
    from island.forcefields.pcff import (
        assign_pcff_parameters,
        define_pcff_model,
        special_pair_policy,
    )
    from island.forcefields.pcff.automatic import graph_system
    from island.workflows.bundle import system_data

    old_bytes = parameters.json_text
    monkeypatch.setitem(
        expanded.PROFILE, "source_sha256", parameters.source.expected_sha256
    )
    m = graph_system(unpack(old_bytes)["charge_record"]["automatic_typing"]["graph"])
    for n, i in enumerate(m.topology.sites):
        m.coordinates.set(i, [n % 3, n // 3, 0.1 * (n % 2)])
    t = type_pcff_atoms(m, parameters.source, profile=expanded.PROFILE_NAME)
    q = assign_automatic_pcff_charges(m, t)
    a = assign_pcff_parameters(m, t, q)
    spec = define_pcff_model(
        a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert spec.payload["model_definition_complete"]
    assert spec.payload["schema"] == "island_pcff_source_model_v1"
    prepared = adopt_forcefield(m, "pcff", spec)
    local = tmp_path / "source.frc"
    local.write_bytes(parameters.source.raw)
    sources = PreparedForceFieldSources(pcff_frc=local)
    path = save_prepared_forcefield(m, prepared, tmp_path / "bundle", sources=sources)
    loaded = load_prepared_forcefield(path, sources=sources)
    assert loaded.prepared.identity == prepared.identity
    assert system_data(loaded.system) == system_data(m)
    bad = unpack(a.json_text)
    bad["assignments"][0]["normalized_values"][0] += 1
    with pytest.raises(PCFFError, match="Contradictory"):
        type(a)(pack(bad), parameters.source).validate_integrity(m)
    bad = unpack(spec.json_text)
    bad["terms"][0]["coefficients"][0] += 1
    with pytest.raises(PCFFError, match="Contradictory"):
        type(spec)(pack(bad), a).validate_integrity(m)
    assert parameters.json_text == old_bytes
    parameters.validate_integrity()


def test_independent_converter_aa_correction_keeps_center(monkeypatch):
    """Synthetic reproduction of the retained n4 converter-role discrepancy."""
    from pathlib import Path

    monkeypatch.syspath_prepend(str(Path("scripts").resolve()))
    from diagnose_pcff_angle_angle import source_aa_coefficients

    raw = b"""!synthetic source rows, not redistributed scientific parameter data
#equivalence cff91
1 1 c3 c c c c c
1 1 n4 n+ n+ n+ n+ n+
#angle-angle cff91
1 1 c n+ c c -1.5
1 1 c c n+ c -4.2
#end
"""
    sections = {
        "Angles": [
            ["1", "1", "1", "2", "3"],
            ["2", "1", "1", "2", "4"],
            ["3", "1", "3", "2", "4"],
        ],
        "Angle Coeffs": [["1", "109"]],
        "AngleAngle Coeffs": [["1", "-1.5", "-4.2", "-4.2", "109", "109", "109"]],
        "Impropers": [["1", "1", "1", "2", "3", "4"]],
    }
    commands, evidence = source_aa_coefficients(
        raw, sections, {1: "c3", 2: "n4", 3: "c3", 4: "c3"}
    )
    assert commands == ["improper_coeff 1 aa -1.5 -1.5 -1.5 109 109 109"]
    assert all(e["source_lines"] == [6] for e in evidence)
