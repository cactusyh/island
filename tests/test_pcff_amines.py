"""Synthetic source software tests; J14 real-source acceptance is separate."""

import builtins
from copy import deepcopy

import pytest
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_domains import domains_source  # noqa: F401
from test_pcff_organic_domains import organic_source  # noqa: F401

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields import PCFFOptions
from island.forcefields.pcff import (
    assign_pcff_source_types,
    expanded,
    load_pcff_automatic_record,
    save_pcff_automatic_record,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.amine_domains import PROFILE_NAME
from island.forcefields.pcff.automatic import PCFFAutomaticTypingResult, chemical_graph
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY
from island.forcefields.pcff.source import PCFFSource, digest

pytest.importorskip("rdkit")


@pytest.fixture
def amine_source(organic_source, monkeypatch):  # noqa: F811
    raw = organic_source.raw.replace(
        b"#equivalence cff91",
        b"1 1 h+ 1.008 H 1 synthetic charged hydrogen\n"
        b"1 1 n4 14.007 N 4 synthetic ammonium\n"
        b"1 1 n+ 14.007 N 4 synthetic ammonium alias\n"
        b"1 1 n_2 14.007 N 3 synthetic urethane\n"
        b"1 1 c_2 12.011 C 3 synthetic carbonyl\n"
        b"1 1 o_2 15.999 O 2 synthetic ester\n#equivalence cff91",
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


@pytest.mark.parametrize(
    "smiles,nt,ht",
    [
        ("CN", "na", "hn"),
        ("CNC", "na", "hn"),
        ("CN(C)C", "na", None),
        ("C[NH3+]", "n4", "h+"),
        ("C[NH2+]C", "n4", "h+"),
        ("C[NH+](C)C", "n4", "h+"),
        ("C[N+](C)(C)C", "n4", None),
    ],
)
def test_local_amine_family_and_explicit_contract(amine_source, smiles, nt, ht):
    s = from_smiles(smiles, random_seed=2026)
    before = deepcopy(s.to_dict())
    t = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    assert t.complete
    n = next(i for i, a in s.topology.sites.items() if a.element == "N")
    hs = [
        j
        for b in s.topology.bonds.values()
        if n in (b.site1, b.site2)
        for j in (b.site1, b.site2)
        if s.topology.sites[j].element == "H"
    ]
    assert t.assignments[n] == nt
    assert all(t.assignments[j] == ht for j in hs)
    assert (
        t.payload["entries"][n]["environment"]["amine_rule"]["formal_charge"]
        == s.topology.sites[n].formal_charge
    )
    assert (
        assign_pcff_source_types(
            s,
            amine_source,
            t.assignments,
            profile=PROFILE_NAME,
            provenance="synthetic expectation",
        ).assignments
        == t.assignments
    )
    if hs:
        bad = t.assignments
        bad[hs[0]] = "hn2" if ht == "hn" else "hn"
        with pytest.raises(PCFFError, match="contradict"):
            assign_pcff_source_types(
                s, amine_source, bad, profile=PROFILE_NAME, provenance="wrong family"
            )
    assert s.to_dict() == before


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)N",
        "COC(=O)N",
        "Nc1ccccc1",
        "c1cc[nH+]cc1",
        "NC(=[NH2+])N",
        "C[NH]",
        "[2H]NC",
        "C1CN1",
        "C[NH2+]C",
        "NNC",
    ],
)
def test_near_misses_keep_historical_decisions(amine_source, smiles):
    s = from_smiles(smiles, random_seed=2026)
    # The charged control is deliberately made incompatible before typing.
    if smiles == "C[NH2+]C":
        next(a for a in s.topology.sites.values() if a.element == "N").mass += 2
    old = type_pcff_atoms(s, amine_source, profile="island_pcff_source_graph_v4")
    new = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    assert old.payload["assignments"] == new.payload["assignments"]
    assert old.payload["diagnostics"] == new.payload["diagnostics"]
    assert not any(
        "amine_rule" in e["environment"] for e in new.payload["entries"].values()
    )


def test_owned_offline_persistence_tampering_and_binding(
    amine_source, tmp_path, monkeypatch
):
    s = from_smiles("CN", random_seed=2026)
    old = type_pcff_atoms(s, amine_source, profile="island_pcff_source_graph_v4")
    t = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    old_text = old.json_text
    graph = chemical_graph(s)
    path = tmp_path / "typing.json"
    save_pcff_automatic_record(t, path)
    real_import = builtins.__import__

    def guarded(name, *a, **kw):
        if name.split(".")[0] in {"rdkit", "openmm", "scipy", "foyer", "parmed"}:
            raise AssertionError("Scientific import during native reconstruction")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", guarded)
    loaded = load_pcff_automatic_record(path, amine_source)
    assert loaded.identity == t.identity
    assert old.json_text == old_text
    payload = loaded.payload
    n = next(i for i, a in s.topology.sites.items() if a.element == "N")
    payload["entries"][n]["environment"]["amine_rule"]["formal_charge"] = 1
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFAutomaticTypingResult(pack(payload), amine_source).validate_integrity()
    s.coordinates.positions = {i: [0.0, 0.0, 0.0] for i in s.topology.sites}
    loaded.validate_integrity(s)
    assert chemical_graph(s) == graph
    s.topology.sites[n].formal_charge = 1
    with pytest.raises(PCFFError, match="graph mismatch"):
        loaded.validate_integrity(s)


def test_profile_is_distinct_from_resolution_policy(tmp_path):
    for policy in (MSI_POLICY, COMPATIBILITY_POLICY):
        PCFFOptions(
            tmp_path / "external.frc",
            (0, 0, 1),
            (0, 0, 1),
            typing_profile=PROFILE_NAME,
            resolution_policy=policy,
        )
    with pytest.raises(Exception, match="fallback policy"):
        PCFFOptions(
            tmp_path / "external.frc", (0, 0, 1), (0, 0, 1), typing_profile=PROFILE_NAME
        )


def test_oligomer_ends_interrepeat_and_remapped_ids(amine_source):
    from dataclasses import replace

    from island import Coordinates, MolecularSystem, Topology
    from island.builders import build_linear_polymer

    s = build_linear_polymer("[*:1]CCN[*:2]", dp=2, random_seed=2026)
    before = deepcopy(s.to_dict())
    t = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    ns = [
        e["environment"]["amine_rule"]
        for e in t.payload["entries"].values()
        if e["type"] == "na"
    ]
    assert sorted(len(e["hydrogen_neighbors"]) for e in ns) == [1, 2]
    remap = {i: 1009 - 7 * j for j, i in enumerate(s.topology.sites)}
    topology = Topology()
    for atom in reversed(list(s.topology.sites.values())):
        topology.add_site(replace(deepcopy(atom), id=remap[atom.id]))
    for bond in reversed(list(s.topology.bonds.values())):
        topology.add_bond(remap[bond.site2], remap[bond.site1], order=bond.order)
    changed = MolecularSystem(topology, Coordinates())
    u = type_pcff_atoms(changed, amine_source, profile=PROFILE_NAME)
    assert u.assignments == {remap[i]: v for i, v in t.assignments.items()}
    s.coordinates.translate([2.0, 3.0, 4.0])
    s.metadata["unrelated"] = "not a typing input"
    assert type_pcff_atoms(s, amine_source, profile=PROFILE_NAME).identity == t.identity
    assert all(s.to_dict()[k] == before[k] for k in ("sites", "bonds"))
