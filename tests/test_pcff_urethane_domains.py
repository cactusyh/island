"""Synthetic-source typing contracts; authentic source evidence is separate."""

import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import pytest
from test_pcff_amines import amine_source  # noqa: F401
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_domains import domains_source  # noqa: F401
from test_pcff_organic_domains import organic_source  # noqa: F401

from island import Coordinates, MolecularSystem, Topology
from island.builders import build_linear_polymer
from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff import assign_pcff_source_types, type_pcff_atoms
from island.forcefields.pcff.automatic import PCFFAutomaticTypingResult
from island.forcefields.pcff.urethane_domains import PROFILE_NAME

pytest.importorskip("rdkit")
CASES = json.loads(
    (
        Path(__file__).parents[1] / "docs/evidence/phase_4j15_declaration.json"
    ).read_text()
)["cases"]


def construct(c):
    return (
        build_linear_polymer(c["psmiles"], dp=c["dp"], random_seed=2026)
        if "psmiles" in c
        else from_smiles(c["smiles"], random_seed=2026)
    )


@pytest.mark.parametrize("c", CASES, ids=[c["name"] for c in CASES])
def test_declared_predicates_and_near_misses(amine_source, c):  # noqa: F811
    from collections import Counter

    s = construct(c)
    before = deepcopy(s.to_dict())
    old = type_pcff_atoms(s, amine_source, profile="island_pcff_source_graph_v5")
    t = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    if c["role"] == "positive":
        assert t.complete, t.payload["diagnostics"]
        counts = Counter(t.assignments.values())
        assert all(counts[k] == n for k, n in c["expected_motif_counts"].items())
        assert (
            assign_pcff_source_types(
                s,
                amine_source,
                t.assignments,
                profile=PROFILE_NAME,
                provenance="synthetic explicit check",
            ).assignments
            == t.assignments
        )
        hs = [
            i
            for i, e in t.payload["entries"].items()
            if e["type"] in ("hn", "hn2") and "nitrogen_family_rule" in e["environment"]
        ]
        if hs:
            bad = t.assignments
            bad[hs[0]] = "hn" if bad[hs[0]] == "hn2" else "hn2"
            with pytest.raises(PCFFError, match="contradict"):
                assign_pcff_source_types(
                    s,
                    amine_source,
                    bad,
                    profile=PROFILE_NAME,
                    provenance="wrong family",
                )
    else:
        assert old.payload["assignments"] == t.payload["assignments"]
        assert old.payload["diagnostics"] == t.payload["diagnostics"]
        assert not any(
            "nitrogen_family_rule" in e["environment"]
            for e in t.payload["entries"].values()
        )
    assert before == s.to_dict()


def test_stable_id_coordinate_independence_and_tampering(amine_source):  # noqa: F811
    c = next(c for c in CASES if c["name"] == "mixed_functional")
    s = construct(c)
    t = type_pcff_atoms(s, amine_source, profile=PROFILE_NAME)
    remap = {i: 1001 - 9 * j for j, i in enumerate(s.topology.sites)}
    top = Topology()
    for a in reversed(list(s.topology.sites.values())):
        top.add_site(replace(deepcopy(a), id=remap[a.id]))
    for b in reversed(list(s.topology.bonds.values())):
        top.add_bond(remap[b.site2], remap[b.site1], order=b.order, aromatic=b.aromatic)
    changed = MolecularSystem(top, Coordinates())
    u = type_pcff_atoms(changed, amine_source, profile=PROFILE_NAME)
    assert u.assignments == {remap[i]: v for i, v in t.assignments.items()}
    s.coordinates.translate([2, 4, 6])
    s.metadata["unrelated"] = "not typing input"
    assert type_pcff_atoms(s, amine_source, profile=PROFILE_NAME).identity == t.identity
    p = t.payload
    i = next(
        i for i, e in p["entries"].items() if "nitrogen_family_rule" in e["environment"]
    )
    p["entries"][i]["environment"]["nitrogen_family_rule"]["id"] = "invented"
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFAutomaticTypingResult(pack(p), amine_source).validate_integrity()
    n = next(i for i, a in s.topology.sites.items() if a.element == "N")
    s.topology.sites[n].formal_charge = 1
    with pytest.raises(PCFFError, match="graph mismatch"):
        t.validate_integrity(s)
