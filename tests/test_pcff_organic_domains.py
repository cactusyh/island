"""Synthetic source software tests. Real-source receipts are separate."""

import builtins
import json
from collections import Counter
from copy import deepcopy
from pathlib import Path

import pytest

pytest.importorskip("rdkit")
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_domains import domains_source  # noqa: F401

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_pcff_source_types,
    expanded,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.automatic import PCFFAutomaticTypingResult, chemical_graph
from island.forcefields.pcff.organic_domains import PROFILE_NAME
from island.forcefields.pcff.source import PCFFSource, digest

CASES = json.loads(
    (Path(__file__).parents[1] / "docs/evidence/phase_4j4_declaration.json").read_text()
)["cases"]


@pytest.fixture
def organic_source(domains_source, monkeypatch):  # noqa: F811
    types = {
        "f": ("F", 18.9984, 1),
        "cl": ("Cl", 35.453, 1),
        "br": ("Br", 79.909, 1),
        "i": ("I", 126.9044, 1),
        "s'": ("S", 32.064, 1),
        "nh+": ("N", 14.007, 3),
        "c+": ("C", 12.011, 3),
        "nr": ("N", 14.007, 3),
        "n3n": ("N", 14.007, 3),
        "n4n": ("N", 14.007, 3),
        "np": ("N", 14.007, 2),
        "n3m": ("N", 14.007, 3),
        "n4m": ("N", 14.007, 3),
        "c3h": ("C", 12.011, 4),
        "c4h": ("C", 12.011, 4),
    }
    existing = {
        r["data"]["type"]
        for sec in domains_source.inventory["sections"]
        if sec["name"] == "atom_types"
        for r in sec["records"]
    }
    lines = "".join(
        f"1 1 {t} {m} {el} {n} synthetic J4 software fixture\n"
        for t, (el, m, n) in types.items()
        if t not in existing
    )
    raw = domains_source.raw.replace(
        b"#equivalence cff91", lines.encode() + b"#equivalence cff91"
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_declared_labels_and_explicit_contract(organic_source, case):
    s = from_smiles(case["smiles"], random_seed=2026)
    before = deepcopy(s)
    t = type_pcff_atoms(s, organic_source, profile=PROFILE_NAME)
    assert t.complete, t.payload["diagnostics"]
    assert dict(Counter(t.assignments.values())) == case["expected_counts"]
    explicit = assign_pcff_source_types(
        s,
        organic_source,
        t.assignments,
        profile=PROFILE_NAME,
        provenance="synthetic controlled labels",
    )
    assert explicit.assignments == t.assignments
    assert chemical_graph(s) == chemical_graph(before)
    wrong = t.payload
    sid = next(iter(wrong["assignments"]))
    wrong["assignments"][sid] = "hw"
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFAutomaticTypingResult(pack(wrong), organic_source).validate_integrity()
    with pytest.raises(PCFFError):
        assign_pcff_source_types(
            s,
            organic_source,
            wrong["assignments"],
            profile=PROFILE_NAME,
            provenance="synthetic incorrect labels",
        )


@pytest.mark.parametrize(
    "smiles,forbidden",
    [
        ("O=O", {"o_1"}),
        ("N#N", {"nt"}),
        ("[2H]Cl", {"h"}),
        ("CS(=O)C", {"s'"}),
        ("O=S=O", {"s'"}),
        ("c1ccncc1", {"nh+"}),
        ("N1CC1", {"n3n"}),
        ("N1CCC1", {"n4n"}),
        ("NC(N)=N", {"nr", "c+"}),
        ("[S-]C", {"s'"}),
        ("[2H]C=S", {"s'"}),
        ("c1cc[n+](C)cc1", {"nh+"}),
    ],
)
def test_near_misses_do_not_absorb_new_rules(organic_source, smiles, forbidden):
    s = from_smiles(smiles, random_seed=2026)
    t = type_pcff_atoms(s, organic_source, profile=PROFILE_NAME)
    assert not forbidden.intersection(t.payload["assignments"].values())


@pytest.mark.parametrize("field,value", [("formal_charge", 1), ("mass", 3.0)])
def test_changed_identity_blocks_hydrogen_halide(organic_source, field, value):
    s = from_smiles("[H]Cl", random_seed=2026)
    h = next(a for a in s.topology.sites.values() if a.element == "H")
    setattr(h, field, value)
    assert not type_pcff_atoms(s, organic_source, profile=PROFILE_NAME).complete


def test_offline_validation_and_frozen_v3(organic_source, monkeypatch):
    s = from_smiles("[H]Cl", random_seed=2026)
    t = type_pcff_atoms(s, organic_source, profile=PROFILE_NAME)
    assert not type_pcff_atoms(
        s, organic_source, profile="island_pcff_source_graph_v3"
    ).complete
    importer = builtins.__import__

    def block(name, *a, **kw):
        if name.split(".")[0] in {"rdkit", "openmm", "scipy", "foyer", "parmed"}:
            raise AssertionError(name)
        return importer(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", block)
    t.validate_integrity(s)
    assert PCFFAutomaticTypingResult(t.json_text, organic_source).identity == t.identity


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
@pytest.mark.parametrize("mutation", ["isotope", "radical", "formal_charge"])
def test_motif_identity_near_miss_cannot_use_explicit_labels(
    organic_source, case, mutation
):
    system = from_smiles(case["smiles"], random_seed=2026)
    valid = type_pcff_atoms(system, organic_source, profile=PROFILE_NAME)
    added = {"h", "s'", "nh+", "c+", "nr", "n3n", "n4n"}
    site = next(
        system.topology.sites[i] for i, t in valid.assignments.items() if t in added
    )
    if mutation == "formal_charge":
        site.formal_charge += 1
    else:
        site.metadata["isotope" if mutation == "isotope" else "radical_electrons"] = 1
    with pytest.raises(PCFFError):
        valid.validate_integrity(system)
    changed = type_pcff_atoms(system, organic_source, profile=PROFILE_NAME)
    assert not changed.complete
    with pytest.raises(PCFFError, match="Explicit chemistry constraints unresolved"):
        assign_pcff_source_types(
            system,
            organic_source,
            valid.assignments,
            profile=PROFILE_NAME,
            provenance="old labels cannot bypass changed chemical identity",
        )
