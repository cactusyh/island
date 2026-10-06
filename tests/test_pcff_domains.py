"""J3 synthetic software contracts; pinned-source acceptance is separate."""

import pytest


def test_rdkit_preserves_declared_isotope_and_radical_identity():
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    water = from_smiles("[2H]O[2H]", random_seed=2026)
    assert [a.metadata.get("isotope", 0) for a in water.topology.sites.values()] == [
        2,
        0,
        2,
    ]
    radical = from_smiles("[CH3]", random_seed=2026)
    assert (
        next(a for a in radical.topology.sites.values() if a.element == "C").metadata[
            "radical_electrons"
        ]
        == 1
    )
    normal = from_smiles("CC", random_seed=2026)
    assert all(
        "isotope" not in a.metadata and "radical_electrons" not in a.metadata
        for a in normal.topology.sites.values()
    )


from copy import deepcopy

from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_fallbacks import fallback_source, make  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_pcff_source_types,
    expanded,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.domains import PROFILE_NAME
from island.forcefields.pcff.source import PCFFSource, digest


@pytest.fixture
def domains_source(synthetic, monkeypatch):  # noqa: F811
    # Invented numerical fixture, not a redistributed PCFF library.
    types = {
        "h": ("H", 1.008, 1),
        "s": ("S", 32.064, 2),
        "hs": ("H", 1.008, 1),
        "dw": ("D", 2.014, 1),
        "o*": ("O", 15.999, 2),
        "hw": ("H", 1.008, 1),
        "n": ("N", 14.007, 3),
        "c_1": ("C", 12.011, 3),
        "c_0": ("C", 12.011, 3),
        "o_1": ("O", 15.999, 1),
        "hn": ("H", 1.008, 1),
        "hn2": ("H", 1.008, 1),
        "npc": ("N", 14.007, 3),
        "c5": ("C", 12.011, 3),
        "cp": ("C", 12.011, 3),
        "nn": ("N", 14.007, 3),
        "n=": ("N", 14.007, 2),
        "n=1": ("N", 14.007, 2),
        "n=2": ("N", 14.007, 2),
        "c=": ("C", 12.011, 3),
        "c=1": ("C", 12.011, 3),
        "c=2": ("C", 12.011, 3),
        "p": ("P", 30.974, 4),
        "hp": ("H", 1.008, 1),
        "o=": ("O", 15.999, 1),
        "o": ("O", 15.999, 2),
        "ho2": ("H", 1.008, 1),
        "na": ("N", 14.007, 3),
    }
    existing = {
        r["data"]["type"]
        for sec in synthetic.inventory["sections"]
        if sec["name"] == "atom_types"
        for r in sec["records"]
    }
    lines = "".join(
        f"1 1 {t} {m} {el} {d} synthetic software fixture\n"
        for t, (el, m, d) in types.items()
        if t not in existing
    )
    raw = synthetic.raw.replace(
        b"#equivalence cff91", lines.encode() + b"#equivalence cff91"
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    return PCFFSource(raw, digest(raw))


@pytest.mark.parametrize(
    "smiles,required",
    [
        ("[H][H]", {"h"}),
        ("S", {"s", "hs"}),
        ("[2H]O[2H]", {"dw", "o*"}),
        ("NC=O", {"n", "c_1", "o_1"}),
        ("Cn1cccc1", {"npc", "c5"}),
        ("Nc1ccccc1", {"nn", "hn2"}),
        ("C=N", {"c=", "n=", "hn"}),
        ("C=NC", {"n=1"}),
        ("CC=NC", {"n=2"}),
        ("[PH3]=O", {"p", "hp", "o="}),
        ("OO", {"o", "ho2"}),
        ("CC(=O)O", {"oh", "ho2"}),
    ],
)
def test_specific_domains_and_explicit_bridge(domains_source, smiles, required):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    s = from_smiles(smiles, random_seed=2026)
    before = deepcopy(s)
    t = type_pcff_atoms(s, domains_source, profile=PROFILE_NAME)
    assert t.complete, t.payload["diagnostics"]
    assert required <= set(t.assignments.values())
    explicit = assign_pcff_source_types(
        s,
        domains_source,
        t.assignments,
        provenance="synthetic labels",
        profile=PROFILE_NAME,
    )
    assert explicit.assignments == t.assignments
    old = type_pcff_atoms(s, domains_source, profile="island_pcff_source_graph_v2")
    assert not old.complete
    s.coordinates.translate([3, 4, 5])
    s.metadata["unrelated"] = "not chemical"
    assert (
        type_pcff_atoms(s, domains_source, profile=PROFILE_NAME).identity == t.identity
    )
    assert (
        s.topology.sites == before.topology.sites
        and s.topology.bonds == before.topology.bonds
    )
    corrupted = unpack(t.json_text)
    corrupted["assignments"][next(iter(t.assignments))] = "fake"
    with pytest.raises(PCFFError):
        type(t)(pack(corrupted), domains_source).validate_integrity()
    changed = deepcopy(s)
    sid = next(iter(changed.topology.sites))
    changed.topology.sites[sid].mass += 0.1
    with pytest.raises(PCFFError):
        t.validate_integrity(changed)


@pytest.mark.parametrize(
    "smiles",
    [
        "O=O",
        "N#N",
        "[H]",
        "[2H][H]",
        "[2H]O",
        "[3H]O[3H]",
        "[PH4+]",
        "P",
        "O=S=O",
        "[Na+]",
        "[Al+3]",
        "[O]O",
    ],
)
def test_similar_unsupported_environments(domains_source, smiles):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    s = from_smiles(smiles, random_seed=2026)
    t = type_pcff_atoms(s, domains_source, profile=PROFILE_NAME)
    assert not t.complete
    with pytest.raises(PCFFError):
        assign_pcff_source_types(
            s,
            domains_source,
            dict.fromkeys(s.topology.sites, "h"),
            provenance="cannot bypass",
            profile=PROFILE_NAME,
        )


def test_automatic_nonbonded_is_explicit_and_family_specific():
    from island.forcefields.pcff.fallbacks import DOMAIN_POLICY, POLICY, resolve

    catalog = {
        "nb": {
            "id": "nb",
            "section": "nonbond(9-6)",
            "namespace": "cff91",
            "types": ["x"],
            "version": "1",
            "normalized_values": [3.0, 0.2],
        }
    }
    eq = {
        "id": "auto",
        "version": "1",
        "data": {"type": "special", "families": {"nonbond": "x", "bond": "wrong"}},
        "section": "auto_equivalence",
        "namespace": "cff91_auto",
    }
    inv = {
        "sections": [
            {"name": "auto_equivalence", "namespace": "cff91_auto", "records": [eq]}
        ]
    }
    assert (
        resolve("nonbond(9-6)", ["special"], catalog, inv, policy=POLICY)["status"]
        == "missing"
    )
    out = resolve("nonbond(9-6)", ["special"], catalog, inv, policy=DOMAIN_POLICY)
    assert out["normalized_values"] == [3.0, 0.2]
    assert out["path"] == "auto_equivalence.nonbond"
    assert out["equivalence_evidence"][0]["record"]["id"] == "auto"
    assert out["selection_policy"] == DOMAIN_POLICY


@pytest.mark.parametrize("smiles", ["S", "[2H]O[2H]", "NC=O", "Cn1cccc1", "C=NC"])
def test_domain_remapping_and_persistence(domains_source, smiles, tmp_path):
    pytest.importorskip("rdkit")
    from island import Coordinates, MolecularSystem, Topology
    from island.chemistry import from_smiles
    from island.forcefields.pcff import (
        load_pcff_automatic_record,
        save_pcff_automatic_record,
    )

    original = from_smiles(smiles, random_seed=2026)
    typed = type_pcff_atoms(original, domains_source, profile=PROFILE_NAME)
    mapping = {i: 700 - 13 * i for i in original.topology.sites}
    topology = Topology()
    for i, atom in reversed(list(original.topology.sites.items())):
        owned = deepcopy(atom)
        owned.id = mapping[i]
        topology.add_site(owned)
    for bond in reversed(list(original.topology.bonds.values())):
        topology.add_bond(
            mapping[bond.site2],
            mapping[bond.site1],
            order=bond.order,
            aromatic=bond.aromatic,
        )
    moved = MolecularSystem(topology, Coordinates())
    result = type_pcff_atoms(moved, domains_source, profile=PROFILE_NAME)
    assert result.assignments == {mapping[i]: v for i, v in typed.assignments.items()}
    save_pcff_automatic_record(result, tmp_path / "typing.json")
    assert (
        load_pcff_automatic_record(
            tmp_path / "typing.json", domains_source, system=moved
        ).identity
        == result.identity
    )
    site = next(iter(moved.topology.sites.values()))
    site.formal_charge = 1
    with pytest.raises(PCFFError):
        result.validate_integrity(moved)
    assert not type_pcff_atoms(moved, domains_source, profile=PROFILE_NAME).complete


def test_bounded_gate_is_not_empty_all_or_duplicate_acceptance():
    import sys
    from pathlib import Path
    from unittest.mock import patch

    with patch.object(
        sys, "path", [str(Path(__file__).parents[1] / "scripts"), *sys.path]
    ):
        from pcff_domain_gates import REQUIRED, numerical_gate
    good = [
        {
            "case": name,
            "passed": True,
            "bundle_hashes": {"manifest.json": "hash"},
            "prepared_identity": "p",
            "model_identity": "m",
            "numerical": [
                {
                    "energy_error": 0.0,
                    "force_max_error": 1e-8,
                    "parameter_rows_and_inventory_match": True,
                }
                for _ in range(2)
            ],
        }
        for name in REQUIRED
    ]
    assert numerical_gate(good, True)
    assert not numerical_gate([], True)
    assert not numerical_gate(good, False)
    assert not numerical_gate([good[0], good[0]], True)
    for key, value in [
        ("error", "after numerical failure"),
        ("bundle_hashes", {}),
        ("passed", False),
        ("numerical", []),
    ]:
        bad = deepcopy(good)
        bad[0][key] = value
        assert not numerical_gate(bad, True)
    bad = deepcopy(good)
    bad[0]["numerical"][0]["force_max_error"] = float("nan")
    assert not numerical_gate(bad, True)


def test_v3_native_search_provenance_and_bundle_roundtrip(fallback_source, tmp_path):  # noqa: F811
    import builtins
    from unittest.mock import patch

    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.forcefields import (
        PreparedForceFieldSources,
        adopt_forcefield,
        load_prepared_forcefield,
        save_prepared_forcefield,
    )
    from island.forcefields.pcff import (
        assign_automatic_pcff_charges,
        assign_pcff_parameters,
        define_pcff_model,
        special_pair_policy,
    )
    from island.forcefields.pcff.fallbacks import DOMAIN_POLICY

    m, _, _, _, _ = make(fallback_source)
    t = type_pcff_atoms(m, fallback_source, profile=PROFILE_NAME)
    q = assign_automatic_pcff_charges(m, t, resolution_policy=DOMAIN_POLICY)
    native = q.payload["native_charge_record"]
    assert native["base_charges"] == dict.fromkeys(m.topology.sites, 0.0)
    assert all(c["searches"][0]["path"] == "direct" for c in native["contributions"])
    bad = unpack(q.json_text)
    bad["native_charge_record"]["contributions"][0]["searches"][0][
        "matched"
    ] = not native["contributions"][0]["searches"][0]["matched"]
    with pytest.raises(PCFFError):
        type(q)(pack(bad), fallback_source).validate_integrity()
    a = assign_pcff_parameters(m, t, q, resolution_policy=DOMAIN_POLICY)
    model = define_pcff_model(
        a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert model.payload["model_definition_complete"]
    path = tmp_path / "synthetic.frc"
    path.write_bytes(fallback_source.raw)
    sources = PreparedForceFieldSources(pcff_frc=path)
    prepared = adopt_forcefield(m, "pcff", model)
    save_prepared_forcefield(m, prepared, tmp_path / "bundle", sources=sources)
    importer = builtins.__import__

    def blocked(name, *args, **kw):
        if name.split(".")[0] in ("openmm", "rdkit", "scipy", "foyer", "parmed"):
            raise AssertionError("optional import")
        return importer(name, *args, **kw)

    with patch("builtins.__import__", blocked):
        loaded = load_prepared_forcefield(tmp_path / "bundle", sources=sources)
        assert loaded.prepared.identity == prepared.identity
        assert pcff_evaluation_identity(
            loaded.prepared.native_result
        ) == pcff_evaluation_identity(model)
    assert loaded.prepared.native_result.identity == model.identity
    bad = unpack(model.json_text)
    bad["compatibility_profile"]["name"] = "wrong"
    with pytest.raises(PCFFError):
        type(model)(pack(bad), a).validate_integrity()


@pytest.mark.parametrize(
    "smiles,element",
    [("S", "S"), ("[2H]O[2H]", "H"), ("[PH3]=O", "O"), ("CC(=O)O", "O")],
)
def test_domain_rules_do_not_erase_invalid_aromatic_identity(
    domains_source, smiles, element
):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    m = from_smiles(smiles, random_seed=2026)
    # Mutate only a chemical flag, not coordinates or input shape.
    atom = next(
        a for a in reversed(list(m.topology.sites.values())) if a.element == element
    )
    atom.metadata["aromatic"] = True
    assert not type_pcff_atoms(m, domains_source, profile=PROFILE_NAME).complete
