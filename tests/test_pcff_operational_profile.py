"""Synthetic profile/source contract tests. Real aromatic acceptance is separate."""

import builtins
import math
from copy import deepcopy
from dataclasses import replace
from itertools import product

import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields import (
    ForceFieldRequest,
    ForceFieldRequestError,
    PCFFOptions,
    PreparedBundleError,
    PreparedForceFieldError,
    PreparedForceFieldSources,
    load_prepared_forcefield,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    PCFFOperationalProfile,
    PCFFOperationalSelection,
    expanded,
    load_pcff_operational_profile,
    model,
    save_pcff_operational_profile,
    source,
)
from island.forcefields.pcff.class2 import FAMILIES
from island.forcefields.pcff.operational_profile import DEFINITION, NAME
from island.workflows import storage


def ring():
    t = Topology()
    coords = Coordinates()
    cs = [11, 23, 37, 41, 59, 67]
    hs = [107, 211, 313, 419, 523, 631]
    for n, (c, h) in enumerate(zip(cs, hs, strict=True)):
        angle = 2 * math.pi * n / 6
        t.add_site(
            AtomSite(
                id=c,
                name="C" + str(c),
                element="C",
                atomic_number=6,
                mass=12.011,
                metadata={"aromatic": True},
            )
        )
        t.add_site(
            AtomSite(id=h, name="H" + str(h), element="H", atomic_number=1, mass=1.008)
        )
        coords.set(c, [1.4 * math.cos(angle), 1.4 * math.sin(angle), 0.01 * n])
        coords.set(h, [2.4 * math.cos(angle), 2.4 * math.sin(angle), 0.01 * n])
    for n, c in enumerate(cs):
        t.add_bond(c, cs[(n + 1) % 6], order=1.5, aromatic=True)
        t.add_bond(c, hs[n], order=1)
    return MolecularSystem(t, coords)


@pytest.fixture
def synthetic_profile_source(tmp_path, monkeypatch):
    raw = """!BIOSYM forcefield 1
#version synthetic.frc 1 synthetic
#atom_types cff91
1 1 cp 12.011 C 3 synthetic aromatic carbon
1 1 hc 1.008 H 1 synthetic hydrogen
1 1 h 1.008 H 1 synthetic parameter hydrogen
#equivalence cff91
1 1 cp cp cp cp cp cp
1 1 hc h h h h h
1 1 h h h h h h
#auto_equivalence cff91_auto
1 1 cp cp cp cp_ cp_ cp_ cp_ cp_ cp_ cp_
1 1 hc h h h_ h_ h_ h_ h_ h_ h_
#bond_increments cff91_auto
1 1 cp cp 0 0
1 1 cp h -0.1 0.1
"""
    for family, (arity, _, units) in FAMILIES.items():
        values = [0] * len(units)
        if family == "quartic_bond":
            values = [1.4, 10, 0, 0]
        elif family == "quartic_angle":
            values = [120, 10, 0, 0]
        elif family == "nonbond(9-6)":
            values = [3, 0.1]
        elif family == "wilson_out_of_plane":
            values = [1, 0]
        raw += f"#{family} cff91\n"
        for labels in product(("cp", "h"), repeat=arity):
            raw += "1 1 " + " ".join(labels) + " " + " ".join(map(str, values)) + "\n"
    raw += "#end\n"
    raw = raw.encode()
    sha = source.digest(raw)
    monkeypatch.setitem(source.PIN, "sha256", sha)
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", sha)
    monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
    monkeypatch.setitem(DEFINITION, "source_sha256", sha)
    monkeypatch.setitem(
        DEFINITION,
        "provenance",
        {
            "repository": "synthetic software fixture",
            "commit": "fixture",
            "source_path": "synthetic.frc",
        },
    )
    path = tmp_path / "synthetic.frc"
    path.write_bytes(raw)
    return path, PCFFOperationalSelection(NAME, sha)


def test_profile_roundtrip_ownership_and_promotion_rejection(tmp_path):
    sel = PCFFOperationalSelection(NAME, DEFINITION["source_sha256"])
    p = sel.profile()
    out = tmp_path / "profile.json"
    save_pcff_operational_profile(p, out)
    assert load_pcff_operational_profile(out).identity == p.identity
    p.payload["authorized_operations"].clear()
    assert p.payload["authorized_operations"]["evaluation"]
    for name, sha in [
        (NAME, "0" * 64),
        ("candidate", DEFINITION["source_sha256"]),
        (NAME, "3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693"),
    ]:
        with pytest.raises(PCFFError):
            PCFFOperationalSelection(name, sha)
    for mutate in ["source", "capabilities", "labels", "policy", "wilson"]:
        data = unpack(p.json_text)
        if mutate == "source":
            data["source_sha256"] = (
                "3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693"
            )
        elif mutate == "capabilities":
            data["authorized_operations"]["model_assembly"] = False
        elif mutate == "labels":
            data["authorized_atom_labels"].append("nh+")
        elif mutate == "policy":
            data["authorized_charge_policy"] = "uniform_normalization"
        else:
            data["wilson"] = "nonzero equilibrium implicitly allowed"
        with pytest.raises(PCFFError):
            save_pcff_operational_profile(
                PCFFOperationalProfile(pack(data)), tmp_path / (mutate + ".json")
            )
    with pytest.raises(PCFFError):
        save_pcff_operational_profile(p, out)


def test_exact_graph_domain_negatives_and_nonmutation():
    s = ring()
    before = deepcopy(s.to_dict())
    p = PCFFOperationalSelection(NAME, DEFINITION["source_sha256"]).profile()
    p.validate_system(s)
    assert s.to_dict() == before
    for change in [
        "isotope",
        "radical",
        "charge",
        "mass",
        "stereo",
        "bond",
        "aromatic",
        "site",
    ]:
        m = deepcopy(s)
        if change == "isotope":
            m.topology.sites[107].metadata["isotope"] = 2
        elif change == "radical":
            m.topology.sites[11].metadata["radical_electrons"] = 1
        elif change == "charge":
            m.topology.sites[11].formal_charge = 1
        elif change == "mass":
            m.topology.sites[107].mass = 2.014
        elif change == "stereo":
            m.topology.sites[11].metadata["chiral_tag"] = "CHI_TETRAHEDRAL_CW"
        elif change == "bond":
            m.topology.bonds[(11, 23)] = replace(m.topology.bonds[(11, 23)], order=1)
        elif change == "aromatic":
            m.topology.sites[11].metadata["aromatic"] = False
        else:
            del m.topology.sites[107]
        with pytest.raises(PCFFError):
            p.validate_system(m)
    moved = deepcopy(s)
    for i in moved.topology.sites:
        moved.coordinates.set(i, [1, 2, 3])
    p.validate_system(moved)


def test_profiled_dispatch_legacy_identity_and_owned_bundle(
    synthetic_profile_source, tmp_path, monkeypatch
):
    path, sel = synthetic_profile_source
    s = ring()
    legacy = prepare_forcefield(
        s,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                path, (0, 0, 1), (0, 0, 1), typing_profile="island_pcff_source_graph_v1"
            ),
        ),
    )
    prepared = prepare_forcefield(
        s,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                path,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v1",
                source_profile=sel,
            ),
        ),
    )
    assert legacy.native_result.identity == prepared.native_result.identity
    assert legacy.metadata["schema"] == "island_prepared_forcefield_v1"
    assert prepared.metadata["schema"] == "island_prepared_forcefield_v2"
    assert legacy.identity != prepared.identity
    assert prepared.operational_profile.payload["source_sha256"] == sel.sha256
    bundle = tmp_path / "bundle"
    inputs = PreparedForceFieldSources(pcff_frc=path)
    save_prepared_forcefield(s, prepared, bundle, sources=inputs)
    importer = builtins.__import__

    def blocked(name, *a, **kw):
        if name.split(".")[0] in {"openmm", "scipy", "rdkit", "foyer", "parmed"}:
            raise AssertionError(name)
        return importer(name, *a, **kw)

    with monkeypatch.context() as m:
        m.setattr(builtins, "__import__", blocked)
        loaded = load_prepared_forcefield(bundle, sources=inputs)
        assert loaded.prepared.identity == prepared.identity
        assert loaded.prepared.native_result.identity == legacy.native_result.identity
    graph = loaded.system
    graph.topology.sites[11].mass += 1
    with pytest.raises(PreparedForceFieldError):
        loaded.prepared.validate_integrity(graph)
    # A recomputed inner and outer checksum cannot authorize source/charge policy changes.
    profile_file = bundle / "pcff-profile.json"
    data = unpack(profile_file.read_text())
    data["authorized_charge_policy"] = "normalized"
    profile_file.write_text(pack(data))
    env = storage.read_json(bundle / "manifest.json")
    env["payload"]["files"]["pcff_profile"]["sha256"] = storage.checksum(
        profile_file.read_bytes()
    )
    env["payload"]["prepared"]["pcff_operational_profile"] = data
    env["sha256"] = storage.checksum(storage.json_bytes(env["payload"]))
    (bundle / "manifest.json").write_bytes(storage.json_bytes(env))
    with pytest.raises(PreparedBundleError):
        load_prepared_forcefield(bundle, sources=inputs)


def test_scope_failure_before_native_execution(synthetic_profile_source, monkeypatch):
    path, sel = synthetic_profile_source
    import island.forcefields.pcff as backend

    def forbidden(*a, **kw):
        raise AssertionError("typing executed before scope rejection")

    monkeypatch.setattr(backend, "type_pcff_atoms", forbidden)
    s = ring()
    s.topology.sites[11].formal_charge = 1
    with pytest.raises(PCFFError):
        prepare_forcefield(
            s,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    path,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v1",
                    source_profile=sel,
                ),
            ),
        )
    with pytest.raises(ForceFieldRequestError):
        PCFFOptions(path, (0, 0, 1), (0, 0, 1), source_profile=sel)


def test_source_model_tampering_and_nonzero_wilson(synthetic_profile_source):
    path, sel = synthetic_profile_source
    s = ring()
    p = prepare_forcefield(
        s,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                path, (0, 0, 1), (0, 0, 1), typing_profile="island_pcff_source_graph_v1"
            ),
        ),
    )
    native = p.native_result
    before = native.json_text
    data = unpack(before)
    for term in data["terms"]:
        if term["family"] == "wilson_out_of_plane":
            term["coefficients"][1] = 0.3
            break
    bad = type(native)(pack(data), native.assignment)
    with pytest.raises(PCFFError):
        sel.profile().validate_native(s, bad)
    data = unpack(before)
    data["terms"][0]["origin"] = "policy_derived_zero"
    data["terms"][0]["source_rows"] = []
    with pytest.raises(PCFFError):
        sel.profile().validate_native(s, type(native)(pack(data), native.assignment))
    assert native.json_text == before


def test_explicit_zero_is_required_source_record_not_missing_fallback(
    synthetic_profile_source, tmp_path, monkeypatch
):
    path, sel = synthetic_profile_source
    good = prepare_forcefield(
        ring(),
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                path,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v1",
                source_profile=sel,
            ),
        ),
    )
    zeros = [
        t
        for t in unpack(good.native_result.json_text)["terms"]
        if t["family"] == "bond-angle"
    ]
    assert len(zeros) == 18 and all(
        t["coefficients"] == [0, 0] and t["origin"] == "source_row" for t in zeros
    )
    raw = path.read_bytes()
    start = raw.index(b"#bond-angle cff91")
    end = raw.index(b"\n#", start + 1)
    absent = raw[:start] + b"#bond-angle cff91\n" + raw[end + 1 :]
    sha = source.digest(absent)
    monkeypatch.setitem(source.PIN, "sha256", sha)
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", sha)
    monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
    monkeypatch.setitem(DEFINITION, "source_sha256", sha)
    badpath = tmp_path / "missing.frc"
    badpath.write_bytes(absent)
    with pytest.raises(PreparedForceFieldError, match="Incomplete PCFF model"):
        prepare_forcefield(
            ring(),
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    badpath,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v1",
                    source_profile=PCFFOperationalSelection(NAME, sha),
                ),
            ),
        )


def test_wrong_source_and_charge_orientation_component_check(
    synthetic_profile_source, tmp_path
):
    path, sel = synthetic_profile_source
    system = ring()
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                path,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v1",
                source_profile=sel,
            ),
        ),
    )
    rows = unpack(prepared.native_result.json_text)["nonbonded"]
    assert {r["site"]: r["charge"] for r in rows} == {
        i: (-0.1 if a.element == "C" else 0.1) for i, a in system.topology.sites.items()
    }
    # Selecting a registered profile does not authorize different source bytes.
    impostor = tmp_path / "candidate.frc"
    impostor.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(PCFFError, match="hash mismatch"):
        prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    impostor,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v1",
                    source_profile=sel,
                ),
            ),
        )
    data = unpack(prepared.native_result.json_text)
    data["source"]["sha256"] = "0" * 64
    with pytest.raises(PCFFError):
        sel.profile().validate_native(
            system,
            type(prepared.native_result)(pack(data), prepared.native_result.assignment),
        )


def test_registered_profile_cannot_normalize_failed_component_charges(
    synthetic_profile_source, tmp_path, monkeypatch
):
    path, _ = synthetic_profile_source
    raw = path.read_bytes().replace(b"cp h -0.1 0.1", b"cp h -0.1 0.2")
    sha = source.digest(raw)
    for table, key in [
        (source.PIN, "sha256"),
        (expanded.PROFILE, "source_sha256"),
        (model.PROFILE, "frc_sha256"),
        (DEFINITION, "source_sha256"),
    ]:
        monkeypatch.setitem(table, key, sha)
    file = tmp_path / "residual.frc"
    file.write_bytes(raw)
    with pytest.raises(PCFFError, match="charge"):
        prepare_forcefield(
            ring(),
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    file,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v1",
                    source_profile=PCFFOperationalSelection(NAME, sha),
                ),
            ),
        )
    assert file.read_bytes() == raw


def test_operational_profile_does_not_authorize_wildcard_fallback_policy(
    synthetic_profile_source,
):
    path, selection = synthetic_profile_source
    with pytest.raises(ForceFieldRequestError, match="resolution options mismatch"):
        PCFFOptions(
            path,
            (0, 0, 1),
            (0, 0, 1),
            typing_profile="island_pcff_source_graph_v1",
            resolution_policy="island_pcff_positional_fallbacks_v2",
            source_profile=selection,
        )
