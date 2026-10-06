"""J9 synthetic software contracts; real pinned-source acceptance is separate."""

import builtins
from dataclasses import replace

import pytest
from test_pcff_operational_profile import (  # noqa: F401
    ring,
    synthetic_profile_source,
)

from island import Coordinates
from island.charge_references.records import pack
from island.exceptions import PCFFError
from island.forcefields import (
    PCFFOptions,
    PreparedForceFieldSources,
    load_prepared_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    PCFFOperationalProfile,
    edit_polymer_topology,
    export_msi2lmp_car_mdf,
    inspect_pcff_profile,
    list_pcff_profiles,
    load_msi2lmp_car_mdf,
    prepare_pcff_polymer,
    reassign_pcff_polymer,
    select_pcff_profile,
)
from island.forcefields.pcff.operational_profile import DEFINITION, NAME
from island.forcefields.pcff.registry import LINKED_BENZENOID
from island.workflows import storage


def linked():
    a = ring()
    b = ring()
    for i, site in b.topology.sites.items():
        a.topology.add_site(replace(site, id=i + 1000))
        a.coordinates.set(i + 1000, b.coordinates.get(i) + [4, 0, 0])
    for bond in b.topology.bonds.values():
        a.topology.add_bond(
            bond.site1 + 1000,
            bond.site2 + 1000,
            order=bond.order,
            aromatic=bond.aromatic,
        )
    for i in (107, 1107):
        a.topology.remove_site(i)
    a.coordinates = Coordinates({i: a.coordinates.get(i) for i in a.topology.sites})
    a.topology.add_bond(11, 1011, order=1)
    for i, site in a.topology.sites.items():
        site.metadata.update(
            chain_id="A", repeat_unit_index=int(i >= 1000), repeat_unit_type="phenylene"
        )
    a.metadata["polymer"] = {
        "sequence": ["P", "P"],
        "coordinate_provenance": "synthetic software geometry",
    }
    return a


def options(path):
    return PCFFOptions(
        path,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        source_profile=select_pcff_profile(
            LINKED_BENZENOID, sha256=DEFINITION["source_sha256"]
        ),
    )


def test_registry_owned_states_and_old_identity():
    entries = list_pcff_profiles()
    assert {p["state"] for p in entries} == {
        "operational",
        "audited",
        "candidate",
        "unsupported",
    }
    entries[0]["definition"]["authorized_atom_labels"].append("nh+")
    assert (
        "nh+" not in inspect_pcff_profile(NAME)["definition"]["authorized_atom_labels"]
    )
    old = select_pcff_profile(NAME, sha256=DEFINITION["source_sha256"]).profile()
    assert (
        old.identity
        == "db14223b8527ca7d4f9cafb3cca2b40610ff6c944cde18108ffa46dcbb38228d"
    )
    for entry in list_pcff_profiles():
        if entry["state"] != "operational":
            with pytest.raises(PCFFError):
                select_pcff_profile(entry["name"], sha256=entry["source_sha256"])
    p = select_pcff_profile(
        LINKED_BENZENOID, sha256=DEFINITION["source_sha256"]
    ).profile()
    for key, value in [
        ("source_sha256", "0" * 64),
        ("authorized_charge_policy", "normalize"),
        ("chemical_rule", "all_polymers"),
        ("unsupported_families", []),
    ]:
        changed = p.payload
        changed[key] = value
        with pytest.raises(PCFFError):
            PCFFOperationalProfile(pack(changed)).validate_integrity()


@pytest.mark.parametrize(
    "change",
    [
        "element",
        "number",
        "charge",
        "isotope",
        "radical",
        "stereo",
        "mass",
        "order",
        "missing_h",
    ],
)
def test_final_graph_domain_boundaries(change):
    s = linked()
    profile = select_pcff_profile(
        LINKED_BENZENOID, sha256=DEFINITION["source_sha256"]
    ).profile()
    profile.validate_system(s)
    site = s.topology.sites[11]
    if change == "element":
        site.element, site.atomic_number = "N", 7
    elif change == "number":
        site.atomic_number = 7
    elif change == "charge":
        site.formal_charge = 1
    elif change == "isotope":
        site.metadata["isotope"] = 13
    elif change == "radical":
        site.metadata["radical_electrons"] = 1
    elif change == "stereo":
        site.metadata["chiral_tag"] = "CHI_TETRAHEDRAL_CW"
    elif change == "mass":
        site.mass += 2
    elif change == "order":
        s.topology.remove_bond(11, 1011)
        s.topology.add_bond(11, 1011, order=2)
    else:
        s.topology.remove_site(211)
        s.coordinates = Coordinates({i: s.coordinates.get(i) for i in s.topology.sites})
    with pytest.raises(PCFFError):
        profile.validate_system(s)


def test_final_graph_strict_diagnostic_update_and_bundle(
    tmp_path,
    synthetic_profile_source,  # noqa: F811
):
    path, _ = synthetic_profile_source
    s = linked()
    before = s.to_dict()
    result = prepare_pcff_polymer(s, options(path))
    assert s.to_dict() == before
    from island.forcefields import adopt_forcefield
    from island.forcefields.pcff import PCFFFinalGraphPreparation

    unprofiled = adopt_forcefield(s, "pcff", result.prepared.native_result)
    with pytest.raises(PCFFError, match="operational profile"):
        PCFFFinalGraphPreparation(result.json_text, unprofiled).validate_integrity(s)
    assert any(r["roles"]["inter_repeat"] for r in result.payload["interactions"])
    diagnostic = prepare_pcff_polymer(s, options(path), mode="diagnostic")
    assert diagnostic["diagnostic_only"] and "prepared" not in diagnostic
    sources = PreparedForceFieldSources(pcff_frc=path)
    save_prepared_forcefield(s, result.prepared, tmp_path / "bundle", sources=sources)
    (tmp_path / "bundle").rename(tmp_path / "relocated")
    loaded = load_prepared_forcefield(tmp_path / "relocated", sources=sources)
    assert loaded.prepared.identity == result.prepared.identity
    assert loaded.system.to_dict() == before
    changed = edit_polymer_topology(
        s, remove_hydrogens=[313, 1313], add_crosslinks=[(37, 1037)]
    )
    with pytest.raises(PCFFError):
        result.prepared.native_result.validate_integrity(changed)
    updated, receipt = reassign_pcff_polymer(s, result, changed, options(path))
    assert (
        updated.prepared.native_result.identity
        != result.prepared.native_result.identity
    )
    assert (
        receipt["locally_retyped_sites"] and receipt["no_cached_interaction_inventory"]
    )
    assert any(r["roles"]["crosslink"] for r in updated.payload["interactions"])
    assert s.to_dict() == before
    bad = updated.payload
    bad["interactions"][0]["roles"]["crosslink"] = not bad["interactions"][0]["roles"][
        "crosslink"
    ]
    with pytest.raises(PCFFError):
        replace(updated, json_text=pack(bad)).validate_integrity(changed)
    # Compatible coordinates do not require reparameterization.
    moved = s.copy()
    moved.coordinates.translate([1, 2, 3])
    result.prepared.validate_integrity(moved)


def test_interop_owned_roundtrip_and_rechecksummed_tamper(tmp_path):
    s = linked()
    # Exact integer metadata and atom roles survive; no external program required.
    s.metadata["integer_map"] = {11: ("end", "crosslink")}
    before = s.to_dict()
    t = {i: "cp" if a.element == "C" else "hc" for i, a in s.topology.sites.items()}
    q = {i: 0.0 for i in t}  # explicitly synthetic software charges
    path = tmp_path / "interop"
    restored = export_msi2lmp_car_mdf(
        s, path, types=t, charges=q, provenance={"origin": "synthetic test"}
    )
    assert restored.to_dict() == before
    assert "/1.5" in (path / "system.mdf").read_text()
    restored.topology.sites[11].mass += 1
    assert load_msi2lmp_car_mdf(path).to_dict() == before
    with pytest.raises(PCFFError):
        export_msi2lmp_car_mdf(s, path, types=t, charges=q)
    mdf = (path / "system.mdf").read_text().replace("/1.5", "/1", 1)
    (path / "system.mdf").write_text(mdf)
    manifest = storage.read_json(path / "manifest.json")
    manifest["files"]["system.mdf"] = storage.checksum(
        (path / "system.mdf").read_bytes()
    )
    (path / "manifest.json").write_bytes(storage.json_bytes(manifest))
    with pytest.raises(PCFFError, match="contradiction"):
        load_msi2lmp_car_mdf(path)
    assert s.to_dict() == before


def test_offline_registry_and_interop_do_not_import_scientific_backends(
    tmp_path, monkeypatch
):
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in ("rdkit", "openmm", "foyer", "scipy", "parmed"):
            raise AssertionError("Unexpected scientific import: " + name)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    s = linked()
    assert list_pcff_profiles()
    profile = select_pcff_profile(
        LINKED_BENZENOID, sha256=DEFINITION["source_sha256"]
    ).profile()
    profile.validate_system(s)
    t = {i: "cp" if a.element == "C" else "hc" for i, a in s.topology.sites.items()}
    export_msi2lmp_car_mdf(s, tmp_path / "files", types=t, charges={i: 0 for i in t})
    assert load_msi2lmp_car_mdf(tmp_path / "files").to_dict() == s.to_dict()


def test_explicit_end_groups_and_default_construction_preserved():
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer, build_polymer_from_sequence
    from island.exceptions import PolymerBuildError

    p = "[*:1]c1ccc([*:2])cc1"
    default = build_linear_polymer(p, dp=2, generate_3d=False)
    assert "end_groups" not in default.metadata["polymer"]
    capped = build_linear_polymer(
        p, dp=2, generate_3d=False, head_end_group="[*]C", tail_end_group="[*]C"
    )
    assert capped.number_of_sites == default.number_of_sites + 6
    assert {
        a.metadata.get("end_group_role") for a in capped.topology.sites.values()
    } == {None, "head", "tail"}
    select_pcff_profile(
        LINKED_BENZENOID, sha256=DEFINITION["source_sha256"]
    ).profile().validate_system(capped)
    mixed = build_polymer_from_sequence(
        {"P": p, "M": "[*:1]c1cc([*:2])ccc1"},
        ["P", "M", "P"],
        generate_3d=False,
        head_end_group="[*]C",
    )
    assert mixed.metadata["polymer"]["sequence"] == ["P", "M", "P"]
    for group in ("[*]C[*]", "C", "[*]=C", "[*]C.C"):
        with pytest.raises(PolymerBuildError):
            build_linear_polymer(p, dp=1, generate_3d=False, head_end_group=group)


def test_rechecksummed_candidate_registration_cannot_promote(tmp_path):
    from island.forcefields.pcff import (
        load_pcff_profile_registration,
        save_pcff_profile_registration,
    )

    path = tmp_path / "registration.json"
    registration = save_pcff_profile_registration(
        "iff_pcff_interface_v1_5_candidate", path
    )
    assert load_pcff_profile_registration(path).identity == registration.identity
    changed = registration.payload
    changed["state"] = "operational"
    changed["capabilities"]["evaluation"] = True
    path.write_text(pack(changed))
    with pytest.raises(PCFFError):
        load_pcff_profile_registration(path)


def test_failed_interop_publication_is_not_a_usable_export(tmp_path, monkeypatch):
    system = linked()
    labels = {
        i: "cp" if a.element == "C" else "hc" for i, a in system.topology.sites.items()
    }
    original = storage.publish

    def failed(path, raw, **kwargs):
        if path.name == "system.mdf":
            raise OSError("synthetic publication failure")
        return original(path, raw, **kwargs)

    monkeypatch.setattr(storage, "publish", failed)
    with pytest.raises(PCFFError, match="publication failure"):
        export_msi2lmp_car_mdf(
            system, tmp_path / "failed", types=labels, charges={i: 0 for i in labels}
        )
    assert not (tmp_path / "failed").exists()
    assert not list(tmp_path.glob(".interop-*"))


def test_intraring_links_and_unsupported_public_preparation_fail_before_source_read(
    tmp_path,
):
    s = linked()
    changed = edit_polymer_topology(
        s, remove_hydrogens=(211, 523), add_crosslinks=((23, 59),)
    )
    with pytest.raises(PCFFError, match="inside one aromatic ring"):
        prepare_pcff_polymer(changed, options(tmp_path / "not-a-source.frc"))
    assert s.number_of_sites == 22


def test_end_group_atom_stereo_preserved_and_bond_stereo_rejected():
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer
    from island.exceptions import PolymerBuildError

    p = "[*:1]CC[*:2]"
    a = build_linear_polymer(p, dp=1, generate_3d=False, head_end_group="[*][C@H](O)C")
    b = build_linear_polymer(p, dp=1, generate_3d=False, head_end_group="[*][C@@H](O)C")

    def label(s):
        return [
            site.metadata["cip_label"]
            for site in s.topology.sites.values()
            if site.metadata.get("end_group_role") == "head"
            and site.metadata.get("cip_label")
        ]

    assert len(label(a)) == len(label(b)) == 1
    assert label(a) != label(b)
    with pytest.raises(PolymerBuildError, match="bond stereochemistry"):
        build_linear_polymer(p, dp=1, generate_3d=False, head_end_group="[*]C/C=C/C")
