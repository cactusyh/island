"""Chemical oracles for termination, stereo transport and public entry semantics."""

import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from island.builders import build_linear_polymer
from island.chemistry import from_rdkit, to_rdkit
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.schema import Branch
from island.forcefields.assignment_core.topology import as_topology, molecule
from island.graph import final_graph
from island.graph.final import _system_from_payload


def smiles(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


@pytest.mark.parametrize("dp", [1, 2, 3, 6])
@pytest.mark.parametrize("repeat", ["*C[C@H](C)*", "*O[Si](C)(C)*", "*[NH]CC*"])
@pytest.mark.parametrize("hydrogens", [False, True])
def test_default_caps_equal_explicit_h_caps(dp, repeat, hydrogens):
    default = build_linear_polymer(
        repeat, dp=dp, generate_3d=False, add_hydrogens=hydrogens
    )
    explicit = build_linear_polymer(
        repeat,
        dp=dp,
        generate_3d=False,
        add_hydrogens=hydrogens,
        head_end_group="*[H]",
        tail_end_group="*[H]",
    )
    mol = to_rdkit(default).mol
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    assert smiles(mol) == smiles(to_rdkit(explicit).mol)
    if not hydrogens:
        assert not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    graph = normalize(PolymerSpec((("A", repeat),), dp))
    assert not graph.blockers
    assert smiles(molecule(graph)) == smiles(mol)


@pytest.mark.parametrize(
    ("repeat", "caps", "expected"),
    [
        (r"*C/C=C\C*", {}, r"C/C=C\CC/C=C\CC/C=C\C"),
        ("*/C=C/*", {}, "C=C/C=C/C=C"),
        (r"*/C=C\*", {}, r"C=C/C=C\C=C"),
        ("*CC*", {"head_end_group": "*/C=C/C"}, "C/C=C/CCCCCC"),
        (
            "*[C@H](*)Cl",
            {"head_end_group": "*F", "tail_end_group": "*Br"},
            "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)Br",
        ),
        ("[*:2]OC[*:1]", {"head_end_group": "*F", "tail_end_group": "*C"}, "FCOCOCOC"),
    ],
)
def test_stereo_matches_independent_capped_smiles(repeat, caps, expected):
    system = build_linear_polymer(repeat, dp=3, generate_3d=False, **caps)
    expected = smiles(Chem.MolFromSmiles(expected))
    assert smiles(to_rdkit(system).mol) == expected
    restored = _system_from_payload(json.loads(json.dumps(system.to_dict())))
    assert smiles(to_rdkit(restored).mol) == expected
    graph = normalize(
        PolymerSpec(
            (("A", repeat),),
            3,
            end_groups=(caps.get("head_end_group"), caps.get("tail_end_group")),
        )
    )
    assert not graph.blockers
    assert smiles(molecule(graph)) == expected


@pytest.mark.parametrize("dp", [1, 2, 3, 6])
@pytest.mark.parametrize("caps", [(None, None), ("*C", "*C"), ("*F", "*Br")])
def test_tacticity_records_only_realizable_centers(dp, caps):
    system = build_linear_polymer(
        "*C[C@H](C)*",
        dp=dp,
        tacticity="syndiotactic",
        generate_3d=False,
        head_end_group=caps[0],
        tail_end_group=caps[1],
    )
    meta = system.metadata["polymer"]
    requested = meta["requested_stereochemical_sequence"]
    actual = [
        a.metadata["cip_label"]
        for a in system.topology.sites.values()
        if a.metadata.get("controllable_stereocenter")
    ]
    assert actual == meta["stereochemical_sequence"]
    for i, (observed, target) in enumerate(zip(actual, requested, strict=True)):
        assert observed == target or (observed is None and i in (0, dp - 1))
    if caps == ("*F", "*Br"):
        assert actual == requested


def test_ez_roundtrip_nonmonotonic_ids_and_identity():
    def system(smi):
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.Compute2DCoords(m)
        return from_rdkit(m, site_ids=list(range(m.GetNumAtoms(), 0, -1))).system

    e, z = system("C/C=C/C"), system(r"C/C=C\C")
    assert smiles(to_rdkit(e).mol) == "C/C=C/C"
    assert smiles(to_rdkit(z).mol) == r"C/C=C\C"
    assert final_graph(e).identity != final_graph(z).identity
    from island.chemistry import system_to_rdkit_graph
    from island.forcefields._emc_graph_input import _from_final

    assert smiles(system_to_rdkit_graph(_from_final(final_graph(e))).mol) == "C/C=C/C"
    from island.forcefields.pcff.automatic import chemical_graph, graph_system

    restored = graph_system(chemical_graph(e))
    assert smiles(system_to_rdkit_graph(restored).mol) == "C/C=C/C"


def test_branch_bond_stereo_survives_dummy_removal():
    graph = normalize(
        PolymerSpec((("A", "*CC*"),), 2, branches=(Branch(2, r"*/C=C\C"),))
    )
    assert not graph.blockers
    assert smiles(molecule(graph)) == smiles(Chem.MolFromSmiles(r"CCC(C)/C=C\C"))
    top = as_topology(graph)
    assert any(b.stereo != "STEREONONE" for b in top.bonds.values())


def test_public_emc_default_uses_completed_graph(monkeypatch):
    from island.forcefields import assignment

    seen = []
    monkeypatch.setattr(
        assignment,
        "assign_emc_graph",
        lambda graph, *a, **kw: seen.append(graph) or graph,
    )
    monkeypatch.setattr(assignment, "assign_emc", lambda *a, **kw: "legacy")
    result = assignment.assign_polymer_parameters(
        "*CCO*", "opls", dp=3, head_end_group="*C"
    )
    assert smiles(to_rdkit(result).mol) == smiles(Chem.MolFromSmiles("CCCOCCOCCO"))
    assert (
        assignment.assign_polymer_parameters("*CC*", "pcff", dp=3, emc_mode="repeat")
        == "legacy"
    )
    with pytest.raises(ValueError, match="emc_mode"):
        assignment.assign_polymer_parameters("*CC*", "pcff", dp=3, nchains=2)
    with pytest.raises(ValueError, match="fixed connectivity"):
        assignment.assign_polymer_parameters(result, "pcff", head_end_group="*F")


@pytest.mark.parametrize("method", ["etkdg", "local_templates"])
def test_double_bond_geometry_is_retained_in_3d(method):
    from island.chemistry._stereo_transfer import verify_bond_geometry

    system = build_linear_polymer(r"*C/C=C\C*", dp=3, coordinate_method=method)
    mol = to_rdkit(system).mol
    verify_bond_geometry(mol)
    Chem.RemoveStereochemistry(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    assert smiles(mol) == smiles(Chem.MolFromSmiles(r"C/C=C\CC/C=C\CC/C=C\C"))


def test_reversed_attachment_labels_preserve_tetrahedral_parity():
    system = build_linear_polymer(
        "[*:2][C@H](F)C[*:1]",
        dp=3,
        generate_3d=False,
        head_end_group="*F",
        tail_end_group="*Br",
    )
    assert smiles(to_rdkit(system).mol) == smiles(
        Chem.MolFromSmiles("FC[C@@H](F)C[C@@H](F)C[C@@H](F)Br")
    )


def test_invalid_stereo_controls_are_rejected():
    from dataclasses import replace

    from island.exceptions import TopologyError

    system = build_linear_polymer("*C/C=C/C*", dp=1, generate_3d=False)
    bond = next(b for b in system.topology.bonds.values() if b.stereo_atoms)
    system.topology.bonds[bond.key] = replace(bond, stereo_atoms=(999, 1000))
    with pytest.raises(TopologyError, match="stereo control"):
        to_rdkit(system)


def test_packing_remaps_ez_controls_with_site_ids():
    from island.packing import PeriodicPackingConfig, pack_multichain_periodic

    system = build_linear_polymer(r"*C/C=C\C*", dp=2)
    graph = final_graph(system)
    result = pack_multichain_periodic(
        [("a", system, graph), ("b", system, graph)],
        config=PeriodicPackingConfig(box_lengths=(40.0, 40.0, 40.0)),
        provenance="J106 E/Z transport regression",
        evidence=("specified cis polymer",),
    )
    mol = to_rdkit(result.system).mol
    fragments = Chem.GetMolFrags(mol, asMols=True)
    assert len(fragments) == 2
    assert all(smiles(part) == smiles(to_rdkit(system).mol) for part in fragments)
    result.graph.validate_integrity(result.system)


def test_wrong_double_bond_geometry_is_rejected():
    from island.chemistry._stereo_transfer import verify_bond_geometry
    from island.exceptions import UnsupportedStereochemistryError

    mol = to_rdkit(build_linear_polymer("*C/C=C/C*", dp=1)).mol
    bond = next(b for b in mol.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOE)
    bond.SetStereo(Chem.BondStereo.STEREOZ)
    with pytest.raises(UnsupportedStereochemistryError, match="double-bond"):
        verify_bond_geometry(mol)
