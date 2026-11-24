"""Whole-graph mapping, native coefficients and actual positive-time integration."""

import json
import os
import re
from importlib.util import find_spec

import numpy as np
import pytest

from island.exceptions import ForceFieldError
from island.forcefields.assignment import (
    assign_emc,
    assign_emc_graph,
    assign_polymer_parameters,
    run_emc_md,
)

pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem import AllChem

from island.builders import build_block_copolymer
from island.chemistry import from_rdkit, from_smiles
from island.core import SimulationBox
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.schema import Branch
from island.graph import final_graph

requires_emc = pytest.mark.skipif(
    not os.environ.get("EMC_ROOT") and find_spec("pyemc") is None, reason="requires EMC"
)


def section(data, name):
    rows = data.split("\n" + name + "\n\n", 1)[1].split("\n\n", 1)[0]
    return [line.split("#")[0].split() for line in rows.splitlines()]


@pytest.mark.parametrize("field", ["pcff", "pcff_ore", "opls/2012/opls-aa"])
@requires_emc
def test_block_all_sites_and_junctions(field):
    system = build_block_copolymer([("*CC*", 3), ("*CCO*", 2), ("*CC*", 1)])
    before = system.to_dict()
    r = assign_polymer_parameters(system, field)
    mapping = {i: v[0] for i, v in r.site_groups.items()}
    assert set(mapping.values()) == set(system.topology.sites) and len(mapping) == 40
    assert r.coordinate_source == "input_3d" and not r.added_hydrogen_ids
    edges = {
        tuple(sorted((mapping[int(a)], mapping[int(b)])))
        for _, _, a, b in section(r.data, "Bonds")
    }
    assert edges == set(system.topology.bonds)
    for row in section(r.data, "Atoms"):
        assert np.allclose(
            [float(i) for i in row[4:7]],
            system.coordinates.get(mapping[int(row[0])]),
            atol=1e-12,
        )
    assert r.site_to_atom_id == {v: k for k, v in mapping.items()}
    assert system.to_dict() == before


@requires_emc
def test_field_gap_at_ether_is_reported():
    ether = build_block_copolymer([("*CC*", 3), ("*CCO*", 2), ("*CC*", 1)])
    with pytest.raises(ForceFieldError, match=r"increment pair"):
        assign_emc_graph(ether, "opls")


def branch_graph():
    return normalize(PolymerSpec((("A", "*CC*"),), 4, branches=(Branch(4, "*CC"),)))


@pytest.mark.parametrize("field", ["pcff", "opls", "opls/2024/opls-ua"])
@requires_emc
def test_branched_graph_and_md(field, tmp_path):
    g = branch_graph()
    r = assign_emc_graph(g, field)
    ua = field.endswith("-ua")
    assert len(r.site_groups) == (10 if ua else 32)
    assert set(r.site_to_atom_id) == set(g.atom_order)
    assert r.coordinate_source == "etkdg" and not r.added_hydrogen_ids
    if ua:
        assert len(r.site_groups[r.site_to_atom_id[4]]) == 2
    else:
        assert all(len(v) == 1 for v in r.site_groups.values())
    bond = [r.site_to_atom_id[4], r.site_to_atom_id[9]]
    assert sorted(bond) in [
        sorted(map(int, row[2:4])) for row in section(r.data, "Bonds")
    ]
    if not (lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS")):
        pytest.skip("set ISLAND_ASSIGNMENT_LAMMPS to run the MD integration")
    run = run_emc_md(r, tmp_path / "md", lammps=lmp, steps=37, sample_every=10)
    assert (run.steps, run.time_fs, run.frames) == (37, 18.5, 5)
    traj = (
        (run.directory / "trajectory.lammpstrj")
        .read_text()
        .split("ITEM: TIMESTEP\n")[1:]
    )
    assert [int(f.splitlines()[0]) for f in traj] == [0, 10, 20, 30, 37]
    assert traj[0].splitlines()[8:] != traj[-1].splitlines()[8:]
    assert (
        json.loads((run.directory / "assignment.json").read_text())["field_sha256"]
        == r.field_sha256
    )
    assert (run.directory / "final.restart").stat().st_size > 0
    with pytest.raises(FileExistsError):
        run_emc_md(r, run.directory, lammps=lmp)


@pytest.mark.parametrize("field", ["pcff", "opls"])
@requires_emc
def test_full_graph_matches_native_chain_coefficients(field):
    g = build_block_copolymer([("*CC*", 4)])
    native = assign_emc("*CC*", field, dp=4)
    actual = assign_emc_graph(g, field)

    def coeffs(params):
        return sorted(
            tuple(row.split()[2:]) + (row.split()[0],)
            for row in params.splitlines()
            if re.match(r"^\w+_coeff\s", row)
        )

    assert coeffs(actual.parameters) == coeffs(native.parameters)
    for label in ["Bonds", "Angles", "Dihedrals"]:
        assert len(section(actual.data, label)) == len(section(native.data, label))


@requires_emc
def test_sparse_ids_box_fragments_and_graph_only():
    m = Chem.AddHs(Chem.MolFromSmiles("CCC.CCO"))
    assert AllChem.EmbedMolecule(m, randomSeed=12) == 0
    xyz = m.GetConformer().GetPositions()
    xyz[list(Chem.GetMolFrags(m)[1])] += 12
    xyz %= 32
    for i, pos in enumerate(xyz):
        m.GetConformer().SetAtomPosition(i, pos)
    ids = [101 + 7 * i for i in range(m.GetNumAtoms())][::-1]
    g = from_rdkit(m, site_ids=ids).system
    g.box = SimulationBox(32, 32, 32)
    r = assign_emc_graph(g, "pcff")
    assert set(r.site_to_atom_id) == set(ids)
    assert re.search(r"(?m)^0 32\s+xlo xhi$", r.data)
    atoms = section(r.data, "Atoms")
    assert len({row[1] for row in atoms}) == 2
    for row in atoms:
        site = r.site_groups[int(row[0])][0]
        shift = (np.array(row[4:7], float) - g.coordinates.get(site)) / 32
        assert np.allclose(shift, np.rint(shift))
    single = from_smiles("CCO", site_ids=[17, 31, 96, 24, 37, 81, 15, 73, 273])
    reread = assign_emc_graph(final_graph(single), "pcff")
    assert set(reread.site_to_atom_id) == set(single.topology.sites)
    assert reread.coordinate_source == "etkdg"


@requires_emc
def test_isotopic_site_and_double_bond_stereo():
    from rdkit.Chem.rdMolTransforms import GetDihedralDeg

    for smiles, trans in [(r"[2H]C/C=C/C", True), (r"[2H]C/C=C\C", False)]:
        mol = Chem.MolFromSmiles(smiles)
        r = assign_emc_graph(mol, "pcff")
        assert any(c in r.input_smiles for c in "/\\") and "[2H]" in r.input_smiles
        atoms = {int(row[0]): row for row in section(r.data, "Atoms")}
        masses = {int(i): float(m) for i, m in section(r.data, "Masses")}
        assert masses[int(atoms[r.site_to_atom_id[0]][2])] == pytest.approx(2.014101778)
        assert set(r.added_hydrogen_ids) == set(range(mol.GetNumAtoms(), len(atoms)))
        assert set(r.input_smiles_site_ids) == set(r.site_to_atom_id)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, tuple(map(float, atoms[r.site_to_atom_id[i]][4:7])))
        check = Chem.Mol(mol)
        check.AddConformer(conf)
        angle = abs(GetDihedralDeg(check.GetConformer(), 1, 2, 3, 4))
        assert (angle > 90) == trans


def test_invalid_graph_options():
    with pytest.raises(ValueError, match="omit dp"):
        assign_polymer_parameters(branch_graph(), "pcff", dp=4)
    with pytest.raises(ForceFieldError, match="terminated"):
        assign_emc_graph(Chem.MolFromSmiles("*CC*"), "opls")


@requires_emc
def test_wrapped_implicit_hydrogens_match_contiguous_molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=19) == 0
    mol = Chem.RemoveHs(mol)
    contiguous = assign_emc_graph(mol, "pcff")
    xyz = mol.GetConformer().GetPositions()
    for i, pos in enumerate(xyz % 32):
        mol.GetConformer().SetAtomPosition(i, pos)
    system = from_rdkit(mol, site_ids=range(len(xyz))).system
    system.box = SimulationBox(32, 32, 32)
    wrapped = assign_emc_graph(system, "pcff")
    assert wrapped.site_groups == contiguous.site_groups
    assert wrapped.added_hydrogen_ids == contiguous.added_hydrogen_ids
    coords = [
        np.array([row[4:7] for row in section(r.data, "Atoms")], dtype=float)
        for r in (wrapped, contiguous)
    ]
    assert np.allclose(coords[0] - coords[0][0], coords[1] - coords[1][0])
