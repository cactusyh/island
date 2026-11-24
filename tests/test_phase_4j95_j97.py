"""Chemical identities, periodic images, native force parity and checkpoint replay."""

import json
import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("scipy")
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from island.epoxy.chemistry import (
    Addition,
    amines,
    epoxides,
    react,
    reference_reactions,
)
from island.epoxy.cure import CureConfig, mixture, run_cure
from island.epoxy.network import bond_images, select, statistics
from island.epoxy.parameters import assign_network, compile_library
from island.epoxy.validation import validate_library


def test_mapped_double_addition():
    _, _, rx = reference_reactions()
    for stage, (pre, e, post) in enumerate(rx, 1):
        assert rdMolDescriptors.CalcMolFormula(pre) == rdMolDescriptors.CalcMolFormula(
            post
        )
        assert Chem.GetFormalCharge(pre) == Chem.GetFormalCharge(post) == 0
        assert pre.GetAtomWithIdx(e.hydrogen).GetNeighbors()[0].GetIdx() == e.nitrogen
        assert post.GetAtomWithIdx(e.hydrogen).GetNeighbors()[0].GetIdx() == e.oxygen
        assert len(amines(post)[e.nitrogen]) == 2 - stage
        assert len(epoxides(post)) == len(epoxides(pre)) - 1
        assert post.GetBondBetweenAtoms(e.carbon, e.nitrogen)
        with pytest.raises(ValueError):
            react(pre, [e, e])
        with pytest.raises(ValueError):
            react(post, [e])
    tertiary = rx[-1][2]
    c, o = next(iter(epoxides(tertiary).items()))
    with pytest.raises(ValueError):
        react(tertiary, [Addition(c, o, rx[-1][1].nitrogen, rx[-1][1].hydrogen)])


def test_periodic_search_and_rewrap():
    cfg = CureConfig(dgeba=2, pacm=1)
    _, mol = mixture(cfg)
    xyz = np.full((mol.GetNumAtoms(), 3), 10.0)
    n = next(iter(amines(mol)))
    c, o = next(iter(epoxides(mol).items()))
    xyz[list(amines(mol))] = [3, 14, 4]
    xyz[n] = [19, 4, 5]
    xyz[c] = [1, 4, 5]
    box = np.full(3, 20.0)

    def chosen(pos):
        return select(mol, pos, box, 3, np.random.default_rng(3), 1)

    e, rec = chosen(xyz)
    assert e == [Addition(c, o, n, amines(mol)[n][0])]
    assert rec["events"][0]["image"] == [-1, 0, 0]
    assert rec["events"][0]["distance_A"] == 2
    shifted = (xyz + [3.3, 4.5, 5]) % box
    ee, rr = chosen(shifted)
    assert ee == e and rr["events"][0]["distance_A"] == pytest.approx(2)
    assert rr["events"][0]["image"] == [0, 0, 0]


def test_periodic_ring_closure():
    mol = Chem.MolFromSmiles("C1CCCCC1")
    for a in mol.GetAtoms():
        a.SetIntProp("monomer", 0)
    box = np.array([12, 12, 12])
    finite = np.array(
        [[11, 4, 4], [0, 4, 4], [1, 5, 4], [1, 6, 4], [0, 6, 4], [11, 5, 4]]
    )
    winding = np.array([[1 + 2 * i, 4, 4] for i in range(6)])
    assert bond_images(mol, finite, box)
    assert statistics(mol, finite, box, 1, 1)["winding_rank"] == 0
    for translation in [0, 3, 11]:
        xyz = (winding + translation) % box
        assert statistics(mol, xyz, box, 1, 1)["winding_rank"] == 1
    finite[1] = [5, 4, 4]
    with pytest.raises(ValueError, match="Ambiguous"):
        bond_images(mol, finite, box)


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    p = os.environ.get("ISLAND_EPOXY_LIBRARY")
    if p:
        return Path(p)
    pytest.importorskip("pyemc")
    out = tmp_path_factory.mktemp("epoxy") / "native"
    compile_library(out)
    return out


def test_network_types_survive_permutation(library):
    lib = json.loads((library / "library.json").read_text())
    m = reference_reactions()[2][-1][2]
    order = np.random.default_rng(82).permutation(m.GetNumAtoms()).tolist()
    permuted = Chem.RenumberAtoms(m, order)
    a, terms = assign_network(m, lib)
    b, other = assign_network(permuted, lib)
    assert b == [a[i] for i in order]

    def keyed(items, ids):
        return sorted((t, tuple(sorted(ids[i] for i in ss))) for t, ss in items)

    for section, items in terms.items():
        # PCFF can express the same improper in different permutations with
        # correspondingly permuted angle-angle coefficients.
        if section == "Impropers":
            assert sorted(tuple(sorted(ss)) for _, ss in items) == sorted(
                tuple(sorted(order[j] for j in ss)) for _, ss in other[section]
            )
            continue
        assert keyed(items, range(len(a))) == keyed(other[section], order)


def test_direct_native_energies_and_forces(library, tmp_path):
    if not (exe := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS")):
        pytest.skip("LAMMPS integration")
    checks = validate_library(library, tmp_path, exe)
    assert len(checks) == 12 and all(x["max_force_error"] < 1e-5 for x in checks)


def test_pause_resume_after_addition(library, tmp_path):
    if not (exe := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS")) or not (
        pack := os.environ.get("ISLAND_PACKMOL")
    ):
        pytest.skip("PACKMOL/LAMMPS integration")
    cfg = CureConfig(
        dgeba=4,
        pacm=2,
        seed=41,
        batch_size=1,
        pre_steps=200,
        batch_steps=100,
        post_steps=100,
        target=0.25,
        cutoff=5.5,
        max_idle=3,
    )
    args = {"lammps": exe, "packmol": pack}
    lib = library / "library.json"
    full = run_cure(tmp_path / "full", lib, cfg, **args)
    out = tmp_path / "split"
    run_cure(out, lib, cfg, stop_after=2, **args)
    resumed = run_cure(out, lib, cfg, **args)
    assert len(resumed) == len(full) and resumed[-1]["statistics"]["crosslinks"] > 0
    assert [r["events"] for r in full] == [r["events"] for r in resumed]
    assert [r["thermo"] for r in full] == [r["thermo"] for r in resumed]
    file = f"stage-{len(full) - 1:04d}/state.npz"
    with np.load(tmp_path / "full" / file) as a, np.load(out / file) as b:
        for k in a.files:
            np.testing.assert_array_equal(a[k], b[k])
    with pytest.raises(ValueError, match="settings"):
        run_cure(out, lib, replace(cfg, seed=3), **args)
    (out / file).write_bytes(b"bad")
    with pytest.raises(ValueError, match="checkpoint"):
        run_cure(out, lib, cfg, **args)
