"""Full-graph periodic GAFF cure, reaction conservation and durable resume."""

import os

import numpy as np
import pytest
from rdkit import Chem

from island.epoxy.chemistry import AMINE, THIOL, epoxides, monomer, react
from island.epoxy.cure import CureConfig, mixture, run_cure
from island.epoxy.gaff import compile_library, parameterize
from island.workflows.recipe import run_recipe


def small_config(**kw):
    return CureConfig(
        dgeba=2,
        pacm=1,
        pre_steps=20,
        batch_steps=20,
        post_steps=20,
        max_idle=1,
        cutoff=4.5,
        target=0.25,
        packing_density=0.3,
        **kw,
    )


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
@pytest.mark.parametrize("template,donor", [(AMINE, "CN"), (THIOL, "CS")])
def test_junction_conservation_and_native_coefficients(
    tmp_path, family, template, donor
):
    from radonpy.core import calc
    from radonpy.ff.gaff import GAFF
    from radonpy.ff.gaff2 import GAFF2

    from island.forcefields.radonpy import RadonPyAssignment
    from island.workflows.recipe_md import radonpy_input

    cfg = small_config()
    lib = compile_library(tmp_path / "lib", family, cfg)
    a, b = monomer("CC1CO1", 42), monomer(donor, 43)
    mol = Chem.CombineMols(a, b)
    Chem.SanitizeMol(mol)
    c, o = next(iter(epoxides(mol).items()))
    d, hs = next((d, h) for d, h in template.donors(mol).items() if h)
    product = react(mol, [template.event_type(c, o, d, hs[0])], template)
    assert [(a.GetAtomicNum(), a.GetFormalCharge()) for a in product.GetAtoms()] == [
        (a.GetAtomicNum(), a.GetFormalCharge()) for a in mol.GetAtoms()
    ]
    xyz = product.GetConformer().GetPositions()
    box = np.ones(3) * 40
    data, params, report = parameterize(product, lib, xyz, box)
    native = Chem.Mol(product)
    engine = {"gaff": GAFF, "gaff2": GAFF2}[family]()
    assert engine.ff_assign(native)
    assert calc.assign_charges(native, charge="gasteiger")
    assignment = RadonPyAssignment(
        native,
        tuple(range(native.GetNumAtoms())),
        family,
        "test",
        "test",
        "gasteiger",
        engine.param.lj_c14,
        engine.param.c_c14,
    )
    assert (data, params) == radonpy_input(assignment, box)
    shifted, shifted_params, shifted_report = parameterize(
        product, lib, xyz + np.array([40, -40, 80]), box
    )
    assert shifted != data and shifted_params == params
    assert shifted_report == report
    assert abs(report["charge"]) < 1e-6


def test_custom_functionality_and_ineligible_amide(tmp_path):
    cfg = small_config(epoxy_smiles="CC1CO1", hardener_smiles="NCCN")
    monomers, mol = mixture(cfg)
    assert len(epoxides(mol)) == 2  # not the old hard-coded four
    assert sum(map(len, cfg.template.donors(monomers[1]).values())) == 4
    with pytest.raises(ValueError, match="donor"):
        small_config(hardener_smiles="CC(=O)N")
    with pytest.raises(ValueError, match="epoxide"):
        small_config(epoxy_smiles="CCO")
    lib = compile_library(tmp_path / "lib", "gaff2", cfg)
    lib["field_sha256"] = {}
    with pytest.raises(ValueError, match="source"):
        parameterize(mol, lib, mol.GetConformer().GetPositions(), [40] * 3)


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
def test_real_custom_cure_resume(tmp_path, family):
    cfg = small_config(epoxy_smiles="CC1CO1", hardener_smiles="NCCN")
    compile_library(tmp_path / "lib", family, cfg)
    library = tmp_path / "lib/library.json"
    kw = {
        "lammps": os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        "packmol": os.environ["ISLAND_PACKMOL"],
    }
    partial = run_cure(tmp_path / "run", library, cfg, stop_after=1, **kw)
    before = (tmp_path / "run/stage-0000/state.npz").read_bytes()
    result = run_cure(tmp_path / "run", library, cfg, **kw)
    assert partial[0] == result[0]
    assert before == (tmp_path / "run/stage-0000/state.npz").read_bytes()
    assert result[-1]["phase"] == "post"
    assert result[-1]["statistics"]["crosslinks"] >= 1
    assert (
        result[-1]["statistics"]["epoxy_conversion"]
        == result[-1]["statistics"]["crosslinks"] / 2
    )
    assert run_cure(tmp_path / "run", library, cfg, **kw) == result
    record_path = tmp_path / "run/stage-0000/record.json"
    raw = record_path.read_text()
    record_path.write_text(raw + " ")
    with pytest.raises(ValueError, match="event record"):
        run_cure(tmp_path / "run", library, cfg, **kw)
    record_path.write_text(raw)
    path = tmp_path / "run/stage-0000/system.params"
    path.write_text(path.read_text() + "# changed\n")
    with pytest.raises(ValueError, match="artifact"):
        run_cure(tmp_path / "run", library, cfg, **kw)


def test_crosslinked_gaff_recipe(tmp_path):
    recipe = {
        "architecture": "crosslinked",
        "forcefield": {"name": "gaff2", "charge": "gasteiger"},
        "density_g_cm3": 0.3,
        "seed": 2026,
        "cure": {
            "dgeba": 2,
            "pacm": 1,
            "epoxy_smiles": "CC1CO1",
            "hardener_smiles": "NCCN",
            "pre_steps": 20,
            "batch_steps": 20,
            "post_steps": 20,
            "max_idle": 1,
            "cutoff": 4.5,
            "target": 0.25,
        },
        "md": {
            "nvt_steps": 20,
            "npt_steps": 20,
            "production_steps": 20,
            "sample_every": 10,
            "segment_steps": 20,
        },
    }
    result = run_recipe(
        recipe,
        tmp_path / "recipe",
        lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        packmol=os.environ["ISLAND_PACKMOL"],
    )
    assert result["status"] == "completed" and result["target_reached"]
    assert result["assignment"]["charge_method"] == "gasteiger"
    assert result["equilibrated"] == "unmeasured"


def test_periodically_winding_graph_has_no_finite_embedding_requirement(tmp_path):
    from island.epoxy.network import statistics
    from island.epoxy.validation import evaluate

    mol = Chem.AddHs(Chem.MolFromSmiles("C1CCCCCCC1"))
    xyz = np.zeros((mol.GetNumAtoms(), 3))
    for i in range(8):
        xyz[i] = [1.5 * i, 6, 6]
        hs = [
            a.GetIdx()
            for a in mol.GetAtomWithIdx(i).GetNeighbors()
            if a.GetAtomicNum() == 1
        ]
        for k, h in enumerate(hs):
            xyz[h] = xyz[i] + [0, (-1) ** k, 0.5]
    for atom in mol.GetAtoms():
        atom.SetIntProp("monomer", 0)
    box = np.array([12, 12, 12])
    assert statistics(mol, xyz, box, 1, 1)["winding_rank"] == 1
    lib = compile_library(tmp_path / "lib", "gaff2", small_config())
    data, params, _ = parameterize(mol, lib, xyz, box)
    shifted, _, _ = parameterize(mol, lib, xyz + box, box)
    first = evaluate(
        tmp_path / "first",
        data,
        params,
        os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        direct_ewald=True,
    )
    second = evaluate(
        tmp_path / "shifted",
        shifted,
        params,
        os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        direct_ewald=True,
    )
    assert np.allclose(first[2], second[2], atol=1e-7)
    assert np.allclose(first[3], second[3], atol=1e-7)


def test_real_thiol_cure(tmp_path):
    cfg = CureConfig(
        dgeba=2,
        petmp=1,
        epoxy_smiles="CC1CO1",
        hardener_smiles="SCCS",
        donor_angle_tolerance=180,
        pre_steps=20,
        batch_steps=20,
        post_steps=20,
        max_idle=1,
        cutoff=4.5,
        target=0.25,
        packing_density=0.3,
    )
    compile_library(tmp_path / "lib", "gaff2", cfg)
    rows = run_cure(
        tmp_path / "run",
        tmp_path / "lib/library.json",
        cfg,
        lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
        packmol=os.environ["ISLAND_PACKMOL"],
    )
    assert rows[-1]["phase"] == "post"
    assert rows[-1]["statistics"]["crosslinks"] >= 1
    assert rows[-1]["statistics"]["thiol_h_conversion"] > 0
