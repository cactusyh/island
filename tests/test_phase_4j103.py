"""Recipe contracts, native export parity and stage-resume integration."""

import json
import os
from copy import deepcopy

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("scipy")
from rdkit import Chem
from rdkit.Chem import AllChem

from island.exceptions import ForceFieldError
from island.forcefields.assignment import assign_polymer_parameters
from island.workflows.recipe import checkpoint, digest, run_recipe
from island.workflows.recipe_config import PolymerRecipe
from island.workflows.recipe_md import radonpy_input


def recipe():
    return {
        "architecture": "homopolymer",
        "polymer": {"repeat_units": {"A": "*CC*"}, "dp": 2},
        "forcefield": {"name": "gaff2", "charge": "gasteiger"},
        "chains": 2,
        "density_g_cm3": 0.2,
        "md": {
            "nvt_steps": 30,
            "npt_steps": 30,
            "production_steps": 35,
            "sample_every": 20,
            "segment_steps": 20,
        },
    }


@pytest.mark.parametrize(
    "change",
    [
        {"chains": True},
        {"density_g_cm3": float("nan")},
        {"schema": 2},
        {"forcefield": {"name": "gaff2"}},
        {"forcefield": {"name": "opls", "profile": "pcff-carboxamide-v1"}},
        {
            "architecture": "crosslinked",
            "polymer": None,
            "chains": 1,
            "cure": {"dgeba": 4, "pacm": 2},
            "forcefield": {"name": "gaff2", "charge": "RESP"},
        },
        {"md": {"timestep_fs": 0}},
        {"polymer": {"repeat_units": {"A": "*CC*"}, "dp": 2, "sequence": ["A", "A"]}},
    ],
)
def test_rejects_ambiguous_or_unsupported_inputs(change):
    with pytest.raises((ValueError, TypeError)):
        PolymerRecipe.from_dict({**recipe(), **change})


def test_blocks_and_branches_use_existing_graph_contract():
    r = recipe()
    r.update(
        architecture="branched",
        polymer={
            "repeat_units": {"A": "*CC*", "B": "*CCO*"},
            "blocks": [["A", 2], ["B", 1]],
            "branches": [{"anchor_site": 2, "psmiles": "*CC"}],
        },
    )
    spec = PolymerRecipe.from_dict(r).polymer_spec()
    assert spec.sequence == ("A", "A", "B")
    assert spec.degree_of_polymerization == 3 and spec.branches[0].anchor_site == 2
    from island.forcefields.assignment_core.topology import molecule, normalize

    assert not normalize(spec).blockers
    assert Chem.GetMolFrags(molecule(normalize(spec)))


def test_checkpoint_failure_and_corruption(tmp_path):
    with pytest.raises(RuntimeError), checkpoint(tmp_path, "input") as (p, exists):
        assert not exists
        (p / "diagnostic.txt").write_text("retained")
        raise RuntimeError("interrupted")
    assert not (tmp_path / "input").exists()
    assert (
        next(tmp_path.glob("failed-input-*/diagnostic.txt")).read_text() == "retained"
    )
    with checkpoint(tmp_path, "input") as (p, exists):
        (p / "source").write_text("parameters")
    with checkpoint(tmp_path, "input") as (p, exists):
        assert exists
    (tmp_path / "input/source").write_text("modified")
    with (
        pytest.raises(ValueError, match="Checkpoint content changed"),
        checkpoint(tmp_path, "input"),
    ):
        pass


def native_engines():
    pytest.importorskip("radonpy")
    if not (exe := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS")) or not (
        packmol := os.environ.get("ISLAND_PACKMOL")
    ):
        pytest.skip("Set the LAMMPS and PACKMOL integration executable paths")
    return {"lammps": exe, "packmol": packmol}


def test_recipe_resume_uses_identical_restart_and_rejects_changed_input(tmp_path):
    engines = native_engines()
    r = recipe()
    direct = tmp_path / "direct"
    full = run_recipe(r, direct, **engines)
    resumed = tmp_path / "resumed"
    assert run_recipe(r, resumed, stop_after="nvt", **engines) == {
        "status": "paused",
        "stage": "nvt",
    }
    initial_hash = digest(resumed / "nvt/final.restart")
    result = run_recipe(r, resumed, **engines)
    assert full == result
    assert digest(resumed / "nvt/final.restart") == initial_hash
    assert (direct / "production/final.dump").read_bytes() == (
        resumed / "production/final.dump"
    ).read_bytes()
    before = (resumed / "production/final.dump").stat().st_mtime_ns
    assert run_recipe(r, resumed, **engines) == result
    assert (resumed / "production/final.dump").stat().st_mtime_ns == before
    assert result["stages"][-1]["frames"] == 3
    changed = deepcopy(r)
    changed["md"]["production_steps"] += 1
    with pytest.raises(ValueError, match="Recipe differs"):
        run_recipe(changed, resumed, **engines)
    (resumed / "production/final.restart").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Checkpoint content changed"):
        run_recipe(r, resumed, **engines)


def test_gaff_graph_preserves_site_ids_and_radonpy_native_energy(tmp_path):
    engines = native_engines()
    from radonpy.sim.lammps import MolToLAMMPSdataBlock

    from island.chemistry import from_rdkit
    from island.epoxy.dynamics import invoke, last_frame

    mol = Chem.AddHs(Chem.MolFromSmiles("CCOc1ccccc1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=41) == 0
    system = from_rdkit(mol, site_ids=range(100, 100 + mol.GetNumAtoms())).system
    assignment = assign_polymer_parameters(system, "gaff2", charge="gasteiger")
    assert assignment.site_ids == tuple(range(100, 100 + mol.GetNumAtoms()))
    data, params = radonpy_input(assignment, [30, 30, 30])
    native = "\n".join(MolToLAMMPSdataBlock(assignment.mol, velocity=False)) + "\n"
    results = []
    for name, source in [("island", data), ("native", native)]:
        out = tmp_path / name
        out.mkdir()
        (out / "system.data").write_text(source)
        styles = "\n".join(l for l in params.splitlines() if "_coeff " not in l)
        read = (
            "read_data system.data\n" + params
            if name == "island"
            else styles + "\nread_data system.data"
        )
        (out / "single.in").write_text(
            "units real\natom_style full\n"
            + read
            + '\nkspace_style ewald 1e-6\nrun 0\nprint "$(pe:%.16g) $(ebond:%.16g) $(eangle:%.16g) $(edihed:%.16g) $(eimp:%.16g)" file energy.txt\nwrite_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g\n'
        )
        invoke(engines["lammps"], ["-in", "single.in"], out, "screen.log")
        results.append(
            (np.loadtxt(out / "energy.txt"), last_frame(out / "final.dump")[-1])
        )
    np.testing.assert_allclose(results[0][0], results[1][0], atol=1e-10, rtol=0)
    np.testing.assert_allclose(results[0][1], results[1][1], atol=1e-10, rtol=0)


def test_missing_opls_terms_preserve_diagnostic_and_packed_input(tmp_path):
    engines = native_engines()
    pytest.importorskip("pyemc")
    r = recipe()
    r["polymer"]["repeat_units"]["A"] = "*CCP(=O)(O)C*"
    r["forcefield"] = {"name": "opls"}
    out = tmp_path / "phosphonate"
    with pytest.raises(ForceFieldError, match="Missing|missing"):
        run_recipe(r, out, **engines)
    status = json.loads((out / "status.json").read_text())
    assert status["status"] == "failed" and status["stage"] == "assignment"
    assert not (out / "result.json").exists() and not (out / "nvt").exists()
    assert (out / "packing/coordinates.npz").is_file()
    assert list(out.glob("failed-assignment-*"))
