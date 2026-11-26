"""Source-bound full-graph GAFF parameterization for periodic reaction stages."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.forcefields.radonpy import assign_radonpy
from island.workflows.recipe_md import radonpy_input


def sources(family):
    import radonpy
    from radonpy.ff import gaff

    db = Path(gaff.__file__).parent / "ff_dat" / (family + ".json")
    return radonpy.__version__, {db.name: sha256(db.read_bytes()).hexdigest()}


def compile_library(directory, family, cfg, charge="gasteiger"):
    if family not in ("gaff", "gaff2") or charge != "gasteiger":
        raise ValueError(
            "Periodic GAFF networks require GAFF/GAFF2 and Gasteiger charges"
        )
    version, hashes = sources(family)
    lib = {
        "backend": "radonpy-graph",
        "field": family,
        "field_sha256": hashes,
        "engine_version": version,
        "charge_method": charge,
        "chemistry": cfg.template.name,
        "epoxy_smiles": cfg.epoxy_smiles,
        "hardener_smiles": cfg.hardener_smiles,
    }
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    (out / "library.json").write_text(json.dumps(lib, indent=2) + "\n")
    return lib


def parameterize(mol, lib, xyz, box, vel=None):
    version, hashes = sources(lib["field"])
    if hashes != lib["field_sha256"] or version != lib["engine_version"]:
        raise ValueError("GAFF source/version differs from the network library")
    if lib["charge_method"] != "gasteiger":
        raise ValueError("Unsupported periodic network charge method")
    xyz, box = np.asarray(xyz), np.asarray(box)
    if (
        xyz.shape != (mol.GetNumAtoms(), 3)
        or box.shape != (3,)
        or not np.isfinite(xyz).all()
        or not np.isfinite(box).all()
        or np.any(box <= 0)
    ):
        raise ValueError("Invalid periodic network coordinates/cell")
    graph = Chem.Mol(mol)
    Chem.SanitizeMol(graph)
    if any(
        a.GetAtomicNum() == 0 or a.GetTotalNumHs() or a.GetNumRadicalElectrons()
        for a in graph.GetAtoms()
    ):
        raise ValueError("Network must be closed-shell with explicit hydrogen sites")
    graph.RemoveAllConformers()
    conf = Chem.Conformer(len(xyz))
    for i, pos in enumerate(xyz):
        conf.SetAtomPosition(i, pos)
    graph.AddConformer(conf)
    # No finite-molecule unwrapping or embedding: graph typing is valid for
    # winding periodic cycles. Coordinates are exported in original site order.
    assignment = assign_radonpy(graph, lib["field"], charge="gasteiger")
    charges = [a.GetDoubleProp("AtomicCharge") for a in assignment.mol.GetAtoms()]
    for ids in Chem.GetMolFrags(graph):
        formal = sum(graph.GetAtomWithIdx(i).GetFormalCharge() for i in ids)
        if abs(sum(charges[i] for i in ids) - formal) > 1e-6:
            raise ValueError("Network component charge is not conserved")
    data, params = radonpy_input(assignment, box)
    if vel is not None:
        vel = np.asarray(vel)
        if vel.shape != xyz.shape or not np.isfinite(vel).all():
            raise ValueError("Invalid network velocities")
        data += (
            "\nVelocities\n\n"
            + "\n".join(
                f"{i + 1} " + " ".join(format(v, ".16g") for v in row)
                for i, row in enumerate(vel)
            )
            + "\n"
        )
    report = {
        "backend": "radonpy-graph",
        "field_sha256": hashes,
        "engine_version": version,
        "charge_method": "gasteiger",
        "atoms": len(charges),
        "charge": sum(charges),
        "atom_types": [a.GetProp("ff_type") for a in assignment.mol.GetAtoms()],
        "charges": charges,
        "lj14_scale": assignment.lj14_scale,
        "coulomb14_scale": assignment.coulomb14_scale,
        "graph_smiles": Chem.MolToSmiles(graph),
        "site_ids": list(assignment.site_ids),
    }
    return data, params, report
