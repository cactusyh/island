"""Build periodic polymer systems and run native MD from one JSON recipe."""

import argparse
import fcntl
import json
import os
import shutil
from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from uuid import uuid4

import numpy as np

from .recipe_config import PolymerRecipe


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def digest(path):
    with path.open("rb") as handle:
        from hashlib import file_digest

        return file_digest(handle, "sha256").hexdigest()


@contextmanager
def checkpoint(root, name):
    """Commit complete stages atomically and verify every retained artifact on reuse."""
    dest = root / name
    if dest.exists():
        files = json.loads((dest / "checkpoint.json").read_text())
        for relative, expected in files.items():
            if digest(dest / relative) != expected:
                raise ValueError(f"Checkpoint content changed: {name}/{relative}")
        yield dest, True
        return
    with TemporaryDirectory(dir=root, prefix="." + name + "-") as scratch:
        work = Path(scratch) / name
        work.mkdir()
        try:
            yield work, False
            files = {
                p.relative_to(work).as_posix(): digest(p)
                for p in sorted(work.rglob("*"))
                if p.is_file()
            }
            write_json(work / "checkpoint.json", files)
            work.rename(dest)
        except Exception:
            work.rename(root / f"failed-{name}-{uuid4().hex[:8]}")
            raise


def _build_chain(recipe, out):
    from rdkit import Chem

    from island.forcefields._emc_graph_input import graph_input
    from island.forcefields.assignment_core.topology import normalize

    graph = normalize(recipe.polymer_spec())
    source = graph_input(graph, recipe.seed)
    conf = Chem.Conformer(source.mol.GetNumAtoms())
    for i, position in enumerate(source.xyz):
        conf.SetAtomPosition(i, position)
    source.mol.RemoveAllConformers()
    source.mol.AddConformer(conf)
    write_json(out / "chain.json", asdict(graph))
    return source.mol


def _pack_chains(recipe, out, packmol):
    from rdkit import Chem

    from island.epoxy.dynamics import pack

    chain = _build_chain(recipe, out)
    xyz, box = pack(
        [chain], [recipe.chains], recipe.density_g_cm3, recipe.seed, packmol, out
    )
    mol = Chem.Mol(chain)
    for _ in range(recipe.chains - 1):
        mol = Chem.CombineMols(mol, chain)
    mol.RemoveAllConformers()
    (out / "molecule.json").write_text(Chem.MolsToJSON([mol]))
    np.savez(out / "coordinates.npz", xyz=xyz, box=box)
    write_json(
        out / "site-map.json",
        [
            {
                "site_id": i + 1,
                "chain": i // chain.GetNumAtoms() + 1,
                "chain_site_id": i % chain.GetNumAtoms() + 1,
            }
            for i in range(len(xyz))
        ],
    )


def _assign(recipe, packed, out, field_dir):
    from rdkit import Chem

    from island.chemistry import from_rdkit
    from island.core import SimulationBox
    from island.forcefields.assignment import assign_polymer_parameters
    from island.forcefields.radonpy import FAMILIES

    from .recipe_md import radonpy_input

    (mol,) = Chem.JSONToMols((packed / "molecule.json").read_text())
    with np.load(packed / "coordinates.npz") as state:
        xyz, box = state["xyz"], state["box"]
    conf = Chem.Conformer(len(xyz))
    for i, position in enumerate(xyz):
        conf.SetAtomPosition(i, position)
    mol.AddConformer(conf)
    system = from_rdkit(mol).system
    system.box = SimulationBox(*box)
    ff = recipe.forcefield
    if ff.name in FAMILIES:
        assignment = assign_polymer_parameters(system, ff.name, charge=ff.charge)
        data, params = radonpy_input(assignment, box)
        charges = [a.GetDoubleProp("AtomicCharge") for a in assignment.mol.GetAtoms()]
        for ids in Chem.GetMolFrags(assignment.mol):
            if (
                abs(
                    sum(charges[i] for i in ids)
                    - sum(
                        assignment.mol.GetAtomWithIdx(i).GetFormalCharge() for i in ids
                    )
                )
                > 1e-4
            ):
                raise ValueError("Assigned charges disagree with a chain formal charge")
        report = {
            "field": assignment.family,
            "engine": assignment.version,
            "field_sha256": {ff.name + ".json": assignment.database_sha256},
            "charge_method": assignment.charge_method,
            "site_groups": {i + 1: [s] for i, s in enumerate(assignment.site_ids)},
            "lj14_scale": assignment.lj14_scale,
            "coulomb14_scale": assignment.coulomb14_scale,
            "atom_types": [a.GetProp("ff_type") for a in assignment.mol.GetAtoms()],
            "charges": charges,
        }
    else:
        assignment = assign_polymer_parameters(
            system, ff.name, seed=recipe.seed, field_dir=field_dir
        )
        data, params = assignment.data, assignment.parameters
        report = asdict(assignment)
        report.pop("data")
        report.pop("parameters")
        report["charge_method"] = "native field increments"
    (out / "system.data").write_text(data)
    (out / "system.params").write_text(params)
    write_json(out / "assignment.json", report)
    shutil.copyfile(packed / "site-map.json", out / "site-map.json")


def _cure_input(root, recipe, field_dir, lammps, packmol):
    from island.epoxy.cure import run_cure
    from island.epoxy.parameters import compile_library
    from island.epoxy.properties import _snapshot
    from island.forcefields.emc import FIELDS

    cfg = recipe.cure_config()
    with checkpoint(root, "reaction-library") as (work, exists):
        if not exists:
            compile_library(
                work / "compiled",
                field=FIELDS[recipe.forcefield.name],
                field_dir=field_dir,
                chemistry="petmp" if cfg.petmp else "pacm",
            )
    library = root / "reaction-library/compiled/library.json"
    library_manifest = json.loads(library.read_text())
    if library_manifest["field_sha256"] != json.loads(
        (root / "field/hashes.json").read_text()
    ):
        raise ValueError("Reaction library differs from the selected field")
    records = run_cure(root / "cure", library, cfg, lammps=lammps, packmol=packmol)
    final = root / "cure" / f"stage-{records[-1]['index']:04d}"
    with checkpoint(root, "assignment") as (work, exists):
        if not exists:
            with np.load(final / "state.npz") as state:
                data = _snapshot(final, state["xyz"], state["box"], state["vel"])
            (work / "system.data").write_text(data)
            shutil.copyfile(final / "system.params", work / "system.params")
            lib = json.loads(library.read_text())
            write_json(
                work / "assignment.json",
                {
                    "field": lib["field"],
                    "field_sha256": lib["field_sha256"],
                    "engine": lib["emc_version"],
                    "library_sha256": digest(library),
                    "charge_method": "native field increments",
                    "cure_endpoint": records[-1],
                },
            )
    return records[-1]


def run_recipe(recipe, directory, *, lammps="lmp", packmol="packmol", stop_after=None):
    """Continue a matching recipe from its last completed, checksummed stage."""
    recipe = PolymerRecipe.from_dict(
        recipe.payload() if isinstance(recipe, PolymerRecipe) else recipe
    )
    if stop_after not in (None, "packing", "assignment", "nvt", "npt"):
        raise ValueError("stop_after must be packing, assignment, nvt or npt")
    if recipe.architecture == "crosslinked" and stop_after == "packing":
        raise ValueError("Crosslinked recipes checkpoint packing within cure")
    out = Path(directory).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _run(recipe, out, lammps, packmol, stop_after)


def _run(recipe, out, lammps, packmol, stop_after):
    from island.forcefields.emc import FIELDS, _field_sources
    from island.forcefields.emc_profiles import install
    from island.forcefields.radonpy import FAMILIES

    from .recipe_md import run_segment

    ff = recipe.forcefield
    field_dir = Path(ff.field_dir).resolve() if ff.field_dir else None
    payload = recipe.payload()
    if field_dir:
        payload["forcefield"]["field_dir"] = str(field_dir)
    manifest = out / "recipe.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != payload:
            raise ValueError("Recipe differs from the saved run")
    else:
        if any(p.name != ".lock" for p in out.iterdir()):
            raise FileExistsError("A new recipe requires an empty output directory")
        write_json(manifest, payload)
    phase = "runtime"
    try:
        executables = {}
        for label, command in (("lammps", lammps), ("packmol", packmol)):
            path = shutil.which(os.fspath(command))
            if path is None:
                raise FileNotFoundError(f"{label} executable not found: {command}")
            executable = Path(path).resolve()
            executables[label] = {"path": str(executable), "sha256": digest(executable)}
        runtime = out / "runtime.json"
        if runtime.exists():
            saved = json.loads(runtime.read_text())
            if any(saved[k]["sha256"] != v["sha256"] for k, v in executables.items()):
                raise ValueError("Simulation executable differs from the saved run")
        else:
            write_json(runtime, executables)
        lammps, packmol = (executables[k]["path"] for k in ("lammps", "packmol"))
        phase = "field"
        write_json(out / "status.json", {"status": "running", "stage": phase})
        with checkpoint(out, "field") as (work, exists):
            if ff.name in FAMILIES:
                from radonpy.ff import gaff

                db = Path(gaff.__file__).parent / "ff_dat" / (ff.name + ".json")
                hashes = {db.name: digest(db)}
            elif ff.profile:
                if not exists:
                    install(ff.profile, work / "installed")
                hashes = json.loads((work / "installed/profile.json").read_text())[
                    "installed"
                ]
            else:
                hashes = _field_sources(FIELDS[ff.name], None, field_dir)[3]
            if exists:
                if json.loads((work / "hashes.json").read_text()) != hashes:
                    raise ValueError("Force-field source changed since this run began")
            else:
                write_json(work / "hashes.json", hashes)
        if ff.profile:
            field_dir = out / "field/installed"
        endpoint = None
        if recipe.architecture == "crosslinked":
            phase = "cure"
            write_json(out / "status.json", {"status": "running", "stage": phase})
            endpoint = _cure_input(out, recipe, field_dir, lammps, packmol)
        else:
            phase = "packing"
            write_json(out / "status.json", {"status": "running", "stage": phase})
            with checkpoint(out, phase) as (work, exists):
                if not exists:
                    _pack_chains(recipe, work, packmol)
            if stop_after == phase:
                return _paused(out, phase)
            phase = "assignment"
            write_json(out / "status.json", {"status": "running", "stage": phase})
            with checkpoint(out, phase) as (work, exists):
                if not exists:
                    _assign(recipe, out / "packing", work, field_dir)
        assignment_report = json.loads((out / "assignment/assignment.json").read_text())
        if endpoint is not None:
            library = out / "reaction-library/compiled/library.json"
            if assignment_report["library_sha256"] != digest(library):
                raise ValueError("Cured assignment differs from its reaction library")
            assignment_report["field_sha256"] = json.loads(library.read_text())[
                "field_sha256"
            ]
        if assignment_report["field_sha256"] != hashes:
            raise ValueError(
                "Assignment parameter source differs from the recipe checkpoint"
            )
        if stop_after == "assignment":
            return _paused(out, "assignment")
        source = out / "assignment"
        outputs = []
        stages = [
            ("nvt", recipe.md.nvt_steps, "nvt"),
            ("npt", recipe.md.npt_steps, "npt"),
            ("production", recipe.md.production_steps, recipe.md.production_ensemble),
        ]
        for phase, steps, ensemble in stages:
            write_json(out / "status.json", {"status": "running", "stage": phase})
            with checkpoint(out, phase) as (work, exists):
                if not exists:
                    run_segment(
                        source,
                        work,
                        recipe.md,
                        steps,
                        ensemble,
                        recipe.seed,
                        lammps,
                        restart=phase != "nvt",
                        velocities=endpoint is not None,
                    )
            source = out / phase
            outputs.append(json.loads((source / "result.json").read_text()))
            if stop_after == phase:
                return _paused(out, phase)
        result = {
            "status": "completed",
            "architecture": recipe.architecture,
            "recipe_sha256": digest(manifest),
            "assignment": assignment_report,
            "stages": outputs,
            "total_md_time_fs": sum(r["time_fs"] for r in outputs),
            "equilibrated": "unmeasured",
            "checkpoint": "production/final.restart",
        }
        if endpoint is not None:
            result["cure"] = endpoint["statistics"]
            result["cure_termination"] = endpoint["reason"]
            result["target_reached"] = (
                endpoint["statistics"]["epoxy_conversion"]
                >= recipe.cure_config().target
            )
        write_json(out / "result.json", result)
        write_json(out / "status.json", {"status": "completed", "stage": "production"})
        return result
    except Exception as exc:
        (out / "result.json").unlink(missing_ok=True)
        write_json(
            out / "status.json",
            {
                "status": "failed",
                "stage": phase,
                "error_type": type(exc).__name__,
                "diagnostic": str(exc),
            },
        )
        raise


def _paused(out, phase):
    result = {"status": "paused", "stage": phase}
    write_json(out / "status.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recipe", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--lammps", default="lmp")
    parser.add_argument("--packmol", default="packmol")
    parser.add_argument("--stop-after", choices=("packing", "assignment", "nvt", "npt"))
    args = parser.parse_args()
    payload = json.loads(args.recipe.read_text())
    recipe = PolymerRecipe.from_dict(payload)
    if recipe.forcefield.field_dir:
        recipe = replace(
            recipe,
            forcefield=replace(
                recipe.forcefield,
                field_dir=str(
                    (
                        args.recipe.resolve().parent / recipe.forcefield.field_dir
                    ).resolve()
                ),
            ),
        )
    result = run_recipe(
        recipe,
        args.directory,
        lammps=args.lammps,
        packmol=args.packmol,
        stop_after=args.stop_after,
    )
    print(
        json.dumps(
            {"status": result["status"], "directory": str(args.directory.resolve())},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
