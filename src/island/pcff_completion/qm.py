"""Explicit QM input adapters and graph-bound deformation templates; no runner."""

import shutil
from itertools import pairwise
from pathlib import Path

import numpy as np

from . import _finite, _verified, digest, inventory, pin, require, signed, write
from .targets import configure

ADAPTERS = {
    "psi4": "psi4",
    "pyscf": "pyscf",
    "orca": "orca",
    "cp2k": "cp2k.psmp",
    "gaussian": "g16",
}


def asset_status(engine, executable=None):
    require(engine in ADAPTERS, "unsupported_qm_adapter")
    path = shutil.which(str(executable or ADAPTERS[engine]))
    return {
        "adapter": engine,
        "status": "available_not_executed" if path else "asset_unavailable",
        "executable": path,
        "executable_sha256": pin(path) if path else None,
        "QM_executed": False,
    }


def deform(system, assessment, scan):
    """Rigid fragment bond/angle/torsion deformations; exact source atom order.

    Values are displacements from the input geometry, never optimized results.
    A bond cut must separate the exact final graph; cyclic scans are rejected.
    """
    p = _verified(system, assessment)
    kind, sites = scan.get("kind"), scan.get("sites", [])
    count = {"bond": 2, "angle": 3, "torsion": 4}.get(kind)
    require(
        count and len(sites) == count and len(set(sites)) == count,
        "scan_sites_required",
    )
    units = "angstrom" if kind == "bond" else "radian"
    require(scan.get("units") == units, "scan_units_required")
    values = _finite(scan.get("displacements", []))
    require(
        values.ndim == 1
        and values.size > 0
        and len(set(values.tolist())) == len(values),
        "unique_scan_displacements_required",
    )
    adjacency = {i: set() for i in p["atom_order"]}
    for b in system.topology.bonds.values():
        adjacency[b.site1].add(b.site2)
        adjacency[b.site2].add(b.site1)
    require(
        all(i in adjacency and j in adjacency[i] for i, j in pairwise(sites)),
        "scan_path_not_in_graph",
    )
    left, right = (sites[0], sites[1]) if kind == "bond" else ((sites[1], sites[2]))
    adjacency[left].remove(right)
    adjacency[right].remove(left)
    group, todo = set(), [right]
    while todo:
        i = todo.pop()
        if i not in group:
            group.add(i)
            todo.extend(adjacency[i] - group)
    require(left not in group, "cyclic_scan_domain_forbidden")
    order = {i: j for j, i in enumerate(p["atom_order"])}
    xyz = np.array([system.coordinates.get(i).tolist() for i in p["atom_order"]])
    origin = xyz[order[right if kind == "torsion" else left]]
    axis = xyz[order[right]] - xyz[order[left]]
    if kind == "angle":
        axis = np.cross(xyz[order[sites[0]]] - origin, axis)
    require(np.linalg.norm(axis) > 1e-12, "degenerate_scan_axis")
    axis = axis / np.linalg.norm(axis)
    selected = [order[i] for i in sorted(group)]
    result = []
    for delta in values:
        x = xyz.copy()
        if kind == "bond":
            require(
                np.linalg.norm(xyz[order[right]] - xyz[order[left]]) + delta > 0,
                "nonpositive_deformed_bond",
            )
            x[selected] += delta * axis
        else:
            v = x[selected] - origin
            x[selected] = origin + (
                v * np.cos(delta)
                + np.cross(axis, v) * np.sin(delta)
                + np.outer(v @ axis, axis) * (1 - np.cos(delta))
            )
        result.append(
            {
                "displacement": float(delta),
                "units": units,
                "atom_order": p["atom_order"],
                "coordinates": x.tolist(),
                "geometry_checksum": digest(x.tolist()),
            }
        )
    return result


def templates(
    system, assessment, target, settings, scan, output, *, engine, executable=None
):
    """Emit XYZ and visibly unfilled engine decks. Never execute a QM program.

    Engine-specific method/dispersion/convergence syntax is deliberately left
    visible for expert completion. The exact supplied plan is retained as JSON;
    no template is claimed executable or a reference output.
    """
    p = _verified(system, assessment)
    t = target.validate(system, assessment)
    plan = configure(target, settings)
    require(engine in ADAPTERS, "unsupported_qm_adapter")
    out = Path(output)
    require(not out.exists(), "output_directory_already_exists")
    geometries = deform(system, assessment, scan)
    out.mkdir(parents=True)
    elements = [system.topology.sites[i].element for i in p["atom_order"]]
    for i, g in enumerate(geometries):
        body = "\n".join(
            e + " " + " ".join(format(v, ".17g") for v in x)
            for e, x in zip(elements, g["coordinates"], strict=True)
        )
        (out / f"geometry-{i:04d}.xyz").write_text(
            f"{len(elements)}\nnot_QM_output\n{body}\n"
        )
        decks = {
            "psi4": f"molecule mol {{\n<charge> <multiplicity>\n{body}\nunits angstrom\n}}\n"
            "set basis <basis>\n<explicit dispersion/convergence>\nenergy('<method>')\n",
            "pyscf": "from pyscf import gto, scf\n"
            f"mol = gto.M(atom={body!r}, unit='Angstrom', basis='<basis>', "
            "charge=<charge>, spin=<2S>)\n<explicit method/dispersion/convergence>\n",
            "orca": f"! <method> <basis> <dispersion> <convergence>\n* xyz <charge> <multiplicity>\n{body}\n*\n",
            "cp2k": "&FORCE_EVAL\n <explicit method/potential/basis/dispersion/convergence>\n"
            f" &SUBSYS\n  &COORD\n{body}\n  &END COORD\n &END SUBSYS\n&END FORCE_EVAL\n",
            "gaussian": f"# <method>/<basis> <dispersion> <convergence>\n\nTemplate only\n\n<charge> <multiplicity>\n{body}\n\n",
        }
        (out / f"input-{i:04d}.{engine}.template").write_text(decks[engine])
    write(out / "plan.json", plan)
    write(out / "geometries.json", geometries)
    record = signed(
        {
            "schema": "island.pcff-qm-input-templates.v1",
            "target_identity": t["identity"],
            "graph_identity": p["graph_identity"],
            "atom_order": p["atom_order"],
            "plan_identity": plan["identity"],
            "scan": scan,
            "files": inventory(out),
            "adapter": asset_status(engine, executable),
            "status": "unfilled_templates_only",
            "reference_data_status": "not_supplied",
            "QM_executed": False,
            "parameter_completeness": "incomplete",
            "execution_authority": False,
        }
    )
    return record
