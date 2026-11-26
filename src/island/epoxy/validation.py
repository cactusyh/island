"""Independent whole-molecule vs local PCFF evaluation with the same coordinates."""

import json
import re
from pathlib import Path

import numpy as np

from island.forcefields.emc_graph import _section

from .chemistry import TEMPLATES, reference_models, reference_reactions
from .dynamics import invoke, last_frame
from .parameters import assign_network, data_text, rows


def evaluate(
    out,
    data,
    params,
    lammps,
    minimize=False,
    *,
    direct_ewald=False,
    ewald_accuracy=1e-7,
    coulomb_table=True,
):
    if not np.isfinite(ewald_accuracy) or not 0 < ewald_accuracy < 1:
        raise ValueError("Ewald accuracy must be finite and in (0,1)")
    out.mkdir(parents=True, exist_ok=False)
    (out / "system.data").write_text(data)
    (out / "system.params").write_text(params)
    script = """units real
atom_style full
boundary p p p
read_data system.data
include system.params
kspace_style pppm 1e-6
thermo_modify norm no
"""
    if direct_ewald:
        script = script.replace(
            "kspace_style pppm 1e-6", f"kspace_style ewald {ewald_accuracy:.16g}"
        )
    if not coulomb_table:
        script += "pair_modify table 0\n"
    if minimize:
        script += "minimize 1e-8 1e-6 2000 10000\n"
    script += 'run 0\nprint "$(ebond:%.15g) $(eangle:%.15g) $(edihed:%.15g) $(eimp:%.15g) $(evdwl:%.15g) $(ecoul:%.15g) $(elong:%.15g)" file energies.txt\n'
    script += "write_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g\n"
    (out / "check.in").write_text(script)
    invoke(lammps, ["-in", "check.in"], out, "screen.log")
    _, xyz, box, _, force = last_frame(out / "final.dump")
    return xyz, box, force, np.loadtxt(out / "energies.txt")


def validate_library(directory, out, lammps, seed=2026):
    from rdkit import Chem

    directory, out = Path(directory), Path(out)
    lib = json.loads((directory / "library.json").read_text())
    template = TEMPLATES[lib.get("chemistry", "pacm")]
    report = []
    for name, m in reference_models(seed, template):
        native = json.loads((directory / (name + ".json")).read_text())
        data = native["data"]
        rawrows = sorted(rows(data, "Atoms"), key=lambda r: int(r[0]))
        source = np.array(
            [native["site_groups"][str(i + 1)][0] for i in range(len(rawrows))]
        )
        inverse = np.argsort(source)
        coords, box, _, _ = evaluate(
            out / name / "relax", data, native["parameters"], lammps, True
        )
        for frame, dx in [
            ("relaxed", 0),
            ("distorted", 0.04),
            ("permuted-distorted", 0.04),
        ]:
            order = np.arange(m.GetNumAtoms())
            if frame == "permuted-distorted":
                order = np.random.default_rng(seed).permutation(order)
            atom, terms = assign_network(Chem.RenumberAtoms(m, order.tolist()), lib)
            xyz = coords + dx * np.random.default_rng(seed).normal(size=coords.shape)
            for i, row in enumerate(rawrows):
                row[4:7] = [str(x) for x in xyz[i]]
            native_data = data.replace(
                _section(data, "Atoms"), "\n".join(" ".join(r) for r in rawrows)
            )
            for axis, side in zip("xyz", box):
                native_data = re.sub(
                    rf"(?m)^.* {axis}lo {axis}hi$",
                    f"0 {side} {axis}lo {axis}hi",
                    native_data,
                )
            _, _, force, energy = evaluate(
                out / name / (frame + "-native"),
                native_data,
                native["parameters"],
                lammps,
            )
            mapped = data_text(
                atom, terms, lib["counts"], xyz[inverse][order], box, [1] * len(xyz)
            )
            _, _, actual, values = evaluate(
                out / name / (frame + "-mapped"), mapped, lib["parameters"], lammps
            )
            diffs = np.abs(energy - values)
            err = float(np.abs(force[inverse][order] - actual).max())
            report.append(
                {
                    "molecule": name,
                    "frame": frame,
                    "energy_component_errors": diffs.tolist(),
                    "max_force_error": err,
                }
            )
            if max(diffs) > 1e-6 or err > 1e-5:
                raise ValueError(f"Native parity failure: {report[-1]}")
    (out / "validation.json").write_text(json.dumps(report, indent=2))
    return report


def opls_coverage(out, seed=2026):
    """Record release-library diagnostics and atom-mapped structures for each stage."""
    from dataclasses import asdict

    from rdkit import Chem

    from island.exceptions import ForceFieldError
    from island.forcefields._emc_graph_input import emc_groups
    from island.forcefields.assignment import assign_emc_graph
    from island.forcefields.emc import _field_sources

    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    d, p, rx = reference_reactions(seed)
    field = "opls/2024/opls-aa"
    _, _, _, hashes = _field_sources(field, None, None)
    report = {"field": field, "field_sha256": hashes, "cases": {}}
    for name, m in [
        ("dgeba", d),
        ("pacm", p),
        ("secondary", rx[0][2]),
        ("tertiary", rx[1][2]),
    ]:
        mol = Chem.Mol(m)
        mol.RemoveAllConformers()
        row = {}
        try:
            a = assign_emc_graph(mol, field, seed=seed)
        except ForceFieldError as exc:
            diagnostic = [
                l.strip()
                for l in str(exc).splitlines()
                if l.startswith(("Warning:", "Error:"))
            ]
            if not diagnostic:
                raise
            row.update(assigned=False, diagnostic=diagnostic)
            groups, order = emc_groups(mol, False, False)
            row["emc_groups"] = groups
            row["emc_traversal_to_source_zero_based"] = order
            missing = []
            for message in diagnostic:
                match = re.search(r"\{group, site\} = \{(\w+), (\d+)\}", message)
                if match:
                    site = order[match[1]][int(match[2])]
                    atom = mol.GetAtomWithIdx(site)
                    missing.append(
                        {
                            "emc_group": match[1],
                            "emc_site": int(match[2]),
                            "source_site": site,
                            "element": atom.GetSymbol(),
                            "neighbors": [a.GetIdx() for a in atom.GetNeighbors()],
                        }
                    )
            row["untyped_sites"] = missing
        else:
            row.update(
                assigned=True, emc_version=a.emc_version, particles=len(a.site_groups)
            )
        if name in ("secondary", "tertiary"):
            row["reaction_zero_based"] = asdict(rx[0 if name == "secondary" else 1][1])
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        row["mapped_product_smiles"] = Chem.MolToSmiles(mol)
        report["cases"][name] = row
    (out / "coverage.json").write_text(json.dumps(report, indent=2))
    return report


def validate_periodic_stage(cure, directory, out, lammps):
    """Compare a winding network before/after translation and atomwise wrapping."""
    from .chemistry import react
    from .cure import bound_cure_inputs, mixture
    from .network import select, statistics

    cure, out = Path(cure), Path(out)
    cfg, lib = bound_cure_inputs(cure, Path(directory) / "library.json")
    _, mol = mixture(cfg)
    records = sorted(cure.glob("stage-*/record.json"))
    for path in records:
        rec = json.loads(path.read_text())
        mol = react(
            mol,
            [cfg.template.event(e) for e in rec["events"]],
            cfg.template,
        )
        if rec["statistics"]["winding_rank"]:
            with np.load(path.parent / "state.npz") as state:
                xyz, box = state["xyz"], state["box"]
            if select(
                mol,
                xyz,
                box,
                cfg.cutoff,
                np.random.default_rng(cfg.seed),
                cfg.batch_size,
                cfg.min_ring,
                template=cfg.template,
                donor_angle_tolerance=cfg.donor_angle_tolerance,
            )[0]:
                break
    else:
        raise ValueError(
            "A winding network with an eligible reaction contact is required for this test"
        )
    with np.load(path.parent / "state.npz") as state:
        xyz, box = state["xyz"], state["box"]
    atoms, terms = assign_network(mol, lib)
    owners = [a.GetIntProp("monomer") + 1 for a in mol.GetAtoms()]
    results = []
    for name, pos in [
        ("original", xyz),
        ("translated", (xyz + box * np.array([0.294, 0.482, 0.717])) % box),
    ]:
        data = data_text(atoms, terms, lib["counts"], pos, box, owners)
        _, _, forces, energies = evaluate(
            out / name, data, lib["parameters"], lammps, direct_ewald=True
        )
        events, info = select(
            mol,
            pos,
            box,
            cfg.cutoff,
            np.random.default_rng(cfg.seed),
            cfg.batch_size,
            cfg.min_ring,
            template=cfg.template,
            donor_angle_tolerance=cfg.donor_angle_tolerance,
        )
        results.append(
            (
                forces,
                energies,
                statistics(mol, pos, box, cfg.dgeba, cfg.hardener_count, cfg.template)[
                    "winding_rank"
                ],
                info["candidate_count"],
                events,
            )
        )
    f, e, rank, n, events = results[0]
    other = results[1]
    ef = float(np.abs(f - other[0]).max())
    ee = float(np.abs(e - other[1]).max())
    if (rank, n, events) != other[2:] or ef > 1e-5 or ee > 1e-6:
        raise ValueError(f"Periodic invariance failure: {ef}, {ee}")
    result = {
        "stage": rec["index"],
        "winding_rank": rank,
        "candidate_count": n,
        "max_force_error": ef,
        "max_energy_error": ee,
    }
    (out / "result.json").write_text(json.dumps(result, indent=2))
    return result
