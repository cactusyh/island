"""Export the same mapped epoxy additions to the LAMMPS REACTION package."""

import json
import re
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .chemistry import Addition, amines, epoxides, react, reference_reactions
from .cure import CureConfig, mixture
from .dynamics import invoke, last_frame
from .network import statistics
from .parameters import assign_network, data_text


def templates(directory, lib):
    out = Path(directory)
    for number, (pre, event, post) in enumerate(reference_reactions()[2], 1):
        changed = set(asdict(event).values())
        keep = set(changed)
        for _ in range(5):
            keep.update(
                n.GetIdx()
                for i in list(keep)
                for n in pre.GetAtomWithIdx(i).GetNeighbors()
            )
        mapping = {i: j for j, i in enumerate(sorted(keep), 1)}
        edges = [
            mapping[i]
            for i in sorted(keep)
            if any(n.GetIdx() not in keep for n in pre.GetAtomWithIdx(i).GetNeighbors())
        ]
        pre_atoms, _ = assign_network(pre, lib)
        post_atoms, _ = assign_network(post, lib)
        charged = [
            mapping[i]
            for i in range(len(pre_atoms))
            if pre_atoms[i][1] != post_atoms[i][1]
        ]
        for label, m in [("pre", pre), ("post", post)]:
            atoms, terms = assign_network(m, lib)
            selected = {
                k: [(t, ss) for t, ss in v if all(s in keep for s in ss)]
                for k, v in terms.items()
            }
            lines = [
                "ISLAND local epoxy addition",
                "",
                f"{len(keep)} atoms",
                "1 fragments",
            ]
            lines.extend(f"{len(v)} {k.lower()}" for k, v in selected.items())
            lines.extend(
                [
                    "",
                    "Fragments",
                    "",
                    "charge " + " ".join(map(str, charged)),
                    "",
                    "Types",
                    "",
                ]
            )
            lines.extend(f"{j} {atoms[i][0]}" for i, j in mapping.items())
            lines.extend(["", "Charges", ""])
            lines.extend(f"{j} {atoms[i][1]}" for i, j in mapping.items())
            lines.extend(["", "Coords", ""])
            lines.extend(
                f"{j} " + " ".join(map(str, m.GetConformer().GetAtomPosition(i)))
                for i, j in mapping.items()
            )
            for k, rows in selected.items():
                if not rows:
                    continue
                lines.extend(["", k, ""])
                for j, (t, ss) in enumerate(rows, 1):
                    lines.append(f"{j} {t} " + " ".join(str(mapping[s]) for s in ss))
            (out / f"{label}{number}.mol").write_text("\n".join(lines) + "\n")
        lines = [
            "ISLAND bijective atom mapping",
            "",
            f"{len(keep)} equivalences",
            f"{len(edges)} edgeIDs",
            "",
            "InitiatorIDs",
            "",
            str(mapping[event.carbon]),
            str(mapping[event.nitrogen]),
            "",
            "EdgeIDs",
            "",
            *map(str, edges),
            "",
            "Equivalences",
            "",
        ]
        lines.extend(f"{j} {j}" for j in mapping.values())
        (out / f"r{number}.map").write_text("\n".join(lines) + "\n")


def continuous(cure, library, out, *, lammps, steps=20000):
    """Start from the packed/equilibrated stage; compare final graph against both mapped reactions."""
    if type(steps) != int or steps < 1:
        raise ValueError("steps must be positive")
    cure, out = Path(cure), Path(out)
    cfg = CureConfig(**json.loads((cure / "config.json").read_text())["config"])
    lib = json.loads(Path(library).read_text())
    if lib["field"] != "pcff":
        raise ValueError(
            "OPLS requires batch cure with explicit hydrogen transfer and relaxation"
        )
    out.mkdir(parents=True, exist_ok=False)
    _, mol = mixture(cfg)
    with np.load(cure / "stage-0000/state.npz") as state:
        xyz, box, vel = [state[k] for k in ("xyz", "box", "vel")]
    atoms, terms = assign_network(mol, lib)
    owners = [a.GetIntProp("monomer") + 1 for a in mol.GetAtoms()]
    (out / "system.data").write_text(
        data_text(atoms, terms, lib["counts"], xyz, box, owners, vel)
    )
    (out / "system.params").write_text(lib["parameters"])
    templates(out, lib)
    script = f"""units real
atom_style full
boundary p p p
read_data system.data extra/bond/per/atom 4 extra/angle/per/atom 16 extra/dihedral/per/atom 48 extra/improper/per/atom 24 extra/special/per/atom 80
include system.params
kspace_style pppm 1e-5
neighbor 2 bin
neigh_modify every 1 delay 0 check yes
molecule p1 pre1.mol
molecule q1 post1.mol
molecule p2 pre2.mol
molecule q2 post2.mol
comm_modify cutoff 14
timestep 0.2
fix rx all bond/react stabilization yes bath 0.008 reset_mol_ids no &
 react first all 40 0 {min(cfg.cutoff, 4.0)} p1 q1 r1.map stabilize_steps 800 custom_charges charge &
 react second all 40 0 {min(cfg.cutoff, 4.0)} p2 q2 r2.map stabilize_steps 800 custom_charges charge
fix ensemble bath_REACT nvt temp {cfg.temperature} {cfg.temperature} 50
thermo 100
thermo_style custom step temp pe press f_rx[1] f_rx[2]
thermo_modify lost error
dump tr all custom 500 traj.dump id x y z vx vy vz fx fy fz
dump_modify tr sort id format float %.16g
run {steps}
write_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g
write_data product.data nocoeff nofix
print "$(f_rx[1]:%.0f) $(f_rx[2]:%.0f)" file counts.txt
"""
    (out / "react.in").write_text(script)
    invoke(lammps, ["-in", "react.in"], out, "screen.log")
    return verify_product(out, cfg, lib, lammps)


def verify_product(out, cfg, lib, lammps):
    _, mol = mixture(cfg)
    owners = [a.GetIntProp("monomer") + 1 for a in mol.GetAtoms()]
    data = (out / "product.data").read_text()
    bonds = re.search(r"(?m)^Bonds[^\n]*\n\n([^\n]+(?:\n[^\n]+)*)", data)
    adj = [set() for _ in mol.GetAtoms()]
    for line in bonds[1].splitlines():
        _, _, a, b = map(int, line.split()[:4])
        adj[a - 1].add(b - 1)
        adj[b - 1].add(a - 1)
    events = []
    for c, o in epoxides(mol).items():
        if o in adj[c]:
            continue
        (n,) = adj[c] & set(amines(mol))
        (h,) = adj[o] & set(amines(mol)[n])
        e = Addition(c, o, n, h)
        mol = react(mol, [e])
        events.append(asdict(e))
    if any(
        {n.GetIdx() for n in a.GetNeighbors()} != adj[a.GetIdx()]
        for a in mol.GetAtoms()
    ):
        raise ValueError(
            "REACTER final connectivity disagrees with atom-mapped additions"
        )
    _, xyz, box, _vel, _ = last_frame(out / "final.dump")
    primary, secondary = np.loadtxt(out / "counts.txt", dtype=int)
    stat = statistics(mol, xyz, box, cfg.dgeba, cfg.pacm)
    if primary != stat["secondary"] + stat["tertiary"] or secondary != stat["tertiary"]:
        raise ValueError("REACTER functional-group accounting mismatch")
    # Compare every final interaction including Class-II improper orientation,
    # charges, energies and forces with the independent topology update path.
    from .validation import evaluate

    atoms, terms = assign_network(mol, lib)
    typed = data_text(atoms, terms, lib["counts"], xyz, box, owners)
    _, _, f, e = evaluate(out / "mapped", typed, lib["parameters"], lammps)
    _, _, rf, re_ = evaluate(out / "native", data, lib["parameters"], lammps)
    report = {
        "events": events,
        "statistics": stat,
        "first": int(primary),
        "second": int(secondary),
        "max_force_error": float(abs(f - rf).max()),
        "max_energy_error": float(abs(e - re_).max()),
    }
    if report["max_force_error"] > 1e-4 or report["max_energy_error"] > 1e-4:
        raise ValueError(f"REACTER force/energy mismatch: {report}")
    (out / "result.json").write_text(json.dumps(report, indent=2))
    return report
