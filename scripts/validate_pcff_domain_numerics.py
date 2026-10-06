"""J3 independent raw-row H2S/D2O oracle. No production resolver builds LAMMPS input."""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from validate_pcff_fallbacks import compare, fmt, raw_row
from validate_prepared_workflow import hashes

from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldSources,
    create_evaluator,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff.domains import PROFILE_NAME
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.workflows.storage import checksum, json_bytes, publish

CASES = {
    "hydrogen_sulfide": ("S", ["h", "s"], ["h", "s", "h"]),
    "heavy_water": ("[2H]O[2H]", ["h*", "o*"], ["h*", "o*", "h*"]),
}


def reference(exe, directory, system, xyz, raw, bondkey, anglekey, charges):
    directory.mkdir()
    ids = sorted(system.topology.sites)
    index = {i: n + 1 for n, i in enumerate(ids)}
    bonds = [tuple(b.key) for b in system.topology.bonds.values()]
    center = next(i for i in ids if sum(i in b for b in bonds) == 2)
    ends = [i for i in ids if i != center]
    assert len(ids) == 3 and len(bonds) == 2
    rows = {
        f: raw_row(raw, f, k)
        for f, k in [
            ("quartic_bond", bondkey),
            ("quartic_angle", anglekey),
            ("bond-bond", anglekey),
            ("bond-angle", anglekey),
        ]
    }
    b, a, bb, ba = [rows[f]["values"] for f in rows]
    assert len(b) == 4 and len(a) == 4 and len(bb) == 1 and len(ba) == 1
    data = "J3 independent raw-source triatomic reference\n\n3 atoms\n2 bonds\n1 angles\n\n3 atom types\n1 bond types\n1 angle types\n\n-50 50 xlo xhi\n-50 50 ylo yhi\n-50 50 zlo zhi\n\nMasses\n\n"
    data += "".join(f"{index[i]} {system.topology.sites[i].mass:.17g}\n" for i in ids)
    data += "\nAtoms # full\n\n" + "".join(
        f"{index[i]} 1 {index[i]} {charges[i]:.17g} " + fmt(xyz[n]) + "\n"
        for n, i in enumerate(ids)
    )
    data += "\nBonds\n\n" + "".join(
        f"{n + 1} 1 {index[i]} {index[j]}\n" for n, (i, j) in enumerate(bonds)
    )
    data += f"\nAngles\n\n1 1 {index[ends[0]]} {index[center]} {index[ends[1]]}\n"
    (directory / "data.lmp").write_text(data)
    script = f"""units real
atom_style full
boundary f f f
bond_style class2
angle_style class2
pair_style lj/class2/coul/cut 30
special_bonds lj 0 0 1 coul 0 0 1
read_data data.lmp
pair_coeff * * 0 1
pair_modify mix sixthpower shift no tail no
bond_coeff 1 {fmt(b)}
angle_coeff 1 {fmt(a)}
angle_coeff 1 bb {fmt(bb + [b[0], b[0]])}
angle_coeff 1 ba {fmt(ba + ba + [b[0], b[0]])}
thermo_style custom step pe ebond eangle evdwl ecoul
run 0
print "$(pe:%.17g) $(ebond:%.17g) $(eangle:%.17g) $(evdwl:%.17g) $(ecoul:%.17g)" file energies.txt
write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g
"""
    # Every pair is graph distance 1 or 2 and excluded. Zero pair coefficients
    # are reference-only: the separately checked native model retains source LJ.
    (directory / "in.lmp").write_text(script)
    cmd = [str(exe), "-in", "in.lmp", "-log", "log.lammps"]
    run = subprocess.run(
        cmd, cwd=directory, capture_output=True, text=True, timeout=180, check=False
    )
    (directory / "stdout.txt").write_text(run.stdout)
    (directory / "stderr.txt").write_text(run.stderr)
    if run.returncode:
        raise RuntimeError(run.stdout[-2000:] + run.stderr)
    e = (
        np.array(list(map(float, (directory / "energies.txt").read_text().split())))
        * 4.184
    )
    lines = (
        (directory / "forces.dump")
        .read_text()
        .split("ITEM: ATOMS id fx fy fz\n")[1]
        .splitlines()
    )
    f = np.array([list(map(float, l.split()[1:])) for l in lines]) * 4.184
    return (
        e,
        f,
        {
            "rows": rows,
            "command": cmd,
            "pair_inventory": {"distance1": 2, "distance2": 1, "remaining": 0},
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    root = args.output
    root.mkdir(parents=True, exist_ok=False)
    raw = args.source.read_bytes()
    declaration = json.loads(
        Path("docs/evidence/phase_4j3_declaration.json").read_text()
    )
    assert checksum(raw) == declaration["source"]["sha256"]
    publish(
        root / "declaration.json",
        json_bytes(
            {
                "parent": checksum(
                    Path("docs/evidence/phase_4j3_declaration.json").read_bytes()
                ),
                "cases": CASES,
                "source": declaration["source"],
                "lammps_sha256": checksum(args.lammps.read_bytes()),
                "tolerances": {"atol": 1e-5, "rtol": 2e-10},
                "pair_oracle": "all 3 pairs excluded by shortest paths; no active nonbonded interaction",
            }
        ),
    )
    outcomes = []
    for name, (smiles, bondkey, anglekey) in CASES.items():
        work = root / name
        work.mkdir()
        out = {"case": name, "passed": False}
        try:
            s = from_smiles(smiles, random_seed=2026)
            prepared = prepare_forcefield(
                s,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        source_path=args.source,
                        typing_profile=PROFILE_NAME,
                        resolution_policy=DOMAIN_POLICY,
                        lj=(0, 0, 1),
                        coulomb=(0, 0, 1),
                    ),
                ),
            )
            evaluator = create_evaluator(s, prepared)
            ids = sorted(s.topology.sites)
            xyz = np.array([s.coordinates.get(i) for i in ids])
            from validate_pcff_expanded import independent_charges

            from island.forcefields.pcff import load_pcff_source, type_pcff_atoms

            t = type_pcff_atoms(s, load_pcff_source(args.source), profile=PROFILE_NAME)
            q, qr = independent_charges(raw, s, t.assignments)
            out["charge_rows"] = qr
            out["numerical"] = []
            for n, x in enumerate(
                (
                    xyz,
                    xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3),
                )
            ):
                r = evaluator.evaluate(dict(zip(ids, x, strict=True)))
                e, f, meta = reference(
                    args.lammps.resolve(),
                    work / f"frame{n}",
                    s,
                    x,
                    raw,
                    bondkey,
                    anglekey,
                    q,
                )
                check = compare(
                    r.potential_energy, np.array([r.forces[i] for i in ids]), e[0], f
                )
                components = dict(r.energy_components)
                np.testing.assert_allclose(
                    components["quartic_bond"], e[1], atol=1e-5, rtol=2e-10
                )
                np.testing.assert_allclose(
                    sum(
                        components.get(k, 0)
                        for k in ["quartic_angle", "bond-bond", "bond-angle"]
                    ),
                    e[2],
                    atol=1e-5,
                    rtol=2e-10,
                )
                # Independently selected source rows and multiplicities, not just energies.
                from collections import Counter

                model_terms = prepared.native_result.payload["terms"]
                assert Counter(t["family"] for t in model_terms) == {
                    "quartic_bond": 2,
                    "quartic_angle": 1,
                    "bond-bond": 1,
                    "bond-angle": 1,
                }
                for term in model_terms:
                    family = term["family"]
                    assert {v["record_id"] for v in term["source_rows"]} == {
                        f"{family}:cff91:{meta['rows'][family]['line']}"
                    }
                center = next(
                    n for n, i in enumerate(ids) if s.topology.sites[i].element != "H"
                )
                ends = [n for n in range(3) if n != center]

                def raw_energy(z, meta=meta, center=center, ends=ends):
                    b, a, bb, ba = [
                        meta["rows"][k]["values"]
                        for k in [
                            "quartic_bond",
                            "quartic_angle",
                            "bond-bond",
                            "bond-angle",
                        ]
                    ]
                    vectors = [z[j] - z[center] for j in ends]
                    lengths = [np.linalg.norm(v) for v in vectors]
                    dr = np.array(lengths) - b[0]
                    theta = np.arccos(np.dot(*vectors) / np.prod(lengths))
                    dt = theta - np.deg2rad(a[0])
                    return 4.184 * (
                        sum(b[k] * sum(dr ** (k + 1)) for k in range(1, 4))
                        + sum(a[k] * dt ** (k + 1) for k in range(1, 4))
                        + bb[0] * np.prod(dr)
                        + ba[0] * sum(dr) * dt
                    )

                fd_errors = []
                for step in (1e-4, 1e-5, 1e-6):
                    fd = np.zeros_like(x)
                    for i in range(3):
                        for j in range(3):
                            hi = x.copy()
                            lo = x.copy()
                            hi[i, j] += step
                            lo[i, j] -= step
                            fd[i, j] = -(raw_energy(hi) - raw_energy(lo)) / (2 * step)
                    fd_errors.append(
                        float(np.max(np.abs(fd - np.array([r.forces[i] for i in ids]))))
                    )
                assert fd_errors[-1] < 1e-5
                out["numerical"].append(
                    {
                        **check,
                        **meta,
                        "components": components,
                        "finite_difference_max_errors": fd_errors,
                        "parameter_rows_and_inventory_match": True,
                    }
                )
            save_prepared_forcefield(
                s,
                prepared,
                work / "bundle",
                sources=PreparedForceFieldSources(pcff_frc=args.source),
            )
            out.update(
                passed=True,
                bundle_hashes=hashes(work / "bundle"),
                prepared_identity=prepared.identity,
                model_identity=prepared.native_result.identity,
            )
        except Exception as error:  # noqa: BLE001 -- retain declared failures
            out["error"] = type(error).__name__ + ": " + str(error)
        publish(work / "outcome.json", json_bytes(out))
        outcomes.append(out)
        print(name, out.get("error", "passed"), flush=True)
    publish(
        root / "report.json",
        json_bytes(
            {
                "cases": outcomes,
                "source_unchanged": checksum(args.source.read_bytes()) == checksum(raw),
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    return 0 if len(outcomes) == 2 and all(r["passed"] for r in outcomes) else 1


if __name__ == "__main__":
    raise SystemExit(main())
