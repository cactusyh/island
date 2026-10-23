"""Read actual LAMMPS text independently of its exporter; no dynamics."""

from math import isclose

from island.analysis.equilibration import require
from island.analysis.npt import FLAGS, _signed


def audit(payload, coordinates_angstrom, directory):
    sections = {}
    current = None
    box = []
    for line in (directory / "system.data").read_text().splitlines()[1:]:
        line = line.split("#")[0].strip()
        if not line:
            continue
        if line.endswith(("xlo xhi", "ylo yhi", "zlo zhi")):
            box.append(list(map(float, line.split()[:2])))
        elif line[0].isalpha():
            current = line
            sections[current] = []
        elif current is not None:
            sections[current].append([float(x) for x in line.split()])
    data, config = payload["numerical_data"], payload["config"]
    require(box == [[0, x] for x in config["box_lengths"]], "Changed box")
    masses = {int(r[0]): r[1] for r in sections["Masses"]}
    lj = {int(r[0]): r[1:] for r in sections["Pair Coeffs"]}
    sites = {s["id"]: s for s in data["sites"]}
    atoms = sections["Atoms"]
    require([int(r[0]) for r in atoms] == sorted(sites), "Changed atom order")
    near = lambda a, b: isclose(
        a, b, rel_tol=2e-15, abs_tol=0
    )  # 17-digit unit conversion only
    for row in atoms:
        sid, _, typ = map(int, row[:3])
        site = sites[sid]
        require(row[3] == site["charge"], "Changed charge")
        require(row[4:] == list(coordinates_angstrom[sid]), "Changed coordinate")
        require(masses[typ] == site["mass"], "Changed mass")
        epsilon, sigma = lj[typ]
        require(
            near(epsilon * 4.184, site["epsilon"]) and sigma == site["sigma"],
            "Changed LJ parameters",
        )
    for section, key in [
        ("Bonds", "bonds"),
        ("Angles", "angles"),
        ("Dihedrals", "torsions"),
        ("Impropers", "impropers"),
    ]:
        require(
            [list(map(int, r[2:])) for r in sections.get(section, [])]
            == [r["sites"] for r in data[key]],
            "Changed " + section,
        )
    commands = (directory / "in.periodic").read_text().splitlines()
    require(
        "units real" in commands and "boundary p p p" in commands,
        "Changed units or boundaries",
    )
    require(
        "pair_style lj/cut/coul/long " + format(config["nonbonded_cutoff"], ".17g")
        in commands,
        "Changed cutoff",
    )
    mixing = "geometric" if data["family"] == "OPLS-AA" else "arithmetic"
    require(
        "pair_modify mix " + mixing + " shift no tail no table 0" in commands,
        "Changed LJ convention",
    )
    require(
        "kspace_style ewald " + format(config["pme_tolerance"], ".17g") in commands,
        "Changed electrostatics",
    )
    lj14, coul14 = config["one_four_scaling"]
    require(
        "special_bonds lj 0 0 "
        + format(lj14, ".17g")
        + " coul 0 0 "
        + format(coul14, ".17g")
        + " angle no dihedral no"
        in commands,
        "Changed exclusions or 1-4 scaling",
    )
    return _signed(
        dict(
            schema="island.engine-input-parity.v1",
            backend_identity=payload["identity"],
            graph_identity=payload["final_graph_identity"],
            assignment_identity=payload["assignment_identity"],
            charge_identity=payload["charge_identity"],
            source=config["source"],
            atom_order="exact",
            coordinates="exact",
            box="exact",
            mass="exact",
            charge="exact",
            lj="17-digit serialization with unit conversion",
            units="real: kcal/mol, Angstrom, elementary charge",
            switching="none",
            dispersion_tail=False,
            exclusions=config["excluded_pair_policy"],
            one_four_scaling=config["one_four_scaling"],
            mixing=mixing,
            electrostatics="Ewald; matching requested tolerance does not guarantee matching truncation error",
            parameter_identity="same bound assignment; bonded energies independently compared by engine",
            passed=True,
            **FLAGS,
        )
    )
