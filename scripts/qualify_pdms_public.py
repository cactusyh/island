"""Independent row/equation/LAMMPS comparison of actual external PDMS graphs.

The raw-row parser and direct Ewald equations below do not import production
parameter parsers. This command never changes scientific readiness.
"""

import argparse
import itertools
import json
import math
import os
import subprocess
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.special import erfc

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.core import SimulationBox
from island.forcefields.pdms_opls import PDMSOPLSOptions, export_lammps
from island.forcefields.research import (
    ForceFieldRequest,
    create_evaluator,
    prepare_forcefield,
)

ENERGY_BOUND = 3e-4
FORCE_BOUND = 3e-4
MM_COULOMB = 1389.35457644382
LMP_COULOMB = 332.06371 * 4.184


def molecule(dp):
    if dp == 0:
        system = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
        anchors = ()
    else:
        system = build_linear_polymer(
            "[*:1]O[Si](C)(C)[*:2]",
            dp=dp,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
            add_hydrogens=True,
            random_seed=2026,
        )
        anchors = (system.metadata["polymer"]["head_site_id"],)
    system.box = SimulationBox(40, 40, 40)
    return system, anchors


def raw_rows(archive):
    with ZipFile(archive) as z:
        names = z.namelist()

        def text(suffix):
            name = next(n for n in names if n.endswith("/param/OPLS/" + suffix))
            return z.read(name).decode()

        settings = {}
        for n, line in enumerate(text("opls_system.in.settings").splitlines(), 1):
            t = line.split("#")[0].split()
            if t:
                key = t[0], int(t[1])
                assert key not in settings
                settings[key] = (n, t[3:] if t[0] == "pair_coeff" else t[2:])
        return settings


class IndependentEnergy:
    """Explicit real + reciprocal + self + background + exclusion correction."""

    def __init__(self, p, length=40, alpha=0.22):
        self.p, self.length, self.alpha = p, length, alpha
        self.ids = {r["site_id"]: i for i, r in enumerate(p["mapping"])}
        self.q = np.array([float(r["charge_e"]) for r in p["mapping"]])
        self.rows = {r["evidence_checksum"]: r for r in p["source"]["rows"]}
        self.images = length * np.array(list(itertools.product(range(-1, 2), repeat=3)))
        h = np.array(
            [v for v in itertools.product(range(-18, 19), repeat=3) if v != (0, 0, 0)]
        )
        self.k = 2 * math.pi / length * h
        k2 = np.sum(self.k**2, axis=1)
        self.weight = 2 * math.pi / length**3 * np.exp(-k2 / (4 * alpha**2)) / k2
        adjacency = {i: set() for i in self.ids}
        for r in p["interactions"]["bond"]:
            a, b = r["sites"]
            adjacency[a].add(b)
            adjacency[b].add(a)
        self.specials = {}
        for i in adjacency:
            reached, frontier = {i}, {i}
            for d in (1, 2, 3):
                frontier = {k for j in frontier for k in adjacency[j]} - reached
                reached.update(frontier)
                for j in frontier:
                    if i < j:
                        self.specials[self.ids[i], self.ids[j]] = 0.5 if d == 3 else 0
        self.pair = {
            self.ids[r["sites"][0]]: list(
                map(float, self.rows[r["row_checksum"]]["values"])
            )
            for r in p["interactions"]["nonbonded"]
        }

    def __call__(self, xyz):
        values = {k: 0.0 for k in ("bond", "angle", "proper_torsion", "lj", "coulomb")}
        for kind in ("bond", "angle", "proper_torsion"):
            for r in self.p["interactions"][kind]:
                ids = [self.ids[i] for i in r["sites"]]
                v = list(map(float, self.rows[r["row_checksum"]]["values"]))
                if kind == "bond":
                    energy = (
                        v[0] * (np.linalg.norm(xyz[ids[0]] - xyz[ids[1]]) - v[1]) ** 2
                    )
                elif kind == "angle":
                    a, b, c = xyz[ids]
                    u, w = a - b, c - b
                    theta = np.arccos(
                        np.clip(
                            np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w), -1, 1
                        )
                    )
                    energy = v[0] * (theta - math.radians(v[1])) ** 2
                else:
                    a, b, c, d = xyz[ids]
                    n1, n2 = np.cross(b - a, c - b), np.cross(c - b, d - c)
                    phi = np.arccos(
                        np.clip(
                            np.dot(n1, n2) / np.linalg.norm(n1) / np.linalg.norm(n2),
                            -1,
                            1,
                        )
                    )
                    energy = sum(
                        k / 2 * (1 + (-1) ** n * np.cos((n + 1) * phi))
                        for n, k in enumerate(v)
                    )
                values[kind] += energy * 4.184
        real, correction = 0.0, 0.0
        for i in range(len(xyz)):
            for j in range(len(xyz)):
                r = np.linalg.norm(xyz[i] - xyz[j] + self.images, axis=1)
                r = r[r > 1e-12]
                real += 0.5 * self.q[i] * self.q[j] * np.sum(erfc(self.alpha * r) / r)
            for j in range(i + 1, len(xyz)):
                delta = xyz[i] - xyz[j]
                delta -= self.length * np.rint(delta / self.length)
                r = np.linalg.norm(delta)
                scale = self.specials.get((i, j), 1)
                correction -= (1 - scale) * self.q[i] * self.q[j] / r
                if r < 11:
                    switch = (
                        1
                        if r <= 9
                        else (121 - r * r) ** 2 * (121 + 2 * r * r - 243) / 40**3
                    )
                    e = math.sqrt(self.pair[i][0] * self.pair[j][0])
                    sigma = math.sqrt(self.pair[i][1] * self.pair[j][1])
                    values["lj"] += (
                        scale
                        * 4
                        * 4.184
                        * e
                        * ((sigma / r) ** 12 - (sigma / r) ** 6)
                        * switch
                    )
        structure = np.exp(1j * (self.k @ xyz.T)) @ self.q
        reciprocal = np.sum(self.weight * np.abs(structure) ** 2)
        self_energy = -self.alpha / math.sqrt(math.pi) * np.dot(self.q, self.q)
        background = (
            -math.pi * np.sum(self.q) ** 2 / (2 * self.alpha**2 * self.length**3)
        )
        # OpenMM electrostatic conversion, kept separate from LAMMPS real units.
        values["coulomb"] = MM_COULOMB * (
            real + reciprocal + self_energy + background + correction
        )
        return values


def lammps_erfc_energy_difference(xyz, charges, alpha):
    """Audit the actual pinned LAMMPS erfc polynomial; do not alter its output."""
    delta = 0.0
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            r = np.linalg.norm(xyz[i] - xyz[j])
            if r < 11:
                z = alpha * r
                t = 1 / (1 + 0.3275911 * z)
                approximate = (
                    t
                    * (
                        0.254829592
                        + t
                        * (
                            -0.284496736
                            + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))
                        )
                    )
                    * np.exp(-z * z)
                )
                delta += (
                    LMP_COULOMB * charges[i] * charges[j] / r * (approximate - erfc(z))
                )
    return float(delta)


def qualify(args):
    assert sha256(args.archive.read_bytes()).hexdigest() == args.archive_sha256
    assert sha256(args.lammps.read_bytes()).hexdigest() == args.lammps_sha256
    source = args.lammps_pair_source.read_bytes()
    assert sha256(source).hexdigest() == args.lammps_pair_source_sha256
    constants = args.lammps_ewald_constants.read_bytes()
    assert sha256(constants).hexdigest() == args.lammps_ewald_constants_sha256
    for expression in (
        b"EWALD_P = 0.3275911",
        b"A1 = 0.254829592",
        b"A2 = -0.284496736",
        b"A3 = 1.421413741",
        b"A4 = -1.453152027",
        b"A5 = 1.061405429",
    ):
        assert expression in constants
    assert b"erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2" in source
    args.output.mkdir(parents=True, exist_ok=False)
    raw = raw_rows(args.archive)
    cases = {}
    for dp in (0, 2, 4):
        system, anchors = molecule(dp)
        xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
        system.coordinates.translate(20 - np.mean(xyz, axis=0))
        options = PDMSOPLSOptions(args.archive, args.archive_sha256, anchors)
        prepared = prepare_forcefield(
            system, ForceFieldRequest("external_pdms_opls", options)
        )
        p = prepared.payload
        prepared.require_executable(system)
        for row in p["source"]["rows"]:
            command = {
                "bond": "bond_coeff",
                "angle": "angle_coeff",
                "proper_torsion": "dihedral_coeff",
                "nonbonded": "pair_coeff",
            }[row["family"]]
            line, vals = raw[command, row["author_type_id"]]
            assert vals == row["values"] and line == row["settings_evidence"]["line"]
        evaluator = create_evaluator(system, prepared)
        try:
            result = evaluator.evaluate()
            xyz = np.array([system.coordinates.get(i) for i in evaluator.ids])
            independent = IndependentEnergy(p)
            energies = independent(xyz)
            component_error = {
                k: abs(energies[k] - result.energy_components[k]) for k in energies
            }
            fd_errors = []
            for element in ("H", "O", "Si", "C"):
                n = next(
                    n for n, a in enumerate(p["mapping"]) if a["element"] == element
                )
                for axis in range(3):
                    h = 5e-5
                    plus, minus = xyz.copy(), xyz.copy()
                    plus[n, axis] += h
                    minus[n, axis] -= h
                    f = -(
                        sum(independent(plus).values())
                        - sum(independent(minus).values())
                    ) / (2 * h)
                    fd_errors.append(abs(f - result.forces[evaluator.ids[n]][axis]))
            work = args.output / f"dp{dp}"
            exports = export_lammps(system, prepared, work)
            script = (
                (work / "system.in").read_text()
                + """pair_modify table 0
kspace_style pppm 1e-10
kspace_modify order 7
thermo 1
thermo_style custom step ebond eangle edihed evdwl ecoul elong pe
thermo_modify format float %.17g
run 0
variable eb equal ebond
variable ea equal eangle
variable ed equal edihed
variable el equal evdwl
variable ec equal ecoul+elong
print "COMPARE ${eb} ${ea} ${ed} ${el} ${ec}"
write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g
"""
            )
            (work / "compare.in").write_text(script)
            run = subprocess.run(
                [
                    str(args.lammps.resolve()),
                    "-in",
                    "compare.in",
                    "-log",
                    "compare.log",
                    "-screen",
                    "compare.screen",
                ],
                cwd=work,
                env={**os.environ, "OMP_NUM_THREADS": "1"},
                capture_output=True,
                text=True,
                timeout=180,
                check=False,
            )
            assert run.returncode == 0, (run.stdout, run.stderr)
            lines = (work / "compare.screen").read_text().splitlines()
            literal = next(
                line
                for line in lines
                if line.startswith("COMPARE ") and "${" not in line
            )
            lmp_e = [float(v) * 4.184 for v in literal.split()[1:]]
            groups = list(energies)
            lmp_component_errors = {
                k: abs(v - result.energy_components[k]) for k, v in zip(groups, lmp_e)
            }
            alpha = float(
                next(
                    line.split("=")[1].strip()
                    for line in lines
                    if "G vector (1/distance) =" in line
                )
            )
            charges = np.array([float(a["charge_e"]) for a in p["mapping"]])
            erfc_difference = lammps_erfc_energy_difference(xyz, charges, alpha)
            unit_constant_difference = energies["coulomb"] * (
                LMP_COULOMB / MM_COULOMB - 1
            )
            expected_lammps_coulomb = (
                energies["coulomb"] + erfc_difference + unit_constant_difference
            )
            algorithm_residual = abs(lmp_e[-1] - expected_lammps_coulomb)
            dump = (work / "forces.dump").read_text().splitlines()
            start = (
                next(n for n, s in enumerate(dump) if s.startswith("ITEM: ATOMS")) + 1
            )
            lmp_f = (
                np.array([list(map(float, s.split()[1:])) for s in dump[start:]])
                * 4.184
            )
            force_error = float(
                np.max(
                    np.abs(lmp_f - np.array([result.forces[i] for i in evaluator.ids]))
                )
            )
            case = {
                "atoms": len(p["mapping"]),
                "raw_rows_reconstructed": 17,
                "independent_component_errors_kj_mol": component_error,
                "independent_force_max_error_kj_mol_A": max(fd_errors),
                "lammps_component_errors_kj_mol": lmp_component_errors,
                "lammps_force_max_error_kj_mol_A": force_error,
                "lammps_exact_equation_energy_parity": max(
                    lmp_component_errors.values()
                )
                <= ENERGY_BOUND,
                "documented_lammps_erfc_energy_difference_kj_mol": erfc_difference,
                "lammps_algorithm_residual_kj_mol": algorithm_residual,
                "lammps_pair_source_sha256": args.lammps_pair_source_sha256,
                "lammps_ewald_constants_sha256": args.lammps_ewald_constants_sha256,
                "export_hashes": exports,
            }
            case["passed"] = (
                max(component_error.values()) <= ENERGY_BOUND
                and max(v for k, v in lmp_component_errors.items() if k != "coulomb")
                <= ENERGY_BOUND
                and algorithm_residual <= ENERGY_BOUND
                and max(fd_errors) <= 1e-4
                and force_error <= FORCE_BOUND
            )
            cases[str(dp)] = case
            (args.output / "qualification.json").write_text(
                json.dumps({"cases": cases}, indent=2) + "\n"
            )
        finally:
            evaluator.close()
    return {
        "cases": cases,
        "passed": all(v["passed"] for v in cases.values()),
        "scientifically_validated": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    parser.add_argument("--lammps-sha256", required=True)
    parser.add_argument("--lammps-pair-source", type=Path, required=True)
    parser.add_argument("--lammps-pair-source-sha256", required=True)
    parser.add_argument("--lammps-ewald-constants", type=Path, required=True)
    parser.add_argument("--lammps-ewald-constants-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = qualify(args)
    (args.output / "qualification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
