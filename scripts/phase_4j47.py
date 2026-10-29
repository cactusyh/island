"""Reproduce J47 static acceptance; nonzero while any required gate is open.

The source-authored run is separate from the historical high-precision numerical
control. No adjusted energy, force, coefficient, charge or tolerance is accepted.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_lammps import (
    ForceFieldRequest,
    OPLSOptions,
    prepare_forcefield,
)
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError, check
from island.pdms import SOLVER, _mm, build_openmm_model, special_pairs
from island.pdms_lammps import (
    ENERGY_BOUND,
    FORCE_BOUND,
    LAMMPSOptions,
    pin,
    reconstruct,
    run_static,
)
from scripts.qualify_pdms_public import (
    MM_COULOMB,
    IndependentEnergy,
    molecule,
    qualify,
    raw_rows,
)

BASE = "9e7d9765d0248e14320b4214db8831771f06e84e"
PAIR_PIN = "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5"
CONSTANTS_PIN = "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c"


def write(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, default=lambda v: v.item()) + "\n"
    )


class Reference:
    """Identical historical model; split force groups only for measurement."""

    def __init__(self, payload, xyz):
        self.mm, self.unit = _mm()
        self.model = build_openmm_model(payload, (40, 40, 40))
        self.model.getForce(4).setReciprocalSpaceForceGroup(5)
        self.model.getForce(5).setForceGroup(6)  # explicit LJ14 bonds
        self.integrator = self.mm.VerletIntegrator(0.0005)
        self.context = self.mm.Context(
            self.model, self.integrator, self.mm.Platform.getPlatformByName("Reference")
        )
        self.xyz = xyz
        self.set(xyz)

    def set(self, xyz, matrix=None):
        if matrix is None:
            matrix = np.eye(3)
        vectors = 4 * matrix.T
        self.context.setPeriodicBoxVectors(*[self.mm.Vec3(*v) for v in vectors])
        self.context.setPositions(xyz * 0.1)

    def energy(self, group=None):
        kw = {} if group is None else {"groups": {group}}
        return float(
            self.context.getState(getEnergy=True, **kw)
            .getPotentialEnergy()
            .value_in_unit(self.unit.kilojoule_per_mole)
        )

    def forces(self, group=None):
        kw = {} if group is None else {"groups": {group}}
        return (
            np.asarray(
                self.context.getState(getForces=True, **kw)
                .getForces(asNumpy=True)
                .value_in_unit(self.unit.kilojoule_per_mole / self.unit.nanometer)
            )
            * 0.1
        )

    def components(self):
        groups = [self.energy(i) for i in range(7)]
        return {
            "bond": groups[0],
            "angle": groups[1],
            "proper_torsion": groups[2],
            "improper": {
                "status": "not_applicable",
                "reason": "no declared-domain improper interactions",
            },
            "lj": groups[3] + groups[6],
            "lj14": groups[6],
            "coulomb_real_corrected": groups[4],
            "coulomb_reciprocal_self_background": groups[5],
            "coulomb": groups[4] + groups[5],
            "total": sum(groups),
        }

    def virial(self):
        values = []
        for a, b in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
            estimates = []
            for h in (1e-4, 5e-5):
                energy = []
                for sign in (1, -1):
                    matrix = np.eye(3)
                    matrix[a, b] += sign * h
                    self.set(self.xyz @ matrix.T, matrix)
                    energy.append(self.energy())
                estimates.append(-(energy[0] - energy[1]) / (2 * h))
            values.append(
                {
                    "values_at_1e-4_and_5e-5": estimates,
                    "convergence_error": abs(estimates[0] - estimates[1]),
                }
            )
        self.set(self.xyz)
        return values

    def close(self):
        del self.context


def special_energies(p, xyz):
    """Independent explicit graph-distance special pairs (never subtract output)."""
    ids = {r["site_id"]: n for n, r in enumerate(p["mapping"])}
    q = [float(r["charge_e"]) for r in p["mapping"]]
    rows = {r["evidence_checksum"]: r for r in p["source"]["rows"]}
    pair = {
        r["sites"][0]: list(map(float, rows[r["row_checksum"]]["values"]))
        for r in p["interactions"]["nonbonded"]
    }
    coul, lj = 0.0, 0.0
    pairs = []
    for (a, b), scale in special_pairs(p).items():
        if not scale:
            continue
        i, j = ids[a], ids[b]
        delta = xyz[i] - xyz[j]
        delta -= 40 * np.rint(delta / 40)
        r = np.linalg.norm(delta)
        coul += scale * MM_COULOMB * q[i] * q[j] / r
        if r < 11:
            switch = (
                1 if r <= 9 else (121 - r * r) ** 2 * (121 + 2 * r * r - 243) / 40**3
            )
            e, s = np.sqrt(np.array(pair[a]) * pair[b])
            lj += scale * 4 * 4.184 * e * ((s / r) ** 12 - (s / r) ** 6) * switch
        pairs.append([a, b])
    return {
        "coulomb14_exact_equation_kj_mol": float(coul),
        "lj14_kj_mol": float(lj),
        "pairs": pairs,
        "partition_note": "Coulomb14 is an independently reconstructed physical partition; LAMMPS ecoul includes exclusion corrections, not a dedicated 14 observable",
    }


def compare(system, prepared, observed):
    backend = prepared.validate_integrity(system)
    p = backend.payload
    xyz = np.asarray([system.coordinates.get(r["site_id"]) for r in p["mapping"]])
    # Prove actual dumped coordinates are identical before any numerical comparison.
    assert np.array_equal(xyz, np.asarray(observed["coordinates_angstrom"]))
    reference = Reference(p, xyz)
    try:
        ref = reference.components()
        lmp = observed["energy_kj_mol"]
        lmp_coul = (
            lmp["coulomb_real_corrected"] + lmp["coulomb_reciprocal_self_background"]
        )
        errors = {
            k: abs(lmp[k] - ref[k])
            for k in ("bond", "angle", "proper_torsion", "lj", "total")
        }
        errors["coulomb"] = abs(lmp_coul - ref["coulomb"])
        force = reference.forces()
        f_error = float(np.max(np.abs(force - observed["forces_kj_mol_angstrom"])))
        virial = reference.virial()
        v_error = [
            abs(v["values_at_1e-4_and_5e-5"][-1] - x)
            for v, x in zip(virial, observed["virial_kj_mol_xx_yy_zz_xy_xz_yz"])
        ]
        fd = []
        for element in ("H", "O", "Si", "C"):
            n = next(n for n, a in enumerate(p["mapping"]) if a["element"] == element)
            for axis in range(3):
                for h in (1e-4, 5e-5):
                    energies = []
                    for sign in (1, -1):
                        x = xyz.copy()
                        x[n, axis] += sign * h
                        reference.set(x)
                        energies.append(reference.energy())
                    fd.append(
                        float(
                            abs(-(energies[0] - energies[1]) / (2 * h) - force[n, axis])
                        )
                    )
        reference.set(xyz)
        independent = IndependentEnergy(p)
        exact = independent(xyz)
        raw_errors = {k: abs(exact[k] - ref[k]) for k in exact}
        independent_force_errors = []
        for element in ("H", "O", "Si", "C"):
            n = next(n for n, a in enumerate(p["mapping"]) if a["element"] == element)
            for axis in range(3):
                plus, minus = xyz.copy(), xyz.copy()
                plus[n, axis] += 5e-5
                minus[n, axis] -= 5e-5
                f = (
                    -(
                        sum(independent(plus).values())
                        - sum(independent(minus).values())
                    )
                    / 1e-4
                )
                independent_force_errors.append(float(abs(f - force[n, axis])))
        independent_trace = []
        for h in (1e-4, 5e-5):
            plus = sum(
                IndependentEnergy(p, length=40 * (1 + h))(xyz * (1 + h)).values()
            )
            minus = sum(
                IndependentEnergy(p, length=40 * (1 - h))(xyz * (1 - h)).values()
            )
            independent_trace.append(float(-(plus - minus) / (2 * h)))
        special = special_energies(p, xyz)
        components = {}
        for k, error in errors.items():
            components[k] = {
                "absolute_error": error,
                "bound": ENERGY_BOUND,
                "passed": error <= ENERGY_BOUND,
                "cause": "precision"
                if k in ("bond", "angle", "proper_torsion", "lj")
                else "long-range solver",
            }
        for k in ("coulomb_real_corrected", "coulomb_reciprocal_self_background"):
            components[k] = {
                "LAMMPS": lmp[k],
                "OpenMM": ref[k],
                "absolute_difference": abs(lmp[k] - ref[k]),
                "cause": "long-range solver",
                "comparable_partition": False,
                "reason": "PPPM and Ewald choose different splitting parameters; only the sum is invariant",
            }
        components["improper"] = ref["improper"]
        components["lj14"] = {
            "independent": special["lj14_kj_mol"],
            "OpenMM": ref["lj14"],
            "absolute_error": abs(special["lj14_kj_mol"] - ref["lj14"]),
            "LAMMPS": "part of evdwl; independent exact source-derived pair partition",
        }
        components["coulomb14"] = special
        components["forces"] = {
            "max_error_kj_mol_angstrom": f_error,
            "bound": FORCE_BOUND,
            "passed": f_error <= FORCE_BOUND,
            "cause": "long-range solver",
        }
        components["virial"] = {
            "OpenMM_affine_energy_derivative": virial,
            "independent_exact_equation_trace_at_1e-4_and_5e-5": independent_trace,
            "absolute_errors_kj_mol": v_error,
            "passed": False,
            "cause": "unresolved",
            "reason": "no inherited virial acceptance tolerance or matched reciprocal solver contract; diagnostic only",
        }
        components["pressure_proxy"] = {
            "status": "scientifically_unqualified",
            "cause": "unresolved",
            "backend_pressure_available": False,
            "NPT_permitted": False,
        }
        return {
            "OpenMM": ref,
            "LAMMPS": lmp,
            "component_comparison": components,
            "OpenMM_solver": SOLVER,
            "LAMMPS_solver": {
                "method": "pppm",
                "tolerance": 1e-4,
                "pair_table": "pinned default",
            },
            "identical_physical_inputs": True,
            "identical_numerical_long_range_settings": False,
            "independent_energy_kj_mol": {k: float(v) for k, v in exact.items()},
            "independent_exact_equation_errors_kj_mol": raw_errors,
            "independent_force_finite_difference": {
                "comparisons": len(independent_force_errors),
                "max_error_kj_mol_angstrom": max(independent_force_errors),
                "bound": 1e-4,
                "passed": max(independent_force_errors) <= 1e-4,
            },
            "OpenMM_finite_difference": {
                "comparisons": len(fd),
                "max_error": max(fd),
                "bound": 1e-4,
                "passed": max(fd) <= 1e-4,
            },
            "acceptance_passed": False,
            "blocked_stages": [
                "identical_long_range_solver_contract",
                "exact_component_energy_parity",
                "qualified_virial_comparison",
                "minimization",
                "short_NVT",
                "fresh_process_MD_continuation",
                "checkpoint_bytes_and_complete_RNG_equality",
            ],
            "MD_executed": False,
            "NPT_executed": False,
        }
    finally:
        reference.close()


def lammps_finite_difference(root, backend, engine, xyz, baseline_forces):
    """Evaluate the unchanged source LAMMPS model at both strict FD steps."""
    from island.pdms_lammps import _execute

    work = Path(root)
    work.mkdir()
    source = work.parent / "lammps"
    for name in ("system.data", "system.in", "coefficients.in"):
        shutil.copyfile(source / name, work / name)
    script = (work / "system.in").read_text() + "thermo 0\nrun 0\n"
    keys = []
    for element in ("H", "O", "Si", "C"):
        n = next(
            n
            for n, a in enumerate(backend.payload["mapping"])
            if a["element"] == element
        )
        for axis in range(3):
            for h in (1e-4, 5e-5):
                keys.append((n, axis, h))
                for sign in (1, -1):
                    coordinate = float(xyz[n, axis]) + sign * h
                    script += (
                        f"set atom {n + 1} {'xyz'[axis]} {coordinate:.17g}\nrun 0\n"
                    )
                    script += f'print "J47FD {len(keys) - 1} {sign} $(pe:%.17g)"\n'
                script += f"set atom {n + 1} {'xyz'[axis]} {float(xyz[n, axis]):.17g}\n"
    (work / "fd.in").write_text(script)
    code = _execute(engine, work, "fd.in")
    check(code == 0, "lammps_finite_difference_process_failed", returncode=code)
    records = {}
    for line in (work / "engine.screen").read_text().splitlines():
        if line.startswith("J47FD ") and "$(" not in line:
            _, index, sign, energy = line.split()
            records[int(index), int(sign)] = float(energy) * 4.184
    values = []
    for index, (n, axis, h) in enumerate(keys):
        force = -(records[index, 1] - records[index, -1]) / (2 * h)
        values.append(
            {
                "atom_index": n,
                "axis": axis,
                "h_angstrom": h,
                "finite_difference_force": force,
                "engine_force": baseline_forces[n][axis],
                "absolute_error": abs(force - baseline_forces[n][axis]),
            }
        )
    maximum = max(v["absolute_error"] for v in values)
    return {
        "comparisons": len(values),
        "bound": 1e-4,
        "max_error": maximum,
        "passed": maximum <= 1e-4,
        "records": values,
    }


def audit_main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    parser.add_argument("--lammps-sha256", required=True)
    parser.add_argument("--lammps-source", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--paper-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    check(
        pin(args.paper) == args.paper_sha256, "paper_hash_mismatch", "asset_unavailable"
    )
    check(
        pin(args.archive) == args.archive_sha256,
        "archive_hash_mismatch",
        "asset_unavailable",
    )
    args.output.mkdir(parents=True, exist_ok=False)
    engine = LAMMPSOptions(args.lammps, args.lammps_sha256)
    raw = raw_rows(args.archive)
    cases = {}
    for dp in (0, 2, 4, 6):
        system, anchors = molecule(dp)
        xyz = np.asarray([system.coordinates.get(i) for i in system.topology.sites])
        system.coordinates.translate(20 - xyz.mean(axis=0))
        request = ForceFieldRequest(
            "opls", OPLSOptions(PROFILE, args.archive, args.archive_sha256, anchors)
        )
        prepared = prepare_forcefield(system, request)
        backend = prepared.validate_integrity(system)
        for row in backend.payload["source"]["rows"]:
            command = {
                "bond": "bond_coeff",
                "angle": "angle_coeff",
                "proper_torsion": "dihedral_coeff",
                "nonbonded": "pair_coeff",
            }[row["family"]]
            n, values = raw[command, row["author_type_id"]]
            assert values == row["values"] and n == row["settings_evidence"]["line"]
        directory = args.output / f"dp{dp}"
        result = run_static(system, prepared, directory, engine=engine)
        assert result["status"] == "static_executed_unqualified", result
        reconstruct(
            directory,
            manifest_sha256=result["manifest_sha256"],
            request=request,
            engine=engine,
        )
        moved = args.output / f"relocated-dp{dp}"
        shutil.copytree(directory, moved)
        # Supply source and executable independently after relocation. Never discover assets.
        code = """import sys
from island.forcefields.pdms_lammps import *
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError, check
from island.pdms_lammps import *
r=ForceFieldRequest("opls",OPLSOptions(PROFILE,sys.argv[3],sys.argv[4],tuple(map(int,sys.argv[7:]))))
s,p,v=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],request=r,engine=LAMMPSOptions(sys.argv[5],sys.argv[6]))
print(p.payload["identity"])
"""
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                str(moved.resolve()),
                result["manifest_sha256"],
                str(args.archive.resolve()),
                args.archive_sha256,
                str(args.lammps.resolve()),
                args.lammps_sha256,
                *map(str, anchors),
            ],
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
        assert proc.stdout.strip() == prepared.payload["identity"]
        shutil.rmtree(moved)
        case = compare(system, prepared, result["observed"])
        xyz = np.array(result["observed"]["coordinates_angstrom"])
        case["LAMMPS_finite_difference"] = lammps_finite_difference(
            directory / "finite-difference-control",
            backend,
            engine,
            xyz,
            result["observed"]["forces_kj_mol_angstrom"],
        )
        shutil.move(
            str(directory / "finite-difference-control"),
            str(args.output / f"finite-difference-dp{dp}"),
        )
        case.update(
            {
                "atoms": len(backend.payload["mapping"]),
                "raw_rows_reconstructed": len(raw),
                "workflow_manifest_sha256": result["manifest_sha256"],
                "relocation_verified": True,
                "fresh_process_reconstruction_verified": True,
            }
        )
        cases[str(dp)] = case
        write(args.output / "cases.json", cases)
        print(
            f"DP{dp}: static executed; total difference {case['component_comparison']['total']['absolute_error']:.8g}; MD blocked",
            flush=True,
        )
    # Reproduce the already declared tighter historical control unchanged. This
    # explains native erfc/unit discrepancies; it cannot grant workflow authority.
    control = qualify(
        SimpleNamespace(
            archive=args.archive,
            archive_sha256=args.archive_sha256,
            lammps=args.lammps,
            lammps_sha256=args.lammps_sha256,
            lammps_pair_source=args.lammps_source
            / "KSPACE/pair_lj_charmm_coul_long.cpp",
            lammps_pair_source_sha256=PAIR_PIN,
            lammps_ewald_constants=args.lammps_source / "ewald_const.h",
            lammps_ewald_constants_sha256=CONSTANTS_PIN,
            output=args.output / "historical-numerical-control",
        )
    )
    write(args.output / "historical-control.json", control)
    from scripts.qualify_pdms_public import LMP_COULOMB

    diagnosis = {
        "native_author_settings_cause": "long-range solver",
        "historical_control_only": True,
        "output_adjustment_for_acceptance": False,
        "constants": {
            "OpenMM_kj_mol_angstrom_e2": MM_COULOMB,
            "LAMMPS_kj_mol_angstrom_e2": LMP_COULOMB,
        },
        "cases": {
            dp: {
                "unit_conversion_difference_kj_mol": cases[dp][
                    "independent_energy_kj_mol"
                ]["coulomb"]
                * (LMP_COULOMB / MM_COULOMB - 1),
                "erfc_polynomial_difference_kj_mol": c[
                    "documented_lammps_erfc_energy_difference_kj_mol"
                ],
                "algorithm_reconstruction_residual_kj_mol": c[
                    "lammps_algorithm_residual_kj_mol"
                ],
                "raw_energy_parity": c["lammps_exact_equation_energy_parity"],
                "causes": ["unit conversion", "precision"],
            }
            for dp, c in control["cases"].items()
        },
    }
    write(args.output / "diagnosis.json", diagnosis)
    summary = {
        "schema": "island.j47-static-acceptance.v1",
        "base_commit": BASE,
        "profile_identity": PROFILE,
        "cases": cases,
        "source_assets": {
            "archive": args.archive_sha256,
            "paper": args.paper_sha256,
            "lammps_pair_source": PAIR_PIN,
            "ewald_constants": CONSTANTS_PIN,
        },
        "exact_static_parity": False,
        "end_to_end_acceptance_gate": 1,
        "historical_control_explains_error": control["passed"],
        "MD_executed": False,
        "NPT_executed": False,
        "limitations": LIMITS,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    write(args.output / "summary.json", {**summary, "identity": digest(summary)})
    return 1


def main(argv=None):
    try:
        return audit_main(argv)
    except PDMSOPLSError as error:
        print(
            json.dumps(
                {
                    "status": error.status,
                    "diagnostics": error.diagnostics,
                    "end_to_end_acceptance_gate": 1,
                    "MD_executed": False,
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                },
                sort_keys=True,
            )
        )
        return 1
    except OSError as error:
        print(
            json.dumps(
                {
                    "status": "asset_unavailable",
                    "diagnostics": [
                        {"stage": "explicit_asset_unavailable", "reason": str(error)}
                    ],
                    "end_to_end_acceptance_gate": 1,
                    "MD_executed": False,
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                },
                sort_keys=True,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
