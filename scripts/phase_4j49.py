"""Raw canonical OPLS electrostatic diagnosis. Exit 1 while source gates fail.

Author runs and numerical controls are physically identical but explicitly NOT
numerically identical. No adjusted energy is used for acceptance. A diagnostic
kernel with libm erfc may demonstrate a backend defect; it cannot silently
replace the source-declared PPPM/table contract or grant MD authority.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from island.electrostatic_audit import COMPONENTS, bind_inputs, finite_difference
from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_registry import (
    OPLSAssets,
    OPLSResolutionError,
    canonical_registry,
    prepare_polymer_forcefield,
    use_opls_assets,
)
from island.forcefields.opls_registry.registry import PAPER_SHA256
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError, check
from island.forcefields.selection.inventory import digest
from island.opls_lammps import LAMMPSOptions, run_md, run_static
from island.pdms import _mm
from island.pdms_lammps import STATIC, _execute, engine_record, parse_static, pin
from scripts.phase_4j47 import Reference, compare, write
from scripts.phase_4j49_workflow import run_control
from scripts.qualify_pdms_public import molecule

BASE = "9e2e5df9cfebf119c0c857c841a401fd98fcb42a"
ALPHA = math.sqrt(-math.log(2e-7)) / 11  # OpenMM Ewald real-space error selection.
MODES = {
    "author": "",
    "no-table": "pair_modify table 0\n",
    "converged-pppm": "pair_modify table 0\nkspace_style pppm 1e-10\nkspace_modify order 7\n",
    "matched-ewald": f"pair_modify table 0\nkspace_style ewald 1e-7\nkspace_modify gewald {ALPHA:.17g} kmax/ewald 20 20 20\n",
}
STEPS = (2e-4, 1e-4, 5e-5)


def packed_components(ref):
    r = ref.components()
    return {
        "bonded": r["bond"] + r["angle"] + r["proper_torsion"],
        "lj": r["lj"],
        "coulomb": r["coulomb"],
        "total": r["total"],
    }


def component_forces(ref):
    return {
        "bonded": ref.forces(0) + ref.forces(1) + ref.forces(2),
        "lj": ref.forces(3) + ref.forces(6),
        "coulomb": ref.forces(4) + ref.forces(5),
        "total": ref.forces(),
    }


def lammps_samples(root, prefix, xyz, probes, engine, *, preparation=None):
    """One real process, raw positive/negative component energies and force dump."""
    root.mkdir()
    script = prefix + STATIC
    keys = []
    for n, axis in probes:
        for h in STEPS:
            for sign in (1, -1):
                displaced = xyz.copy()
                displaced[n, axis] += sign * h
                keys.append(displaced.tobytes())
                script += (
                    f"set atom {n + 1} {'xyz'[axis]} {displaced[n, axis]:.17g}\nrun 0\n"
                )
                script += (
                    f'print "J49FD {len(keys) - 1} $(ebond:%.17g) $(eangle:%.17g) '
                    '$(edihed:%.17g) $(evdwl:%.17g) $(ecoul:%.17g) $(elong:%.17g) $(pe:%.17g)"\n'
                )
            script += f"set atom {n + 1} {'xyz'[axis]} {xyz[n, axis]:.17g}\n"
    (root / "audit.in").write_text(script)
    code = _execute(engine, root, "audit.in")
    check(code == 0, "electrostatic_diagnostic_process_failed", returncode=code)
    rows = {}
    for line in (root / "engine.screen").read_text().splitlines():
        if line.startswith("J49FD ") and "$(" not in line:
            fields = line.split()
            e = np.array(list(map(float, fields[2:]))) * 4.184
            rows[keys[int(fields[1])]] = {
                "bonded": float(sum(e[:3])),
                "lj": float(e[3]),
                "coulomb": float(e[4] + e[5]),
                "total": float(e[6]),
            }
    check(len(rows) == len(keys), "missing_raw_fd_samples")
    if preparation is not None:
        observed = parse_static(root, preparation)
    else:
        lines = (root / "engine.screen").read_text().splitlines()
        i = next(
            i
            for i, line in enumerate(lines)
            if line.strip().startswith("Step ") and "E_bond" in line
        )
        v = np.array(list(map(float, lines[i + 1].split())))
        dump = (root / "state.dump").read_text().splitlines()
        j = next(i for i, line in enumerate(dump) if line.startswith("ITEM: ATOMS"))
        atoms = np.array([list(map(float, line.split())) for line in dump[j + 1 :]])
        observed = {
            "forces_kj_mol_angstrom": (atoms[:, 7:] * 4.184).tolist(),
            "coordinates_angstrom": atoms[:, 4:7].tolist(),
            "energy_kj_mol": dict(
                zip(
                    (
                        "bond",
                        "angle",
                        "proper_torsion",
                        "improper",
                        "lj",
                        "coulomb_real_corrected",
                        "coulomb_reciprocal_self_background",
                        "total",
                    ),
                    v[1:9] * 4.184,
                )
            ),
            "virial_kj_mol_xx_yy_zz_xy_xz_yz": (
                v[9:] * 40**3 * 4.184 / 68568.415
            ).tolist(),
        }
    check(
        np.array_equal(xyz, observed["coordinates_angstrom"]),
        "coordinate_order_changed",
    )
    return rows, observed


def analyze_mode(root, prefix, xyz, probes, engine, ref, preparation=None):
    samples, observed = lammps_samples(
        root, prefix, xyz, probes, engine, preparation=preparation
    )
    forces = component_forces(ref)
    raw_total_force = np.array(observed["forces_kj_mol_angstrom"])
    # Bonded/LJ analytical equations have independently passed J47/J48. LAMMPS
    # has no native separate Coulomb force observable for this combined pair style.
    # Label the subtraction explicitly; do not call it a direct engine measurement.
    lmpforces = {
        **forces,
        "total": raw_total_force,
        "coulomb": raw_total_force - forces["bonded"] - forces["lj"],
    }
    fd = finite_difference(
        lambda x: samples[x.tobytes()], xyz, lmpforces, probes=probes
    )
    reference = packed_components(ref)
    e = observed["energy_kj_mol"]
    raw = {
        "bonded": e["bond"] + e["angle"] + e["proper_torsion"],
        "lj": e["lj"],
        "coulomb": e["coulomb_real_corrected"]
        + e["coulomb_reciprocal_self_background"],
        "total": e["total"],
    }
    return {
        "LAMMPS": observed,
        "OpenMM_energy_kj_mol": reference,
        "OpenMM_force_components": {k: v.tolist() for k, v in forces.items()},
        "LAMMPS_force_components": {k: v.tolist() for k, v in lmpforces.items()},
        "force_component_measurement": {
            "total": "raw LAMMPS dump",
            "coulomb": "LAMMPS total minus independently verified OpenMM bonded and LJ; not a native separate observable",
            "bonded_and_lj": "OpenMM analytical reference; raw LAMMPS component energies differentiated independently",
        },
        "energy_error_kj_mol": {k: abs(raw[k] - reference[k]) for k in COMPONENTS},
        "force_error_kj_mol_angstrom": {
            k: float(np.max(np.abs(lmpforces[k] - forces[k]))) for k in COMPONENTS
        },
        "energy_bound": 3e-4,
        "force_bound": 3e-4,
        "finite_difference": fd,
        "numerical_input_addendum": prefix.split("include coefficients.in\n")[-1]
        if preparation
        else prefix,
        "OpenMM_coulomb_partition": {
            k: v for k, v in ref.components().items() if k.startswith("coulomb")
        },
        "solver_partitions_comparable": "only in matched-ewald control (same alpha); reciprocal truncations are separately converged",
        "diagnostic_only": True,
        "MD_authority": False,
    }


class ChargedReference:
    """Isolated charged backend control, no chemical profile or invented LJ rows."""

    def __init__(self, xyz, charges):
        mm, unit = _mm()
        self.unit = unit
        model = mm.System()
        coul = mm.NonbondedForce()
        coul.setNonbondedMethod(coul.Ewald)
        coul.setCutoffDistance(1.1)
        coul.setEwaldErrorTolerance(1e-7)
        coul.setUseDispersionCorrection(False)
        for q in charges:
            model.addParticle(1)
            coul.addParticle(
                q, 1, 0
            )  # pure Coulomb control; LJ is absent, not a filled missing row
        coul.setForceGroup(4)
        coul.setReciprocalSpaceForceGroup(5)
        model.addForce(coul)
        model.setDefaultPeriodicBoxVectors(
            mm.Vec3(4, 0, 0), mm.Vec3(0, 4, 0), mm.Vec3(0, 0, 4)
        )
        self.integrator = mm.VerletIntegrator(0.0005)
        self.context = mm.Context(
            model, self.integrator, mm.Platform.getPlatformByName("Reference")
        )
        self.set(xyz)

    def set(self, xyz):
        self.context.setPositions(xyz * 0.1)

    def energy(self, group=None):
        groups = {} if group is None else {"groups": {group}}
        return float(
            self.context.getState(getEnergy=True, **groups)
            .getPotentialEnergy()
            .value_in_unit(self.unit.kilojoule_per_mole)
        )

    def forces(self, group=None):
        groups = {} if group is None else {"groups": {group}}
        return (
            np.asarray(
                self.context.getState(getForces=True, **groups)
                .getForces(asNumpy=True)
                .value_in_unit(self.unit.kilojoule_per_mole / self.unit.nanometer)
            )
            * 0.1
        )

    def components(self):
        return {
            "bond": self.energy(0),
            "angle": self.energy(1),
            "proper_torsion": self.energy(2),
            "lj": self.energy(3),
            "coulomb": self.energy(),
            "total": self.energy(),
            "coulomb_real_corrected": self.energy(4),
            "coulomb_reciprocal_self_background": self.energy(5),
        }

    def close(self):
        del self.context


def charged_prefix(xyz, charges):
    atoms = "\n".join(
        f"{i + 1} 1 1 {q:.17g} {x:.17g} {y:.17g} {z:.17g}"
        for i, (q, (x, y, z)) in enumerate(zip(charges, xyz))
    )
    data = f"J49 pure-Coulomb algorithm control\n\n{len(xyz)} atoms\n1 atom types\n\n0 40 xlo xhi\n0 40 ylo yhi\n0 40 zlo zhi\n\nMasses\n\n1 1\n\nAtoms # full\n\n{atoms}\n"
    prefix = "units real\natom_style full\nboundary p p p\npair_style coul/long 11\nkspace_style pppm 1e-4\nread_data ../system.data\npair_coeff * *\n"
    return data, prefix


def process_relocation(root, result, args, identity):
    moved = root.parent / (root.name + "-relocated")
    shutil.copytree(root, moved)
    archive = moved.parent / "relocated-source.zip"
    shutil.copyfile(args.archive, archive)
    code = """import sys
from island.opls_lammps import reconstruct,LAMMPSOptions
from island.forcefields.opls_registry import OPLSAssets
from island.electrostatic_audit import bind_inputs
s,p,r=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(sys.argv[3],sys.argv[4]),engine=LAMMPSOptions(sys.argv[5],sys.argv[6]))
print(bind_inputs(s,p).identity)
"""
    identities = []
    for seed in (1, 99397):
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                str(moved.resolve()),
                result["manifest_sha256"],
                str(archive.resolve()),
                args.archive_sha256,
                str(args.lammps.resolve()),
                args.lammps_sha256,
            ],
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
        )
        identities.append(proc.stdout.strip())
    check(identities == [identity] * 2, "fresh_process_binding_changed")
    shutil.rmtree(moved)
    archive.unlink()
    return {
        "verified": True,
        "process_hash_seeds": [1, 99397],
        "binding_identities": identities,
    }


def audit(args):
    check(pin(args.paper) == PAPER_SHA256, "paper_hash_mismatch", "asset_unavailable")
    args.output.mkdir(parents=True, exist_ok=False)
    engine = LAMMPSOptions(args.lammps, args.lammps_sha256)
    exact = (
        LAMMPSOptions(args.exact_kernel, args.exact_kernel_sha256)
        if args.exact_kernel
        else None
    )
    write(
        args.output / "kspace-inventory.json",
        {
            "author": {"path": str(args.lammps.resolve()), **engine_record(engine)},
            "exact_erfc_diagnostic": {
                "path": str(args.exact_kernel.resolve()),
                **engine_record(exact),
            }
            if exact
            else None,
            "author_input": "lj/charmm/coul/long 9 11; PPPM 1e-4; default table; p p p",
            "diagnostic_precision": "double; exact-erfc is a MODIFIED executable, never author execution authority",
        },
    )
    write(args.output / "canonical-registry.json", canonical_registry())
    assets = OPLSAssets(args.archive, args.archive_sha256)
    cases = {}
    systems = []
    for dp in (0, 2, 4, 6):
        system, _ = molecule(dp)
        xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
        system.coordinates.translate(20 - xyz.mean(axis=0))
        with use_opls_assets(assets):
            prepared = prepare_polymer_forcefield(system, forcefield="opls")
        binding = bind_inputs(system, prepared)
        backend = prepared.validate_integrity(system)
        source = backend.validate_integrity(system)
        directory = args.output / f"dp{dp}"
        result = run_static(system, prepared, directory, engine=engine)
        check(result["status"] == "static_executed_unqualified", "static_failed")
        case = compare(system, backend, result["observed"])
        xyz = np.array(result["observed"]["coordinates_angstrom"])
        probes = [
            (
                next(
                    n
                    for n, a in enumerate(source.payload["mapping"])
                    if a["element"] == el
                ),
                axis,
            )
            for el in ("H", "O", "Si", "C")
            for axis in range(3)
        ]
        reference = Reference(source.payload, xyz)
        try:
            forces = component_forces(reference)

            def evaluate(x, reference=reference):
                reference.set(x)
                return packed_components(reference)

            mmfd = finite_difference(evaluate, xyz, forces, probes=probes)
            reference.set(xyz)
            modes = {}
            base = (directory / "lammps/system.in").read_text()
            for name, addendum in {
                **MODES,
                **({"exact-erfc-ewald": MODES["matched-ewald"]} if exact else {}),
            }.items():
                work = args.output / f"diagnostic-dp{dp}-{name}"
                # Coefficient and data bytes are copied, never regenerated/modified by a numerical control.
                prefix = base.replace(
                    "read_data system.data", f"read_data ../dp{dp}/lammps/system.data"
                ).replace(
                    "include coefficients.in",
                    f"include ../dp{dp}/lammps/coefficients.in",
                )
                modes[name] = analyze_mode(
                    work,
                    prefix + addendum,
                    xyz,
                    probes,
                    exact if name == "exact-erfc-ewald" else engine,
                    reference,
                    source,
                )
            write(args.output / f"components-dp{dp}.json", modes)
            write(args.output / f"finite-difference-openmm-dp{dp}.json", mmfd)
        finally:
            reference.close()
        write(
            args.output / f"input-binding-dp{dp}.json",
            {**binding.payload, "identity": binding.identity},
        )
        relocation = process_relocation(directory, result, args, binding.identity)
        try:
            run_md(system, prepared, args.output / f"forbidden-md-dp{dp}")
        except OPLSResolutionError as error:
            check(
                isinstance(error, OPLSResolutionError)
                and error.status == "scientifically_unqualified",
                "MD_gate_unclassified",
            )
        else:
            raise AssertionError("Unqualified source reached MD")
        check(
            not (args.output / f"forbidden-md-dp{dp}").exists(),
            "blocked_MD_created_artifact",
        )
        systems.append((dp, system, prepared, modes.get("exact-erfc-ewald")))
        cases[str(dp)] = {
            **case,
            "input_binding_identity": binding.identity,
            "public_forcefield": "opls",
            "user_profile_argument": None,
            "user_charge_method_argument": None,
            "selected_profile": prepared.payload["selected_profile"],
            "source_identity": binding.payload["source_identity"],
            "charge_origin": prepared.payload["charge_origin"],
            "domain": prepared.payload["domain"],
            "workflow_manifest_sha256": result["manifest_sha256"],
            "relocation": relocation,
            "OpenMM_FD_passed": mmfd["passed"],
            "mode_summary": {
                k: {
                    "energy_errors": v["energy_error_kj_mol"],
                    "force_errors": v["force_error_kj_mol_angstrom"],
                    "FD_passed": v["finite_difference"]["passed"],
                    "FD_maxima": v["finite_difference"]["max_error_by_component"],
                }
                for k, v in modes.items()
            },
            "minimization": "blocked_by_source_numerical_gates",
            "NVT": "blocked_by_source_numerical_gates",
            "continuation_checkpoint_RNG": "blocked; no MD artifacts created",
        }
        write(args.output / "cases.json", cases)
        print(
            f"DP{dp}: raw author static reproduced; component FD and explicit numerical controls recorded; MD blocked",
            flush=True,
        )
    controls = {}
    for name, xyz, charges in (
        (
            "charged-dimer",
            np.array([[14.0, 19.0, 20.0], [17.17, 20.23, 21.09]]),
            (1.0, -1.0),
        ),
        (
            "periodic-charged-dimer",
            np.array([[0.41, 19.0, 20.0], [38.17, 20.23, 21.09]]),
            (1.0, -1.0),
        ),
    ):
        work = args.output / name
        work.mkdir()
        data, prefix = charged_prefix(xyz, charges)
        (work / "system.data").write_text(data)
        ref = ChargedReference(xyz, charges)
        try:
            probes = list(np.ndindex(xyz.shape))
            f = component_forces(ref)

            def energy(x, ref=ref):
                ref.set(x)
                return packed_components(ref)

            mmfd = finite_difference(energy, xyz, f, probes=probes)
            ref.set(xyz)
            modes = {
                mode: analyze_mode(
                    work / mode,
                    prefix + addendum,
                    xyz,
                    probes,
                    exact if mode == "exact-erfc-ewald" else engine,
                    ref,
                )
                for mode, addendum in {
                    **MODES,
                    **({"exact-erfc-ewald": MODES["matched-ewald"]} if exact else {}),
                }.items()
            }
            controls[name] = {
                "coordinates_angstrom": xyz.tolist(),
                "charges_e": charges,
                "chemical_profile": None,
                "pure_coulomb_control": True,
                "OpenMM_FD": mmfd,
                "modes": modes,
            }
        finally:
            ref.close()
    write(args.output / "charged-controls.json", controls)
    diagnostic_workflows = {}
    if exact:
        controls_pass = all(
            c["modes"]["exact-erfc-ewald"]["finite_difference"]["passed"]
            and all(
                v <= 3e-4
                for v in c["modes"]["exact-erfc-ewald"]["energy_error_kj_mol"].values()
            )
            and all(
                v <= 3e-4
                for v in c["modes"]["exact-erfc-ewald"][
                    "force_error_kj_mol_angstrom"
                ].values()
            )
            for c in controls.values()
        )
        check(
            controls_pass,
            "generic_backend_control_failed",
            "scientifically_unqualified",
        )
        for dp, system, prepared, control in systems:
            diagnostic_workflows[str(dp)] = run_control(
                system,
                prepared,
                args.output / f"workflow-dp{dp}",
                engine=exact,
                static_control=control,
                solver=MODES["matched-ewald"],
            )
        write(args.output / "diagnostic-workflows.json", diagnostic_workflows)
    summary = {
        "schema": "island.j49-electrostatic-audit.v1",
        "base_commit": BASE,
        "profile_identity": PROFILE,
        "cases": cases,
        "new_physical_parameters": [],
        "new_profiles": [],
        "end_to_end_acceptance_gate": 1,
        "scientific_acceptance_gate": 1,
        "MD_executed": False,
        "author_profile_MD_executed": False,
        "diagnostic_MD_executed": bool(diagnostic_workflows),
        "diagnostic_workflows": diagnostic_workflows,
        "bounded_diagnostic_acceptance_gate": 0
        if diagnostic_workflows
        and all(
            c["bounded_diagnostic_workflow_passed"]
            for c in diagnostic_workflows.values()
        )
        else 1,
        "NPT_executed": False,
        "source_assets": {"archive": args.archive_sha256, "paper": PAPER_SHA256},
        "readiness": {
            "production_validated": False,
            "simulation_readiness": "not_established",
        },
        "limitations": LIMITS,
        "diagnostic_changes": [
            "pair-table off",
            "PPPM accuracy convergence",
            "matched Ewald alpha; converged reciprocal cutoff",
            "explicit separate exact-erfc diagnostic executable",
        ],
        "source_authority_unchanged": True,
        "corrections_applied_to_author_outputs": False,
        "missing_conventions": [
            "authorized common numerical solver/table contract",
            "qualified pressure/virial bound",
        ],
        "author_profile_downstream_blocked": [
            "LAMMPS minimization",
            "LAMMPS NVT",
            "MD continuation",
            "checkpoint and RNG equality",
            "NPT",
        ],
    }
    write(args.output / "summary.json", {**summary, "identity": digest(summary)})
    return 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "paper", "lammps", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    for name in ("archive-sha256", "lammps-sha256"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--exact-kernel", type=Path)
    p.add_argument("--exact-kernel-sha256")
    try:
        return audit(p.parse_args(argv))
    except (OPLSResolutionError, PDMSOPLSError, OSError) as error:
        print(
            json.dumps(
                {
                    "status": getattr(error, "status", "asset_unavailable"),
                    "diagnostics": getattr(
                        error,
                        "diagnostics",
                        [{"stage": "explicit_asset_unavailable", "reason": str(error)}],
                    ),
                    "end_to_end_acceptance_gate": 1,
                    "author_profile_MD_executed": False,
                },
                sort_keys=True,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
