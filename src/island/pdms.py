"""Explicit local-research PDMS workflow; no native or production promotion.

OpenMM Ewald implements the *declared* author DFT/NPA physical model with a
separately recorded tighter numerical solver. NPT/pressure remain unavailable.
"""

import json
from copy import deepcopy
from hashlib import sha256
from math import isfinite, pi, sqrt
from pathlib import Path

import numpy as np

from island.evaluation.models import EvaluationResult, fingerprint
from island.forcefields.external_pdms import LIMITS
from island.forcefields.pdms_opls import (
    PDMSOPLSError,
    check,
    export_lammps,
    save_preparation,
)

GROUPS = ("bond", "angle", "proper_torsion", "lj", "coulomb")
SOLVER = {
    "method": "Ewald",
    "error_tolerance": 1e-7,
    "platform": "Reference",
    "author_original": "LAMMPS PPPM 1e-4; not the same numerical solver",
}
FD_BOUND = 1e-4  # Existing strict bound, in kJ/(mol angstrom).


def _mm():
    try:
        import openmm
        from openmm import unit

        return openmm, unit
    except ImportError as error:
        raise PDMSOPLSError("asset_unavailable", "openmm_not_installed") from error


def box_lengths(system):
    check(system.box is not None, "explicit_periodic_box_required", "typing_blocked")
    check(
        tuple(system.box.periodic) == (True, True, True),
        "three_periodic_dimensions_required",
        "typing_blocked",
    )
    lengths = tuple(float(x) for x in system.box.lengths)
    check(
        all(isfinite(x) and x > 22 for x in lengths),
        "box_must_exceed_twice_11A_cutoff",
        "typing_blocked",
    )
    return lengths


def special_pairs(p):
    adjacent = {r["site_id"]: set() for r in p["mapping"]}
    for row in p["interactions"]["bond"]:
        a, b = row["sites"]
        adjacent[a].add(b)
        adjacent[b].add(a)
    result = {}
    for i in adjacent:
        visited, frontier = {i}, {i}
        for distance in (1, 2, 3):
            frontier = {j for k in frontier for j in adjacent[k]} - visited
            visited.update(frontier)
            for j in frontier:
                if i < j:
                    result[i, j] = 0.5 if distance == 3 else 0.0
    return result


def build_openmm_model(p, lengths):
    """No generic XML defaults, Class-II substitution, or charge manipulation."""
    mm, _ = _mm()
    model = mm.System()
    rows = {r["evidence_checksum"]: r for r in p["source"]["rows"]}
    atoms = p["mapping"]
    index = {r["site_id"]: n for n, r in enumerate(atoms)}
    for atom in atoms:
        model.addParticle(float(atom["mass_dalton"]))
    model.setDefaultPeriodicBoxVectors(
        *[
            mm.Vec3(*(lengths[i] * 0.1 if j == i else 0 for j in range(3)))
            for i in range(3)
        ]
    )
    bond, angle = mm.HarmonicBondForce(), mm.HarmonicAngleForce()
    proper = mm.CustomTorsionForce(
        "0.5*(k1*(1+cos(theta))+k2*(1-cos(2*theta))+k3*(1+cos(3*theta))+k4*(1-cos(4*theta)))"
    )
    for n in range(1, 5):
        proper.addPerTorsionParameter("k" + str(n))
    for kind, force in (("bond", bond), ("angle", angle), ("proper_torsion", proper)):
        force.setUsesPeriodicBoundaryConditions(True)
        for r in p["interactions"][kind]:
            values = list(map(float, rows[r["row_checksum"]]["values"]))
            ids = [index[i] for i in r["sites"]]
            if kind == "bond":
                force.addBond(*ids, values[1] * 0.1, values[0] * 836.8)
            elif kind == "angle":
                force.addAngle(*ids, values[1] * pi / 180, values[0] * 8.368)
            else:
                force.addTorsion(*ids, [v * 4.184 for v in values])
    switch = "select(step(r-0.9),(1.21-r*r)^2*(1.21+2*r*r-2.43)/0.4^3,1)"
    lj = mm.CustomNonbondedForce(
        "4*sqrt(e1*e2)*((sqrt(s1*s2)/r)^12-(sqrt(s1*s2)/r)^6)*" + switch
    )
    lj.addPerParticleParameter("s")
    lj.addPerParticleParameter("e")
    lj.setNonbondedMethod(lj.CutoffPeriodic)
    lj.setCutoffDistance(1.1)
    lj.setUseLongRangeCorrection(False)
    coul = mm.NonbondedForce()
    coul.setNonbondedMethod(coul.Ewald)
    coul.setCutoffDistance(1.1)
    coul.setEwaldErrorTolerance(SOLVER["error_tolerance"])
    coul.setUseDispersionCorrection(False)
    coul.setExceptionsUsePeriodicBoundaryConditions(True)
    pair = {
        r["sites"][0]: list(map(float, rows[r["row_checksum"]]["values"]))
        for r in p["interactions"]["nonbonded"]
    }
    for atom in atoms:
        epsilon, sigma = pair[atom["site_id"]]
        lj.addParticle([sigma * 0.1, epsilon * 4.184])
        coul.addParticle(float(atom["charge_e"]), 1, 0)
    lj14 = mm.CustomBondForce("step(1.1-r)*4*e*((s/r)^12-(s/r)^6)*" + switch)
    lj14.addPerBondParameter("s")
    lj14.addPerBondParameter("e")
    lj14.setUsesPeriodicBoundaryConditions(True)
    for (a, b), scale in special_pairs(p).items():
        i, j = index[a], index[b]
        lj.addExclusion(i, j)
        coul.addException(
            i,
            j,
            scale * float(atoms[i]["charge_e"]) * float(atoms[j]["charge_e"]),
            1,
            0,
        )
        if scale:
            lj14.addBond(
                i,
                j,
                [
                    sqrt(pair[a][1] * pair[b][1]) * 0.1,
                    scale * sqrt(pair[a][0] * pair[b][0]) * 4.184,
                ],
            )
    for group, force in enumerate((bond, angle, proper, lj, coul)):
        force.setForceGroup(group)
        model.addForce(force)
    lj14.setForceGroup(3)
    model.addForce(lj14)
    return model


class PDMSOPLSEvaluator:
    def __init__(self, system, prepared):
        prepared.require_executable(system)
        self.system = deepcopy(system)
        self.prepared = prepared
        self.lengths = box_lengths(system)
        self.ids = [r["site_id"] for r in prepared.payload["mapping"]]
        mm, _ = _mm()
        self.model = build_openmm_model(prepared.payload, self.lengths)
        self.integrator = mm.VerletIntegrator(0.0005)
        self.context = mm.Context(
            self.model, self.integrator, mm.Platform.getPlatformByName("Reference")
        )
        self.model_identity = fingerprint(
            {
                "preparation": prepared.payload["identity"],
                "box": self.lengths,
                "solver": SOLVER,
            }
        )

    def validate_system(self, system):
        self.prepared.require_executable(system)
        check(
            box_lengths(system) == self.lengths, "stale_periodic_box", "typing_blocked"
        )

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        self.validate_system(self.system)
        check(coordinate_unit == "angstrom", "coordinate_unit_must_be_angstrom")
        xyz = (
            coordinates
            if coordinates is not None
            else {i: self.system.coordinates.get(i) for i in self.ids}
        )
        check(
            set(xyz) == set(self.ids), "coordinate_atom_ids_mismatch", "typing_blocked"
        )
        array = np.array([xyz[i] for i in self.ids], dtype=float)
        check(
            array.shape == (len(self.ids), 3) and np.isfinite(array).all(),
            "invalid_coordinates",
        )
        self.context.setPositions(array * 0.1)
        return self._read(self.context, array)

    def _read(self, context, xyz):
        mm, unit = _mm()
        components = {
            name: context.getState(getEnergy=True, groups={g})
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
            for g, name in enumerate(GROUPS)
        }
        forces = (
            context.getState(getForces=True)
            .getForces(asNumpy=True)
            .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            * 0.1
        )
        coord = fingerprint(
            {"ids": self.ids, "xyz": np.asarray(xyz).tolist(), "box": self.lengths}
        )
        return EvaluationResult(
            sum(components.values()),
            components,
            {i: tuple(map(float, forces[n])) for n, i in enumerate(self.ids)},
            coord,
            self.prepared.payload["identity"],
            self.model_identity,
            fingerprint({"coordinates": coord, "model": self.model_identity}),
            "OpenMM external PDMS OPLS-compatible local research",
            mm.__version__,
            "Reference",
            {
                "Ewald_tolerance": SOLVER["error_tolerance"],
                "native_opls_coverage": False,
                "pressure_available": False,
                "scientifically_validated": False,
            },
        )

    def close(self):
        if self.context is not None:
            del self.context
            self.context = None


def finite_difference(evaluator, coordinates=None, *, all_components=True):
    """Strict per-run check, never a tolerance-widening retry."""
    xyz = (
        {
            i: np.array(evaluator.system.coordinates.get(i), dtype=float)
            for i in evaluator.ids
        }
        if coordinates is None
        else {i: np.array(coordinates[i], dtype=float) for i in evaluator.ids}
    )
    baseline = evaluator.evaluate(xyz)
    selected = (
        evaluator.ids
        if all_components
        else list(
            {
                next(
                    a["site_id"]
                    for a in evaluator.prepared.payload["mapping"]
                    if a["element"] == e
                )
                for e in ("H", "O", "Si", "C")
            }
        )
    )
    comparisons = []
    for i in selected:
        for axis in range(3):
            for h in (1e-4, 5e-5):
                plus, minus = deepcopy(xyz), deepcopy(xyz)
                plus[i][axis] += h
                minus[i][axis] -= h
                force = -(
                    evaluator.evaluate(plus).potential_energy
                    - evaluator.evaluate(minus).potential_energy
                ) / (2 * h)
                comparisons.append(abs(force - baseline.forces[i][axis]))
    evaluator.evaluate(xyz)
    maximum = max(comparisons)
    return {
        "passed": maximum <= FD_BOUND,
        "comparisons": len(comparisons),
        "max_error_kj_mol_angstrom": maximum,
        "bound": FD_BOUND,
        "all_components": all_components,
    }


def run_pdms_md(
    system,
    prepared,
    directory,
    *,
    steps=1000,
    temperature_kelvin=300.0,
    timestep_fs=0.5,
    seed=2026,
    report_interval=100,
    minimize_tolerance_kj_mol_nm=1.0,
    minimize_max_iterations=10000,
):
    """Minimize, check strict FD, run short NVT, retain trajectory/checkpoint.

    This explicit research workflow has no NPT, pressure or equilibration claim.
    The supplied graph and source are rechecked before any backend execution.
    """
    from island.workflows.bundle import system_data

    prepared.require_executable(system)
    check(type(steps) is int and 0 < steps <= 100000, "short_nvt_step_limit")
    check(type(seed) is int and 0 < seed < 2**31, "positive_seed_required")
    check(
        type(report_interval) is int and report_interval > 0,
        "positive_report_interval_required",
    )
    check(
        isfinite(timestep_fs) and 0 < timestep_fs <= 0.5,
        "timestep_must_be_at_most_0_5fs",
    )
    check(
        isfinite(temperature_kelvin) and temperature_kelvin > 0,
        "positive_temperature_required",
    )
    check(
        isfinite(minimize_tolerance_kj_mol_nm) and minimize_tolerance_kj_mol_nm > 0,
        "positive_minimization_tolerance_required",
    )
    check(
        type(minimize_max_iterations) is int and minimize_max_iterations > 0,
        "positive_minimization_budget_required",
    )
    root = Path(directory)
    check(not root.exists(), "output_directory_exists")
    evaluator = PDMSOPLSEvaluator(system, prepared)
    mm, unit = _mm()
    root.mkdir(parents=True)
    try:
        save_preparation(prepared, root / "preparation.json")
        (root / "system.json").write_text(pack_system(system_data(system)))
        exports = export_lammps(system, prepared, root / "lammps")
        xyz = np.array([system.coordinates.get(i) for i in evaluator.ids])
        evaluator.context.setPositions(xyz * 0.1)
        initial = evaluator._read(evaluator.context, xyz)
        (root / "progress.json").write_text(
            json.dumps({"stage": "minimization"}) + "\n"
        )
        mm.LocalEnergyMinimizer.minimize(
            evaluator.context,
            minimize_tolerance_kj_mol_nm / sqrt(3),
            minimize_max_iterations,
        )
        state = evaluator.context.getState(getPositions=True, getForces=True)
        xyz = (
            np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            * 10
        )
        fmax_nm = float(
            np.max(
                np.linalg.norm(
                    state.getForces(asNumpy=True).value_in_unit(
                        unit.kilojoule_per_mole / unit.nanometer
                    ),
                    axis=1,
                )
            )
        )
        minimum = evaluator._read(evaluator.context, xyz)
        # OpenMM uses RMS force for convergence, report both RMS and max explicitly.
        rms = float(
            np.sqrt(
                np.mean(np.sum(np.array(list(minimum.forces.values())) ** 2, axis=1))
            )
            * 10
        )
        check(
            rms <= minimize_tolerance_kj_mol_nm,
            "minimization_not_converged",
            "scientifically_unqualified",
            rms_force_kj_mol_nm=rms,
        )
        coords = {i: xyz[n] for n, i in enumerate(evaluator.ids)}
        (root / "progress.json").write_text(
            json.dumps({"stage": "finite_difference"}) + "\n"
        )
        fd = finite_difference(evaluator, coords, all_components=True)
        (root / "finite-difference.json").write_text(json.dumps(fd, indent=2) + "\n")
        check(
            fd["passed"], "finite_difference_failed", "scientifically_unqualified", **fd
        )
        prepared.require_executable(system)
        evaluator.close()
        (root / "progress.json").write_text(json.dumps({"stage": "NVT"}) + "\n")
        integrator = mm.LangevinMiddleIntegrator(
            temperature_kelvin, 1.0, timestep_fs * 0.001
        )
        integrator.setRandomNumberSeed(seed)
        context = mm.Context(
            evaluator.model, integrator, mm.Platform.getPlatformByName("Reference")
        )
        context.setPositions(xyz * 0.1)
        context.setVelocitiesToTemperature(temperature_kelvin, seed)
        frames, thermo = [], []

        def capture(step):
            s = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            x = (
                np.array(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
                * 10
            )
            v = (
                np.array(
                    s.getVelocities(asNumpy=True).value_in_unit(
                        unit.nanometer / unit.picosecond
                    )
                )
                * 10
            )
            pe = float(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            ke = float(s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
            check(
                np.isfinite(x).all() and np.isfinite(v).all() and isfinite(pe + ke),
                "nonfinite_md_state",
                "scientifically_unqualified",
            )
            frames.append(
                {
                    "step": step,
                    "coordinates_angstrom": x.tolist(),
                    "velocities_angstrom_ps": v.tolist(),
                }
            )
            thermo.append(
                {
                    "step": step,
                    "time_ps": step * timestep_fs * 0.001,
                    "potential_kj_mol": pe,
                    "kinetic_kj_mol": ke,
                    "temperature_kelvin": 2 * ke / (3 * len(x) * 0.00831446261815324),
                }
            )

        capture(0)
        done = 0
        while done < steps:
            count = min(report_interval, steps - done)
            prepared.require_executable(system)
            integrator.step(count)
            done += count
            capture(done)
        checkpoint = context.createCheckpoint()
        (root / "checkpoint.chk").write_bytes(checkpoint)
        (root / "openmm-system.xml").write_text(
            mm.XmlSerializer.serialize(evaluator.model)
        )
        (root / "integrator.xml").write_text(mm.XmlSerializer.serialize(integrator))
        (root / "trajectory.json").write_text(json.dumps(frames) + "\n")
        result = {
            "status": "completed_local_research",
            "profile_identity": prepared.payload["profile_identity"],
            "preparation_identity": prepared.payload["identity"],
            "solver": SOLVER,
            "steps": steps,
            "timestep_fs": timestep_fs,
            "seed": seed,
            "initial_energy_kj_mol": initial.potential_energy,
            "minimized_energy_kj_mol": minimum.potential_energy,
            "minimization_rms_force_kj_mol_nm": rms,
            "minimization_max_force_kj_mol_nm": fmax_nm,
            "finite_difference": fd,
            "thermodynamics": thermo,
            "lammps_export_hashes": exports,
            "limitations": LIMITS,
            "scientific_status": "scientifically_unqualified",
            "NPT_run": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "openmm_version": mm.__version__,
        }
        (root / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        (root / "progress.json").write_text(json.dumps({"stage": "completed"}) + "\n")
        files = {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in root.iterdir()
            if p.is_file()
        }
        (root / "manifest.json").write_text(
            json.dumps({"files": files}, indent=2) + "\n"
        )
        context = None
        return result
    finally:
        evaluator.close()


def pack_system(value):
    from island.charge_references.records import pack

    return pack(value)


def restore_pdms_system(payload):
    """Keep the historical nonperiodic workflow deserializer unchanged."""
    from island.core import SimulationBox
    from island.workflows.bundle import system_from

    value = deepcopy(payload)
    box = value["box"]
    check(
        isinstance(box, dict) and set(box) == {"lengths", "periodic"},
        "saved_periodic_box_required",
        "typing_blocked",
    )
    value["box"] = None
    system = system_from(value)
    system.box = SimulationBox(*box["lengths"], periodic=tuple(box["periodic"]))
    box_lengths(system)
    return system


def resume_pdms_md(
    directory, output, *, options, manifest_sha256, steps=1000, report_interval=100
):
    """Resume an explicit, independently pinned local NVT checkpoint, offline."""
    from island.charge_references.records import unpack
    from island.forcefields.pdms_opls import load_preparation

    root, target = Path(directory), Path(output)
    try:
        raw = (root / "manifest.json").read_bytes()
    except OSError as error:
        raise PDMSOPLSError(
            "asset_unavailable", "checkpoint_manifest_unavailable", reason=str(error)
        ) from error
    check(sha256(raw).hexdigest() == manifest_sha256, "manifest_hash_mismatch")
    try:
        files = json.loads(raw)["files"]
        check(
            type(files) is dict
            and all(type(n) is str and type(p) is str for n, p in files.items()),
            "invalid_checkpoint_manifest",
        )
    except (ValueError, KeyError, TypeError) as error:
        raise PDMSOPLSError(
            "parameter_incomplete", "invalid_checkpoint_manifest", reason=str(error)
        ) from error
    for name, pin in files.items():
        check(Path(name).name == name, "invalid_manifest_path")
        check((root / name).is_file(), "checkpoint_asset_missing", file=name)
        check(
            sha256((root / name).read_bytes()).hexdigest() == pin,
            "checkpoint_asset_hash_mismatch",
            file=name,
        )
    check(not target.exists(), "output_directory_exists")
    check(type(steps) is int and 0 < steps <= 100000, "short_nvt_step_limit")
    check(
        type(report_interval) is int and report_interval > 0,
        "positive_report_interval_required",
    )
    system = restore_pdms_system(unpack((root / "system.json").read_text()))
    prepared = load_preparation(
        root / "preparation.json", files["preparation.json"], system, options
    )
    previous = json.loads((root / "result.json").read_text())
    mm, unit = _mm()
    check(
        mm.__version__ == previous["openmm_version"],
        "checkpoint_backend_version_mismatch",
    )
    model = build_openmm_model(prepared.payload, box_lengths(system))
    check(
        mm.XmlSerializer.serialize(model) == (root / "openmm-system.xml").read_text(),
        "stale_backend_model",
    )
    integrator = mm.XmlSerializer.deserialize((root / "integrator.xml").read_text())
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.loadCheckpoint((root / "checkpoint.chk").read_bytes())
        # No velocity reinitialization or reseeding after checkpoint loading.
        frames = []
        done = 0
        while done < steps:
            prepared.require_executable(system)
            count = min(report_interval, steps - done)
            integrator.step(count)
            done += count
            s = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            x = (
                np.array(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
                * 10
            )
            v = (
                np.array(
                    s.getVelocities(asNumpy=True).value_in_unit(
                        unit.nanometer / unit.picosecond
                    )
                )
                * 10
            )
            pe = float(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            ke = float(s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
            check(
                np.isfinite(x).all() and np.isfinite(v).all() and isfinite(pe + ke),
                "nonfinite_md_state",
            )
            frames.append(
                {
                    "step": context.getStepCount(),
                    "coordinates_angstrom": x.tolist(),
                    "velocities_angstrom_ps": v.tolist(),
                    "potential_kj_mol": pe,
                    "kinetic_kj_mol": ke,
                    "time_ps": context.getStepCount() * previous["timestep_fs"] * 0.001,
                    "temperature_kelvin": 2 * ke / (3 * len(x) * 0.00831446261815324),
                }
            )
        target.mkdir(parents=True)
        for name in (
            "preparation.json",
            "system.json",
            "openmm-system.xml",
            "integrator.xml",
        ):
            (target / name).write_bytes((root / name).read_bytes())
        (target / "checkpoint.chk").write_bytes(context.createCheckpoint())
        (target / "trajectory.json").write_text(json.dumps(frames) + "\n")
        report = {
            **previous,
            "steps": context.getStepCount(),
            "continued_steps": steps,
            "parent_manifest_sha256": manifest_sha256,
            "last_potential_kj_mol": frames[-1]["potential_kj_mol"],
            "thermodynamics": [
                {
                    k: v
                    for k, v in frame.items()
                    if k not in ("coordinates_angstrom", "velocities_angstrom_ps")
                }
                for frame in frames
            ],
        }
        (target / "result.json").write_text(json.dumps(report, indent=2) + "\n")
        pins = {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in target.iterdir()
            if p.is_file()
        }
        (target / "manifest.json").write_text(
            json.dumps({"files": pins}, indent=2) + "\n"
        )
        return report
    finally:
        del context
