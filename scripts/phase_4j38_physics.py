"""Pinned Foyer.apply/ParmEd independent periodic OPLS energy and force check."""

from pathlib import Path

import numpy as np

from scripts.phase_4j32_real_fixture import XML
from scripts.phase_4j38_contract import require
from scripts.validate_oplsaa_energy import reference, reference_value


def check(backend, frames):
    import openmm as mm

    from island._periodic_openmm import build
    from island.periodic import periodic_chemical_view

    model, witness = reference(periodic_chemical_view(backend.system), XML)
    p = backend.payload
    cfg = p["config"]
    lengths = cfg["box_lengths"]
    model.setDefaultPeriodicBoxVectors(
        *[
            mm.Vec3(*(lengths[j] / 10 if i == j else 0 for j in range(3)))
            for i in range(3)
        ]
    )
    for f in model.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setNonbondedMethod(mm.NonbondedForce.Ewald)
            f.setCutoffDistance(cfg["nonbonded_cutoff"] / 10)
            f.setEwaldErrorTolerance(cfg["pme_tolerance"])
            f.setExceptionsUsePeriodicBoundaryConditions(True)
        elif isinstance(f, mm.CustomNonbondedForce):
            f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            f.setCutoffDistance(cfg["nonbonded_cutoff"] / 10)
        else:
            f.setUsesPeriodicBoundaryConditions(True)
    production = build(p["numerical_data"], cfg, mm)
    for f in production.getForces():
        f.setForceGroup({"bond": 0, "angle": 1, "rb_proper": 2}.get(f.getName(), 3))
    snapshots = []
    for i in (0, len(frames) - 1):
        f = frames[i]
        box = [mm.Vec3(*v) for v in f["box_vectors_nm"]]
        model.setDefaultPeriodicBoxVectors(*box)
        production.setDefaultPeriodicBoxVectors(*box)
        a = reference_value(model, np.array(f["positions_nm"]) * 10)
        b = reference_value(production, np.array(f["positions_nm"]) * 10)
        error = abs(a[0] - b[0])
        force = float(np.max(abs(a[1] - b[1])))
        components = {k: abs(a[2][k] - b[2][k]) for k in a[2]}
        require(
            error <= 1e-5 and force <= 1e-5 and max(components.values()) <= 1e-5,
            f"Independent periodic reference mismatch: {error}, {force}, {components}",
        )
        snapshots.append(
            {
                "frame_index": i,
                "energy_error_kj_mol": error,
                "force_error_kj_mol_angstrom": force,
                "component_errors_kj_mol": components,
            }
        )
    return {
        "passed": True,
        "reference": "pinned Foyer.apply -> ParmEd geometric mixing periodic Ewald vs ISLAND; same OpenMM Reference engine",
        "independent_engine": False,
        "lammps_kspace_executed": False,
        "source": str(Path(XML)),
        "witness": witness,
        "snapshots": snapshots,
        "energy_tolerance_kj_mol": 1e-5,
        "force_tolerance_kj_mol_angstrom": 1e-5,
    }
