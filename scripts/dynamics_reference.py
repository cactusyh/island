"""Independent original-prmtop OpenMM CustomIntegrator reference, validation only."""

import numpy as np

# Declared before acceptance runs. Coordinates A, velocities A/ps, E kJ/mol,
# forces kJ/(mol*A). These are trajectory comparisons, not physical validation.
TOLERANCES = {
    "position_atol": 1e-8,
    "velocity_atol": 1e-6,
    "energy_atol": 2e-7,
    "energy_rtol": 2e-10,
    "kinetic_atol": 1e-8,
    "force_atol": 2e-6,
    "force_rtol": 2e-9,
}


def reference_trajectory(
    prmtop, mapping, coordinates, velocities, masses, timestep_fs, checkpoints
):
    import openmm as mm
    from openmm import app, unit

    model = app.AmberPrmtopFile(str(prmtop)).createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=None,
        removeCMMotion=False,
        hydrogenMass=None,
    )
    source_masses = {}
    for index, site in mapping.items():
        source_masses[site] = model.getParticleMass(index).value_in_unit(unit.dalton)
        # Both paths explicitly use the authoritative graph masses. No redistribution.
        model.setParticleMass(index, masses[site] * unit.dalton)
    for force in model.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.setUseDispersionCorrection(False)
            force.setUseSwitchingFunction(False)
    assert model.getNumConstraints() == 0
    assert all(
        type(force).__name__
        in (
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        )
        for force in model.getForces()
    )
    integrator = mm.CustomIntegrator(timestep_fs * 0.001)
    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.addComputePerDof("x", "x+dt*v")
    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.setKineticEnergyExpression("m*v*v/2")
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(
            np.array([coordinates[mapping[i]] for i in range(len(mapping))])
            * 0.1
            * unit.nanometer
        )
        context.setVelocities(
            np.array([velocities[mapping[i]] for i in range(len(mapping))])
            * 0.1
            * unit.nanometer
            / unit.picosecond
        )
        reports = []
        previous = 0
        for step in checkpoints:
            if step > previous:
                integrator.step(step - previous)
            previous = step
            state = context.getState(
                getEnergy=True, getForces=True, getPositions=True, getVelocities=True
            )
            x = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
            v = (
                state.getVelocities(asNumpy=True).value_in_unit(
                    unit.nanometer / unit.picosecond
                )
                * 10
            )
            f = (
                state.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoule_per_mole / unit.nanometer
                )
                * 0.1
            )
            reports.append(
                {
                    "step": step,
                    "coordinates": {mapping[i]: x[i] for i in range(len(mapping))},
                    "velocities": {mapping[i]: v[i] for i in range(len(mapping))},
                    "forces": {mapping[i]: f[i] for i in range(len(mapping))},
                    "potential": state.getPotentialEnergy().value_in_unit(
                        unit.kilojoule_per_mole
                    ),
                    "kinetic": state.getKineticEnergy().value_in_unit(
                        unit.kilojoule_per_mole
                    ),
                }
            )
        return reports, source_masses
    finally:
        del context, integrator


def compare(result, references):
    result.validate_integrity()
    maximum = {
        "position": 0.0,
        "velocity": 0.0,
        "potential": 0.0,
        "kinetic": 0.0,
        "force": 0.0,
    }
    for frame, ref in zip(result.frames, references, strict=True):
        assert frame.step == ref["step"]
        for name, actual, expected, atol, rtol in (
            (
                "potential",
                frame.potential_energy,
                ref["potential"],
                TOLERANCES["energy_atol"],
                TOLERANCES["energy_rtol"],
            ),
            (
                "kinetic",
                frame.kinetic_energy,
                ref["kinetic"],
                TOLERANCES["kinetic_atol"],
                TOLERANCES["energy_rtol"],
            ),
        ):
            np.testing.assert_allclose(actual, expected, atol=atol, rtol=rtol)
            maximum[name] = max(maximum[name], abs(actual - expected))
        for site in frame.coordinates:
            for name, actual, expected, atol, rtol in (
                (
                    "position",
                    frame.coordinates[site],
                    ref["coordinates"][site],
                    TOLERANCES["position_atol"],
                    0,
                ),
                (
                    "velocity",
                    frame.velocities[site],
                    ref["velocities"][site],
                    TOLERANCES["velocity_atol"],
                    0,
                ),
                (
                    "force",
                    frame.evaluation.forces[site],
                    ref["forces"][site],
                    TOLERANCES["force_atol"],
                    TOLERANCES["force_rtol"],
                ),
            ):
                np.testing.assert_allclose(actual, expected, atol=atol, rtol=rtol)
                maximum[name] = max(
                    maximum[name], float(np.max(np.abs(np.array(actual) - expected)))
                )
    return maximum
