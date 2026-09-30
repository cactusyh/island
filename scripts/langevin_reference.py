"""Independent source-prmtop BAOAB reference with prescribed Gaussian increments.

Validation only: no production noise replay API. Full-step velocities are explicit.
"""

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
    prmtop, mapping, coordinates, velocities, masses, options, checkpoints, noise
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
    integrator = mm.CustomIntegrator(options.timestep_fs * 0.001)
    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.addComputePerDof("x", "x+0.5*dt*v")
    integrator.addGlobalVariable(
        "c", np.exp(-options.friction_per_ps * options.timestep_fs * 0.001)
    )
    integrator.addGlobalVariable(
        "q", -np.expm1(-2 * options.friction_per_ps * options.timestep_fs * 0.001)
    )
    integrator.addGlobalVariable("kT", 0.00831446261815324 * options.temperature_kelvin)
    integrator.addPerDofVariable("noise", 0)
    integrator.addComputePerDof("v", "c*v+sqrt(q*kT/m)*noise")
    integrator.addComputePerDof("x", "x+0.5*dt*v")
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
            # Each prescribed row uses ascending stable-ID, xyz ordering.
            ordered = {site: i for i, site in enumerate(sorted(masses))}
            for index in range(previous, step):
                integrator.setPerDofVariableByName(
                    "noise",
                    [
                        mm.Vec3(*noise[index, ordered[mapping[i]]])
                        for i in range(len(mapping))
                    ],
                )
                integrator.step(1)
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
