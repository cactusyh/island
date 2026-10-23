"""Independent periodic AmberPrmtopFile energy/force/component cross-check.

Both paths use OpenMM Reference; this checks parameter import/construction,
not an independent integrator or LAMMPS KSPACE validation.
"""

import argparse
import json
from pathlib import Path


def check(case):
    import numpy as np
    import openmm as mm
    from openmm import app, unit

    from island._periodic_openmm import build
    from island.periodic_md import _backend

    case = Path(case)
    p = json.loads((case / 'uninterrupted/backend.json').read_text());_backend(p)
    cfg = p['config'];data = p['numerical_data']
    ids = sorted(s['id'] for s in data['sites'])
    mapping = p['parameter_assignment']['native_payload']['record']['mapping']
    if [mapping[str(i)] for i in range(len(ids))] != ids:
        raise ValueError('Reference comparison requires verified atom ordering')
    source = app.AmberPrmtopFile(str(case / 'native-components/result.prmtop'))
    lengths = cfg['box_lengths']
    vectors = [mm.Vec3(*(lengths[j] / 10 if i == j else 0 for j in range(3))) for i in range(3)]
    source.topology.setPeriodicBoxVectors(vectors * unit.nanometer)
    reference = source.createSystem(nonbondedMethod=app.Ewald, nonbondedCutoff=cfg['nonbonded_cutoff'] / 10 * unit.nanometer,
        constraints=None, rigidWater=False, removeCMMotion=False, ewaldErrorTolerance=cfg['pme_tolerance'])
    production = build(data, cfg, mm)
    groups = {'HarmonicBondForce': 0, 'HarmonicAngleForce': 1, 'PeriodicTorsionForce': 2, 'NonbondedForce': 3,
              'CustomNonbondedForce': 3, 'CustomBondForce': 3}
    for system in (reference, production):
        for force in system.getForces():
            force.setForceGroup(groups[type(force).__name__])
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
                force.setUseSwitchingFunction(False)
                force.setExceptionsUsePeriodicBoundaryConditions(True)
            if isinstance(force, (mm.HarmonicBondForce, mm.HarmonicAngleForce, mm.PeriodicTorsionForce)):
                force.setUsesPeriodicBoundaryConditions(True)
    frames = json.loads((case / 'uninterrupted/frames.json').read_text())['frames']
    snapshots = []
    contexts = []
    try:
        for system in (reference, production):
            integrator = mm.VerletIntegrator(.001)
            contexts.append((mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference')), integrator))
        for index in (0, len(frames) - 1):
            frame = frames[index]
            measurements = []
            for context, _ in contexts:
                context.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in frame['box_vectors_nm']])
                context.setPositions(np.array(frame['positions_nm']) * unit.nanometer)
                state = context.getState(getEnergy=True, getForces=True)
                measurements.append((state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                    state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer) * .1,
                    [context.getState(getEnergy=True, groups=1 << g).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) for g in range(4)]))
            a, b = measurements
            energy_error = abs(a[0] - b[0]);force_error = float(np.max(abs(a[1] - b[1])))
            component_errors = {name: abs(x - y) for name, x, y in zip(('bond', 'angle', 'torsion', 'nonbonded'), a[2], b[2], strict=True)}
            if energy_error > 1e-6 or force_error > 1e-5 or max(component_errors.values()) > 1e-6:
                raise ValueError(f'Periodic reference mismatch: {energy_error}, {force_error}, {component_errors}')
            snapshots.append({'frame_index': index, 'energy_error_kj_mol': energy_error,
                'force_error_kj_mol_angstrom': force_error, 'component_errors_kj_mol': component_errors,
                'reference_components_kj_mol': dict(zip(('bond', 'angle', 'torsion', 'nonbonded'), a[2], strict=True))})
    finally:
        contexts.clear()
    return {'passed': True, 'reference': 'OpenMM AmberPrmtopFile Ewald vs ISLAND periodic construction; Reference platform',
            'independent_engine': False, 'lammps_kspace_executed': False, 'snapshots': snapshots,
            'energy_tolerance_kj_mol': 1e-6, 'force_tolerance_kj_mol_angstrom': 1e-5}


def main():
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('--case', type=Path, required=True);p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    try:
        result = check(args.case)
    except Exception as error:  # noqa: BLE001 -- fail with preserved diagnostics
        result = {'passed': False, 'diagnostic': str(error)}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False);stream.write('\n')
    return int(not result['passed'])


if __name__ == '__main__':
    raise SystemExit(main())
