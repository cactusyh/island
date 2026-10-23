"""Independent read-only source-equation check; never a charge provider.

Reconstructs constrained ESP/RESP normal equations from archived ESP/grid and
optimized geometry, without importing resp or ISLAND's production extractor.
The saved grid has 10 decimal places, so reconstruction discrepancies are
reported as diagnostics, never used to repair the authoritative native vector.
"""

import json
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


def check(directory):
    directory = Path(directory)
    data = json.loads((directory / 'input.json').read_text())
    results = []
    for index, component in enumerate(data['components']):
        root = directory / f'component-{index}'
        geometry = np.load(root / 'optimized-geometry-bohr.npy', allow_pickle=False)
        # The worker defines units angstrom. Native driver converts its geometry
        # using this published constant, and grid.dat is in angstrom.
        conversion = 0.52917721092
        points = np.loadtxt(root / '1_j36_grid.dat')
        esp = np.loadtxt(root / '1_j36_grid_esp.dat')
        inverse = conversion / np.linalg.norm(points[:, None, :] - geometry[None, :, :] * conversion, axis=2)
        n = len(geometry)
        matrix = np.zeros((n + 1, n + 1))
        matrix[:n, :n] = inverse.T @ inverse
        matrix[:n, n] = matrix[n, :n] = 1
        rhs = np.r_[inverse.T @ esp, component['formal_charge']]
        solution = np.linalg.solve(matrix, rhs)
        options = data['resp_options']
        converged = False
        for iteration in range(options['MAX_IT']):
            restrained = matrix.copy()
            for j, atom in enumerate(component['atoms']):
                if atom['element'] != 'H' or not options['IHFREE']:
                    restrained[j, j] += options['RESP_A'] / np.sqrt(solution[j] ** 2 + options['RESP_B'] ** 2)
            current = np.linalg.solve(restrained, rhs)
            delta = float(np.max(abs(current[:n] - solution[:n])))
            solution = current
            if delta <= options['TOLER']:
                converged = True
                break
        native = np.load(root / 'native-resp-charges.npy', allow_pickle=False)
        with localcontext() as context:
            context.prec = 200
            binary_total = str(sum((Decimal.from_float(float(q)) for q in native), Decimal(0)))
        results.append({'component': index, 'grid_points': len(points), 'iterations': iteration + 1,
                        'fit_iteration_converged': converged, 'last_iteration_max_change_e': delta,
                        'independent_reconstruction_max_difference_e': float(np.max(abs(solution[:n] - native))),
                        'native_decimal_total_e': str(sum((Decimal(repr(float(q))) for q in native), Decimal(0))),
                        'native_binary_exact_total_e': binary_total,
                        'formal_total_e': component['formal_charge']})
    return {'schema': 'island.j36-independent-equations.v1', 'source': 'CDS RESP espfit.py fit/iterate/restraint, Bayly et al. equations 10,12-14',
            'charge_provider': False, 'normalization': False, 'redistribution': False, 'components': results}
