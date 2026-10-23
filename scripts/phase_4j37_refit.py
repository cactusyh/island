"""Independent external CDS RESP refit of retained raw ESP, no ISLAND imports.

This replays the fitting calculation, not geometry optimization or Psi4 SCF.
Native float64 output must reproduce the authoritative charge bits.
"""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    import numpy as np
    import psi4
    import resp
    from resp import espfit

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    root = args.directory
    data = json.loads((root / 'input.json').read_text())
    result = {'psi4_version': psi4.__version__, 'resp_module_version': resp.__version__,
              'quantum_scf_reexecuted': False, 'native_resp_refit_executed': True, 'components': [], 'passed': False}
    try:
        if psi4.__version__ != '1.9.1' or resp.__version__ != '0.8':
            raise ValueError('Wrong pinned external runtime')
        for index, component in enumerate(data['components']):
            directory = root / f'component-{index}'
            fit = np.load(directory / 'fit-input.npz', allow_pickle=False)
            geometry = np.load(directory / 'geometry-bohr.npy', allow_pickle=False)
            grid = np.load(directory / 'grid-bohr.npy', allow_pickle=False)
            esp = np.load(directory / 'esp-native.npy', allow_pickle=False)
            if not np.array_equal(fit['esp'].view(np.uint64), esp.view(np.uint64)):
                raise ValueError('Fitting ESP input differs from raw native ESP')
            # Reconstruct inverse distance paths independently of original loop.
            inverse = 1 / np.linalg.norm(grid[:, None, :] - geometry[None, :, :], axis=2)
            inverse_error = float(np.max(abs(inverse - fit['invr'])))
            if inverse_error > 1e-13:
                raise ValueError('Fitting geometry/grid association mismatch')
            options = json.loads((directory / 'fit-options.json').read_text())
            fitted, labels, note = espfit.fit(options, {'natoms': len(component['atoms']), 'symbols': [a['element'] for a in component['atoms']],
                'mol_charge': component['formal_charge'], 'invr': [fit['invr']], 'esp_values': [esp]})
            expected = np.load(directory / 'charges-native.npy', allow_pickle=False)
            bit_exact = np.array_equal(fitted[-1].view(np.uint64), expected.view(np.uint64))
            if not bit_exact or labels != ['ESP', 'RESP'] or note:
                raise ValueError('External RESP refit does not reproduce retained charge bits')
            result['components'].append({'atom_ids': [a['atom_id'] for a in component['atoms']], 'bit_exact': bit_exact,
                'inverse_distance_max_error_bohr_inverse': inverse_error,
                'charges_sha256': hashlib.sha256((directory / 'charges-native.npy').read_bytes()).hexdigest()})
        result['passed'] = True
    except Exception as error:  # noqa: BLE001 -- retain raw replay diagnostics
        result['diagnostic'] = str(error)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False);stream.write('\n')
    return int(not result['passed'])


if __name__ == '__main__':
    raise SystemExit(main())
