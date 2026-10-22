"""Independent native Psi4/RESP runner. No ISLAND charge or typing imports.

The authoritative charge vector is the untouched RESP API ndarray. Text tables
in results.out are not parsed. A failed exact charge-total gate is preserved;
this script never normalizes, redistributes, rounds or repairs fitted charges.
"""

import argparse
import hashlib
import json
import os
import traceback
from decimal import Decimal
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def execute(input_path, output):
    import numpy as np
    import psi4
    import resp

    data = json.loads(Path(input_path).read_text())
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write(output / 'input.json', data)
    rows, components = [], []
    result = {'schema': 'island.j36-native-quantum-result.v1',
              'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
              'psi4_version': psi4.__version__, 'resp_module_version': resp.__version__,
              'final_graph_identity': data['final_graph_identity'], 'input_sha256': sha(input_path),
              'method': data['method'], 'basis': data['basis'],
              'normalization': False, 'redistribution': False,
              'status': 'failed', 'acceptance_met': False, 'stage': 'geometry_optimization'}
    psi4.set_num_threads(int(os.environ.get('ISLAND_J36_QM_THREADS', '4')))
    psi4.set_memory('2 GB')
    old = Path.cwd()
    try:
        for index, component in enumerate(data['components']):
            directory = output / f'component-{index}'
            directory.mkdir()
            os.chdir(directory)
            psi4.core.clean()
            psi4.core.clean_options()
            psi4.core.set_output_file(str(directory / 'psi4.out'), False)
            atoms = component['atoms']
            xyz = np.array([a['position_angstrom'] for a in atoms])
            centered = xyz - xyz.mean(axis=0)
            text = f"{component['formal_charge']} {component['multiplicity']}\n" + '\n'.join(f"{a['element']} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g}" for a, r in zip(atoms, centered, strict=True)) + '\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n'
            Path('geometry-input.dat').write_text(text)
            molecule = psi4.geometry(text)
            molecule.set_name('j36')
            psi4.set_options({**data['geometry_optimization'], 'basis': data['optimization_basis']})
            energy = psi4.optimize(data['optimization_method'], molecule=molecule)
            geometry = molecule.geometry().np.copy()
            np.save('optimized-geometry-bohr.npy', geometry, allow_pickle=False)
            write('optimized-geometry.json', {'atom_ids': [a['atom_id'] for a in atoms], 'elements': [a['element'] for a in atoms], 'units': 'bohr', 'positions': geometry.tolist(), 'optimization_energy_hartree': float(energy)})
            result['stage'] = 'resp_fit'
            options = data['resp_options']
            write('resp-options.json', options)
            # qf[-1] is the full-precision native RESP vector (qf[0] is ESP).
            qf = resp.resp([molecule], options)
            np.save('native-resp-charges.npy', qf[-1], allow_pickle=False)
            np.save('native-esp-charges.npy', qf[0], allow_pickle=False)
            if 'Charge fitting did not converge' in Path('results.out').read_text():
                raise ValueError('Native RESP fitting did not converge')
            values = qf[-1].tolist()
            if not all(type(q) is float and np.isfinite(q) for q in values):
                raise ValueError('Native RESP returned nonfinite/nonfloat charges')
            current = [{'atom_id': atom['atom_id'], 'element': atom['element'], 'charge_e': value, 'decimal_e': repr(value)} for atom, value in zip(atoms, values, strict=True)]
            write('charge-vector.json', current)
            rows.extend(current)
            total = sum((Decimal(r['decimal_e']) for r in current), Decimal(0))
            components.append({'atom_ids': [a['atom_id'] for a in atoms], 'formal_total_e': component['formal_charge'], 'provided_total_e': str(total), 'exact_charge_conservation': total == component['formal_charge'], 'optimized_geometry_sha256': sha('optimized-geometry-bohr.npy'), 'native_charge_array_sha256': sha('native-resp-charges.npy')})
            write(output / 'partial-result.json', {**result, 'rows': rows, 'components': components})
        result['stage'] = 'component_total_validation'
        result['rows'] = sorted(rows, key=lambda r: r['atom_id'])
        result['components'] = components
        result['charge_generation_executed'] = True
        if not all(c['exact_charge_conservation'] for c in components):
            result['classification'] = 'charge_conservation_failed'
            result['diagnostic'] = 'Unchanged full-precision native RESP vector fails exact decimal component-total conservation; no normalization, redistribution or rounding permitted.'
        else:
            result.update(status='charges_generated', classification='native_resp_charge_generation_complete', acceptance_met=True)
    except Exception as error:  # noqa: BLE001 -- archive diagnostics, never call a failed fit acceptance
        result.update(classification='quantum_execution_failure', diagnostic=str(error), exception_type=type(error).__name__)
        (output / 'traceback.txt').write_text(traceback.format_exc())
    finally:
        psi4.core.flush_outfile()
        os.chdir(old)
        result['artifact_sha256'] = {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
        write(output / 'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = execute(args.input, args.output)
    print(result['classification'])
    return 0 if result['acceptance_met'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
