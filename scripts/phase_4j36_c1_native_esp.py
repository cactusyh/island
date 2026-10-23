"""Independent native Psi4 ESP -> CDS RESP execution, without ISLAND imports.

Uses retained optimized geometry, recomputes HF/6-31G* wavefunctions and ESP,
and fits unchanged native float64 arrays. No normalization or redistribution.
"""

import argparse
import hashlib
import json
import os
import traceback
from fractions import Fraction
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    import numpy as np
    import psi4
    import resp
    from resp import espfit, vdw_surface
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    original, output = args.native.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    data = json.loads((original / 'input.json').read_text())
    write(output / 'input.json', data)
    result = {'schema': 'island.native-esp-resp.v1', 'charge_origin': 'resp_psi4', 'charge_method': 'RESP',
              'psi4_version': psi4.__version__, 'resp_module_version': resp.__version__,
              'method': data['method'], 'basis': data['basis'], 'final_graph_identity': data['final_graph_identity'],
              'input_sha256': sha(original / 'input.json'), 'normalization': False, 'redistribution': False,
              'charge_conservation_passed': False, 'acceptance_met': False, 'components': [],
              'geometry_optimization': 'retained J36 native optimized geometry, not reoptimized',
              'stage': 'native_esp', 'classification': 'execution_failure', 'dynamics_executed': False}
    old = Path.cwd()
    psi4.set_num_threads(4)
    psi4.set_memory('2 GB')
    try:
        for index, component in enumerate(data['components']):
            root = output / f'component-{index}'
            root.mkdir()
            os.chdir(root)
            psi4.core.clean()
            psi4.core.clean_options()
            psi4.core.set_output_file(str(root / 'psi4.out'), False)
            atoms = component['atoms']
            geometry = np.load(original / f'component-{index}/optimized-geometry-bohr.npy', allow_pickle=False)
            np.save('geometry-bohr.npy', geometry, allow_pickle=False)
            text = f"{component['formal_charge']} {component['multiplicity']}\n" + '\n'.join(f"{a['element']} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g}" for a, r in zip(atoms, geometry, strict=True)) + '\nunits bohr\nsymmetry c1\nno_com\nno_reorient\n'
            Path('geometry-input.dat').write_text(text)
            molecule = psi4.geometry(text)
            psi4.set_options({'basis': data['basis'], **data['resp_options']['PSI4_OPTIONS'], 'reference': 'rhf'})
            energy, wfn = psi4.energy(data['method'], molecule=molecule, return_wfn=True)
            if not np.array_equal(molecule.geometry().np, geometry):
                raise ValueError('Native molecule geometry changed on import')
            symbols = [molecule.symbol(i) for i in range(len(atoms))]
            if symbols != [a['element'] for a in atoms]:
                raise ValueError('Native molecule atom elements/order changed')
            conversion = 0.52917721092  # Pinned CDS RESP driver constant.
            angstrom = geometry * conversion
            shells = [vdw_surface(angstrom, symbols, scale, data['resp_options']['VDW_POINT_DENSITY'], {})[0] for scale in data['resp_options']['VDW_SCALE_FACTORS']]
            points = np.concatenate(shells)
            points_bohr = points / conversion
            esp = psi4.core.ESPPropCalc(wfn).compute_esp_over_grid_in_memory(psi4.core.Matrix.from_array(points_bohr)).np.copy()
            np.save('grid-bohr.npy', points_bohr, allow_pickle=False)
            np.save('esp-native.npy', esp, allow_pickle=False)
            # Independent official output path verifies the API's units/semantics.
            # Its rounded values are never used in fitting or as charges.
            np.savetxt('grid.dat', points_bohr, fmt='%.17g')
            psi4.oeprop(wfn, 'GRID_ESP')
            rendered = np.loadtxt('grid_esp.dat')
            difference = float(np.max(abs(esp - rendered)))
            if difference > 5.1e-11:
                raise ValueError('Native ESP API does not match official GRID_ESP units/semantics')
            inverse = np.empty((len(points), len(atoms)))
            for i in range(len(points)):
                for j in range(len(atoms)):
                    inverse[i, j] = 1 / np.linalg.norm(points[i] - angstrom[j])
            inverse *= conversion
            options = {**data['resp_options'], 'CONSTRAINT_CHARGE': [], 'CONSTRAINT_GROUP': [], 'WEIGHT': [1], 'RESTRAINT': True}
            write('fit-options.json', options)
            np.savez('fit-input.npz', invr=inverse, esp=esp)
            result['stage'] = 'native_resp_fit'
            fitted, labels, note = espfit.fit(options, {'natoms': len(atoms), 'symbols': symbols, 'mol_charge': component['formal_charge'], 'invr': [inverse], 'esp_values': [esp]})
            if 'did not converge' in note:
                raise ValueError(note)
            charges = fitted[-1]
            if charges.dtype != np.float64 or not np.isfinite(charges).all():
                raise ValueError('Invalid native fit output')
            np.save('charges-native.npy', charges, allow_pickle=False)
            rows = [{'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': float(q), 'ieee754_hex': float(q).hex()} for a, q in zip(atoms, charges, strict=True)]
            write('charges-native.json', rows)
            total = sum((Fraction.from_float(float(q)) for q in charges), Fraction(0))
            result['components'].append({'atom_ids': [a['atom_id'] for a in atoms], 'elements': symbols,
                                         'geometry_sha256': sha('geometry-bohr.npy'), 'source_geometry_sha256': sha(original / f'component-{index}/optimized-geometry-bohr.npy'),
                                         'charge_array_sha256': sha('charges-native.npy'), 'formal_charge': component['formal_charge'], 'multiplicity': component['multiplicity'],
                                         'total_fraction': [total.numerator, total.denominator], 'exact_conservation': total == component['formal_charge'],
                                         'native_esp_vs_official_text_max_difference_au': difference,
                                         'scf_energy_hartree': float(energy), 'fit_labels': labels, 'fit_note': note})
        result['stage'] = 'native_component_and_system_total_validation'
        result['charge_conservation_passed'] = all(c['exact_conservation'] for c in result['components'])
        result['classification'] = 'native_charge_totals_exact' if result['charge_conservation_passed'] else 'authoritative_native_precision_insufficient'
    except Exception as error:  # noqa: BLE001 -- preserve failed native calculations
        result.update(diagnostic=str(error), exception_type=type(error).__name__)
        (output / 'traceback.txt').write_text(traceback.format_exc())
    finally:
        psi4.core.flush_outfile()
        os.chdir(old)
        result['artifacts'] = {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
        write(output / 'result.json', result)
    print(result['classification'])
    return 0 if result['charge_conservation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
