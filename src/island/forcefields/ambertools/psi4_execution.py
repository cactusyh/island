"""Path-independent native Psi4/RESP execution evidence and exact charge gate.

A readable failed execution record is not an accepted provided-charge vector.
Native arrays, not printed tables, remain authoritative, with no charge repair.
"""

import json
from decimal import Decimal
from pathlib import Path

import numpy as np

from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.resp import bound_file, file_hash

SCHEMA = 'island.psi4-resp-execution.v1'


def signed(body):
    return {**body, 'identity': reference_digest(body)}


def reconstruct(directory, graph, *, expected_input, runtime_identity, worker_sha256):
    directory = Path(directory)
    data = json.loads((directory / 'input.json').read_text())
    require(reference_digest(data) == reference_digest(expected_input), 'Stale geometry/mapping, method, basis, formal charge or graph input')
    require(data['final_graph_identity'] == graph.identity and data['graph'] == graph.payload, 'Wrong quantum final graph')
    raw = json.loads((directory / 'result.json').read_text())
    require(raw.get('schema') == 'island.j36-native-quantum-result.v1', 'Unsupported native result schema')
    require(raw.get('charge_origin') == 'resp_psi4' and raw.get('charge_method') == 'RESP' and raw.get('quantum_software') == 'Psi4', 'External OPLS-AA is not native Psi4/RESP')
    require(raw['input_sha256'] == file_hash(directory / 'input.json'), 'Stale quantum input file')
    require(raw['final_graph_identity'] == graph.identity and raw['method'] == data['method'] and raw['basis'] == data['basis'], 'Stale graph or changed quantum method/basis')
    require(raw['normalization'] is False and raw['redistribution'] is False, 'Charge repair prohibited')
    require(raw['psi4_version'] == '1.9.1' and raw['resp_module_version'] == '0.8', 'Wrong Psi4/RESP version')
    require(file_hash(directory / 'worker-source.py') == worker_sha256, 'Changed executed quantum provider source')
    required = {'input.json'}
    if raw.get('charge_generation_executed'):
        for i in range(len(data['components'])):
            required.update(f'component-{i}/{name}' for name in ('geometry-input.dat', 'optimized-geometry-bohr.npy', 'optimized-geometry.json', 'resp-options.json', 'native-resp-charges.npy', 'native-esp-charges.npy', 'charge-vector.json', 'psi4.out', 'results.out', '1_j36_grid.dat', '1_j36_grid_esp.dat'))
    require(required <= set(raw['artifact_sha256']), 'Missing quantum evidence inventory')
    for name, expected in raw['artifact_sha256'].items():
        bound_file(directory, name, expected)
    artifacts = {**raw['artifact_sha256'], 'result.json': file_hash(directory / 'result.json'), 'worker-source.py': worker_sha256}
    base = {'schema': SCHEMA, 'final_graph_identity': graph.identity, 'runtime_identity': runtime_identity,
            'input_identity': reference_digest(data), 'native_result_identity': reference_digest(raw),
            'worker_sha256': worker_sha256, 'artifact_sha256': artifacts,
            'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
            'psi4_version': raw['psi4_version'], 'resp_module_version': raw['resp_module_version'],
            'method': data['method'], 'basis': data['basis'], 'optimization_method': data['optimization_method'],
            'optimization_basis': data['optimization_basis'], 'geometry_optimization': data['geometry_optimization'],
            'normalization': False, 'redistribution': False, 'native_am1bcc_success': False,
            'production_validated': False, 'simulation_readiness': 'not_established',
            'acceptance_met': False, 'charge_generation_executed': raw.get('charge_generation_executed', False)}
    if not base['charge_generation_executed']:
        require(raw['acceptance_met'] is False and raw['classification'] == 'quantum_execution_failure', 'Contradictory incomplete quantum result')
        return signed({**base, 'status': 'failed', 'classification': raw['classification'], 'stage': raw['stage'], 'diagnostic': raw['diagnostic']})
    rows, components, optimized = [], [], []
    require(len(raw['components']) == len(data['components']), 'Incomplete native component inventory')
    for index, source in enumerate(data['components']):
        root = directory / f'component-{index}'
        require(json.loads((root / 'resp-options.json').read_text()) == data['resp_options'], 'Changed native fitting policy')
        atoms = source['atoms']
        ids = [a['atom_id'] for a in atoms]
        require(all(type(i) is int for i in ids) and ids == sorted(ids), 'Reordered/boolean atom mapping')
        values = np.load(root / 'native-resp-charges.npy', allow_pickle=False)
        geometry = np.load(root / 'optimized-geometry-bohr.npy', allow_pickle=False)
        require(values.dtype == np.float64 and values.shape == (len(ids),) and np.isfinite(values).all(), 'Invalid native RESP ndarray')
        require(geometry.dtype == np.float64 and geometry.shape == (len(ids), 3) and np.isfinite(geometry).all(), 'Invalid optimized geometry')
        geo = json.loads((root / 'optimized-geometry.json').read_text())
        require(geo['atom_ids'] == ids and geo['elements'] == [a['element'] for a in atoms] and geo['units'] == 'bohr' and geo['positions'] == geometry.tolist(), 'Stale optimized geometry/atom mapping')
        current = [{'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': float(q), 'decimal_e': repr(float(q))} for a, q in zip(atoms, values, strict=True)]
        require(current == json.loads((root / 'charge-vector.json').read_text()), 'Rounded/truncated/reordered charge file differs from native ndarray')
        total = sum((Decimal(row['decimal_e']) for row in current), Decimal(0))
        comp = {'atom_ids': ids, 'formal_total_e': source['formal_charge'], 'provided_total_e': str(total),
                'exact_charge_conservation': total == source['formal_charge'],
                'optimized_geometry_sha256': file_hash(root / 'optimized-geometry-bohr.npy'),
                'native_charge_array_sha256': file_hash(root / 'native-resp-charges.npy')}
        require(comp == raw['components'][index], 'Changed component totals or geometry/charge identity')
        require('Charge fitting did not converge' not in (root / 'results.out').read_text(), 'Unconverged native RESP fit')
        rows.extend(current)
        components.append(comp)
        optimized.append({'atom_ids': ids, 'geometry_sha256': comp['optimized_geometry_sha256']})
    rows.sort(key=lambda r: r['atom_id'])
    require(rows == raw['rows'] and [r['atom_id'] for r in rows] == [s['id'] for s in graph.payload['sites']], 'Wrong native row coverage/order')
    exact = all(c['exact_charge_conservation'] for c in components)
    require(raw['acceptance_met'] is exact, 'Contradictory native exact-charge outcome')
    classification = 'native_resp_charge_generation_complete' if exact else 'charge_conservation_failed'
    require(raw['classification'] == classification, 'Wrong native failure classification')
    return signed({**base, 'status': 'charges_validated' if exact else 'failed', 'classification': classification,
                   'stage': 'component_total_validation', 'acceptance_met': exact, 'components': components,
                   'rows': rows, 'units': 'elementary_charge', 'atom_mapping': [{'atom_id': r['atom_id'], 'element': r['element']} for r in rows],
                   'optimized_geometry_identity': reference_digest(optimized)})


def load_execution(directory, graph, *, expected_identity, expected_input, runtime_identity, worker_sha256):
    saved = json.loads((Path(directory) / 'execution-record.json').read_text())
    rebuilt = reconstruct(directory, graph, expected_input=expected_input, runtime_identity=runtime_identity, worker_sha256=worker_sha256)
    require(saved == rebuilt and rebuilt['identity'] == expected_identity, 'Changed Psi4/RESP execution record or provenance')
    return rebuilt
