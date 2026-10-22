"""Real failed quantum evidence plus software negatives; no synthetic acceptance."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from island.forcefields.ambertools.provided import ProvidedChargeError
from island.forcefields.ambertools.psi4_execution import load_execution, reconstruct
from island.forcefields.ambertools.psi4_runtime import (
    Psi4RuntimeError,
    probe,
    validate_runtime,
)
from island.forcefields.ambertools.resp import file_hash
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input, validate_input
from scripts.phase_4j36_reference import check
from scripts.validate_phase_4j36 import BASE, WORKER, evaluate, signed
from scripts.validate_phase_4j36_receipt import validate_declarations, validate_matrix

EVIDENCE = Path('docs/evidence/phase_4j36')


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


@pytest.fixture(scope='module')
def native(tmp_path_factory):
    root = tmp_path_factory.mktemp('actual-native-quantum')
    with tarfile.open(EVIDENCE / 'artifacts/native-quantum.tar.gz') as tar:
        tar.extractall(root, filter='data')
    return root


@pytest.fixture
def case(native, tmp_path):
    root = tmp_path / 'native'
    shutil.copytree(native / 'quantum-linear', root)
    graph = build_fixture('linear', 'GAFF2').fixture.packed.graph
    runtime = json.loads((EVIDENCE / 'environment.json').read_text())
    kwargs = {'expected_input': build_input('linear'), 'runtime_identity': runtime['identity'], 'worker_sha256': file_hash(WORKER)}
    return root, graph, kwargs


def test_real_native_vectors_fail_unchanged(native):
    totals = {'linear': ['1.06E-16', '-2.65E-17'], 'network': ['3.307E-17']}
    for kind, expected_totals in totals.items():
        root = native / ('quantum-' + kind)
        graph = build_fixture(kind, 'GAFF2').fixture.packed.graph
        runtime = json.loads((EVIDENCE / 'environment.json').read_text())
        result = reconstruct(root, graph, expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))
        assert result['charge_generation_executed'] is True
        assert result['acceptance_met'] is False
        assert result['classification'] == 'charge_conservation_failed'
        assert [c['provided_total_e'] for c in result['components']] == expected_totals
        assert result['normalization'] is result['redistribution'] is False
        reference = check(root)
        assert reference['charge_provider'] is False
        assert all(c['fit_iteration_converged'] for c in reference['components'])
        assert all(c['independent_reconstruction_max_difference_e'] < 2e-11 for c in reference['components'])


def test_relocation_and_fresh_process(case, tmp_path):
    root, graph, kwargs = case
    record = reconstruct(root, graph, **kwargs)
    write(root / 'execution-record.json', record)
    target = tmp_path / 'elsewhere'
    shutil.copytree(root, target)
    assert load_execution(target, graph, expected_identity=record['identity'], **kwargs) == record
    code = '''import json,sys
from island.forcefields.ambertools.psi4_execution import load_execution
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input
r=load_execution(sys.argv[1],build_fixture('linear','GAFF2').fixture.packed.graph,expected_identity=sys.argv[2],expected_input=build_input('linear'),runtime_identity=sys.argv[3],worker_sha256=sys.argv[4]);print(r['identity'])'''
    child = subprocess.run([sys.executable, '-c', code, str(target), record['identity'], kwargs['runtime_identity'], kwargs['worker_sha256']], env={**os.environ, 'PYTHONPATH': '.:src'}, text=True, capture_output=True, check=False)
    assert child.returncode == 0, child.stderr
    assert child.stdout.strip() == record['identity']


@pytest.mark.parametrize('field,value', [('method', 'B3LYP'), ('basis', 'STO-3G'), ('final_graph_identity', 'stale'), ('charge_origin', 'external_oplsaa'), ('charge_method', 'oplsaa_provided'), ('quantum_software', 'OPLS-AA'), ('psi4_version', '1.8'), ('normalization', True), ('redistribution', True), ('acceptance_met', True)])
def test_resigned_native_metadata_rejected(case, field, value):
    root, graph, kwargs = case
    raw = json.loads((root / 'result.json').read_text())
    raw[field] = value
    write(root / 'result.json', raw)
    with pytest.raises(ProvidedChargeError):
        reconstruct(root, graph, **kwargs)


@pytest.mark.parametrize('change', ['order', 'round', 'truncate', 'boolean', 'geometry', 'options', 'source', 'esp', 'missing_inventory', 'wrong_mapping'])
def test_intermediates_rejected_even_after_hash_update(case, change):
    root, graph, kwargs = case
    raw = json.loads((root / 'result.json').read_text())
    path = root / 'component-0/charge-vector.json'
    rows = json.loads(path.read_text())
    if change == 'order':
        rows.reverse()
    elif change == 'round':
        rows[0]['charge_e'] = round(rows[0]['charge_e'], 3)
    elif change == 'truncate':
        rows.pop()
    elif change == 'boolean':
        rows[0]['charge_e'] = True
    elif change == 'wrong_mapping':
        rows[0]['atom_id'] = 9999
    elif change == 'geometry':
        path = root / 'component-0/optimized-geometry.json'
        geo = json.loads(path.read_text())
        geo['positions'][0][0] += 1
        write(path, geo)
    elif change == 'options':
        path = root / 'component-0/resp-options.json'
        options = json.loads(path.read_text())
        options['RESP_A'] = 0.5
        write(path, options)
    elif change == 'source':
        path = root / 'worker-source.py'
        path.write_text(path.read_text() + '# tampered\n')
    elif change == 'esp':
        path = root / 'component-0/1_j36_grid_esp.dat'
        path.write_text('changed')
        # Hash mutation must remain bound to the trusted execution identity.
        saved = json.loads((root / 'execution-record.json').read_text())
        with pytest.raises(ProvidedChargeError):
            load_execution(root, graph, expected_identity=saved['identity'], **kwargs)
        return
    else:
        del raw['artifact_sha256']['component-0/psi4.out']
    if change in ('order', 'round', 'truncate', 'boolean', 'wrong_mapping'):
        write(path, rows)
    if str(path.relative_to(root)) in raw['artifact_sha256']:
        raw['artifact_sha256'][str(path.relative_to(root))] = file_hash(path)
    write(root / 'result.json', raw)
    with pytest.raises(ProvidedChargeError):
        reconstruct(root, graph, **kwargs)


@pytest.mark.parametrize('change', ['basis', 'method', 'geometry', 'order', 'formal_charge', 'multiplicity', 'graph'])
def test_fresh_input_policy_rejects_mutation(change):
    data = build_input('linear')
    if change in ('basis', 'method'):
        data[change] = 'wrong'
    elif change == 'geometry':
        data['components'][0]['atoms'][0]['position_angstrom'][0] += 1
    elif change == 'order':
        data['components'][0]['atoms'].reverse()
    elif change in ('formal_charge', 'multiplicity'):
        data['components'][0][change] += 1
    else:
        data['graph']['sites'][0]['element'] = 'N'
    with pytest.raises(ProvidedChargeError):
        validate_input(data)


def test_missing_psi4_no_path_fallback(tmp_path):
    with pytest.raises(Psi4RuntimeError, match='Missing declared Psi4 Python') as failure:
        probe(tmp_path / 'absent')
    assert failure.value.classification == 'asset_unavailable'


@pytest.mark.parametrize('change', ['version', 'source', 'manifest', 'missing'])
def test_runtime_mutations_rejected(tmp_path, change):
    prefix = tmp_path / 'runtime'
    (prefix / 'bin').mkdir(parents=True)
    for name in ('python', 'psi4'):
        (prefix / 'bin' / name).write_text('software test only')
    body = {'schema': 'island.psi4-resp-runtime.v1', 'prefix': str(prefix), 'files': {f'bin/{n}': file_hash(prefix / 'bin' / n) for n in ('python', 'psi4')}, 'psi4_version': '1.9.1', 'resp_package_version': '1.0', 'resp_module_version': '0.8'}
    record = signed(body)
    expected = record['identity']
    if change == 'version':
        body['psi4_version'] = '0.0'
        record = signed(body)
        expected = record['identity']
    elif change == 'source':
        (prefix / 'bin/psi4').write_text('changed')
    elif change == 'manifest':
        record['prefix'] = 'elsewhere'
    else:
        (prefix / 'bin/psi4').unlink()
    with pytest.raises(Psi4RuntimeError):
        validate_runtime(record, expected_identity=expected, execute_probe=False)


def test_cli_missing_runtime_preserves_nonzero_diagnostic(tmp_path):
    record = json.loads((EVIDENCE / 'environment.json').read_text())
    record['prefix'] = str(tmp_path / 'missing')
    record = signed({k: v for k, v in record.items() if k != 'identity'})
    manifest = tmp_path / 'manifest.json'
    write(manifest, record)
    output = tmp_path / 'output'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j36', '--runtime', str(manifest), '--runtime-identity', record['identity'], '--output', str(output)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src'})
    assert child.returncode == 1, child.stderr
    diagnostic = json.loads((output / 'result.json').read_text())
    assert diagnostic['classification'] == 'asset_unavailable'
    assert diagnostic['acceptance_met'] is False and diagnostic['cases'] == {}


@pytest.mark.parametrize('change', ['acceptance', 'missing', 'dynamics', 'origin', 'case_swap', 'readiness'])
def test_resigned_matrix_cannot_claim_acceptance(native, tmp_path, change):
    runtime = json.loads((EVIDENCE / 'environment.json').read_text())
    rows = {}
    for kind in ('linear', 'network'):
        shutil.copytree(native / ('quantum-' + kind), tmp_path / kind)
        rows[kind] = evaluate(tmp_path / kind, kind, runtime)
    body = {'schema': 'island.j36-matrix.v1', 'base_commit': BASE, 'runtime_identity': runtime['identity'], 'runtime_validation': 'complete', 'acceptance_met': False, 'production_validated': False, 'simulation_readiness': 'not_established', 'native_am1bcc_success': False, 'external_oplsaa_is_resp': False, 'cases': rows}
    assert validate_matrix(signed(body), tmp_path, runtime)['acceptance_met'] is False
    mutated = deepcopy(body)
    if change == 'acceptance':
        mutated['acceptance_met'] = True
    elif change == 'missing':
        del mutated['cases']['linear']
    elif change == 'dynamics':
        mutated['cases']['linear']['dynamics_executed'] = True
    elif change == 'origin':
        mutated['external_oplsaa_is_resp'] = True
    elif change == 'readiness':
        mutated['production_validated'] = True
    else:
        mutated['cases']['network'] = deepcopy(mutated['cases']['linear'])
    with pytest.raises(ProvidedChargeError):
        validate_matrix(signed(mutated), tmp_path, runtime)


def test_truncated_native_array_rejected(case):
    root, graph, kwargs = case
    path = root / 'component-0/native-resp-charges.npy'
    np.save(path, np.load(path, allow_pickle=False)[:-1], allow_pickle=False)
    raw = json.loads((root / 'result.json').read_text())
    raw['artifact_sha256'][str(path.relative_to(root))] = file_hash(path)
    write(root / 'result.json', raw)
    with pytest.raises(ProvidedChargeError, match='ndarray'):
        reconstruct(root, graph, **kwargs)


@pytest.mark.parametrize('change', ['origin', 'method', 'basis', 'case', 'dynamics'])
def test_receipt_declarations_bound_to_execution(native, tmp_path, change):
    runtime = json.loads((EVIDENCE / 'environment.json').read_text())
    rows, executions = {}, {}
    for kind in ('linear', 'network'):
        shutil.copytree(native / ('quantum-' + kind), tmp_path / kind)
        rows[kind] = evaluate(tmp_path / kind, kind, runtime)
        executions[kind] = json.loads((tmp_path / kind / 'execution-record.json').read_text())
    receipt = {'provider': {'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
                            'psi4_version': '1.9.1', 'resp_package_version': '1.0', 'resp_module_version': '0.8',
                            'method': 'HF', 'basis': '6-31G*', 'optimization_method': 'HF', 'optimization_basis': 'STO-3G', 'fit_stages': 1},
               'cases': {kind: {**row, 'final_graph_identity': executions[kind]['final_graph_identity'], 'components': executions[kind]['components']} for kind, row in rows.items()},
               'not_executed': ['RESP-based GAFF2 preparation', 'minimization', 'NVT', 'NPT', 'split/uninterrupted MD comparison', 'MD checkpoint relocation and fresh-process continuation']}
    validate_declarations(receipt, {'cases': rows}, executions)
    receipt = deepcopy(receipt)
    if change == 'origin':
        receipt['provider']['charge_origin'] = 'external_oplsaa'
    elif change in ('method', 'basis'):
        receipt['provider'][change] = 'wrong'
    elif change == 'case':
        receipt['cases']['linear']['components'][0]['provided_total_e'] = '0'
    else:
        receipt['not_executed'].remove('NPT')
    with pytest.raises(ProvidedChargeError):
        validate_declarations(receipt, {'cases': rows}, executions)
