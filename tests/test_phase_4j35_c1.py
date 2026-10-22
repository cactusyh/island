"""Software contract fixtures only: none of these vectors are quantum evidence."""

import json
import os
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest

from island.forcefields.ambertools.provided import ProvidedChargeError
from island.forcefields.ambertools.resp import RESPCharges, file_hash
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j35_c1 import signed


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False) + '\n')


@pytest.fixture
def evidence(tmp_path):
    """Synthetic serialization fixture. Reviewer trust is supplied only by test."""
    graph = build_fixture('linear', 'GAFF2').fixture.packed.graph
    ref = json.loads(Path('docs/evidence/phase_4j35/references/linear.json').read_text())
    roles = ('psi4_input', 'psi4_output', 'esp_grid', 'resp_input', 'resp_output', 'charge_vector', 'psi4_executable', 'resp_implementation_source', 'environment_manifest')
    assets = {}
    for role in roles:
        path = tmp_path / (role + '.json')
        write(path, ref['rows'] if role == 'charge_vector' else {'software_fixture': role})
        assets[role] = {'path': path.name, 'sha256': file_hash(path)}
    ref.update(method='RESP/Psi4', provenance='SYNTHETIC SOFTWARE TEST ONLY',
               source_sha256=assets['psi4_output']['sha256'], generator_sha256=assets['resp_implementation_source']['sha256'])
    ref = signed({k: v for k, v in ref.items() if k != 'identity'})
    p = {'schema': 'island.resp-psi4-result.v1', 'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4', 'psi4_version': 'synthetic-1', 'method': 'HF', 'basis': '6-31G*', 'geometry_optimization': {'executed': True, 'settings': {'maxiter': 100}, 'input_geometry_identity': 'test-only'}, 'formal_charge': 0, 'multiplicity': 1, 'atom_mapping': [{'atom_id': s['id'], 'element': s['element']} for s in graph.payload['sites']], 'normalization': False, 'redistribution': False, 'resp_implementation': 'test-only', 'resp_version': 'synthetic-1', 'provenance': ref['provenance'], 'final_graph_identity': graph.identity, 'evidence_references': ['synthetic-test-only'], 'reference': ref, 'native_am1bcc_success': False}
    review = {'schema': 'island.resp-independent-verification.v1', 'actual_quantum_and_fit_executed': True, 'source_kind': 'Psi4_ESP_RESP_fit', 'reviewer': 'SYNTHETIC TEST NOT SCIENTIFIC EVIDENCE', 'evidence_references': p['evidence_references'], 'reference_identity': ref['identity'], 'settings': {k: p[k] for k in ('psi4_version', 'method', 'basis', 'geometry_optimization', 'formal_charge', 'multiplicity', 'atom_mapping', 'resp_implementation', 'resp_version')}, 'asset_hashes': {k: v['sha256'] for k, v in assets.items()}}
    write(tmp_path / 'review.json', review)
    assets['independent_verification'] = {'path': 'review.json', 'sha256': file_hash(tmp_path / 'review.json')}
    p = signed({**p, 'assets': assets})
    write(tmp_path / 'resp-result.json', p)
    return tmp_path, graph, p


def test_synthetic_roundtrip_and_relocation(evidence, tmp_path):
    import shutil
    root, graph, p = evidence
    record = RESPCharges.load(root, graph, expected_identity=p['identity'])
    provided = record.provided(root, graph, expected_identity=p['identity'])
    assert provided.payload['rows'] == p['reference']['rows']
    relocated = tmp_path / 'relocated'
    relocated.mkdir()
    for path in root.iterdir():
        if path.is_file():
            shutil.copy2(path, relocated / path.name)
    assert RESPCharges.load(relocated, graph, expected_identity=p['identity']).payload == record.payload
    code = 'from island.forcefields.ambertools.resp import RESPCharges; from scripts.phase_4j34_real_fixture import build_fixture; import sys; r=RESPCharges.load(sys.argv[1],build_fixture("linear","GAFF2").fixture.packed.graph,expected_identity=sys.argv[2]); print(r.payload["identity"])'
    child = subprocess.run([sys.executable, '-c', code, str(relocated), p['identity']], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == 0, child.stderr
    assert child.stdout.strip() == p['identity']


@pytest.mark.parametrize('field,value', [('charge_origin', 'external_oplsaa'), ('charge_method', 'oplsaa_provided'), ('charge_method', 'AM1-BCC'), ('quantum_software', 'OPLS-AA'), ('psi4_version', 'changed'), ('basis', 'changed'), ('method', 'changed'), ('final_graph_identity', 'stale'), ('provenance', 'changed'), ('evidence_references', []), ('normalization', True), ('multiplicity', True), ('formal_charge', 1)])
def test_resigned_metadata_rejected(evidence, field, value):
    root, graph, p = evidence
    mutated = deepcopy(p)
    mutated[field] = value
    mutated = signed({k: v for k, v in mutated.items() if k != 'identity'})
    write(root / 'resp-result.json', mutated)
    # External trusted identity cannot be replaced by a self-signed assertion.
    with pytest.raises(ProvidedChargeError):
        RESPCharges.load(root, graph, expected_identity=p['identity'])
    # Even a caller replacing that anchor cannot skip independent review.
    with pytest.raises(ProvidedChargeError):
        RESPCharges.load(root, graph, expected_identity=mutated['identity'])


@pytest.mark.parametrize('role', ['psi4_input', 'psi4_output', 'esp_grid', 'resp_input', 'resp_output', 'charge_vector', 'psi4_executable', 'resp_implementation_source', 'environment_manifest', 'independent_verification'])
def test_mutated_intermediates_rejected(evidence, role):
    root, graph, p = evidence
    path = root / p['assets'][role]['path']
    path.write_text(path.read_text() + 'mutation')
    with pytest.raises(ProvidedChargeError):
        RESPCharges.load(root, graph, expected_identity=p['identity'])


@pytest.mark.parametrize('change', ['reorder', 'round', 'truncate', 'boolean', 'nonfinite', 'stale_mapping'])
def test_vector_and_mapping_rejected(evidence, change):
    root, graph, p = evidence
    rows = deepcopy(p['reference']['rows'])
    if change == 'reorder':
        rows.reverse()
    elif change == 'truncate':
        rows.pop()
    elif change == 'boolean':
        rows[0]['charge_e'] = True
    elif change == 'nonfinite':
        rows[0]['charge_e'] = 'NaN'
    elif change == 'round':
        rows[0]['charge_e'] = round(rows[0]['charge_e'], 1)
    else:
        p['atom_mapping'].reverse()
    write(root / p['assets']['charge_vector']['path'], rows)
    p['assets']['charge_vector']['sha256'] = file_hash(root / p['assets']['charge_vector']['path'])
    p = signed({k: v for k, v in p.items() if k != 'identity'})
    write(root / 'resp-result.json', p)
    with pytest.raises(ProvidedChargeError):
        RESPCharges.load(root, graph, expected_identity=p['identity'])


def test_real_opls_reference_cannot_enter_resp(tmp_path):
    ref = json.loads(Path('docs/evidence/phase_4j35/references/linear.json').read_text())
    write(tmp_path / 'resp-result.json', {**ref, 'schema': 'island.resp-psi4-result.v1', 'charge_origin': 'external_oplsaa'})
    graph = build_fixture('linear', 'GAFF2').fixture.packed.graph
    with pytest.raises(ProvidedChargeError, match='Only actual RESP'):
        RESPCharges.load(tmp_path, graph, expected_identity=ref['identity'])


def test_acceptance_cli_preserves_blocked_diagnostic(tmp_path):
    out = tmp_path / 'new'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j35_c1', '--output', str(out)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src'})
    assert child.returncode == 1, child.stderr
    result = json.loads((out / 'result.json').read_text())
    assert result['acceptance_met'] is False
    assert all(not row['dynamics_executed'] for row in result['cases'].values())
    assert all(row['charge_origin'] == 'external_oplsaa' for row in result['external_controls'].values())
    assert 'BLOCKED' in child.stdout


@pytest.fixture(scope='module')
def blocked_diagnostic():
    from scripts.validate_phase_4j35_c1 import diagnose
    return diagnose()


@pytest.mark.parametrize('mutation', ['opls_as_resp', 'success', 'dynamics', 'native', 'controls_missing', 'case_missing', 'audit', 'readiness'])
def test_resigned_diagnostic_mutations_rejected(blocked_diagnostic, mutation):
    from scripts.validate_phase_4j35_c1_receipt import validate_diagnostic
    p = deepcopy(blocked_diagnostic)
    if mutation == 'opls_as_resp':
        p['external_controls']['linear']['charge_method'] = 'RESP'
        p['external_controls']['linear']['charge_origin'] = 'resp_psi4'
    elif mutation == 'success':
        p['acceptance_met'] = True
    elif mutation == 'dynamics':
        p['cases']['linear']['dynamics_executed'] = True
    elif mutation == 'native':
        p['native_am1bcc']['acceptance_met'] = True
    elif mutation == 'controls_missing':
        p['external_controls'].pop('network')
    elif mutation == 'case_missing':
        p['cases'].pop('linear')
    elif mutation == 'audit':
        p['environment_audit'] = []
    else:
        p['production_validated'] = True
    p = signed({k: v for k, v in p.items() if k != 'identity'})
    with pytest.raises(ProvidedChargeError):
        validate_diagnostic(p)


def test_diagnostic_readability_is_not_acceptance(blocked_diagnostic):
    from scripts.validate_phase_4j35_c1_receipt import validate_diagnostic
    validate_diagnostic(blocked_diagnostic)
    assert blocked_diagnostic['acceptance_met'] is False


@pytest.fixture
def software_receipt(tmp_path, blocked_diagnostic):
    """Minimal software receipt fixture, never archived as acceptance evidence."""
    from scripts.validate_phase_4j35_c1 import BASE
    root = tmp_path / 'software-receipt'
    root.mkdir()
    write(root / 'result.json', blocked_diagnostic)
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json'], text=True).splitlines()
    write(root / 'historical-hashes.json', {n: file_hash(n) for n in names})
    impl = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    write(root / 'command.json', {'exit_code': 0, 'implementation_commit': impl})
    (root / 'software.xml').write_text('<testsuites><testsuite tests="1" failures="0" errors="0" skipped="0"><testcase name="synthetic_schema_only"/></testsuite></testsuites>')
    from scripts.validate_phase_4j34_c2_receipt import junit_counts
    files = ['src/island/forcefields/ambertools/resp.py', 'scripts/validate_phase_4j35_c1.py', 'scripts/validate_phase_4j35_c1_receipt.py', 'tests/test_phase_4j35_c1.py']
    p = {'schema': 'island.j35-c1-receipt.v1', 'base_commit': BASE, 'classification': 'asset_unavailable', 'acceptance_met': False, 'production_validated': False, 'simulation_readiness': 'not_established', 'implementation_commit': impl, 'implementation_files': {n: file_hash(n) for n in files}, 'artifacts': {path.name: file_hash(path) for path in root.iterdir()}, 'tests': {n: {'command': 'command.json', 'junit': 'software.xml', 'counts': junit_counts(root / 'software.xml')} for n in ('focused', 'j22-j35', 'full_supported')}}
    write(root / 'receipt.json', signed(p))
    return root


@pytest.mark.parametrize('diagnostic_only,expected_exit', [(True, 0), (False, 1)])
def test_receipt_cli_does_not_claim_acceptance(software_receipt, diagnostic_only, expected_exit):
    args = [sys.executable, '-m', 'scripts.validate_phase_4j35_c1_receipt', str(software_receipt / 'receipt.json')]
    if diagnostic_only:
        args.append('--diagnostic-only')
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    assert result.returncode == expected_exit, result.stdout + result.stderr
    assert 'VALID BLOCKED DIAGNOSTIC' in result.stdout


def test_full_receipt_rejects_resigned_opls_as_resp(software_receipt):
    root = software_receipt
    result = json.loads((root / 'result.json').read_text())
    result['external_controls']['linear'].update(charge_origin='resp_psi4', charge_method='RESP', resp_psi4=True)
    write(root / 'result.json', signed({k: v for k, v in result.items() if k != 'identity'}))
    receipt = json.loads((root / 'receipt.json').read_text())
    receipt['artifacts']['result.json'] = file_hash(root / 'result.json')
    write(root / 'receipt.json', signed({k: v for k, v in receipt.items() if k != 'identity'}))
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j35_c1_receipt', str(root / 'receipt.json'), '--diagnostic-only'], capture_output=True, text=True, check=False)
    assert child.returncode == 1
    assert 'External OPLS-AA is not RESP/Psi4' in child.stdout


def test_full_receipt_rejects_top_level_false_resp_claim(software_receipt):
    root = software_receipt
    receipt = json.loads((root / 'receipt.json').read_text())
    receipt.update(charge_origin='external_oplsaa', charge_method='RESP', quantum_software='Psi4')
    write(root / 'receipt.json', signed({k: v for k, v in receipt.items() if k != 'identity'}))
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j35_c1_receipt', str(root / 'receipt.json'), '--diagnostic-only'], capture_output=True, text=True, check=False)
    assert child.returncode == 1
    assert 'cannot claim a generated quantum charge' in child.stdout
