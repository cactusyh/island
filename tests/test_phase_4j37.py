"""Raw external RESP provenance mutation controls; never scientific evidence."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.exceptions import AnalysisError
from scripts.phase_4j37_reference import audit, digest, read, residual
from scripts.validate_phase_4j37 import validate_bundle


@pytest.fixture(scope='module')
def raw(tmp_path_factory):
    root = tmp_path_factory.mktemp('j37-raw')
    for name in ('docs/evidence/phase_4j36/artifacts/native-quantum.tar.gz', 'docs/evidence/phase_4j36_c1/artifacts/native-memory.tar.gz'):
        with tarfile.open(name) as archive:
            archive.extractall(root, filter='data')
    return root


def inputs(kind):
    return read(f'docs/evidence/phase_4j36_c2/memory-{kind}-charges.json'), read('docs/evidence/phase_4j36/environment.json')


@pytest.mark.parametrize('kind', ['linear', 'network'])
def test_independent_raw_provenance_and_residual(raw, kind):
    charges, runtime = inputs(kind)
    result = audit(raw, kind, charges, runtime)
    assert result['raw_provenance_verified'] and result['import_bit_exact']
    assert result['system']['residual_fraction'] == charges['policy']['system']['residual_fraction']
    assert result['system']['absolute_residual_e'] > 0


@pytest.mark.parametrize('mutation', ['vector', 'order', 'component', 'source', 'executable', 'origin', 'policy', 'policy_text', 'basis', 'boolean'])
def test_consistently_resigned_charge_mutations(raw, mutation):
    charges, runtime = inputs('linear')
    if mutation == 'vector':charges['rows'][0]['charge_e'] += 1e-6
    elif mutation == 'order':charges['rows'].reverse()
    elif mutation == 'component':charges['source_audit']['components'][0]['formal_total_e'] = 1
    elif mutation == 'source':runtime['files']['lib/python3.11/site-packages/resp/espfit.py'] = '0' * 64
    elif mutation == 'executable':runtime['files']['bin/psi4'] = '0' * 64
    elif mutation == 'origin':charges['charge_origin'] = 'external_oplsaa'
    elif mutation == 'policy':charges['policy']['components'][0]['bound_fraction'] = [1, 1]
    elif mutation == 'policy_text':charges['policy']['justification'] = 'arbitrary tolerance'
    elif mutation == 'basis':charges['source_audit']['basis'] = 'STO-3G'
    else:charges['rows'][0]['charge_e'] = True
    # Re-sign all nested identities: self-consistency is not provenance.
    def resign(value):
        if isinstance(value, list):return [resign(v) for v in value]
        if isinstance(value, dict):
            value = {k: resign(v) for k, v in value.items()}
            if 'identity' in value:value['identity'] = digest({k: v for k, v in value.items() if k != 'identity'})
        return value
    charges, runtime = resign(charges), resign(runtime)
    with pytest.raises(ValueError):
        audit(raw, 'linear', charges, runtime)


@pytest.mark.parametrize('filename', ['charges-native.npy', 'geometry-input.dat', 'fit-input.npz', 'psi4.out'])
def test_mutated_raw_artifacts(raw, tmp_path, filename):
    copy = tmp_path / 'raw';shutil.copytree(raw, copy)
    path = copy / 'native-esp-linear/component-0' / filename
    path.write_bytes(path.read_bytes()[:-1])
    charges, runtime = inputs('linear')
    with pytest.raises(ValueError, match='Mutated raw'):
        audit(copy, 'linear', charges, runtime)


def test_compensated_residual_is_not_charge_repair():
    source = [.1, -.1 + .003]
    saved = deepcopy(source)
    assert not residual(source, 0)['accepted']
    assert source == saved


def test_missing_raw_cli_fails_with_diagnostic(tmp_path):
    out = tmp_path / 'diagnostic.json'
    argv = [sys.executable, '-m', 'scripts.phase_4j37_reference', '--root', str(tmp_path / 'absent'), '--fixture', 'linear', '--charges', 'docs/evidence/phase_4j36_c2/memory-linear-charges.json', '--runtime', 'docs/evidence/phase_4j36/environment.json', '--output', str(out)]
    child = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert child.returncode != 0 and read(out)['status'] == 'blocked'
    assert read(out)['diagnostic']


def test_audit_inventory_cannot_claim_acceptance_and_relocates(tmp_path):
    output = tmp_path / 'inventory'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j37', '--output', str(output), '--audit-only'], capture_output=True, text=True, check=False)
    assert child.returncode == 1
    matrix = read(output / 'result.json')
    assert matrix['reviews'].keys() == {'linear', 'network'}
    assert not matrix['acceptance_met'] and matrix['cases'] == {}
    assert matrix['provenance_relocation'] and matrix['provenance_fresh_process']
    assert validate_bundle(output, 'linear') == validate_bundle(output / 'relocated-provenance', 'linear')
    path = output / 'linear-review.json';review = read(path);review['system']['residual_fraction'] = [0, 1]
    review['identity'] = digest({k: v for k, v in review.items() if k != 'identity'});path.write_text(json.dumps(review))
    with pytest.raises(ValueError, match='changed'):
        validate_bundle(output, 'linear')


def test_real_fresh_gaff2_matrix_subprocess(tmp_path):
    python = Path('../island-validation/phase4g1-env/bin/python')
    assert python.is_file(), 'Required declared execution environment missing'
    output = tmp_path / 'actual-molecular-workflows'
    env = {**os.environ, 'PYTHONPATH': '.:src', 'OPENMM_CPU_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    child = subprocess.run([str(python), '-m', 'scripts.validate_phase_4j37', '--output', str(output)], env=env, capture_output=True, text=True, check=False)
    assert child.returncode == 0, child.stdout + child.stderr + (output / 'result.json').read_text()
    result = read(output / 'result.json')
    assert result['acceptance_met']
    for kind in ('linear', 'network'):
        row = result['cases'][kind]
        assert row['comparison_passed'] and row['fresh_process_continuation'] and row['periodic_reference_passed']
        assert read(output / (kind + '-refit.json'))['passed']
        assert read(output / (kind + '-physics.json'))['passed']


@pytest.fixture
def synthetic_contract():
    """Schema controls only; this synthetic record is never acceptance evidence."""
    from scripts.validate_phase_4j37 import BASE
    matrix = read('docs/evidence/phase_4j36_c2/outcomes.json')
    matrix.update(schema='island.j37-outcomes.v1', base_commit=BASE, charge_origin='RESP/Psi4', native_am1bcc_success=False,
                  provenance_relocation=True, provenance_fresh_process=True, reviews={'linear': 'a' * 64, 'network': 'b' * 64})
    for kind, row in matrix['cases'].items():
        row.update(worker_exit_code=0, periodic_reference_passed=True, independent_review_identity=matrix['reviews'][kind])
    matrix['identity'] = digest({k: v for k, v in matrix.items() if k != 'identity'})
    receipt = {'schema': 'island.j37-receipt.v1', 'base_commit': BASE, 'acceptance_met': True, 'charge_origin': 'RESP/Psi4', 'native_am1bcc_success': False,
               'production_validated': False, 'simulation_readiness': 'not_established', 'outcomes_identity': matrix['identity'],
               'implementation_commit': 'c' * 40, 'tests': {name: {} for name in ('focused', 'j22-j36', 'full_supported')}}
    receipt['identity'] = digest(receipt)
    return receipt, matrix


@pytest.mark.parametrize('mutation', ['all_blocked', 'missing_case', 'origin', 'readiness', 'source_review', 'reference_failure', 'missing_child', 'fixture_key', 'checkpoint_claim'])
def test_resigned_receipt_contract_mutations(synthetic_contract, mutation):
    from scripts.validate_phase_4j37_receipt import validate_structure
    receipt, matrix = synthetic_contract
    if mutation == 'all_blocked':
        for row in matrix['cases'].values():row.update(status='blocked', acceptance_met=False, npt_executed=False)
    elif mutation == 'missing_case':del matrix['cases']['network']
    elif mutation == 'origin':receipt['charge_origin'] = 'external_oplsaa'
    elif mutation == 'readiness':receipt['production_validated'] = True
    elif mutation == 'source_review':matrix['cases']['linear']['independent_review_identity'] = '0' * 64
    elif mutation == 'reference_failure':matrix['cases']['network']['periodic_reference_passed'] = False
    elif mutation == 'missing_child':matrix['provenance_fresh_process'] = False
    elif mutation == 'fixture_key':matrix['cases']['linear']['case'] = 'gaff2/network'
    else:matrix['cases']['linear']['comparison_passed'] = False
    matrix['identity'] = digest({k: v for k, v in matrix.items() if k != 'identity'})
    receipt['outcomes_identity'] = matrix['identity']
    receipt['identity'] = digest({k: v for k, v in receipt.items() if k != 'identity'})
    with pytest.raises((ValueError, AnalysisError)):
        validate_structure(receipt, matrix)


def test_synthetic_structure_is_not_full_receipt_acceptance(synthetic_contract, tmp_path):
    from scripts.validate_phase_4j37_receipt import validate, validate_structure
    receipt, matrix = synthetic_contract
    validate_structure(receipt, matrix)
    # A readable outcome contract without raw/workflow/artifact evidence cannot
    # pass the real receipt entry point.
    with pytest.raises((ValueError, KeyError)):
        validate(receipt, tmp_path)
