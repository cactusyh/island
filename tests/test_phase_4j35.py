"""Provided-charge software negative controls and real, external-source MD."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.exceptions import IslandError
from island.forcefields.ambertools.models import digest
from island.forcefields.ambertools.provided import (
    GAFF2ProvidedCharges,
    ProvidedChargeError,
)
from island.forcefields.ambertools.provided_preparation import validate_backend_record
from scripts.phase_4j32_real_fixture import FOYER
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j35 import CASES, validate_case


def child(*args):
    return subprocess.run([str(FOYER), *map(str, args)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src', 'OPENMM_CPU_THREADS': '1', 'OMP_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1'})


@pytest.fixture(scope='module', params=['linear', 'network'])
def supplied(request):
    kind = request.param
    reference = json.loads(Path(f'docs/evidence/phase_4j35/references/{kind}.json').read_text())
    graph = build_fixture(kind, 'GAFF2').fixture.packed.graph
    record = GAFF2ProvidedCharges.from_reference(graph, reference, expected_reference_identity=reference['identity'])
    return graph, reference, record


def test_valid_source_vector_unchanged(supplied):
    graph, reference, record = supplied
    charges = record.validate(graph, expected_reference_identity=reference['identity'])
    assert list(charges) == [r['atom_id'] for r in reference['rows']]
    assert list(charges.values()) == [r['charge_e'] for r in reference['rows']]
    assert record.payload['native_am1bcc_success'] is False
    assert record.payload['diagnostic'] == 'provided_charge_validated'
    copy = record.payload
    copy['rows'].clear()
    assert record.payload['rows']


@pytest.mark.parametrize('mutation', ['missing', 'extra', 'duplicate', 'reordered', 'bool_id', 'bool', 'nan', 'infinity', 'units', 'component', 'nonneutral', 'graph', 'source', 'method', 'rounded', 'normalized', 'evidence', 'native_success'])
def test_resigned_vector_mutations_rejected(supplied, mutation):
    graph, reference, record = supplied
    p = record.payload
    if mutation == 'missing':
        p['rows'].pop()
    elif mutation == 'extra':
        p['rows'].append({**p['rows'][0], 'atom_id': 999999})
    elif mutation == 'duplicate':
        p['rows'][1]['atom_id'] = p['rows'][0]['atom_id']
    elif mutation == 'reordered':
        p['rows'].reverse()
    elif mutation == 'bool_id':
        p['rows'][0]['atom_id'] = True
    elif mutation in ('bool', 'nan', 'infinity'):
        p['rows'][0]['charge_e'] = {'bool': True, 'nan': float('nan'), 'infinity': float('inf')}[mutation]
    elif mutation == 'units':
        p['units'] = 'coulomb'
    elif mutation == 'component':
        p['reference']['components'][0]['formal_total_e'] = 1
    elif mutation == 'nonneutral':
        p['rows'][0]['charge_e'] += .001
    elif mutation == 'graph':
        p['final_graph_identity'] = '0' * 64
    elif mutation == 'source':
        p['source_identity'] = '0' * 64
    elif mutation == 'method':
        p['method_identity'] = digest('am1bcc')
    elif mutation == 'rounded':
        for row in p['rows']:
            row['charge_e'] = round(row['charge_e'], 1)
    elif mutation == 'normalized':
        p['rows'][0]['charge_e'] += .000001
        p['rows'][1]['charge_e'] -= .000001
    elif mutation == 'evidence':
        p['evidence_references'] = []
    else:
        p['native_am1bcc_success'] = True
    p['identity'] = digest({k: v for k, v in p.items() if k != 'identity'})
    with pytest.raises(ProvidedChargeError):
        GAFF2ProvidedCharges(json.dumps(p)).validate(graph, expected_reference_identity=reference['identity'])


@pytest.mark.parametrize('mutation', ['source', 'method', 'precision', 'component', 'reorder', 'wrong_anchor'])
def test_self_signed_reference_is_not_source_verification(supplied, mutation):
    from island.forcefields.ambertools.provided import reference_digest

    graph, reference, _ = supplied
    changed = deepcopy(reference)
    if mutation == 'precision':
        changed['rows'][0]['decimal_e'] = '0.1234567890123456'
    elif mutation == 'component':
        changed['components'][0]['provided_total_e'] = '1'
    elif mutation == 'reorder':
        changed['rows'].reverse()
    else:
        changed[mutation] = 'mutated software negative control'
    changed['identity'] = reference_digest({k: v for k, v in changed.items() if k != 'identity'})
    with pytest.raises(ProvidedChargeError):
        GAFF2ProvidedCharges.from_reference(graph, changed, expected_reference_identity=reference['identity'])


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    configured = os.environ.get('ISLAND_J35_RUN')
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp('j35') / 'runs'
    runtime = Path('docs/evidence/phase_4j34/runtime.json')
    p = child('-m', 'scripts.validate_phase_4j35', '--runtime', runtime, '--runtime-identity', json.loads(runtime.read_text())['identity'], '--output', root)
    assert p.returncode == 0, p.stdout + p.stderr
    return root


@pytest.mark.parametrize('name', CASES)
def test_real_acceptance_outcomes(runs, name):
    row = json.loads((runs / name.replace('/', '-') / 'result.json').read_text())
    validate_case(name, row)
    if name.startswith('gaff2/'):
        assert row['native_am1bcc_success'] is False
        changed = {**row, 'native_am1bcc_success': True}
        with pytest.raises(Exception, match='Wrong provided charge'):
            validate_case(name, changed)


@pytest.mark.parametrize('kind', ['linear', 'network'])
def test_real_continuation_and_independent_reconstruction(runs, kind, tmp_path):
    root = runs / ('gaff2-' + kind)
    comparison = compare_continuations(root / 'uninterrupted', root / 'resumed', root / 'relocated')
    assert comparison['passed'] is True
    assert comparison['checkpoint_hashes'] == comparison['resumed_checkpoint_hashes']
    assert comparison['rng_state_hashes'] == comparison['resumed_rng_state_hashes']
    assert load_continuation(root / 'split') == load_continuation(root / 'relocated')
    reference = tmp_path / 'reference.json'
    p = child('-m', 'scripts.phase_4j35_charge_reference', '--fixture', kind, '--output', reference)
    assert p.returncode == 0, p.stderr
    assert json.loads(reference.read_text()) == json.loads((root / 'independent-charge-reference.json').read_text())
    p = child('-c', 'import sys;from island.analysis.npt_continuation import load_continuation;load_continuation(sys.argv[1])', root / 'resumed')
    assert p.returncode == 0, p.stderr
    backend = json.loads((root / 'uninterrupted/backend.json').read_text())
    record = validate_backend_record(backend)
    assert record['native_am1bcc_success'] is False
    assert all(c['record']['charge_outcome']['qm_run'] is False for c in record['components'])


@pytest.mark.parametrize('mutation', ['charge', 'mapping', 'graph', 'source', 'method', 'membership', 'coordinates'])
def test_offline_preparation_mutations(runs, mutation):
    payload = json.loads((runs / 'gaff2-linear/uninterrupted/backend.json').read_text())
    record = validate_backend_record(payload)
    if mutation == 'charge':
        payload['numerical_data']['sites'][0]['charge'] += .001
    elif mutation == 'mapping':
        a, b = list(record['mapping'])[:2]
        record['mapping'][a], record['mapping'][b] = record['mapping'][b], record['mapping'][a]
    elif mutation == 'graph':
        record['packed_graph_raw'] = '{}'
    elif mutation == 'source':
        record['force_field_data']['sha256'] = '0' * 64
    elif mutation == 'method':
        record['charge_method'] = 'am1bcc'
    elif mutation == 'membership':
        record['components'][0]['atom_ids'].pop()
    else:
        record['input_coordinate_identity'] = '0' * 64
    payload['parameter_assignment']['native_payload']['record'] = record
    payload['typed_graph']['evidence'][0] = digest(record)
    with pytest.raises((ValueError, KeyError, IslandError)):
        validate_backend_record(payload)


@pytest.mark.parametrize('failure', ['relocation', 'child'])
def test_damaged_relocated_bundle_fails_in_child(runs, tmp_path, failure):
    root = tmp_path / 'relocated'
    shutil.copytree(runs / 'gaff2-linear/relocated', root)
    (root / ('restart.json' if failure == 'relocation' else 'frames.json')).write_text('{}')
    p = child('-c', 'import sys;from island.analysis.npt_continuation import load_continuation;load_continuation(sys.argv[1])', root)
    assert p.returncode != 0 and p.stderr


def test_missing_runtime_diagnostic_cli(tmp_path):
    root = tmp_path / 'missing-runtime'
    p = child('-m', 'scripts.validate_phase_4j35', '--runtime', tmp_path / 'absent', '--runtime-identity', '0' * 64, '--output', root, '--case', 'gaff2/linear')
    assert p.returncode != 0
    row = json.loads((root / 'result.json').read_text())
    assert row['status'] == 'asset_unavailable' and not row['npt_executed']
    assert (root / 'failure.log').is_file()


def test_native_acceptance_stays_nonzero(tmp_path):
    p = child('-m', 'scripts.validate_phase_4j34_c2', '--output', tmp_path / 'native')
    assert p.returncode != 0


@pytest.mark.parametrize('mutation', ['bool', 'nan', 'infinity', 'units', 'total', 'component', 'precision', 'reorder', 'missing', 'duplicate'])
def test_reference_semantics_even_with_new_anchor(supplied, mutation):
    from island.forcefields.ambertools.provided import reference_digest

    graph, reference, _ = supplied
    r = deepcopy(reference)
    if mutation in ('bool', 'nan', 'infinity'):
        r['rows'][0]['charge_e'] = {'bool': True, 'nan': float('nan'), 'infinity': float('inf')}[mutation]
    elif mutation == 'units':
        r['units'] = 'coulomb'
    elif mutation == 'total':
        r['rows'][0].update(charge_e=1.0, decimal_e='1.0')
    elif mutation == 'component':
        r['components'][0]['provided_total_e'] = '1'
    elif mutation == 'precision':
        r['rows'][0]['decimal_e'] = '0.1234567890123456'
    elif mutation == 'reorder':
        r['rows'].reverse()
    elif mutation == 'missing':
        r['rows'].pop()
    else:
        r['rows'][1]['atom_id'] = r['rows'][0]['atom_id']
    if mutation in ('nan', 'infinity'):
        # Canonical independent evidence cannot even serialize nonfinite values.
        with pytest.raises(ValueError):
            reference_digest(r)
        return
    r['identity'] = reference_digest({k: v for k, v in r.items() if k != 'identity'})
    with pytest.raises(ProvidedChargeError):
        GAFF2ProvidedCharges.from_reference(graph, r, expected_reference_identity=r['identity'])


def test_high_precision_native_serialization_control(tmp_path):
    """Synthetic precision control, never a physical charge-source receipt."""
    code = '''
import json,sys
from decimal import Decimal
from pathlib import Path
from island.forcefields.ambertools.provided import GAFF2ProvidedCharges,reference_digest
from island.forcefields.ambertools.provided_preparation import parameterize_provided
from scripts.phase_4j34_real_fixture import build_fixture
r=json.load(open('docs/evidence/phase_4j35/references/network.json'))
r['provenance']='Synthetic serialization control only; no scientific charge evidence'
r['method']='software_precision_control'
for row,delta in zip(r['rows'][:2], (Decimal('0.000000123456789'),Decimal('-0.000000123456789'))):
 row['decimal_e']=str(Decimal(row['decimal_e'])+delta)
 row['charge_e']=float(row['decimal_e'])
for component in r['components']:
 component['provided_total_e']=str(sum((Decimal(row['decimal_e']) for row in r['rows'] if row['atom_id'] in component['atom_ids']),Decimal(0)))
r['identity']=reference_digest({k:v for k,v in r.items() if k!='identity'})
f=build_fixture('network','GAFF2').fixture
p=GAFF2ProvidedCharges.from_reference(f.packed.graph,r,expected_reference_identity=r['identity'])
runtime=json.load(open('docs/evidence/phase_4j34/runtime.json'))
n=parameterize_provided(f.packed.system,f.packed.graph,p,Path(sys.argv[1]),runtime=runtime,runtime_identity=runtime['identity'],expected_reference_identity=r['identity'])
for row in r['rows']:
 assert abs(n.imported_result.charge_result.assignments[row['atom_id']].charge-row['charge_e'])<=1e-14
'''
    from island.forcefields.ambertools.runtime import declared_environment

    runtime = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
    p = subprocess.run([str(FOYER), '-c', code, str(tmp_path / 'precision')], capture_output=True, text=True, check=False, env=declared_environment(runtime))
    assert p.returncode == 0, p.stdout + p.stderr


@pytest.mark.parametrize('stage', ['provided_charge_source_verification', 'preparation'])
def test_failure_classification_preserves_diagnostics(tmp_path, monkeypatch, stage):
    from types import SimpleNamespace

    from scripts import validate_phase_4j35 as driver

    runtime = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
    monkeypatch.setattr(driver, 'validate_runtime', lambda *a, **k: None)
    monkeypatch.setattr(driver.sys, 'executable', str(FOYER))
    if stage == 'provided_charge_source_verification':
        monkeypatch.setattr(driver.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=1, stdout='', stderr='injected missing source'))
    else:
        def failed(*args, **kwargs):
            raise OSError('Injected parameterization failure; software control')
        monkeypatch.setattr(driver, 'prepare', failed)
    row = driver.worker(tmp_path / 'case', 'gaff2/linear', runtime, runtime['identity'], {**os.environ, 'PYTHONPATH': '.:src'})
    assert row['acceptance_met'] is False
    assert row['npt_executed'] is False
    assert row['charge_diagnostic'] == ('provided_charge_source_unverified' if stage == 'provided_charge_source_verification' else 'parameterization_incomplete')
    assert row['diagnostic']['stage'] == stage
    assert (tmp_path / 'case/failure.log').is_file()
