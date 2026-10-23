"""Native evidence audits and software negative controls; no synthetic MD success."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
import pytest

from island.forcefields.ambertools.provided import ProvidedChargeError
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_precision import (
    PrecisionRESPCharges,
    audit_in_memory,
    audit_native,
    exact_decimal,
    load_audit,
    load_in_memory_audit,
    signed,
)
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input
from scripts.validate_phase_4j36_c1 import PARENT, WORKER
from scripts.validate_phase_4j36_c1_receipt import validate_outcomes


@pytest.fixture(scope='module')
def original(tmp_path_factory):
    root = tmp_path_factory.mktemp('original-native-outputs')
    with tarfile.open(PARENT / 'artifacts/native-quantum.tar.gz') as archive:
        archive.extractall(root, filter='data')
    return root


@pytest.fixture
def case(original, tmp_path):
    root = tmp_path / 'linear'
    shutil.copytree(original / 'quantum-linear', root)
    parent = json.loads((PARENT / 'receipt.json').read_text())
    runtime = json.loads((PARENT / 'environment.json').read_text())
    kwargs = {'expected_identity': parent['cases']['linear']['execution_identity'], 'expected_input': build_input('linear'), 'runtime_identity': runtime['identity'], 'worker_sha256': file_hash(WORKER)}
    graph = build_fixture('linear', 'GAFF2').fixture.packed.graph
    return root, graph, runtime, kwargs


def test_authoritative_residual_not_decimal_rounding(case):
    root, graph, _, kwargs = case
    normal = audit_native(root, graph, **kwargs)
    with localcontext() as context:
        context.prec = 3
        low = audit_native(root, graph, **kwargs)
    assert normal == low
    assert normal['charge_conservation_passed'] is False
    assert normal['classification'] == 'authoritative_native_precision_insufficient'
    for c in normal['components']:
        assert Fraction(*c['native_total_fraction']) != c['formal_total_e']
        assert c['json_roundtrip_bit_exact'] and c['input_coordinate_mapping_exact']
    for r in normal['rows']:
        assert Decimal(r['decimal_e']) == Decimal.from_float(r['charge_e'])
        assert float.fromhex(r['ieee754_hex']) == r['charge_e']


@pytest.mark.parametrize('value', [0.0, 0.1, -0.3, 5e-324, 1.7976931348623157e308])
def test_exact_binary_expansion_independent_of_decimal_context(value):
    with localcontext() as ctx:
        ctx.prec = 2
        fraction = Fraction.from_float(value)
        assert Decimal(exact_decimal(fraction)) == Decimal.from_float(value)


def test_strict_loader_blocks_before_review_or_publication(case, monkeypatch):
    root, graph, runtime, kwargs = case
    record = audit_native(root, graph, **kwargs)
    (root / 'precision-audit.json').write_text(json.dumps(record))
    monkeypatch.setattr('island.forcefields.ambertools.resp_precision.validate_runtime', lambda *a, **kw: runtime)
    def prohibited(*args, **kw):
        raise AssertionError('Failed native charges must not reach reviewed publication')
    monkeypatch.setattr('island.forcefields.ambertools.resp_precision.RESPCharges.load', prohibited)
    with pytest.raises(ProvidedChargeError) as error:
        PrecisionRESPCharges.load(root, graph, runtime=runtime, expected_runtime_identity=runtime['identity'],
                                  expected_audit_identity=record['identity'], review_directory=root / 'never-created',
                                  expected_review_identity='0' * 64, **kwargs)
    assert error.value.classification == 'authoritative_native_precision_insufficient'


def test_relocation_and_fresh_process(case, tmp_path):
    root, graph, _, kwargs = case
    record = audit_native(root, graph, **kwargs)
    (root / 'precision-audit.json').write_text(json.dumps(record))
    target = tmp_path / 'moved'
    shutil.copytree(root, target)
    assert load_audit(target, graph, expected_audit_identity=record['identity'], **kwargs) == record
    script = '''import sys,json
from island.forcefields.ambertools.resp_precision import load_audit
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input
from scripts.validate_phase_4j36_c1 import PARENT,WORKER
from island.forcefields.ambertools.resp import file_hash
p=json.loads((PARENT/'receipt.json').read_text());r=load_audit(sys.argv[1],build_fixture('linear','GAFF2').fixture.packed.graph,expected_audit_identity=sys.argv[2],expected_identity=p['cases']['linear']['execution_identity'],expected_input=build_input('linear'),runtime_identity=p['runtime_identity'],worker_sha256=file_hash(WORKER));print(r['identity'])'''
    child = subprocess.run([sys.executable, '-c', script, str(target), record['identity']], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == 0, child.stderr
    assert child.stdout.strip() == record['identity']


@pytest.mark.parametrize('change', ['ids', 'elements', 'coordinates', 'charge', 'source', 'origin', 'formal', 'method'])
def test_stale_or_mutated_original_rejected(case, change):
    root, graph, _, kwargs = case
    if change in ('ids', 'elements', 'coordinates', 'formal', 'method'):
        data = deepcopy(kwargs['expected_input'])
        if change == 'ids':
            data['components'][0]['atoms'].reverse()
        elif change == 'elements':
            data['components'][0]['atoms'][0]['element'] = 'N'
        elif change == 'coordinates':
            data['components'][0]['atoms'][0]['position_angstrom'][0] += 1
        elif change == 'formal':
            data['components'][0]['formal_charge'] = 1
        else:
            data['basis'] = 'STO-3G'
        kwargs['expected_input'] = data
    elif change == 'source':
        (root / 'worker-source.py').write_text('changed')
    elif change == 'origin':
        p = json.loads((root / 'result.json').read_text());p['charge_origin'] = 'external_oplsaa'
        (root / 'result.json').write_text(json.dumps(p))
    else:
        p = root / 'component-0/native-resp-charges.npy'
        values = np.load(p, allow_pickle=False);values[0] += .01
        np.save(p, values, allow_pickle=False)
    with pytest.raises(ProvidedChargeError):
        audit_native(root, graph, **kwargs)


@pytest.fixture
def memory_case(original, tmp_path):
    root = tmp_path / 'memory'
    with tarfile.open('docs/evidence/phase_4j36_c1/artifacts/native-memory.tar.gz') as archive:
        archive.extractall(tmp_path, filter='data')
    shutil.copytree(tmp_path / 'native-esp-linear', root)
    parent = json.loads((PARENT / 'receipt.json').read_text())
    graph = build_fixture('linear', 'GAFF2').fixture.packed.graph
    kwargs = {'native_directory': original / 'quantum-linear', 'in_memory_worker_sha256': file_hash('scripts/phase_4j36_c1_native_esp.py'),
              'expected_identity': parent['cases']['linear']['execution_identity'], 'expected_input': build_input('linear'),
              'runtime_identity': parent['runtime_identity'], 'worker_sha256': file_hash(WORKER)}
    return root, graph, kwargs


def test_native_memory_lossless_but_total_still_fails(memory_case):
    root, graph, kwargs = memory_case
    record = audit_in_memory(root, graph, **kwargs)
    assert record['charge_conservation_passed'] is False
    assert record['serialization_loss_detected'] is False
    assert record['source_significand_bits'] == 53
    (root / 'precision-audit.json').write_text(json.dumps(record))
    assert load_in_memory_audit(root, graph, expected_audit_identity=record['identity'], **kwargs) == record


@pytest.mark.parametrize('change', ['round', 'reorder', 'nonfinite', 'boolean', 'geometry', 'esp', 'basis', 'origin', 'worker'])
def test_native_memory_mutations_rejected(memory_case, change):
    root, graph, kwargs = memory_case
    raw = json.loads((root / 'result.json').read_text())
    if change in ('basis', 'origin'):
        raw['basis' if change == 'basis' else 'charge_origin'] = 'changed'
    elif change == 'worker':
        (root / 'worker-source.py').write_text('changed')
    else:
        path = root / ('component-0/' + {'round': 'charges-native.json', 'reorder': 'charges-native.json', 'boolean': 'charges-native.json', 'nonfinite': 'charges-native.npy', 'geometry': 'geometry-bohr.npy', 'esp': 'esp-native.npy'}[change])
        if path.suffix == '.json':
            rows = json.loads(path.read_text())
            if change == 'round': rows[0]['charge_e'] = round(rows[0]['charge_e'], 3)
            elif change == 'reorder': rows.reverse()
            else: rows[0]['charge_e'] = True
            path.write_text(json.dumps(rows))
        else:
            array = np.load(path, allow_pickle=False)
            array.flat[0] = float('nan') if change == 'nonfinite' else array.flat[0] + 1
            np.save(path, array, allow_pickle=False)
        raw['artifacts'][str(path.relative_to(root))] = file_hash(path)
    (root / 'result.json').write_text(json.dumps(raw))
    with pytest.raises(ProvidedChargeError):
        audit_in_memory(root, graph, **kwargs)


def test_cli_nonzero_preserves_diagnostic(tmp_path):
    target = tmp_path / 'run'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j36_c1', '--offline', '--output', str(target)], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == 1, child.stderr
    result = json.loads((target / 'result.json').read_text())
    audits = {k: json.loads((target / k / 'precision-audit.json').read_text()) for k in ('linear', 'network')}
    validate_outcomes(result, audits, result['runtime_identity'], {k: json.loads((target / ('memory-' + k) / 'precision-audit.json').read_text()) for k in ('linear', 'network')})
    assert all(c['classification'] == 'authoritative_native_precision_insufficient' for c in result['cases'].values())
    result['cases']['linear']['dynamics_executed'] = True
    with pytest.raises(ProvidedChargeError):
        validate_outcomes(signed({k: v for k, v in result.items() if k != 'identity'}), audits, result['runtime_identity'], {k: json.loads((target / ('memory-' + k) / 'precision-audit.json').read_text()) for k in ('linear', 'network')})
