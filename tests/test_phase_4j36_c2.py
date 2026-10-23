"""Finite-precision RESP policy controls; fixtures use genuine native outputs."""

import json
import math
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.forcefields.ambertools.provided import (
    GAFF2ProvidedCharges,
    ProvidedChargeError,
)
from island.forcefields.ambertools.resp_policy import (
    FiniteRESPCharges,
    create,
    make_policy,
    measure,
)
from island.forcefields.ambertools.resp_precision import signed
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_c2_reference import check, parse_npy


@pytest.fixture(scope='module')
def sources(tmp_path_factory):
    root = tmp_path_factory.mktemp('real-resp-native')
    for name in ('docs/evidence/phase_4j36/artifacts/native-quantum.tar.gz', 'docs/evidence/phase_4j36_c1/artifacts/native-memory.tar.gz'):
        with tarfile.open(name) as archive:
            archive.extractall(root, filter='data')
    return root


def case(sources, label):
    kind = label.removeprefix('memory-')
    graph = build_fixture(kind, 'GAFF2').fixture.packed.graph
    audit = json.loads(Path('docs/evidence/phase_4j36_c1', label + '.json').read_text())
    runtime = json.loads(Path('docs/evidence/phase_4j36/environment.json').read_text())
    native = sources / (('native-esp-' if label.startswith('memory-') else 'quantum-') + kind)
    independent = check(native, audit)
    return graph, audit, runtime, independent


@pytest.mark.parametrize('label', ['linear', 'network', 'memory-linear', 'memory-network'])
def test_real_source_vectors_conserve_under_declared_policy(sources, label):
    graph, audit, runtime, independent = case(sources, label)
    charges = create(audit, graph, runtime, independent)
    expected = {r['atom_id']: r['charge_e'] for r in audit['rows']}
    assert charges.validate(graph, expected_reference_identity=charges.identity) == expected
    assert GAFF2ProvidedCharges(charges.raw).validate(graph, expected_reference_identity=charges.identity) == expected
    policy = charges.payload['policy']
    assert policy['accepted']
    assert all(c['within_source_half_ulp_envelope'] for c in policy['components'])
    assert all(c['relative_residual_formal'] is None for c in policy['components'])
    assert charges.payload['rows'] == audit['rows']  # no normalization, no rewriting
    assert all(c['absolute_residual_e'] != '0' for c in policy['components'])


def test_compensated_and_exact_source_analysis():
    result = measure([1.0, 1e-16, -1.0], 0)
    assert result['naive_sum_e'] == 0
    assert result['compensated_fsum_e'] == 1e-16
    assert float(result['decimal_from_binary_total_e']) == result['compensated_fsum_e']
    assert result['accepted']
    assert not measure([.1, -.1 + .003], 0)['accepted']
    assert not measure([.1, -.1 + 1e-8], 0)['accepted']


@pytest.mark.parametrize('value', [True, float('nan'), float('inf'), -float('inf')])
def test_invalid_values(value):
    with pytest.raises(ProvidedChargeError):
        measure([value], 0)


@pytest.mark.parametrize('mutation', ['vector', 'ids', 'element', 'component', 'formal', 'origin', 'source_hash', 'executable', 'absolute', 'relative', 'sum_method', 'normalization', 'independent', 'precision'])
def test_persisted_mutations_invalidated(sources, tmp_path, mutation):
    graph, audit, runtime, independent = case(sources, 'memory-linear')
    good = create(audit, graph, runtime, independent)
    bad = deepcopy(good.payload)
    if mutation == 'vector':bad['rows'][0]['charge_e'] += .001
    elif mutation == 'ids':bad['rows'].reverse()
    elif mutation == 'element':bad['rows'][0]['element'] = 'N'
    elif mutation == 'component':bad['source_audit']['components'][0]['atom_ids'].reverse()
    elif mutation == 'formal':bad['source_audit']['components'][0]['formal_total_e'] = 1
    elif mutation == 'origin':bad['charge_origin'] = 'external_oplsaa'
    elif mutation == 'source_hash':bad['runtime']['files']['lib/python3.11/site-packages/resp/espfit.py'] = '0' * 64
    elif mutation == 'executable':bad['runtime']['files']['bin/psi4'] = '0' * 64
    elif mutation == 'absolute':bad['policy']['components'][0]['absolute_tolerance_e'] = '.1'
    elif mutation == 'relative':bad['policy']['components'][0]['relative_tolerance'] = .1
    elif mutation == 'sum_method':bad['policy']['summation_method'] = 'naive-only'
    elif mutation == 'normalization':bad['normalization'] = True
    elif mutation == 'independent':bad['independent']['system_residual_fraction'] = [0, 1]
    else:bad['source_audit']['source_significand_bits'] = 24
    bad = signed({k: v for k, v in bad.items() if k != 'identity'})
    path = tmp_path / 'bad.json';path.write_text(json.dumps(bad))
    with pytest.raises(ProvidedChargeError):
        FiniteRESPCharges.load(path, graph, expected_identity=good.identity)
    # Re-signing the outer record does not legitimize an altered nested policy.
    with pytest.raises(ProvidedChargeError):
        FiniteRESPCharges.load(path, graph, expected_identity=bad['identity'])


def test_independent_parser_truncation_rejected(sources, tmp_path):
    original = sources / 'native-esp-linear/component-0/charges-native.npy'
    path = tmp_path / 'truncated.npy';path.write_bytes(original.read_bytes()[:-1])
    with pytest.raises(ValueError, match='Truncated'):
        parse_npy(path)


def test_wrong_graph_and_source_rejected(sources):
    graph, audit, runtime, independent = case(sources, 'linear')
    audit['final_graph_identity'] = '0' * 64
    with pytest.raises(ProvidedChargeError):
        create(audit, graph, runtime, independent)
    audit['charge_origin'] = 'external_oplsaa'
    with pytest.raises(ProvidedChargeError):
        make_policy(audit, runtime)


def test_relocation_and_fresh_process(sources, tmp_path):
    graph, audit, runtime, independent = case(sources, 'memory-network')
    good = create(audit, graph, runtime, independent)
    path = tmp_path / 'original.json';good.save(path)
    moved = tmp_path / 'new/path/resp.json';moved.parent.mkdir(parents=True);shutil.copyfile(path, moved)
    code = '''import sys
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from scripts.phase_4j34_real_fixture import build_fixture
r=FiniteRESPCharges.load(sys.argv[1],build_fixture('network','GAFF2').fixture.packed.graph,expected_identity=sys.argv[2]);print(r.identity)'''
    child = subprocess.run([sys.executable, '-c', code, str(moved), good.identity], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == 0, child.stderr
    assert child.stdout.strip() == good.identity
    assert moved.read_bytes() == path.read_bytes()
    assert math.isfinite(float(good.payload['policy']['system']['absolute_residual_e']))


@pytest.mark.parametrize('kind', ['linear', 'network'])
def test_real_gaff2_preparation_and_continuation(sources, tmp_path, kind):
    """Regression: linear AmberParm addition previously lost native charge bits."""
    from island.forcefields.ambertools.runtime import declared_environment
    graph, audit, quantum, independent = case(sources, 'memory-' + kind)
    charges = create(audit, graph, quantum, independent)
    path = tmp_path / 'charges.json';charges.save(path)
    runtime = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
    output = tmp_path / 'execution'
    command = ['../island-validation/phase4g1-env/bin/python', '-m', 'scripts.validate_phase_4j36_c2', '--output', str(output), '--worker', kind, '--charges', str(path), '--expected', charges.identity]
    child = subprocess.run(command, capture_output=True, text=True, env=declared_environment(runtime), check=False, timeout=240)
    assert child.returncode == 0, child.stdout + child.stderr + (output / 'result.json').read_text()
    row = json.loads((output / 'result.json').read_text())
    comparison = json.loads((output / 'comparison.json').read_text())
    assert row['acceptance_met'] and row['charge_origin'] == 'resp_psi4'
    assert comparison['passed'] and comparison['checkpoint_hashes'] == comparison['resumed_checkpoint_hashes']
    assert comparison['rng_state_hashes'] == comparison['resumed_rng_state_hashes']
    assert json.loads((output / 'resp-charges.json').read_text())['rows'] == audit['rows']


def test_workflow_outcome_hash_domain():
    from island.analysis.npt import _signed
    from scripts.validate_phase_4j36_c2_receipt import verify_outcome_identity
    body = {'scope': 'software hash-domain control; not scientific acceptance'}
    actual = _signed(body)
    verify_outcome_identity(actual, actual['identity'])
    # Charge-record canonical hashes intentionally differ (no terminal newline).
    wrong = signed(body)
    with pytest.raises(ProvidedChargeError, match='outcome binding'):
        verify_outcome_identity(wrong, wrong['identity'])
    actual['scope'] = 'mutated'
    with pytest.raises(ProvidedChargeError, match='outcome binding'):
        verify_outcome_identity(actual, actual['identity'])
