"""Exact-graph source controls. Mutations are software negative controls, not MD."""

import json
import os
import shutil
import subprocess
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.forcefields.ambertools.provided import reference_digest
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from scripts.phase_4j38_c3_input import build_input, validate_input
from scripts.phase_4j38_c3_reference import audit
from scripts.phase_4j38_fixture import build
from scripts.validate_phase_4j38_c3_receipt import validate_raw

ROOT = Path(__file__).resolve().parents[1]
PYTHON = ROOT.parent / 'island-validation/phase4g1-env/bin/python'
ENV = {**os.environ, 'PYTHONPATH': '.:src', 'OPENMM_CPU_THREADS': '1', 'OMP_NUM_THREADS': '1'}


def read(path):
    return json.loads(Path(path).read_text())


def resign(value):
    value['identity'] = reference_digest({k: v for k, v in value.items() if k != 'identity'})
    return value


def cli(module, *args):
    return subprocess.run([str(PYTHON), '-m', module, *map(str, args)], env=ENV, cwd=ROOT, capture_output=True, text=True, check=False)


@pytest.fixture(scope='module')
def raw(tmp_path_factory):
    root = tmp_path_factory.mktemp('c3-native-source')
    with tarfile.open(ROOT / 'docs/evidence/phase_4j38_c3/fresh-source.tar.gz') as archive:
        archive.extractall(root, filter='data')
    return root / 'raw'


@pytest.mark.parametrize('kind,count', [('pe', 52), ('aromatic', 44)])
def test_fresh_exact_graph_charge_conservation(raw, kind, count):
    data = read(raw / (kind + '-charges.json'))
    graph = build(kind)[0].packed.graph
    assert data['final_graph_identity'] == graph.identity
    vector = FiniteRESPCharges(json.dumps(data)).validate(graph, expected_reference_identity=data['identity'])
    assert len(vector) == count
    assert data['charge_origin'] == 'resp_psi4' and not data['native_am1bcc_success']
    assert data['policy']['accepted'] and all(c['accepted'] for c in data['policy']['components'])
    review = audit(raw, kind, data, read(raw / 'runtime.json'), read(raw / 'anchors.json'))
    assert review == read(raw / (kind + '-review.json'))
    assert review['import_bit_exact'] and review['raw_provenance_verified']


@pytest.mark.parametrize('kind', ['pe', 'aromatic'])
@pytest.mark.parametrize('mutation', ['graph', 'ids', 'coordinates', 'formal_charge', 'component', 'basis', 'method'])
def test_wrong_quantum_input_rejected(kind, mutation):
    data = build_input(kind)
    if mutation == 'graph':
        data['final_graph_identity'] = build('peo')[0].packed.graph.identity
    elif mutation == 'ids':
        data['components'][0]['atoms'].reverse()
    elif mutation == 'coordinates':
        data['components'][0]['atoms'][0]['position_angstrom'][0] += .01
    elif mutation == 'formal_charge':
        data['components'][0]['formal_charge'] = 1
    elif mutation == 'component':
        data['components'].reverse()
    else:
        data[mutation] = 'unsupported'
    with pytest.raises(ValueError):
        validate_input(data)


@pytest.mark.parametrize('mutation', ['vector', 'order', 'component', 'source', 'executable', 'origin', 'tolerance', 'boolean', 'rounding', 'normalize'])
def test_resigned_import_mutations_rejected(raw, mutation):
    data = deepcopy(read(raw / 'pe-charges.json'))
    if mutation == 'vector':
        data['rows'][0]['charge_e'] += .001
    elif mutation == 'order':
        data['rows'].reverse()
    elif mutation == 'component':
        data['source_audit']['components'][0]['formal_total_e'] = 1
        resign(data['source_audit'])
    elif mutation == 'source':
        data['source_audit']['worker_sha256'] = '0' * 64
        resign(data['source_audit'])
    elif mutation == 'executable':
        data['runtime']['files']['bin/python'] = '0' * 64
        resign(data['runtime'])
    elif mutation == 'origin':
        data['charge_origin'] = 'external_oplsaa'
    elif mutation == 'tolerance':
        data['policy']['system']['allowed_residual_e'] = 1.0
        resign(data['policy'])
    elif mutation == 'boolean':
        data['rows'][0]['charge_e'] = True
    elif mutation == 'rounding':
        for r in data['rows']:
            r['charge_e'] = round(r['charge_e'], 3)
    elif mutation == 'normalize':
        data['normalization'] = True
    resign(data)
    with pytest.raises((ValueError, AttributeError)):
        FiniteRESPCharges(json.dumps(data)).validate(build('pe')[0].packed.graph, expected_reference_identity=data['identity'])


@pytest.mark.parametrize('name', [
    'native-esp-pe/component-0/charges-native.npy',
    'native-esp-aromatic/component-1/geometry-bohr.npy',
    'quantum-pe/component-0/geometry-input.dat',
    'native-esp-pe/component-0/psi4.out',
    'native-esp-pe/result.json',
])
def test_changed_raw_files_rejected(raw, tmp_path, name):
    root = tmp_path / 'raw'; shutil.copytree(raw, root)
    p = root / name; p.write_bytes(p.read_bytes() + b'NEGATIVE CONTROL')
    with pytest.raises(ValueError):
        validate_raw(root)


def test_resigned_raw_anchor_mutation_rejected(raw, tmp_path):
    root = tmp_path / 'raw'; shutil.copytree(raw, root)
    data = read(root / 'anchors.json')
    data['cases']['pe']['charges_sha256'] = '0' * 64
    (root / 'anchors.json').write_text(json.dumps(data))
    with pytest.raises(ValueError, match='anchors'):
        validate_raw(root)


def test_relocation_independent_fresh_process(raw, tmp_path):
    root = tmp_path / 'moved'; shutil.copytree(raw, root)
    for kind in ('pe', 'aromatic'):
        output = tmp_path / (kind + '.json')
        result = cli('scripts.phase_4j38_c3_reference', '--root', root, '--fixture', kind,
                     '--charges', root / (kind + '-charges.json'), '--runtime', root / 'runtime.json',
                     '--anchors', root / 'anchors.json', '--output', output)
        assert result.returncode == 0, result.stdout + result.stderr
        assert read(output) == read(raw / (kind + '-review.json'))


def test_worker_rejects_borrowed_charge_graph_and_preserves_diagnostic(raw, tmp_path):
    root = tmp_path / 'raw'; shutil.copytree(raw, root)
    shutil.copyfile(root / 'aromatic-charges.json', root / 'pe-charges.json')
    output = tmp_path / 'case'
    result = cli('scripts.validate_phase_4j38_c3', '--worker', 'pe', '--raw', root, '--output', output)
    assert result.returncode == 1
    row = read(output / 'result.json')
    assert not row['execution_met'] and not row['software_acceptance_met']
    assert row['diagnostic']['stage'] == 'charge_source_audit'
    assert (output / 'failure.log').is_file() and not (output / 'uninterrupted').exists()


def test_missing_runtime_keeps_both_cells_blocked(tmp_path, monkeypatch):
    from island.forcefields.ambertools.psi4_runtime import Psi4RuntimeError
    from scripts import validate_phase_4j38_c3 as driver
    def missing(*args, **kwargs):
        raise Psi4RuntimeError('Missing declared Psi4 Python: negative-control path', 'asset_unavailable')
    monkeypatch.setattr(driver, 'validate_runtime', missing)
    monkeypatch.setattr(driver, 'recheck', lambda output: shutil.copytree(ROOT / 'docs/evidence/phase_4j38_c2/run', output))
    output = tmp_path / 'output'
    value = driver.run(output)
    assert value['new_md_executions'] == 0 and value['acceptance_met'] is False
    assert all(value['cases']['gaff2/' + k]['charge_source_status'] == 'asset_unavailable' for k in ('pe', 'aromatic'))
    assert (output / 'failure.log').is_file()


def test_native_aromatic_source_remains_blocked():
    from scripts.phase_4j38_c2 import SourceIncompatibleError, require_native_source
    decision = read(ROOT / 'docs/evidence/phase_4j38_c2/run/native-source-decision.json')
    assert [c['residual_e'] for c in decision['components']] == ['-0.230'] * 2
    with pytest.raises(SourceIncompatibleError):
        require_native_source(decision)


@pytest.fixture(scope='module')
def runs(tmp_path_factory, raw):
    configured = os.environ.get('ISLAND_J38_C3_RUN')
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp('c3-workflows')
    archive_path = ROOT / 'docs/evidence/phase_4j38_c3/artifacts.tar.gz'
    if archive_path.exists():
        with tarfile.open(archive_path) as archive:
            archive.extractall(root, filter='data')
        return root / 'run'
    result = cli('scripts.validate_phase_4j38_c3', '--output', root / 'run', '--from-reviewed', raw)
    assert result.returncode == 1, result.stdout + result.stderr
    return root / 'run'


@pytest.mark.parametrize('kind', ['pe', 'aromatic'])
def test_actual_gaff2_workflow_checkpoint_and_rng_equality(runs, kind):
    from island.analysis.npt_continuation import compare_continuations
    case = runs / ('gaff2-' + kind)
    row = read(case / 'result.json')
    assert row['software_acceptance_met'] and row['frame_count'] == 321 and row['density_blocks'] == 4
    comparison = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
    assert comparison['passed']
    assert comparison['checkpoint_hashes'] == comparison['resumed_checkpoint_hashes']
    assert comparison['rng_state_hashes'] == comparison['resumed_rng_state_hashes']
    assert row['box_volume_changed'] and not row['acceptance_met']
    assert row['gates']['thermodynamic_equilibration']['status'] == 'insufficient_sampling'


@pytest.mark.parametrize('field', ['positions_nm', 'box_vectors_nm', 'velocities_nm_ps', 'potential_energy_kj_mol', 'measured_temperature_kelvin', 'stage', 'checkpoint_base64'])
def test_gaff2_mutated_continuation_rejected(runs, tmp_path, field):
    from island.analysis.npt import _body, _signed
    from island.analysis.npt_continuation import load_continuation
    from island.exceptions import AnalysisError
    destination = tmp_path / 'resumed'
    shutil.copytree(runs / 'gaff2-aromatic/resumed', destination)
    record = read(destination / 'frames.json')
    frame = record['frames'][-1]
    if field in ('positions_nm', 'box_vectors_nm', 'velocities_nm_ps'):
        frame[field][0][0] += .001
    elif field in ('stage', 'checkpoint_base64'):
        frame[field] = 'software-negative-control'
    else:
        frame[field] += 1
    (destination / 'frames.json').write_text(json.dumps(_signed(_body(record))))
    with pytest.raises((ValueError, AnalysisError)):
        load_continuation(destination)


def test_failed_fresh_process_continuation_nonzero_preserves_checkpoint(runs, tmp_path):
    source = runs / 'gaff2-pe/relocated'
    result = cli('scripts.validate_phase_4j33_c1', '--resume', source, '--destination', tmp_path / 'bad',
                 '--restart', '0' * 64, '--parent', '0' * 64, '--config', '0' * 64)
    assert result.returncode != 0 and 'Stale checkpoint' in result.stderr
    assert (source / 'frames.json').is_file() and not (tmp_path / 'bad').exists()


def test_matrix_preserves_native_failure_and_separates_scientific_gates(runs):
    from scripts.validate_phase_4j38_c3 import matrix
    previous = read(runs / 'c2-recheck/matrix.json')
    outcomes = read(runs / 'outcomes.json')['cases']
    value = matrix(previous, outcomes)
    assert value == read(runs / 'matrix.json')
    assert value['new_md_executions'] == 2 and value['unmet_charge_cells'] == ['oplsaa/aromatic']
    assert value['cases']['oplsaa/aromatic'] == previous['cases']['oplsaa/aromatic']
    for key in ('gaff2/pe', 'gaff2/aromatic'):
        row = value['cases'][key]
        assert row['stages']['charge_completeness']['status'] == 'passed'
        assert row['stages']['numerical_reference']['status'] == 'passed'
        assert row['stages']['equilibration']['status'] == 'insufficient_sampling'
        assert row['stages']['scientific_validation']['status'] == 'blocked'



def test_missing_reviewed_source_cli_is_asset_unavailable_before_md(tmp_path):
    output = tmp_path / 'missing-reviewed-output'
    result = cli('scripts.validate_phase_4j38_c3', '--output', output, '--from-reviewed', tmp_path / 'no-such-reviewed-bundle')
    assert result.returncode == 1, result.stdout + result.stderr
    value = read(output / 'matrix.json')
    for kind in ('pe', 'aromatic'):
        row = value['cases']['gaff2/' + kind]
        assert row['charge_source_status'] == 'asset_unavailable'
        assert row['stages']['charge_completeness']['status'] == 'asset_unavailable'
        assert row['execution_origin'] == 'not_executed' and not row['new_dynamics_executed']
        assert not (output / ('gaff2-' + kind)).exists()
    assert (output / 'failure.log').is_file()


def test_missing_charge_file_worker_preserves_asset_diagnostic(raw, tmp_path):
    root = tmp_path / 'raw'; shutil.copytree(raw, root)
    (root / 'pe-charges.json').unlink()
    output = tmp_path / 'case'
    result = cli('scripts.validate_phase_4j38_c3', '--output', output, '--worker', 'pe', '--raw', root)
    assert result.returncode == 1
    row = read(output / 'result.json')
    assert row['classification'] == 'asset_unavailable'
    assert row['diagnostic']['stage'] == 'charge_source_audit'
    assert row['gates']['source_completeness']['status'] == 'asset_unavailable'
    assert not row['execution_met'] and not (output / 'uninterrupted').exists()
