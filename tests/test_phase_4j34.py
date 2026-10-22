"""J34 runtime negative controls and real native GAFF2 continuation evidence."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.exceptions import AmberToolsUnavailableError, AnalysisError
from island.forcefields.ambertools import runtime as api
from scripts.phase_4j32_real_fixture import FOYER
from scripts.validate_phase_4j34 import CASES, validate_case

RUNTIME = Path('docs/evidence/phase_4j34/runtime.json')
PREPARATION = Path('docs/evidence/phase_4j34/linear-native-preparation.tar.gz')
PREPARATION_SHA = '6f3191258aa760c7adde1b3902809e2b79860b01c79b181d6f23faa97cd401fa'


def child(*args):
    return subprocess.run([str(FOYER), *map(str, args)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src', 'OPENMM_CPU_THREADS': '1', 'OMP_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1'})


@pytest.fixture
def toy_runtime(tmp_path, monkeypatch):
    # Synthetic schema fixture: these /bin/true copies are never scientific or
    # native force-field acceptance evidence. Actual GAFF2 runs are tested below.
    home = tmp_path / 'toy'
    (home / 'bin').mkdir(parents=True)
    (home / 'conda-meta').mkdir()
    for name in (*api.TOOLS, 'python'):
        shutil.copyfile('/bin/true', home / 'bin' / name)
        (home / 'bin' / name).chmod(0o755)
    for name in api.SOURCES:
        p = home / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('gaff2.dat software schema fixture\n')
    files = [str(p.relative_to(home)) for p in home.rglob('*') if p.is_file()]
    (home / 'conda-meta/ambertools-toy.json').write_text(json.dumps({'name': 'ambertools', 'version': 'toy', 'build': 'software-control', 'files': files}))
    original = api._run

    def probe(argv, env, **kwargs):
        if Path(argv[0]) == home / 'bin/python':
            return subprocess.CompletedProcess(argv, 0, '[]' if argv[1] == '-c' else 'No broken requirements found.', '')
        return original(argv, env, **kwargs)

    monkeypatch.setattr(api, '_run', probe)
    record = api.capture_runtime(home)
    return record, api.declared_environment(record)


@pytest.mark.parametrize('name', api.TOOLS)
def test_missing_required_executable(toy_runtime, name):
    record, env = toy_runtime
    Path(record['amberhome'], 'bin', name).unlink()
    with pytest.raises(AmberToolsUnavailableError, match='Missing required executable'):
        api.validate_runtime(record, expected_identity=record['identity'], environment=env)


def test_wrong_path_resolution(toy_runtime, tmp_path):
    record, env = toy_runtime
    shadow = tmp_path / 'shadow'
    shadow.mkdir()
    shutil.copyfile('/bin/true', shadow / 'antechamber')
    (shadow / 'antechamber').chmod(0o755)
    env['PATH'] = str(shadow) + os.pathsep + env['PATH']
    with pytest.raises(AmberToolsUnavailableError, match='Wrong PATH executable resolution'):
        api.validate_runtime(record, expected_identity=record['identity'], environment=env)


@pytest.mark.parametrize('mutation', ['signature', 'resigned_manifest', 'version', 'source', 'package', 'executable', 'library', 'dependencies'])
def test_manifest_environment_and_source_mutations(toy_runtime, monkeypatch, mutation):
    record, env = toy_runtime
    changed = deepcopy(record)
    if mutation in ('signature', 'resigned_manifest', 'version'):
        changed['package']['version'] = 'wrong-version'
        if mutation != 'signature':
            changed = api.signed({k: v for k, v in changed.items() if k != 'identity'})
    elif mutation == 'source':
        Path(record['amberhome'], api.SOURCES[1]).write_text('mutated software negative control')
    elif mutation == 'package':
        Path(record['amberhome'], 'conda-meta/ambertools-toy.json').write_text('{}')
    elif mutation == 'executable':
        Path(record['amberhome'], 'bin/antechamber').write_bytes(b'corrupted executable')
    elif mutation == 'library':
        changed['libraries'][next(iter(changed['libraries']))] = '0' * 64
        changed = api.signed({k: v for k, v in changed.items() if k != 'identity'})
    else:
        original = api._run

        def missing(argv, env, **kwargs):
            if list(argv[1:]) == ['-m', 'pip', 'check']:
                return subprocess.CompletedProcess(argv, 1, 'missing required dependency', '')
            return original(argv, env, **kwargs)

        monkeypatch.setattr(api, '_run', missing)
    pin = changed['identity'] if mutation in ('version', 'library') else record['identity']
    with pytest.raises(AmberToolsUnavailableError):
        api.validate_runtime(changed, expected_identity=pin, environment=env)


def test_missing_runtime_cli_preserves_diagnostic(tmp_path):
    root = tmp_path / 'unavailable'
    p = child('-m', 'scripts.validate_phase_4j34', '--runtime', tmp_path / 'missing.json', '--runtime-identity', '0' * 64, '--output', root, '--case', 'gaff2/network')
    assert p.returncode != 0
    row = json.loads((root / 'result.json').read_text())
    assert row['status'] == 'asset_unavailable' and row['npt_executed'] is False
    assert row['diagnostic']['stage'] == 'runtime_validation'
    assert (root / 'failure.log').is_file()


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    configured = os.environ.get('ISLAND_J34_RUN')
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp('j34') / 'runs'
    runtime = json.loads(RUNTIME.read_text())
    p = child('-m', 'scripts.validate_phase_4j34', '--runtime', RUNTIME, '--runtime-identity', runtime['identity'], '--linear-preparation', PREPARATION, '--linear-preparation-sha256', PREPARATION_SHA, '--output', root)
    assert p.returncode == 0, p.stdout + p.stderr
    return root


@pytest.mark.parametrize('name', CASES)
def test_frozen_matrix_outcomes(runs, name):
    row = json.loads((runs / name.replace('/', '-') / 'result.json').read_text())
    validate_case(name, row)
    if name != 'pcff/network':
        changed = deepcopy(row)
        changed['comparison_passed'] = False
        with pytest.raises(Exception, match='Missing NPT evidence'):
            validate_case(name, changed)


@pytest.mark.parametrize('kind', ['linear', 'network'])
def test_native_gaff2_exact_continuation_and_reconstruction(runs, kind):
    case = runs / ('gaff2-' + kind)
    comparison = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
    assert comparison['passed'] is True
    assert comparison['checkpoint_hashes'] == comparison['resumed_checkpoint_hashes']
    assert comparison['rng_state_hashes'] == comparison['resumed_rng_state_hashes']
    assert load_continuation(case / 'split') == load_continuation(case / 'relocated')
    independent = json.loads((case / 'independent-assignment.json').read_text())
    assert independent['parmchk2_estimation_lines'] == []
    assert independent['max_force_error_kj_mol_angstrom'] <= 1e-5
    p = child('-c', 'import sys;from island.analysis.npt_continuation import load_continuation;load_continuation(sys.argv[1])', case / 'resumed')
    assert p.returncode == 0, p.stderr


def test_actual_runtime_matches_frozen_manifest():
    record = json.loads(RUNTIME.read_text())
    api.validate_runtime(record, expected_identity=record['identity'], environment=api.declared_environment(record))


@pytest.mark.parametrize('mutation', ['charge', 'membership', 'images', 'source', 'coordinates'])
def test_resigned_component_lineage_mutations_rejected(runs, mutation):
    from island.exceptions import InvalidAmberImportResultError
    from island.forcefields.ambertools.models import digest
    from island.forcefields.ambertools.multimolecule import validate_backend_record

    payload = json.loads((runs / 'gaff2-linear/uninterrupted/backend.json').read_text())
    record = validate_backend_record(payload)
    entry = record['components'][0]
    if mutation == 'charge':
        entry['charges'][str(entry['atom_ids'][0])] += 0.01
    elif mutation == 'membership':
        entry['atom_ids'].pop()
    elif mutation == 'images':
        entry['image_shifts'][str(entry['atom_ids'][0])][0] += 1
    elif mutation == 'source':
        record['force_field_data']['sha256'] = '0' * 64
    else:
        record['input_coordinate_identity'] = '0' * 64
    payload['parameter_assignment']['native_payload']['record'] = record
    payload['typed_graph']['evidence'][0] = digest(record)
    with pytest.raises(InvalidAmberImportResultError):
        validate_backend_record(payload)


@pytest.fixture(scope='module')
def receipt_files(runs, tmp_path_factory):
    import tarfile

    from island.analysis.npt import FLAGS
    from island.analysis.npt_continuation import checksum

    root = tmp_path_factory.mktemp('j34-receipt')
    shutil.copyfile(RUNTIME, root / 'runtime.json')
    shutil.copyfile(PREPARATION, root / PREPARATION.name)
    with tarfile.open(root / 'workflows.tar.gz', 'w:gz') as archive:
        for path in sorted(runs.iterdir()):
            archive.add(path, arcname=path.name)
    matrix = json.loads((runs / 'matrix.json').read_text())
    receipt = api.signed({
        'schema': 'island.j34-receipt.v1',
        'base_commit': '3866e5366092816df2238e96ef35a779ed2dcf5b',
        'implementation_commit': '1' * 40,
        'runtime_identity': json.loads(RUNTIME.read_text())['identity'],
        'cases': matrix['cases'], 'acceptance_met': True,
        'workflow_archive': 'workflows.tar.gz',
        'linear_preparation_archive': PREPARATION.name,
        'bundle_files': {str(p.relative_to(runs)): checksum(p) for p in runs.rglob('*') if p.is_file()},
        'artifacts': {p.name: checksum(p) for p in root.iterdir()}, **FLAGS,
    })
    return root, receipt


def test_real_receipt_acceptance(receipt_files):
    from scripts.validate_phase_4j34_receipt import validate

    root, receipt = receipt_files
    validate(receipt, root, check_runtime=False)


@pytest.mark.parametrize('mutation', ['missing_case', 'blocked_success', 'wrong_case', 'runtime', 'readiness', 'artifact', 'archive_unbound'])
def test_resigned_receipt_mutations_rejected(receipt_files, mutation):
    from scripts.validate_phase_4j34_receipt import validate

    root, original = receipt_files
    receipt = deepcopy(original)
    if mutation == 'missing_case':
        del receipt['cases']['gaff2/linear']
    elif mutation == 'blocked_success':
        receipt['cases']['gaff2/linear']['status'] = 'asset_unavailable'
    elif mutation == 'wrong_case':
        receipt['cases']['gaff2/network']['case'] = 'oplsaa/network'
    elif mutation == 'runtime':
        receipt['runtime_identity'] = '0' * 64
    elif mutation == 'readiness':
        receipt['production_validated'] = True
    elif mutation == 'artifact':
        receipt['artifacts']['workflows.tar.gz'] = '0' * 64
    else:
        del receipt['artifacts']['workflows.tar.gz']
    receipt = api.signed({k: v for k, v in receipt.items() if k != 'identity'})
    with pytest.raises(AnalysisError):
        validate(receipt, root, check_runtime=False)


def test_changed_frozen_preparation_rejected(tmp_path):
    from scripts.phase_4j34_assignment import prepare_linear
    from scripts.phase_4j34_real_fixture import build_fixture

    runtime = json.loads(RUNTIME.read_text())
    root = tmp_path / 'preparation'
    root.mkdir()
    with pytest.raises(AmberToolsUnavailableError, match='archive changed'):
        prepare_linear(build_fixture('linear', 'GAFF2').fixture, root, runtime, runtime['identity'], PREPARATION, '0' * 64)
