"""Software negative controls over independently executed J33-C1 trajectories."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.npt import _body, _signed
from island.analysis.npt_continuation import (
    checksum,
    compare_continuations,
    load_continuation,
    relocate_continuation,
)
from island.exceptions import AnalysisError
from scripts.phase_4j32_real_fixture import FOYER
from scripts.validate_phase_4j33_c1 import validate_case


def child(*args):
    return subprocess.run([str(FOYER), *map(str, args)], env={**os.environ, 'PYTHONPATH': '.:src', 'OPENMM_CPU_THREADS': '1', 'OMP_NUM_THREADS': '1'}, capture_output=True, text=True, check=False)


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    configured = os.environ.get('ISLAND_J33_C1_RUN')
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp('j33-c1') / 'runs'
    process = child('-m', 'scripts.validate_phase_4j33_c1', '--output', root)
    assert process.returncode in (0, 1), process.stderr
    assert (root / 'matrix.json').exists(), process.stdout + process.stderr
    for name in ('oplsaa/linear', 'oplsaa/network', 'pcff/network', 'gaff2/network'):
        validate_case(name, json.loads((root / name.replace('/', '-') / 'result.json').read_text()))
    return root


@pytest.mark.parametrize('kind', ['linear', 'network'])
def test_real_independent_split_relocated_child_resume(runs, kind):
    case = runs / ('oplsaa-' + kind)
    result = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
    assert result['passed']
    assert result['tolerances'] == {'absolute': 0, 'relative': 0, 'checkpoint_and_rng': 'byte_exact'}
    assert result['uninterrupted_trajectory_identity'] == result['resumed_trajectory_identity']
    assert result['checkpoint_hashes'] == result['resumed_checkpoint_hashes']
    assert result['rng_state_hashes'] == result['resumed_rng_state_hashes']
    assert json.loads((case / 'child-command.json').read_text())['exit_code'] == 0
    original = load_continuation(case / 'split')[1]
    moved = load_continuation(case / 'relocated')[1]
    assert original == moved


@pytest.mark.parametrize('field', ['positions_nm', 'box_vectors_nm', 'velocities_nm_ps', 'potential_energy_kj_mol', 'kinetic_energy_kj_mol', 'measured_temperature_kelvin', 'checkpoint_base64', 'stage'])
def test_changed_retained_frames_rejected(runs, tmp_path, field):
    dest = tmp_path / 'resumed'
    shutil.copytree(runs / 'oplsaa-linear/resumed', dest)
    record = json.loads((dest / 'frames.json').read_text())
    frame = record['frames'][-1]
    if field in ('positions_nm', 'box_vectors_nm', 'velocities_nm_ps'):
        frame[field][0][0] += 0.001
    elif field in ('stage', 'checkpoint_base64'):
        frame[field] = 'mutated-software-negative-control'
    else:
        frame[field] += 1
    (dest / 'frames.json').write_text(json.dumps(_signed(_body(record))))
    with pytest.raises((AnalysisError, ValueError)):
        load_continuation(dest)


@pytest.mark.parametrize('mutation', ['density', 'rng_state', 'checkpoint', 'parent', 'configuration', 'source', 'assignment', 'pressure', 'barostat', 'lineage'])
def test_consistently_resigned_mutations_cannot_prove_acceptance(runs, tmp_path, mutation):
    dest = tmp_path / 'resumed'
    shutil.copytree(runs / 'oplsaa-linear/resumed', dest)
    restart_path = dest / 'restart.json'
    restart = json.loads(restart_path.read_text())
    if mutation == 'density':
        path = dest / 'analysis.json'
        analysis = json.loads(path.read_text())
        analysis['density']['mean'] += 1
        path.write_text(json.dumps(_signed(_body(analysis))))
    elif mutation == 'lineage':
        path = dest / 'continuation.json'
        record = json.loads(path.read_text())
        record['parent_restart_identity'] = '0' * 64
        path.write_text(json.dumps(_signed(_body(record))))
    elif mutation in ('pressure', 'barostat'):
        path = dest / 'system.xml'
        raw = path.read_text()
        raw = raw.replace('pressure="1"', 'pressure="2"') if mutation == 'pressure' else raw.replace('MonteCarloBarostat', 'CMMotionRemover')
        path.write_text(raw)
        # Actual model/config checks run in the fresh child below.
    else:
        key = {'rng_state': 'rng_state_sha256', 'checkpoint': 'checkpoint_identity', 'parent': 'parent_workflow_identity', 'configuration': 'npt_config_identity', 'source': 'source_identity', 'assignment': 'assignment_identity'}[mutation]
        restart[key] = '0' * 64
    restart['files'] = {name: checksum(dest / name) for name in restart['files']}
    restart_path.write_text(json.dumps(_signed(_body(restart))))
    if mutation in ('pressure', 'barostat'):
        # Restart the split bundle instead; a completed bundle is not resumable.
        split = tmp_path / 'split'
        shutil.copytree(runs / 'oplsaa-linear/relocated', split)
        shutil.copyfile(dest / 'system.xml', split / 'system.xml')
        record = json.loads((split / 'restart.json').read_text())
        record['files']['system.xml'] = checksum(split / 'system.xml')
        record = _signed(_body(record))
        (split / 'restart.json').write_text(json.dumps(record))
        p = child('-m', 'scripts.validate_phase_4j33_c1', '--resume', split, '--destination', tmp_path / 'bad', '--restart', record['identity'], '--parent', record['parent_workflow_identity'], '--config', record['npt_config_identity'])
        assert p.returncode != 0, p.stdout + p.stderr
    else:
        with pytest.raises((AnalysisError, ValueError)):
            load_continuation(dest)


def test_stale_expected_checkpoint_config_parent_and_relocation(runs, tmp_path):
    source = runs / 'oplsaa-linear/relocated'
    for key in ('expected_restart', 'expected_parent', 'expected_config'):
        with pytest.raises(AnalysisError):
            load_continuation(source, **{key: '0' * 64})
    destination = tmp_path / 'already-exists'
    destination.mkdir()
    with pytest.raises(AnalysisError, match='destination exists'):
        relocate_continuation(source, destination)
    assert (source / 'restart.json').exists()


def test_child_failure_preserves_split_diagnostics(runs, tmp_path):
    source = runs / 'oplsaa-linear/relocated'
    p = child('-m', 'scripts.validate_phase_4j33_c1', '--resume', source, '--destination', tmp_path / 'new', '--restart', '0' * 64, '--parent', '0' * 64, '--config', '0' * 64)
    assert p.returncode != 0
    assert 'Stale checkpoint' in p.stderr
    assert (source / 'frames.json').exists()


@pytest.mark.parametrize('name', ['oplsaa/linear', 'oplsaa/network', 'pcff/network', 'gaff2/network'])
def test_case_outcomes_cannot_be_promoted(runs, name):
    row = json.loads((runs / name.replace('/', '-') / 'result.json').read_text())
    validate_case(name, row)
    changed = deepcopy(row)
    if name.startswith('oplsaa'):
        changed['comparison_passed'] = False
    elif name.startswith('pcff'):
        changed.update(status='complete', npt_executed=True)
    else:
        changed.update(status='complete', acceptance_met=True)
    with pytest.raises(AnalysisError):
        validate_case(name, changed)


def test_unexpected_pcff_success_actual_cli_fails(tmp_path):
    output = tmp_path / 'pcff-success'
    code = '''
import sys
import scripts.validate_phase_4j33_c1 as driver
driver.prepare = lambda *args: object()
sys.argv = ['driver', '--output', sys.argv[1], '--case', 'pcff/network']
raise SystemExit(driver.main())
'''
    process = child('-c', code, output)
    assert process.returncode != 0
    record = json.loads((output / 'result.json').read_text())
    assert record['acceptance_met'] is False and record['npt_executed'] is False
    assert 'Unexpected PCFF preparation success' in (output / 'failure.log').read_text()


def test_fresh_offline_reconstruction_forbids_openmm(runs):
    code = '''
import builtins,sys
old=builtins.__import__
def guard(name,*args,**kw):
 if name.split('.')[0] in {'openmm','foyer','rdkit','parmed'}: raise RuntimeError('offline forbidden '+name)
 return old(name,*args,**kw)
builtins.__import__=guard
from island.analysis.npt_continuation import load_continuation
load_continuation(sys.argv[1])
'''
    p = child('-c', code, runs / 'oplsaa-network/resumed')
    assert p.returncode == 0, p.stdout + p.stderr


@pytest.mark.parametrize('field', ['positions_nm', 'box_vectors_nm', 'volume_nm3', 'velocities_nm_ps', 'potential_energy_kj_mol', 'kinetic_energy_kj_mol', 'measured_temperature_kelvin', 'checkpoint_sha256', 'checkpoint_base64', 'stage'])
def test_comparator_reports_detailed_mismatches(runs, monkeypatch, field):
    # Software negative control of comparison itself; integrity checks above
    # separately reject mutated on-disk records. Never positive MD evidence.
    import island.analysis.npt_continuation as module

    case = runs / 'oplsaa-linear'
    full, restart = load_continuation(case / 'uninterrupted')
    altered = deepcopy(full)
    altered['frames'][-1][field] = 'injected-software-negative-control'
    original = module.load_continuation

    def injected(path, **kwargs):
        if Path(path) == case / 'uninterrupted':
            return altered, restart
        return original(path, **kwargs)

    monkeypatch.setattr(module, 'load_continuation', injected)
    comparison = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
    assert comparison['passed'] is False
    assert {'frame_index': len(full['frames']) - 1, 'step': full['frames'][-1]['step'], 'field': field} in comparison['mismatches']


@pytest.mark.parametrize('field,value', [('stage', 'execution'), ('message', 'unrelated typing failure'), ('origin', 'unrelated.module'), ('exception_type', 'ValueError')])
def test_pcff_wrong_rejection_contract(runs, field, value):
    row = json.loads((runs / 'pcff-network/result.json').read_text())
    row['diagnostic'][field] = value
    with pytest.raises(AnalysisError):
        validate_case('pcff/network', row)


@pytest.fixture(scope='module')
def corrective_receipt(runs, tmp_path_factory):
    import tarfile

    directory = tmp_path_factory.mktemp('corrective-receipt')
    archive = directory / 'continuation.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for path in sorted(runs.rglob('*')):
            if path.is_file():
                tar.add(path, arcname=str(path.relative_to(runs)))
    matrix = json.loads((runs / 'matrix.json').read_text())
    receipt = _signed({
        'schema': 'island.j33-c1-receipt.v1',
        'base_commit': '6f34bf370746ff01fc0491490c0b3dda4a696ec0',
        'implementation_commit': '50048cfd45d542db8489498cb6ebfada7de37879',
        'historical_evidence_preserved': True,
        'artifacts': {'continuation.tar.gz': checksum(archive)},
        'continuation_archive': 'continuation.tar.gz',
        'bundle_files': {str(p.relative_to(runs)): checksum(p) for p in runs.rglob('*') if p.is_file()},
        'cases': matrix['cases'], 'acceptance_met': False,
        'opls_pcff_contract_met': True,
        'production_validated': False, 'simulation_readiness': 'not_established',
        'software_schema_fixture': True,
    })
    return receipt, directory


def test_full_receipt_entrypoint_and_asset_gate(corrective_receipt):
    from scripts.validate_phase_4j33_c1_receipt import validate

    receipt, directory = corrective_receipt
    validate(receipt, directory)
    path = directory / 'receipt.json'
    path.write_text(json.dumps(receipt))
    process = child('-m', 'scripts.validate_phase_4j33_c1_receipt', path)
    assert process.returncode == 0, process.stdout + process.stderr
    process = child('-m', 'scripts.validate_phase_4j33_c1_receipt', path, '--require-all')
    assert process.returncode != 0 and 'All-family acceptance unmet' in process.stdout


@pytest.mark.parametrize('mutation', ['missing_case', 'blocked_success', 'case_identity', 'claimed_acceptance', 'artifact_hash', 'inventory_hash', 'readiness'])
def test_resigned_corrective_receipt_mutations(corrective_receipt, mutation):
    from scripts.validate_phase_4j33_c1_receipt import validate

    receipt, directory = corrective_receipt
    changed = deepcopy(receipt)
    if mutation == 'missing_case':
        del changed['cases']['oplsaa/network']
    elif mutation == 'blocked_success':
        changed['cases']['oplsaa/network'].update(status='blocked', acceptance_met=False)
    elif mutation == 'case_identity':
        changed['cases']['pcff/network']['case'] = 'oplsaa/network'
    elif mutation == 'claimed_acceptance':
        changed['acceptance_met'] = True
    elif mutation == 'artifact_hash':
        changed['artifacts']['continuation.tar.gz'] = '0' * 64
    elif mutation == 'inventory_hash':
        changed['bundle_files']['matrix.json'] = '0' * 64
    else:
        changed['production_validated'] = True
    with pytest.raises(AnalysisError):
        validate(_signed(_body(changed)), directory)
