"""J34 all-family finite-run periodic NPT acceptance with declared GAFF2 assets."""

import argparse
import json
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

from island.analysis.npt import FLAGS, _signed, _write_json, run_npt_workflow
from island.analysis.npt_continuation import (
    compare_continuations,
    relocate_continuation,
    require,
)
from island.exceptions import AmberToolsUnavailableError
from island.forcefields.ambertools.runtime import declared_environment, validate_runtime
from scripts.phase_4j32_assignment import independent_check, prepare
from scripts.phase_4j34_assignment import prepare_linear
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j33_c1 import validate_case as validate_c1_case
from scripts.validate_phase_4j33_c1 import worker as pcff_worker

CASES = ('oplsaa/linear', 'oplsaa/network', 'gaff2/linear', 'gaff2/network', 'pcff/network')
SUCCESS = CASES[:-1]


def validate_case(name, row):
    require(name in CASES and row.get('case') == name, 'Wrong J34 case identity')
    require(all(row.get(k) == v and type(row.get(k)) is type(v) for k, v in FLAGS.items()), 'Changed readiness flags')
    if name == 'pcff/network':
        validate_c1_case(name, row)
        return
    require(row.get('status') == 'complete' and row.get('acceptance_met') is True, 'Required NPT success unmet: ' + name)
    for field in ('npt_executed', 'comparison_passed', 'relocation', 'fresh_process_continuation', 'box_volume_changed'):
        require(row.get(field) is True, 'Missing NPT evidence: ' + field)
    require(row.get('density_blocks') == 4 and row.get('equilibrium_status') == 'not_established', 'Missing finite-run density evidence')
    for key in ('final_graph_identity', 'assignment_identity', 'comparison_identity', 'density_frame_identity', 'runtime_binding_identity'):
        require(isinstance(row.get(key), str) and len(row[key]) == 64, 'Missing J34 binding: ' + key)
    if name.startswith('gaff2/'):
        require(isinstance(row.get('runtime_identity'), str) and len(row['runtime_identity']) == 64, 'Missing GAFF2 runtime binding')


def worker(output, name, runtime, expected_runtime, environment, retained=None, retained_sha256=None):
    output = Path(output)
    if name == 'pcff/network':
        return pcff_worker(output, name)
    output.mkdir(parents=True, exist_ok=False)
    row = {'case': name, 'status': 'failed', 'acceptance_met': False, 'npt_executed': False, **FLAGS}
    stage = 'runtime_validation'
    try:
        family, kind = name.split('/')
        if family == 'gaff2':
            validate_runtime(runtime, expected_identity=expected_runtime, environment=environment)
            row['runtime_identity'] = runtime['identity']
        stage = 'fixture'
        fixture = build_fixture(kind, 'GAFF2' if family == 'gaff2' else 'OPLS-AA')
        row['final_graph_identity'] = fixture.fixture.packed.graph.identity
        _write_json(output / 'fixture.json', fixture.fixture.evidence())
        stage = 'preparation'
        backend = prepare_linear(fixture.fixture, output, runtime, expected_runtime, retained, retained_sha256) if name == 'gaff2/linear' else prepare(fixture.fixture, family, output)
        _write_json(output / 'independent-assignment.json', independent_check(backend))
        config = replace(fixture.npt, source_identity=backend.payload['config']['source']['sha256'])
        binding = _signed({'schema': 'island.j34-runtime-workflow-binding.v1', 'runtime_identity': runtime['identity'] if family == 'gaff2' else None, 'family': config.family, 'final_graph_identity': backend.payload['final_graph_identity'], 'assignment_identity': backend.payload['assignment_identity'], 'source_identity': config.source_identity, 'backend_identity': backend.payload['identity'], 'npt_config_identity': config.identity, **FLAGS})
        _write_json(output / 'runtime-binding.json', binding)
        stage = 'uninterrupted_execution'
        full = output / 'uninterrupted'
        workflow = run_npt_workflow(backend, config, full, continuation=True)
        row['npt_executed'] = True
        stage = 'split_execution'
        split = output / 'split'
        boundary = config.npt_steps // 2
        run_npt_workflow(backend, config, split, stop_npt_step=boundary)
        stage = 'relocation'
        relocated = output / 'relocated'
        restart = relocate_continuation(split, relocated)
        stage = 'fresh_process_continuation'
        resumed = output / 'resumed'
        command = [sys.executable, '-m', 'scripts.validate_phase_4j33_c1', '--resume', str(relocated), '--destination', str(resumed), '--restart', restart['identity'], '--parent', restart['parent_workflow_identity'], '--config', config.identity]
        child = subprocess.run(command, capture_output=True, text=True, check=False, env=environment)
        _write_json(output / 'child-command.json', {'argv': command, 'exit_code': child.returncode})
        (output / 'child.log').write_text(child.stdout + child.stderr)
        require(child.returncode == 0, 'Fresh-process continuation failed; see child.log')
        stage = 'comparison'
        comparison = compare_continuations(full, resumed, relocated)
        _write_json(output / 'comparison.json', comparison)
        require(comparison['passed'], 'Uninterrupted/split mismatch: ' + str(comparison['mismatches']))
        volumes = [f['volume_nm3'] for f in workflow['frames'] if f['stage'] == 'npt']
        row.update(status='complete', acceptance_met=True, comparison_passed=True, relocation=True, fresh_process_continuation=True, box_volume_changed=len(set(volumes)) > 1, density_blocks=workflow['analysis']['density']['blocks'], equilibrium_status='not_established', density_frame_identity=comparison['density_frame_identity'], comparison_identity=comparison['identity'], split_boundary_npt_step=boundary, assignment_identity=backend.payload['assignment_identity'], runtime_binding_identity=binding['identity'])
        validate_case(name, row)
    except Exception as error:  # noqa: BLE001 -- retain every failed acceptance diagnostic
        (output / 'failure.log').write_text(traceback.format_exc())
        row.update(status='asset_unavailable' if isinstance(error, AmberToolsUnavailableError) else 'failed', acceptance_met=False, diagnostic={'stage': stage, 'exception_type': type(error).__name__, 'message': str(error)})
    _write_json(output / 'result.json', _signed(row))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--runtime-identity', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--linear-preparation', type=Path, help='Explicit frozen native preparation for this exact J34 graph; never a fallback')
    parser.add_argument('--linear-preparation-sha256')
    parser.add_argument('--activated-child', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        runtime = json.loads(args.runtime.read_text())
    except (OSError, ValueError):
        runtime = {'amberhome': '/missing-declared-gaff2-runtime'}
    environment = declared_environment(runtime)
    if not args.activated_child:
        # Activation is explicit and scoped to a child; historical validators and
        # the caller's environment keep their original executable discovery.
        return subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j34', *sys.argv[1:], '--activated-child'], env=environment, check=False).returncode
    if args.case:
        return int(not worker(args.output, args.case, runtime, args.runtime_identity, environment, args.linear_preparation, args.linear_preparation_sha256)['acceptance_met'])
    args.output.mkdir(parents=True, exist_ok=False)
    _write_json(args.output / 'runtime.json', runtime)
    matrix = {'schema': 'island.j34-acceptance.v1', 'cases': {}, 'acceptance_met': False, 'runtime_identity': args.runtime_identity, **FLAGS}
    for name in CASES:
        matrix['cases'][name] = worker(args.output / name.replace('/', '-'), name, runtime, args.runtime_identity, environment, args.linear_preparation, args.linear_preparation_sha256)
        _write_json(args.output / 'matrix.json', _signed(matrix))
    matrix['acceptance_met'] = all(row['acceptance_met'] for row in matrix['cases'].values())
    _write_json(args.output / 'matrix.json', _signed(matrix))
    return int(not matrix['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
