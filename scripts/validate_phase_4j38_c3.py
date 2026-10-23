"""Fresh exact-graph RESP/Psi4 generation and GAFF2 benchmark correction.

The driver exits nonzero while any benchmark cell remains blocked. Individual
software execution success is distinct from scientific validation.
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

from island.analysis.npt import (
    FLAGS,
    NPTConfiguration,
    _signed,
    _write_json,
    run_npt_workflow,
)
from island.analysis.npt_continuation import relocate_continuation
from island.forcefields.ambertools.provided import require
from island.forcefields.ambertools.psi4_runtime import environment, validate_runtime
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from island.forcefields.ambertools.runtime import declared_environment
from scripts.phase_4j32_assignment import independent_check
from scripts.phase_4j38_c3_charges import charges, write
from scripts.phase_4j38_c3_input import build_input
from scripts.phase_4j38_c3_reference import audit, sha
from scripts.phase_4j38_contract import base_row, gate
from scripts.phase_4j38_fixture import build
from scripts.validate_phase_4j36_c2 import prepare
from scripts.validate_phase_4j38 import finish, input_provenance, postprocess
from scripts.validate_phase_4j38_c2 import run as recheck

BASE = '3d298cccd7f430340c82d8e069f387a352e4f19b'
KINDS = ('pe', 'aromatic')
QUANTUM = Path('docs/evidence/phase_4j36/environment.json')


def command(argv, stem, *, env=None):
    started = time.time()
    with stem.with_suffix('.log').open('w') as log:
        child = subprocess.run(argv, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    write(stem.with_suffix('.command.json'), {
        'argv': argv, 'exit_code': child.returncode, 'elapsed_seconds': time.time() - started,
        'implementation_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    })
    return child.returncode


def worker(output, kind, raw):
    """No dynamics before independent exact-graph source and policy validation."""
    output.mkdir(parents=True, exist_ok=False)
    row = base_row('gaff2/' + kind)
    stage = 'charge_source_audit'
    try:
        fixture, evidence = build(kind)
        _write_json(output / 'fixture.json', evidence)
        row['input_provenance'] = input_provenance(fixture, evidence)
        runtime = json.loads(QUANTUM.read_text())
        require(json.loads((raw / 'anchors.json').read_text()) == json.loads(Path('scripts/phase_4j38_c3_anchors.json').read_text()), 'Unreviewed exact-graph source anchors; review fresh quantum results before MD')
        saved = json.loads((raw / (kind + '-charges.json')).read_text())
        review = audit(raw, kind, saved, runtime, json.loads((raw / 'anchors.json').read_text()))
        require(review == json.loads((raw / (kind + '-review.json')).read_text()), 'Changed independent source review')
        provided = FiniteRESPCharges.load(raw / (kind + '-charges.json'), fixture.packed.graph, expected_identity=review['charge_identity'])
        provided.save(output / 'resp-charges.json')
        write(output / 'raw-review.json', review)
        row.update(provided_charge_identity=provided.identity, independent_review_identity=review['identity'], native_am1bcc_success=False)
        row['gates']['source_completeness'] = gate('passed', 'Fresh exact-graph native RESP/Psi4 arrays independently audited; unchanged J36-C2 binary64 conservation policy')
        stage = 'gaff2_preparation'
        amber = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
        backend = prepare(fixture, output, amber, provided)
        write(output / 'independent-assignment.json', independent_check(backend))
        q = provided.validate(fixture.packed.graph, expected_reference_identity=provided.identity)
        error = max(abs(s['charge'] - q[s['id']]) for s in backend.payload['numerical_data']['sites'])
        require(error <= 1e-14, 'Native PRMTOP charge roundtrip exceeds retained serialization bound')
        write(output / 'charge-roundtrip.json', {'authoritative_vector_unchanged': True, 'max_native_prmtop_roundtrip_error_e': error, 'existing_serialization_bound_e': 1e-14})
        config = NPTConfiguration(family='GAFF2', source_identity=backend.payload['config']['source']['sha256'],
            final_graph_identity=fixture.packed.graph.identity, packing_identity=fixture.packed.plan.identity,
            box_lengths_nm=tuple(x / 10 for x in fixture.packed.system.box.lengths), nvt_steps=64, npt_steps=256)
        stage = 'uninterrupted'
        run_npt_workflow(backend, config, output / 'uninterrupted', continuation=True)
        stage = 'independent_reference'
        from scripts.phase_4j37_physics import check
        write(output / 'physics.json', check(output))
        row['gates']['numerical_reference'] = gate('passed', 'Independent AmberPrmtopFile periodic Ewald energy/force/components; same OpenMM Reference engine, no LAMMPS KSPACE claim')
        stage = 'split'
        run_npt_workflow(backend, config, output / 'split', stop_npt_step=config.npt_steps // 2)
        stage = 'relocation'
        restart = relocate_continuation(output / 'split', output / 'relocated')
        stage = 'child_resume'
        argv = [sys.executable, '-m', 'scripts.validate_phase_4j33_c1', '--resume', str(output / 'relocated'),
                '--destination', str(output / 'resumed'), '--restart', restart['identity'],
                '--parent', restart['parent_workflow_identity'], '--config', config.identity]
        require(command(argv, output / 'resume') == 0, 'Fresh-process resume failed; see resume.log')
        stage = 'analysis'
        row = postprocess(output, row, fixture, evidence)
    except Exception as error:  # noqa: BLE001 -- preserve partial artifacts and fail the cell
        (output / 'failure.log').write_text(traceback.format_exc())
        row['diagnostic'] = {'stage': stage, 'exception': type(error).__name__, 'message': str(error)}
        if stage == 'charge_source_audit' and isinstance(error, FileNotFoundError):
            row['classification'] = 'asset_unavailable'
            row['gates']['source_completeness'] = gate('asset_unavailable', str(error))
        row['gates']['software_execution'] = gate('failed', str(error))
    return finish(row, output)


def matrix(previous, outcomes):
    """Additive update: native OPLS decision and retained cells remain frozen."""
    from copy import deepcopy
    value = deepcopy(previous)
    value.pop('identity')
    value.update(schema='island.j38-c3-matrix.v1', base_commit=BASE,
                 c2_matrix_identity=previous['identity'], new_md_executions=0)
    for kind, result in outcomes.items():
        row = value['cases']['gaff2/' + kind]
        row['c2_row_identity'] = previous['cases']['gaff2/' + kind].get('identity')
        passed = result['execution_met']
        source = result['gates']['source_completeness']['status'] == 'passed'
        row.update(charge_source_status='complete' if source else result.get('classification', 'source_incomplete'),
                   execution_origin='fresh_j38_c3' if passed else 'not_executed', new_dynamics_executed=passed,
                   c3_result_identity=result['identity'], provided_charge_identity=result.get('provided_charge_identity'),
                   independent_review_identity=result.get('independent_review_identity'))
        row['source_artifacts']['status'] = 'present' if source else 'asset_unavailable'
        row['charge_paths']['resp_psi4_exact_graph'].update(status=row['charge_source_status'], exact_graph_verified=source,
                                                         review_identity=result.get('independent_review_identity'))
        for name in ('typing', 'charge_completeness', 'parameter_completeness', 'executable_backend', 'numerical_reference'):
            okay = source if name == 'charge_completeness' else passed
            status = 'passed' if okay else ('asset_unavailable' if name == 'charge_completeness' and row['charge_source_status'] == 'asset_unavailable' else 'not_executed')
            row['stages'][name] = gate(status,
                                      'Fresh exact-graph C3 evidence' if okay else str(result.get('diagnostic', 'Charge gate blocks execution')))
        row['stages']['equilibration'] = gate('insufficient_sampling' if passed else 'not_executed', 'One short software acceptance trajectory; no thermodynamic equilibrium established')
        value['new_md_executions'] += int(passed)
    value['unmet_charge_cells'] = [key for key, row in value['cases'].items() if row['charge_source_status'] not in ('complete', 'retained_verified')]
    value['acceptance_met'] = False
    return _signed(value)


def run(output, *, reviewed=None):
    output.mkdir(parents=True, exist_ok=False)
    recheck(output / 'c2-recheck')
    raw = output / 'raw'; raw.mkdir()
    runtime = json.loads(QUANTUM.read_text())
    outcomes = {}
    try:
        validate_runtime(runtime, expected_identity=runtime['identity'])
        write(raw / 'runtime.json', runtime)
        if reviewed is not None:
            shutil.rmtree(raw)
            shutil.copytree(reviewed, raw)
            anchors = json.loads((raw / 'anchors.json').read_text())
            require(anchors == json.loads(Path('scripts/phase_4j38_c3_anchors.json').read_text()), 'Unreviewed exact-graph source anchors')
        else:
            for kind in KINDS:
                write(raw / (kind + '-input.json'), build_input(kind))
                for worker_name, args, target in (
                    ('phase_4j36_quantum_reference.py', ['--input', str(raw / (kind + '-input.json'))], 'quantum-' + kind),
                    ('phase_4j38_c3_native_esp.py', ['--native', str(raw / ('quantum-' + kind))], 'native-esp-' + kind),
                ):
                    script = Path('scripts') / worker_name
                    argv = [runtime['prefix'] + '/bin/python', str(script.resolve()), *args, '--output', str(raw / target)]
                    command(argv, raw / (target + '-command'), env=environment(runtime['prefix']))
                    if (raw / target).is_dir():
                        shutil.copyfile(script, raw / target / 'worker-source.py')
                charges(raw, kind, runtime)
            anchors = {'schema': 'island.j38-c3-raw-anchors.v1', 'runtime_sha256': sha(QUANTUM), 'cases': {}}
            for kind in KINDS:
                anchors['cases'][kind] = {'files': {str(f.relative_to(raw)): sha(f) for prefix in ('quantum-', 'native-esp-') for f in sorted((raw / (prefix + kind)).rglob('*')) if f.is_file()}, 'charges_sha256': sha(raw / (kind + '-charges.json'))}
            write(raw / 'anchors.json', anchors)
        for kind in KINDS:
            review = audit(raw, kind, json.loads((raw / (kind + '-charges.json')).read_text()), runtime, anchors)
            write(raw / (kind + '-review.json'), review)
            amber = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
            argv = [sys.executable, '-m', 'scripts.validate_phase_4j38_c3', '--output', str(output / ('gaff2-' + kind)), '--worker', kind, '--raw', str(raw)]
            command(argv, output / (kind + '-workflow'), env=declared_environment(amber))
            outcomes[kind] = json.loads((output / ('gaff2-' + kind) / 'result.json').read_text())
    except Exception as error:  # noqa: BLE001 -- never invent a charge or execute MD after a source failure
        (output / 'failure.log').write_text(traceback.format_exc())
        for kind in KINDS:
            if kind not in outcomes:
                row = base_row('gaff2/' + kind)
                row['classification'] = getattr(error, 'classification', 'asset_unavailable' if isinstance(error, FileNotFoundError) else 'source_incomplete')
                row['diagnostic'] = {'stage': 'exact_graph_quantum_source', 'message': str(error)}
                row['gates']['source_completeness'] = gate('asset_unavailable' if row['classification'] == 'asset_unavailable' else 'failed', str(error))
                row.update(execution_met=False, software_acceptance_met=False)
                outcomes[kind] = _signed(row)
    write(output / 'outcomes.json', _signed({'schema': 'island.j38-c3-outcomes.v1', 'cases': outcomes, 'acceptance_met': False, **FLAGS}))
    result = matrix(json.loads((output / 'c2-recheck/matrix.json').read_text()), outcomes)
    write(output / 'matrix.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--from-reviewed', type=Path, help='Replay the independently audited C3 raw bundle; never a historical different-graph charge vector')
    parser.add_argument('--worker', choices=KINDS, help=argparse.SUPPRESS)
    parser.add_argument('--raw', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        return int(not worker(args.output.resolve(), args.worker, args.raw.resolve())['software_acceptance_met'])
    return int(not run(args.output.resolve(), reviewed=args.from_reviewed)['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
