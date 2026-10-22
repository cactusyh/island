"""J33-C1 real independent split/uninterrupted NPT acceptance and diagnostics."""

import argparse
import os
import shutil
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

from island.analysis.npt import FLAGS, _signed, _write_json, run_npt_workflow
from island.analysis.npt_continuation import (
    compare_continuations,
    load_continuation,
    relocate_continuation,
    require,
    resume_continuation,
)
from island.exceptions import AmberToolsUnavailableError
from scripts.phase_4j31_contract import PCFF_MESSAGE, typing_evidence
from scripts.phase_4j32_assignment import independent_check, prepare
from scripts.phase_4j32_real_fixture import AMBER, FRC
from scripts.phase_4j33_real_fixture import build_j33_fixture
from scripts.validate_phase_4j32 import pcff_contract

CASES = ('oplsaa/linear', 'oplsaa/network', 'pcff/network', 'gaff2/network')


def validate_case(name, row):
    require(name in CASES and row['case'] == name, 'Wrong case inventory/identity')
    require(all(row.get(k) == v for k, v in FLAGS.items()), 'Readiness flags changed')
    if name.startswith('oplsaa/'):
        require(row.get('status') == 'complete' and row.get('acceptance_met') is True, 'Mandatory OPLS continuation unmet')
        for key in ('npt_executed', 'box_volume_changed', 'fresh_process_continuation', 'relocation', 'comparison_passed'):
            require(row.get(key) is True, 'Missing continuation acceptance: ' + key)
        require(row.get('density_blocks') == 4 and row.get('equilibrium_status') == 'not_established', 'Missing four-block finite-run diagnostics')
        require(row.get('density_frame_identity') and row.get('comparison_identity'), 'Missing density/comparison identities')
    elif name == 'pcff/network':
        diagnostic = row.get('diagnostic', {})
        require(row.get('status') == 'blocked' and row.get('acceptance_met') is True and row.get('npt_executed') is False, 'PCFF expected rejection unmet')
        require(diagnostic.get('stage') == 'typing' and diagnostic.get('classification') == 'expected_typing_rejection', 'Wrong PCFF failure stage/class')
        fixture = build_j33_fixture('network', 'PCFF').fixture
        require(diagnostic.get('exception_type') == 'PCFFError' and diagnostic.get('message') == PCFF_MESSAGE and diagnostic.get('origin') == 'island.forcefields.pcff.automatic.bridge_data', 'Wrong PCFF exception/origin')
        require(row.get('final_graph_identity') == fixture.packed.graph.identity and diagnostic.get('typing_evidence') == typing_evidence(fixture, FRC), 'PCFF source/domain/case contract unmet')
    else:
        require(row.get('status') == 'asset_unavailable' and row.get('acceptance_met') is False and row.get('npt_executed') is False and row.get('missing_executables'), 'GAFF2 asset-unavailable control unmet')


def worker(output, name):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    family, kind = name.split('/')
    row = {'case': name, 'status': 'failed', 'acceptance_met': False, 'npt_executed': False, **FLAGS}
    stage = 'fixture'
    try:
        fixture = build_j33_fixture(kind, 'OPLS-AA' if family == 'oplsaa' else family.upper())
        row['final_graph_identity'] = fixture.fixture.packed.graph.identity
        stage = 'preparation'
        if family == 'pcff':
            evidence = typing_evidence(fixture.fixture, FRC)
            try:
                prepare(fixture.fixture, family, output)
            except Exception as error:  # noqa: BLE001 -- frozen precise contract below
                require(pcff_contract(fixture.fixture, evidence, error), 'Unexpected PCFF preparation failure: ' + str(error))
                row.update(status='blocked', acceptance_met=True, diagnostic={'stage': 'typing', 'classification': 'expected_typing_rejection', 'exception_type': type(error).__name__, 'message': str(error), 'origin': 'island.forcefields.pcff.automatic.bridge_data', 'typing_evidence': evidence})
                (output / 'expected-typing-rejection.log').write_text(traceback.format_exc())
            else:
                raise ValueError('Unexpected PCFF preparation success; no NPT is permitted')
        else:
            if family == 'gaff2':
                missing = [str(AMBER / 'bin' / name) for name in ('antechamber', 'parmchk2', 'tleap', 'sqm') if not os.access(AMBER / 'bin' / name, os.X_OK) or (name != 'sqm' and shutil.which(name) is None)]
                if missing:
                    row['missing_executables'] = missing
                    raise AmberToolsUnavailableError('Pinned AmberTools executable(s) missing: ' + ', '.join(missing))
            backend = prepare(fixture.fixture, family, output)
            _write_json(output / 'assignment.json', independent_check(backend))
            config = replace(fixture.npt, source_identity=backend.payload['config']['source']['sha256'])
            stage = 'uninterrupted_execution'
            full = output / 'uninterrupted'
            uninterrupted = run_npt_workflow(backend, config, full, continuation=True)
            row['npt_executed'] = True
            stage = 'split_execution'
            split = output / 'split'
            boundary = config.npt_steps // 2
            run_npt_workflow(backend, config, split, stop_npt_step=boundary)
            stage = 'relocation'
            moved = output / 'relocated'
            restart = relocate_continuation(split, moved)
            stage = 'fresh_process_continuation'
            resumed = output / 'resumed'
            command = [sys.executable, '-m', 'scripts.validate_phase_4j33_c1', '--resume', str(moved), '--destination', str(resumed), '--restart', restart['identity'], '--parent', restart['parent_workflow_identity'], '--config', config.identity]
            child = subprocess.run(command, text=True, capture_output=True, check=False)
            _write_json(output / 'child-command.json', {'argv': command, 'exit_code': child.returncode})
            (output / 'child.log').write_text(child.stdout + child.stderr)
            require(child.returncode == 0, 'Fresh-process continuation failed; see child.log')
            stage = 'comparison'
            comparison = compare_continuations(full, resumed, moved)
            _write_json(output / 'comparison.json', comparison)
            require(comparison['passed'], 'Split/uninterrupted mismatch: ' + str(comparison['mismatches']))
            analysis = uninterrupted['analysis']
            npt = [f['volume_nm3'] for f in uninterrupted['frames'] if f['stage'] == 'npt']
            row.update(status='complete', acceptance_met=True, box_volume_changed=len(set(npt)) > 1, fresh_process_continuation=True, relocation=True, comparison_passed=True, density_blocks=analysis['density']['blocks'], equilibrium_status='not_established', density_frame_identity=analysis['frame_content_identity'], comparison_identity=comparison['identity'], split_boundary_npt_step=boundary, assignment_identity=backend.payload['assignment_identity'])
        validate_case(name, row)
    except Exception as error:  # noqa: BLE001 -- preserve all diagnostic artifacts
        (output / 'failure.log').write_text(traceback.format_exc())
        row.update(status='asset_unavailable' if isinstance(error, AmberToolsUnavailableError) else 'failed', acceptance_met=False, diagnostic={'stage': stage, 'exception_type': type(error).__name__, 'message': str(error)})
    _write_json(output / 'result.json', _signed(row))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--destination', type=Path)
    parser.add_argument('--restart')
    parser.add_argument('--parent')
    parser.add_argument('--config')
    args = parser.parse_args()
    if args.resume:
        try:
            resume_continuation(args.resume, args.destination, expected_restart=args.restart, expected_parent=args.parent, expected_config=args.config)
            load_continuation(args.destination)
        except Exception:  # noqa: BLE001 -- diagnostic exit without suppressing traceback
            traceback.print_exc()
            return 1
        return 0
    require(args.output is not None, '--output required')
    if args.case:
        return int(not worker(args.output, args.case)['acceptance_met'])
    args.output.mkdir(parents=True, exist_ok=False)
    report = {'schema': 'island.j33-c1-matrix.v1', 'cases': {}, 'acceptance_met': False, **FLAGS}
    for name in CASES:
        report['cases'][name] = worker(args.output / name.replace('/', '-'), name)
        _write_json(args.output / 'matrix.json', _signed(report))
    report['acceptance_met'] = all(r['acceptance_met'] for r in report['cases'].values())
    report['opls_pcff_contract_met'] = all(report['cases'][name]['acceptance_met'] for name in CASES if not name.startswith('gaff2'))
    _write_json(args.output / 'matrix.json', _signed(report))
    return int(not report['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
