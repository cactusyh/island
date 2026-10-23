"""RESP-only finite precision audit followed by fresh GAFF2 acceptance."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tarfile
import traceback
from dataclasses import replace
from pathlib import Path

from island.analysis.npt import FLAGS, _signed, _write_json, run_npt_workflow
from island.analysis.npt_continuation import (
    compare_continuations,
    relocate_continuation,
    require,
)
from island.forcefields import adopt_forcefield
from island.forcefields.ambertools.provided_preparation import parameterize_provided
from island.forcefields.ambertools.psi4_runtime import (
    validate_runtime as validate_quantum,
)
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges, create
from island.forcefields.ambertools.resp_precision import audit_in_memory
from island.forcefields.ambertools.runtime import declared_environment, validate_runtime
from island.periodic import (
    PeriodicForceFieldConfig,
    periodic_chemical_view,
    prepare_periodic_forcefield,
)
from scripts.phase_4j32_assignment import independent_check
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input
from scripts.validate_phase_4j35 import validate_case
from scripts.validate_phase_4j36_c1 import audit_case

BASE = '2a62eeb0b26e8225a0d46db00ba02f44f1f386ba'
PARENT = Path('docs/evidence/phase_4j36_c1')


def charge_inputs(output):
    quantum = json.loads(Path('docs/evidence/phase_4j36/environment.json').read_text())
    parent = json.loads(Path('docs/evidence/phase_4j36/receipt.json').read_text())
    validate_quantum(quantum, expected_identity=parent['runtime_identity'])
    native = output / 'quantum-evidence';native.mkdir()
    for name in ('docs/evidence/phase_4j36/artifacts/native-quantum.tar.gz', str(PARENT / 'artifacts/native-memory.tar.gz')):
        with tarfile.open(name) as archive:
            archive.extractall(native, filter='data')
    policies = {}
    for kind in ('linear', 'network'):
        graph = build_fixture(kind, 'GAFF2').fixture.packed.graph
        original = audit_case(native / ('quantum-' + kind), kind, quantum, parent)
        memory = audit_in_memory(native / ('native-esp-' + kind), graph, native_directory=native / ('quantum-' + kind),
                                 in_memory_worker_sha256=file_hash('scripts/phase_4j36_c1_native_esp.py'), expected_identity=parent['cases'][kind]['execution_identity'],
                                 expected_input=build_input(kind), runtime_identity=quantum['identity'], worker_sha256=file_hash('scripts/phase_4j36_quantum_reference.py'))
        for label, audit, directory in ((kind, original, native / ('quantum-' + kind)), ('memory-' + kind, memory, native / ('native-esp-' + kind))):
            require(audit == json.loads((PARENT / (label + '.json')).read_text()), 'Native evidence differs from reviewed C1 anchor')
            auditpath = output / (label + '-audit.json');_write_json(auditpath, audit)
            reference = output / (label + '-independent.json')
            command = [sys.executable, '-m', 'scripts.phase_4j36_c2_reference', '--directory', str(directory), '--audit', str(auditpath), '--output', str(reference)]
            child = subprocess.run(command, capture_output=True, text=True, check=False)
            _write_json(output / (label + '-reference-command.json'), {'argv': command, 'exit_code': child.returncode})
            (output / (label + '-reference.log')).write_text(child.stdout + child.stderr)
            require(child.returncode == 0, 'Independent raw charge parsing failed')
            charges = create(audit, graph, quantum, json.loads(reference.read_text()))
            charges.validate(graph, expected_reference_identity=charges.identity)
            charges.save(output / (label + '-charges.json'))
            policies[label] = charges.identity
    return policies


def prepare(fixture, root, runtime, provided):
    native = parameterize_provided(fixture.packed.system, fixture.packed.graph, provided, root / 'native-components', runtime=runtime, runtime_identity=runtime['identity'], expected_reference_identity=provided.identity)
    prepared = adopt_forcefield(periodic_chemical_view(fixture.packed.system), 'gaff2', native)
    imported = native.imported_result
    config = PeriodicForceFieldConfig(family='GAFF2', source={'family': 'GAFF2', 'sha256': imported.source_sha256, 'source': imported.source}, box_lengths=fixture.packed.system.box.lengths, nonbonded_cutoff=5.0, lj_mixing_rule='lorentz_berthelot_12_6', one_four_scaling=(imported.nonbonded_policy.lj_scale_14, imported.nonbonded_policy.coulomb_scale_14), electrostatics_method='ewald', pme_tolerance=1e-5,
        provenance='J36-C2 authoritative RESP/Psi4 charges; finite-binary64 policy, unchanged vector; no AM1-BCC or OPLS-AA substitution', evidence=(fixture.packed.graph.identity, fixture.packed.plan.identity, runtime['identity'], provided.identity))
    backend = prepare_periodic_forcefield(fixture.packed.system, fixture.packed.graph, prepared, config)
    shutil.copytree(root / 'native-components', root / 'source-assets/amber')
    return backend


def worker(output, kind, charge_path, expected):
    output.mkdir(parents=True, exist_ok=False)
    row = {'case': 'gaff2/' + kind, 'status': 'failed', 'acceptance_met': False, 'npt_executed': False, **FLAGS}
    stage = 'runtime'
    try:
        runtime = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
        validate_runtime(runtime, expected_identity=runtime['identity'], environment=os.environ)
        fixture = build_fixture(kind, 'GAFF2')
        stage = 'charge_validation'
        provided = FiniteRESPCharges.load(charge_path, fixture.fixture.packed.graph, expected_identity=expected)
        provided.save(output / 'resp-charges.json')
        row.update(runtime_identity=runtime['identity'], final_graph_identity=fixture.fixture.packed.graph.identity,
                   charge_diagnostic='provided_charge_validated', native_am1bcc_success=False,
                   provided_charge_identity=provided.identity, charge_reference_identity=provided.identity,
                   charge_source_identity=provided.payload['source_audit']['identity'], charge_origin='resp_psi4', charge_method='RESP/Psi4',
                   policy_identity=provided.payload['policy']['identity'])
        stage = 'fresh_gaff2_preparation'
        backend = prepare(fixture.fixture, output, runtime, provided)
        _write_json(output / 'independent-assignment.json', independent_check(backend))
        q = provided.validate(fixture.fixture.packed.graph, expected_reference_identity=expected)
        errors = [abs(s['charge'] - q[s['id']]) for s in backend.payload['numerical_data']['sites']]
        _write_json(output / 'charge-roundtrip.json', {'authoritative_vector_unchanged': True, 'policy_identity': row['policy_identity'], 'max_native_prmtop_roundtrip_error_e': max(errors), 'existing_serialization_bound_e': 1e-14})
        require(max(errors) <= 1e-14, 'Parameter charge serialization changed authoritative values beyond existing bound')
        config = replace(fixture.npt, source_identity=backend.payload['config']['source']['sha256'])
        stage = 'uninterrupted'
        full = output / 'uninterrupted';workflow = run_npt_workflow(backend, config, full, continuation=True)
        row['npt_executed'] = True
        stage = 'split'
        split = output / 'split';boundary = config.npt_steps // 2
        run_npt_workflow(backend, config, split, stop_npt_step=boundary)
        stage = 'relocation'
        relocated = output / 'relocated';restart = relocate_continuation(split, relocated)
        stage = 'child_resume'
        resumed = output / 'resumed'
        command = [sys.executable, '-m', 'scripts.validate_phase_4j33_c1', '--resume', str(relocated), '--destination', str(resumed), '--restart', restart['identity'], '--parent', restart['parent_workflow_identity'], '--config', config.identity]
        child = subprocess.run(command, capture_output=True, text=True, check=False)
        _write_json(output / 'child-command.json', {'argv': command, 'exit_code': child.returncode});(output / 'child.log').write_text(child.stdout + child.stderr)
        require(child.returncode == 0, 'Fresh-process continuation failed')
        stage = 'comparison'
        comparison = compare_continuations(full, resumed, relocated);_write_json(output / 'comparison.json', comparison)
        require(comparison['passed'], 'Split/uninterrupted mismatch: ' + str(comparison['mismatches']))
        row.update(status='complete', acceptance_met=True, comparison_passed=True, relocation=True, fresh_process_continuation=True,
                   box_volume_changed=len({f['volume_nm3'] for f in workflow['frames'] if f['stage'] == 'npt'}) > 1,
                   density_blocks=workflow['analysis']['density']['blocks'], equilibrium_status='not_established',
                   density_frame_identity=comparison['density_frame_identity'], comparison_identity=comparison['identity'],
                   split_boundary_npt_step=boundary, assignment_identity=backend.payload['assignment_identity'], runtime_binding_identity=provided.identity)
        validate_case('gaff2/' + kind, row)
    except Exception as error:  # noqa: BLE001 -- retain failures, never domain success
        (output / 'failure.log').write_text(traceback.format_exc())
        row.update(status='failed', acceptance_met=False, diagnostic={'stage': stage, 'type': type(error).__name__, 'message': str(error)})
    _write_json(output / 'result.json', _signed(row))
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--audit-only', action='store_true')
    p.add_argument('--worker', choices=('linear', 'network'), help=argparse.SUPPRESS)
    p.add_argument('--charges', type=Path, help=argparse.SUPPRESS)
    p.add_argument('--expected', help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.worker:
        return int(not worker(args.output, args.worker, args.charges, args.expected)['acceptance_met'])
    args.output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'island.j36-c2-outcomes.v1', 'base_commit': BASE, 'acceptance_met': False, 'cases': {}, **FLAGS}
    try:
        identities = charge_inputs(args.output)
        result['charge_identities'] = identities
        if args.audit_only:
            result['scope'] = 'charge_policy_only_no_workflow_acceptance'
        else:
            runtime = json.loads(Path('docs/evidence/phase_4j34/runtime.json').read_text())
            for kind in ('linear', 'network'):
                command = [sys.executable, '-m', 'scripts.validate_phase_4j36_c2', '--output', str(args.output / kind), '--worker', kind, '--charges', str(args.output / ('memory-' + kind + '-charges.json')), '--expected', identities['memory-' + kind]]
                child = subprocess.run(command, env=declared_environment(runtime), capture_output=True, text=True, check=False)
                _write_json(args.output / (kind + '-command.json'), {'argv': command, 'exit_code': child.returncode})
                (args.output / (kind + '.log')).write_text(child.stdout + child.stderr)
                casepath = args.output / kind / 'result.json'
                result['cases'][kind] = json.loads(casepath.read_text()) if casepath.is_file() else {'acceptance_met': False, 'diagnostic': 'Worker failed without result'}
                _write_json(args.output / 'result.json', _signed(result))
            result['acceptance_met'] = all(result['cases'][k]['acceptance_met'] for k in ('linear', 'network'))
    except Exception as error:  # noqa: BLE001 -- preserve partial matrix and gate failure
        (args.output / 'failure.log').write_text(traceback.format_exc());result['diagnostic'] = str(error)
    _write_json(args.output / 'result.json', _signed(result))
    return int(not result['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
