"""Execute fresh native Psi4/RESP; never substitute historical provided charges.

Each output directory must be new. Reconstruction mode validates genuine native
artifacts without rerunning quantum chemistry. Neither mode treats a failed
exact-charge gate as workflow acceptance.
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest
from island.forcefields.ambertools.psi4_execution import reconstruct
from island.forcefields.ambertools.psi4_runtime import environment, validate_runtime
from island.forcefields.ambertools.resp import file_hash
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input

BASE = '5bd1498745bab65a57ca1d2b752c7c1638aa5a86'
WORKER = Path(__file__).with_name('phase_4j36_quantum_reference.py').resolve()


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def signed(data):
    return {**data, 'identity': reference_digest(data)}


def evaluate(directory, kind, runtime):
    record = reconstruct(directory, build_fixture(kind, 'GAFF2').fixture.packed.graph,
                         expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))
    write(Path(directory) / 'execution-record.json', record)
    # Publishing a RESP provider additionally requires independent review under
    # J35-C1. A quantum vector alone never asserts downstream MD acceptance.
    return {'status': record['status'], 'classification': record['classification'],
            'acceptance_met': False, 'charge_generation_executed': record['charge_generation_executed'],
            'charge_conservation_passed': record['acceptance_met'], 'execution_identity': record['identity'],
            'dynamics_executed': False, 'downstream': 'not_executed',
            'unmet_gate': 'exact_component_charge_conservation' if not record['acceptance_met'] else 'reviewed_RESP_provider_and_downstream_workflow'}


def run(output, runtime, expected_runtime, *, native_results=None):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'island.j36-matrix.v1', 'base_commit': BASE, 'runtime_identity': expected_runtime,
              'acceptance_met': False, 'cases': {}, 'production_validated': False,
              'simulation_readiness': 'not_established', 'native_am1bcc_success': False,
              'external_oplsaa_is_resp': False}
    try:
        validate_runtime(runtime, expected_identity=expected_runtime)
        result['runtime_validation'] = 'complete'
        for kind in ('linear', 'network'):
            root = output / kind
            try:
                if native_results:
                    shutil.copytree(Path(native_results) / ('quantum-' + kind), root)
                else:
                    data = build_input(kind)
                    input_path = output / (kind + '-input.json')
                    write(input_path, data)
                    scratch = output / (kind + '-scratch')
                    scratch.mkdir()
                    argv = [str(Path(runtime['prefix']) / 'bin/python'), str(WORKER), '--input', str(input_path), '--output', str(root)]
                    start = time.monotonic()
                    before = file_hash(WORKER)
                    with (output / (kind + '.log')).open('w') as log:
                        process = subprocess.run(argv, env=environment(runtime['prefix']), cwd=scratch, stdout=log, stderr=subprocess.STDOUT, check=False)
                    write(output / (kind + '-command.json'), {'argv': argv, 'exit_code': process.returncode, 'elapsed_seconds': time.monotonic() - start, 'worker_sha256': before})
                    if before != file_hash(WORKER):
                        raise ValueError('Executed worker source changed during execution')
                    shutil.copy2(WORKER, root / 'worker-source.py')
                result['cases'][kind] = evaluate(root, kind, runtime)
            except Exception as error:  # noqa: BLE001 -- preserve diagnostics and evaluate other safe cases
                result['cases'][kind] = {'status': 'failed', 'classification': getattr(error, 'classification', 'execution_integrity_failure'),
                                         'acceptance_met': False, 'diagnostic': str(error), 'exception_type': type(error).__name__,
                                         'dynamics_executed': False, 'downstream': 'not_executed'}
            write(output / 'result.json', signed(result))
    except Exception as error:  # noqa: BLE001 -- runtime failure is not a charge-domain conclusion
        result.update(classification=getattr(error, 'classification', 'runtime_execution_failure'), diagnostic=str(error), exception_type=type(error).__name__)
    write(output / 'result.json', signed(result))
    return signed(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--runtime-identity', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reconstruct-native', type=Path)
    args = parser.parse_args()
    result = run(args.output, json.loads(args.runtime.read_text()), args.runtime_identity, native_results=args.reconstruct_native)
    print(json.dumps({'acceptance_met': result['acceptance_met'], 'cases': result['cases'], 'diagnostic': result.get('diagnostic')}))
    return 0 if result['acceptance_met'] else 1


if __name__ == '__main__':
    sys.exit(main())
