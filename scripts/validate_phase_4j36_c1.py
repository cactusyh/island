"""Read-only precision audit of genuine retained native Psi4/RESP outputs.

Reconstruction is not fresh quantum generation and cannot establish MD acceptance.
"""

import argparse
import json
import shutil
import tarfile
import tempfile
from pathlib import Path

from island.forcefields.ambertools.psi4_runtime import validate_runtime
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_precision import (
    audit_in_memory,
    audit_native,
    signed,
)
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input

BASE = 'a53dcd4bc1222dbb3ded245aa52df2f134726f7c'
PARENT = Path('docs/evidence/phase_4j36')
WORKER = Path('scripts/phase_4j36_quantum_reference.py')


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')


def audit_case(root, kind, runtime, parent):
    return audit_native(root, build_fixture(kind, 'GAFF2').fixture.packed.graph,
                        expected_identity=parent['cases'][kind]['execution_identity'],
                        expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))


def run(output, *, offline=False):
    output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'island.j36-c1-outcomes.v1', 'base_commit': BASE,
              'scope': 'read-only reconstruction of fresh J36 native quantum evidence; no new quantum or MD run',
              'acceptance_met': False, 'dynamics_executed': False, 'cases': {}, 'native_esp_cases': {},
              'production_validated': False, 'simulation_readiness': 'not_established',
              'native_am1bcc_success': False, 'external_oplsaa_is_resp': False}
    try:
        runtime = json.loads((PARENT / 'environment.json').read_text())
        parent = json.loads((PARENT / 'receipt.json').read_text())
        if not offline:
            validate_runtime(runtime, expected_identity=parent['runtime_identity'])
        result.update(runtime_identity=runtime['identity'], runtime_validation='offline_record_binding_only' if offline else 'live_hashes_and_probe_complete')
        with tempfile.TemporaryDirectory(prefix='j36-c1-native-') as tmp:
            with tarfile.open(PARENT / 'artifacts/native-quantum.tar.gz') as archive:
                archive.extractall(tmp, filter='data')
            with tarfile.open('docs/evidence/phase_4j36_c1/artifacts/native-memory.tar.gz') as archive:
                archive.extractall(tmp, filter='data')
            for kind in ('linear', 'network'):
                try:
                    target = output / kind
                    shutil.copytree(Path(tmp) / ('quantum-' + kind), target)
                    audit = audit_case(target, kind, runtime, parent)
                    write(target / 'precision-audit.json', audit)
                    memory_root = output / ('memory-' + kind)
                    shutil.copytree(Path(tmp) / ('native-esp-' + kind), memory_root)
                    memory = audit_in_memory(memory_root, build_fixture(kind, 'GAFF2').fixture.packed.graph, native_directory=target,
                                             in_memory_worker_sha256=file_hash('scripts/phase_4j36_c1_native_esp.py'),
                                             expected_identity=parent['cases'][kind]['execution_identity'], expected_input=build_input(kind),
                                             runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))
                    write(memory_root / 'precision-audit.json', memory)
                    result['native_esp_cases'][kind] = {'audit_identity': memory['identity'], 'classification': memory['classification'],
                                                       'charge_conservation_passed': memory['charge_conservation_passed'], 'dynamics_executed': False}
                    result['cases'][kind] = {'audit_identity': audit['identity'], 'classification': audit['classification'],
                                              'stage': audit['stage'], 'charge_conservation_passed': audit['charge_conservation_passed'],
                                              'acceptance_met': False, 'dynamics_executed': False, 'downstream': 'not_executed'}
                except Exception as error:  # noqa: BLE001 -- diagnostic failure is never acceptance
                    result['cases'][kind] = {'classification': getattr(error, 'classification', 'integrity_failure'),
                                              'diagnostic': str(error), 'acceptance_met': False, 'dynamics_executed': False}
                write(output / 'result.json', signed(result))
    except Exception as error:  # noqa: BLE001 -- distinguish runtime assets from source precision
        result.update(classification=getattr(error, 'classification', 'integrity_failure'), diagnostic=str(error))
    write(output / 'result.json', signed(result))
    return signed(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--offline', action='store_true', help='Artifact reconstruction only; no claim of present runtime availability')
    args = parser.parse_args()
    record = run(args.output, offline=args.offline)
    print(json.dumps(record))
    return 0 if record['acceptance_met'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
