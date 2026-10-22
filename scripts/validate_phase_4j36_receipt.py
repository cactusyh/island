"""J36 evidence integrity; diagnostic validity never implies RESP acceptance."""

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.psi4_execution import load_execution
from island.forcefields.ambertools.psi4_runtime import SCHEMA as RUNTIME_SCHEMA
from island.forcefields.ambertools.resp import file_hash
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_input import build_input
from scripts.validate_phase_4j36 import BASE, WORKER, evaluate


def signed(record):
    require(record['identity'] == reference_digest({k: v for k, v in record.items() if k != 'identity'}), 'Changed evidence identity')


def validate_matrix(record, directory, runtime):
    signed(record)
    require(record['schema'] == 'island.j36-matrix.v1' and record['base_commit'] == BASE, 'Wrong J36 schema/base')
    require(record['runtime_identity'] == runtime['identity'] and record['runtime_validation'] == 'complete', 'Unverified runtime')
    require(record['production_validated'] is False and record['simulation_readiness'] == 'not_established', 'Changed readiness')
    require(record['acceptance_met'] is False and record['native_am1bcc_success'] is False and record['external_oplsaa_is_resp'] is False, 'Unproven quantum/MD acceptance')
    require(set(record['cases']) == {'linear', 'network'}, 'Missing/extra quantum cases')
    for kind, row in record['cases'].items():
        root = Path(directory) / kind
        saved = json.loads((root / 'execution-record.json').read_text())
        load_execution(root, build_fixture(kind, 'GAFF2').fixture.packed.graph, expected_identity=saved['identity'],
                       expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))
        require(row == evaluate(root, kind, runtime), 'Contradictory case/graph/execution/outcome')
        require(row['classification'] == 'charge_conservation_failed' and row['charge_generation_executed'] is True and row['charge_conservation_passed'] is False, 'Frozen J36 failed gate cannot be relabeled')
    return record


def validate_declarations(receipt, matrix, executions):
    expected_provider = {'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
                         'psi4_version': '1.9.1', 'resp_package_version': '1.0', 'resp_module_version': '0.8',
                         'method': 'HF', 'basis': '6-31G*', 'optimization_method': 'HF', 'optimization_basis': 'STO-3G', 'fit_stages': 1}
    require(receipt['provider'] == expected_provider, 'Changed receipt provider provenance/method/basis')
    require(set(executions) == set(matrix['cases']) == {'linear', 'network'}, 'Incomplete receipt execution inventory')
    expected_cases = {kind: {**row, 'final_graph_identity': executions[kind]['final_graph_identity'],
                             'components': executions[kind]['components']} for kind, row in matrix['cases'].items()}
    require(receipt['cases'] == expected_cases, 'Receipt cases contradict native execution evidence')
    require(receipt['not_executed'] == ['RESP-based GAFF2 preparation', 'minimization', 'NVT', 'NPT',
                                      'split/uninterrupted MD comparison', 'MD checkpoint relocation and fresh-process continuation'], 'Unexecuted downstream stages changed')


def validate(receipt, directory):
    directory = Path(directory).resolve()
    signed(receipt)
    require(receipt['schema'] == 'island.j36-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong J36 receipt/base')
    require(receipt['production_validated'] is False and receipt['simulation_readiness'] == 'not_established', 'Changed readiness')
    require(receipt['acceptance_met'] is False and receipt['classification'] == 'charge_conservation_failed', 'RESP acceptance unmet')
    require(receipt['dynamics_executed'] is False and receipt['native_am1bcc_success'] is False, 'Unexecuted dynamics/native charge claim')
    for name, expected in receipt['artifacts'].items():
        path = (directory / name).resolve()
        require(path.is_relative_to(directory) and file_hash(path) == expected, 'Changed/unsafe artifact: ' + name)
    require({'environment.json', 'environment-explicit.txt', 'artifacts/native-quantum.tar.gz', 'artifacts/reconstruction.tar.gz', 'historical-hashes.json'} <= set(receipt['artifacts']), 'Incomplete J36 evidence inventory')
    runtime = json.loads((directory / 'environment.json').read_text())
    signed(runtime)
    require(runtime['schema'] == RUNTIME_SCHEMA and runtime['identity'] == receipt['runtime_identity'], 'Changed runtime provenance')
    require(runtime['psi4_version'] == '1.9.1' and runtime['resp_package_version'] == '1.0' and runtime['resp_module_version'] == '0.8', 'Wrong native runtime')
    with tempfile.TemporaryDirectory(prefix='j36-receipt-') as tmp:
        with tarfile.open(directory / 'artifacts/reconstruction.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        matrix = json.loads((Path(tmp) / 'result.json').read_text())
        validate_matrix(matrix, tmp, runtime)
        require(matrix['identity'] == receipt['matrix_identity'], 'Changed matrix binding')
    with tempfile.TemporaryDirectory(prefix='j36-native-receipt-') as tmp:
        with tarfile.open(directory / 'artifacts/native-quantum.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        executions = {}
        for kind, row in matrix['cases'].items():
            root = Path(tmp) / ('quantum-' + kind)
            executions[kind] = load_execution(root, build_fixture(kind, 'GAFF2').fixture.packed.graph, expected_identity=row['execution_identity'],
                           expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash(WORKER))
    validate_declarations(receipt, matrix, executions)
    inventory = json.loads((directory / 'historical-hashes.json').read_text())
    require(receipt['historical_file_count'] == len(inventory), 'Changed preservation count')
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json'], text=True).splitlines()
    require(set(names) == set(inventory), 'Incomplete historical evidence inventory')
    for name, expected in inventory.items():
        require(file_hash(name) == expected and hashlib.sha256(subprocess.check_output(['git', 'show', BASE + ':' + name])).hexdigest() == expected, 'Historical bytes changed: ' + name)
    require(set(receipt['tests']) == {'focused', 'j22-j35', 'full_supported'}, 'Missing verification suites')
    from scripts.validate_phase_4j34_c2_receipt import junit_counts
    for gate in receipt['tests'].values():
        require({gate['junit'], gate['command']} <= set(receipt['artifacts']), 'Unbound test evidence')
        counts = junit_counts(directory / gate['junit'])
        command = json.loads((directory / gate['command']).read_text())
        require(counts == gate['counts'] and counts['failures'] == counts['errors'] == 0 and counts['passed'] > 0, 'Changed or failed test evidence')
        require(command['exit_code'] == 0 and command['implementation_commit'] == receipt['implementation_commit'], 'Wrong tested implementation')
    required = {'src/island/forcefields/ambertools/psi4_runtime.py', 'src/island/forcefields/ambertools/psi4_execution.py', 'scripts/validate_phase_4j36.py', 'scripts/validate_phase_4j36_receipt.py', 'scripts/phase_4j36_quantum_reference.py', 'scripts/phase_4j36_reference.py', 'scripts/phase_4j36_input.py', 'tests/test_phase_4j36.py'}
    require(required <= set(receipt['implementation_files']), 'Missing implementation bindings')
    for name, expected in receipt['implementation_files'].items():
        require(file_hash(name) == expected and hashlib.sha256(subprocess.check_output(['git', 'show', receipt['implementation_commit'] + ':' + name])).hexdigest() == expected, 'Tested implementation changed: ' + name)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()), args.receipt.parent)
    except Exception as error:  # noqa: BLE001 -- explicit diagnostic gate
        print('FAIL:', error)
        return 1
    print('VALID FAILED-ACCEPTANCE EVIDENCE: exact native RESP component totals fail; GAFF2 MD unexecuted')
    return 0 if args.diagnostic_only else 1


if __name__ == '__main__':
    raise SystemExit(main())
