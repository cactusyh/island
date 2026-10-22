"""Validate blocked C2 evidence; default nonzero means fresh acceptance is unmet."""

import argparse
import json
import subprocess
from pathlib import Path
from xml.etree import ElementTree

from island.forcefields.ambertools.charge_audit import require
from island.forcefields.ambertools.native_precision import RUNTIME_ID, validate_study
from island.forcefields.ambertools.runtime import (
    checksum,
    declared_environment,
    digest,
    validate_runtime,
)
from scripts.validate_phase_4j34_c2 import RUNTIME, STUDY, reconstruct

BASE = '03682038e067a3ce54132adb4beaba84626b8794'
NOT_EXECUTED = 'not_executed_no_verified_native_high_precision_charge_vector'
DOWNSTREAM = ('fresh_gaff2_preparation', 'minimization', 'nvt', 'npt', 'split_uninterrupted', 'relocation', 'fresh_process_continuation')


def validate_outcome(record, study, runtime):
    require(record == reconstruct(study, runtime, Path(study['native_archive'])), 'Native outcome does not reconstruct')
    return record


def junit_counts(path):
    root = ElementTree.parse(path).getroot()
    suites = list(root) if root.tag == 'testsuites' else [root]
    result = {k: sum(int(s.get(k, 0)) for s in suites) for k in ('tests', 'failures', 'errors', 'skipped')}
    result['passed'] = result['tests'] - result['failures'] - result['errors'] - result['skipped']
    return result


def validate(receipt, directory, *, check_runtime=True):
    directory = Path(directory)
    require(receipt.get('identity') == digest({k: v for k, v in receipt.items() if k != 'identity'}), 'Changed receipt identity')
    require(receipt['schema'] == 'island.j34-c2-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong schema/base')
    require(receipt['status'] == 'blocked' and receipt['acceptance_met'] is False, 'Unestablished fresh acceptance')
    require(receipt['production_validated'] is False and receipt['simulation_readiness'] == 'not_established', 'Readiness changed')
    require(receipt['runtime_identity'] == RUNTIME_ID, 'Wrong runtime identity')
    for name, sha in receipt['artifacts'].items():
        path = directory / name
        require(path.resolve().is_relative_to(directory.resolve()) and checksum(path) == sha, 'Changed receipt artifact: ' + name)
    study = validate_study(json.loads(STUDY.read_text()), STUDY)
    require(receipt['study_identity'] == study['identity'], 'Wrong study binding')
    runtime = json.loads(RUNTIME.read_text())
    if check_runtime:
        validate_runtime(runtime, expected_identity=RUNTIME_ID, environment=declared_environment(runtime))
    outcome = json.loads((directory / receipt['outcome']).read_text())
    require(receipt['outcome'] in receipt['artifacts'], 'Unbound native outcome')
    validate_outcome(outcome, study, runtime)
    require(receipt['outcome_identity'] == outcome['identity'], 'Wrong native outcome binding')
    require(receipt['downstream_gates'] == dict.fromkeys(DOWNSTREAM, NOT_EXECUTED), 'False downstream gate claims')
    for path, sha in json.loads((STUDY.parent / 'artifacts/historical-hashes.json').read_text()).items():
        require(checksum(path) == sha, 'Historical evidence changed: ' + path)
    for path, sha in receipt['additional_protected_files'].items():
        require(checksum(path) == sha, 'Protected manifest changed: ' + path)
    commit = receipt['implementation_commit']
    require(len(commit) == 40 and all(c in '0123456789abcdef' for c in commit), 'Missing implementation SHA')
    for path, sha in receipt['implementation_files'].items():
        require(checksum(path) == sha, 'Implementation changed: ' + path)
        # Bind the actual committed implementation, not a documentation commit.
        import hashlib
        content = subprocess.check_output(['git', 'show', commit + ':' + path])
        require(hashlib.sha256(content).hexdigest() == sha, 'Tested implementation commit differs: ' + path)
    require(set(receipt['tests']) == {'focused', 'j22-j34', 'full_supported'}, 'Incomplete test gates')
    for name, gate in receipt['tests'].items():
        require(gate['junit'] in receipt['artifacts'], 'Unbound JUnit gate')
        counts = junit_counts(directory / gate['junit'])
        require(counts == gate['counts'] and counts['failures'] == counts['errors'] == 0 and counts['passed'] > 0, 'Changed or failed test counts')
        command = json.loads((directory / gate['command']).read_text())
        require(gate['command'] in receipt['artifacts'] and command['exit_code'] == 0 and command['implementation_commit'] == commit, 'Wrong test command/commit')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    try:
        receipt = validate(json.loads(args.receipt.read_text()), args.receipt.parent, check_runtime=not args.offline)
    except Exception as error:  # noqa: BLE001 -- preserve validator failure reason
        print('INVALID:', error)
        return 1
    print('VALID BLOCKED DIAGNOSTIC: native high-precision AM1-BCC unavailable; fresh acceptance unmet')
    return 0 if args.diagnostic_only else int(not receipt['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
