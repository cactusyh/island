"""Validate C1 diagnostic integrity; exact-conservation acceptance remains blocked."""

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_precision import audit_in_memory
from scripts.validate_phase_4j36_c1 import BASE, PARENT, audit_case


def verify_signed(record):
    require(record['identity'] == reference_digest({k: v for k, v in record.items() if k != 'identity'}), 'Changed precision evidence identity')


def validate_outcomes(record, audits, runtime_identity, memory_audits):
    verify_signed(record)
    require(record['schema'] == 'island.j36-c1-outcomes.v1' and record['base_commit'] == BASE, 'Wrong precision schema/base')
    require(record['runtime_identity'] == runtime_identity and record['runtime_validation'] in ('live_hashes_and_probe_complete', 'offline_record_binding_only'), 'Wrong runtime binding')
    require(record['acceptance_met'] is False and record['dynamics_executed'] is False and record['native_am1bcc_success'] is False and record['external_oplsaa_is_resp'] is False, 'Unproven workflow or charge origin')
    require(record['production_validated'] is False and record['simulation_readiness'] == 'not_established', 'Changed readiness')
    require(set(record['cases']) == set(audits) == {'linear', 'network'}, 'Missing/extra precision cases')
    require(set(record['native_esp_cases']) == set(memory_audits) == {'linear', 'network'}, 'Missing native in-memory ESP cases')
    for kind, memory in memory_audits.items():
        verify_signed(memory)
        require(memory['charge_conservation_passed'] is False and memory['classification'] == 'authoritative_native_precision_insufficient', 'Unverified native ESP charge success')
        require(record['native_esp_cases'][kind] == {'audit_identity': memory['identity'], 'classification': memory['classification'], 'charge_conservation_passed': False, 'dynamics_executed': False}, 'Native ESP outcome mismatch')
    for kind, audit in audits.items():
        verify_signed(audit)
        require(audit['classification'] == 'authoritative_native_precision_insufficient' and audit['charge_conservation_passed'] is False, 'Native precision failure relabeled')
        require(audit['serialization_loss_detected'] is False and audit['normalization'] is False and audit['redistribution'] is False and audit['charge_tolerance_e'] == '0', 'Changed source precision/charge policy')
        require(record['cases'][kind] == {'audit_identity': audit['identity'], 'classification': audit['classification'], 'stage': audit['stage'],
                                          'charge_conservation_passed': False, 'acceptance_met': False, 'dynamics_executed': False, 'downstream': 'not_executed'}, 'Contradictory outcome/provenance')


def validate(receipt, directory):
    directory = Path(directory).resolve()
    verify_signed(receipt)
    require(receipt['schema'] == 'island.j36-c1-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong corrective receipt/base')
    require(receipt['acceptance_met'] is False and receipt['classification'] == 'authoritative_native_precision_insufficient', 'Unproven RESP acceptance')
    require(receipt['production_validated'] is False and receipt['simulation_readiness'] == 'not_established' and receipt['dynamics_executed'] is False, 'Readiness or dynamics changed')
    require(receipt['charge_origin'] == 'resp_psi4' and receipt['native_am1bcc_success'] is False, 'Charge origin changed')
    for name, expected in receipt['artifacts'].items():
        path = (directory / name).resolve()
        require(path.is_relative_to(directory) and file_hash(path) == expected, 'Changed/unsafe corrective artifact: ' + name)
    require({'result.json', 'source-audit.json', 'historical-hashes.json', 'linear.json', 'network.json'} <= set(receipt['artifacts']), 'Missing precision evidence')
    runtime = json.loads((PARENT / 'environment.json').read_text())
    parent = json.loads((PARENT / 'receipt.json').read_text())
    audits, memory_audits = {}, {}
    with tempfile.TemporaryDirectory(prefix='j36-c1-receipt-') as tmp:
        with tarfile.open(PARENT / 'artifacts/native-quantum.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        with tarfile.open(directory / 'artifacts/native-memory.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        for kind in ('linear', 'network'):
            audits[kind] = audit_case(Path(tmp) / ('quantum-' + kind), kind, runtime, parent)
            require(audits[kind] == json.loads((directory / (kind + '.json')).read_text()), 'Changed precision audit/atom mapping/source values')
            from scripts.phase_4j34_real_fixture import build_fixture
            from scripts.phase_4j36_input import build_input
            memory_audits[kind] = audit_in_memory(Path(tmp) / ('native-esp-' + kind), build_fixture(kind, 'GAFF2').fixture.packed.graph,
                                                native_directory=Path(tmp) / ('quantum-' + kind), in_memory_worker_sha256=file_hash('scripts/phase_4j36_c1_native_esp.py'),
                                                expected_identity=parent['cases'][kind]['execution_identity'], expected_input=build_input(kind),
                                                runtime_identity=runtime['identity'], worker_sha256=file_hash('scripts/phase_4j36_quantum_reference.py'))
            require(memory_audits[kind] == json.loads((directory / ('memory-' + kind + '.json')).read_text()), 'Changed native ESP audit')
    result = json.loads((directory / 'result.json').read_text())
    validate_outcomes(result, audits, runtime['identity'], memory_audits)
    require(receipt['outcomes_identity'] == result['identity'], 'Changed corrective outcome binding')
    source = json.loads((directory / 'source-audit.json').read_text())
    verify_signed(source)
    require(receipt['source_audit_identity'] == source['identity'] and source['runtime_identity'] == runtime['identity'], 'Changed source audit binding')
    require(source['higher_precision_native_export_identified'] is False and source['native_solver_dtype'] == 'float64', 'Unverified higher precision claim')
    for name, expected in source['source_files'].items():
        require(runtime['files'].get(name) == expected, 'Changed pinned native source fingerprint')
    require({'lib/python3.11/site-packages/resp/driver.py', 'lib/python3.11/site-packages/resp/espfit.py'} <= set(source['source_files']), 'Incomplete precision source audit')
    hashes = json.loads((directory / 'historical-hashes.json').read_text())
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json'], text=True).splitlines()
    require(set(hashes) == set(names), 'Incomplete historical inventory')
    for name, expected in hashes.items():
        require(file_hash(name) == expected and hashlib.sha256(subprocess.check_output(['git', 'show', BASE + ':' + name])).hexdigest() == expected, 'Historical evidence changed: ' + name)
    required = {'src/island/forcefields/ambertools/resp_precision.py', 'scripts/validate_phase_4j36_c1.py', 'scripts/validate_phase_4j36_c1_receipt.py', 'tests/test_phase_4j36_c1.py', 'scripts/phase_4j36_c1_native_esp.py', 'scripts/phase_4j36_c1_source_audit.py'}
    require(required <= set(receipt['implementation_files']), 'Missing implementation binding')
    for name, expected in receipt['implementation_files'].items():
        require(file_hash(name) == expected and hashlib.sha256(subprocess.check_output(['git', 'show', receipt['implementation_commit'] + ':' + name])).hexdigest() == expected, 'Tested implementation changed')
    require(set(receipt['tests']) == {'focused', 'j22-j36', 'full_supported'}, 'Missing verification suites')
    from scripts.validate_phase_4j34_c2_receipt import junit_counts
    for gate in receipt['tests'].values():
        require({gate['junit'], gate['command']} <= set(receipt['artifacts']), 'Unbound test evidence')
        counts = junit_counts(directory / gate['junit'])
        command = json.loads((directory / gate['command']).read_text())
        require(counts == gate['counts'] and counts['errors'] == counts['failures'] == 0 and counts['passed'] > 0, 'Failed/changed test evidence')
        require(command['exit_code'] == 0 and command['implementation_commit'] == receipt['implementation_commit'], 'Wrong tested implementation')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()), args.receipt.parent)
    except Exception as error:  # noqa: BLE001 -- preserve nonzero integrity failures
        print('FAIL:', error)
        return 1
    print('VALID BLOCKED DIAGNOSTIC: native float64 RESP residuals persist; no charge repair or RESP MD')
    return 0 if args.diagnostic_only else 1


if __name__ == '__main__':
    raise SystemExit(main())
