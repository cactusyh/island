"""Validate additive J35-C1 diagnostics; readable blocked evidence is not acceptance."""

import argparse
import json
import subprocess
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.resp import file_hash
from scripts.validate_phase_4j35_c1 import BASE


def signed(p):
    require(p.get('identity') == reference_digest({k: v for k, v in p.items() if k != 'identity'}), 'Changed corrective evidence identity')


def validate(p, directory):
    signed(p)
    require(p.get('schema') == 'island.j35-c1-receipt.v1' and p.get('base_commit') == BASE, 'Wrong corrective receipt/base')
    require(p.get('production_validated') is False and p.get('simulation_readiness') == 'not_established', 'Readiness changed')
    require(p.get('acceptance_met') is False and p.get('classification') == 'asset_unavailable', 'RESP acceptance cannot be inferred from external charges')
    directory = Path(directory).resolve()
    for name, sha in p['artifacts'].items():
        path = (directory / name).resolve()
        require(path.is_relative_to(directory) and file_hash(path) == sha, 'Changed corrective artifact: ' + name)
    require({'result.json', 'historical-hashes.json'} <= set(p['artifacts']), 'Missing diagnostic/history binding')
    result = json.loads((directory / 'result.json').read_text())
    signed(result)
    validate_diagnostic(result)
    import hashlib
    inventory = json.loads((directory / 'historical-hashes.json').read_text())
    tree = subprocess.check_output(['git', 'ls-tree', '-r', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json'], text=True)
    entries = {line.split('\t', 1)[1]: line.split('\t', 1)[0].split()[2] for line in tree.splitlines()}
    require(set(inventory) == set(entries), 'Incomplete historical inventory')
    for name, sha in inventory.items():
        data = Path(name).read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        require(file_hash(name) == sha and blob == entries[name], 'Historical evidence changed: ' + name)
    required = {'src/island/forcefields/ambertools/resp.py', 'scripts/validate_phase_4j35_c1.py', 'scripts/validate_phase_4j35_c1_receipt.py', 'tests/test_phase_4j35_c1.py'}
    require(required <= set(p['implementation_files']), 'Incomplete implementation binding')
    require(set(p['tests']) == {'focused', 'j22-j35', 'full_supported'}, 'Missing test suites')
    from scripts.validate_phase_4j34_c2_receipt import junit_counts
    for gate in p['tests'].values():
        require(gate['junit'] in p['artifacts'] and gate['command'] in p['artifacts'], 'Unbound test evidence')
        counts = junit_counts(directory / gate['junit'])
        command = json.loads((directory / gate['command']).read_text())
        require(counts == gate['counts'] and counts['errors'] == counts['failures'] == 0 and counts['passed'] > 0, 'Failed/changed test evidence')
        require(command['exit_code'] == 0 and command['implementation_commit'] == p['implementation_commit'], 'Wrong tested implementation')
    for name, sha in p['implementation_files'].items():
        require(file_hash(name) == sha and hashlib.sha256(subprocess.check_output(['git', 'show', p['implementation_commit'] + ':' + name])).hexdigest() == sha, 'Tested implementation changed')
    return p


def validate_diagnostic(result):
    signed(result)
    require(isinstance(result.get('environment_audit'), list) and len(result['environment_audit']) == 4, 'Missing runtime audit')
    from scripts.validate_phase_4j35_c1 import CODE, PYTHONS
    for row, python in zip(result['environment_audit'], PYTHONS, strict=True):
        require(row['argv'] == [python, '-c', CODE] and row['exit_code'] == 0, 'Incomplete environment audit')
        audit = json.loads(row['stdout'])
        require(audit['psi4_module'] == audit['resp_module'] == 'None' and audit['psi4_executable'] is None, 'Runtime audit contradicts missing assets')
    require(result['schema'] == 'island.j35-c1-diagnostic.v1' and result['base_commit'] == BASE, 'Wrong diagnostic schema/base')
    require(result['acceptance_met'] is False and result['classification'] == 'asset_unavailable', 'Contradictory quantum acceptance')
    require(result['charge_origin'] == 'resp_psi4' and result['charge_method'] == 'RESP' and result['quantum_software'] == 'Psi4', 'Wrong requested provider')
    require(result['production_validated'] is False and result['simulation_readiness'] == 'not_established', 'Diagnostic readiness changed')
    require(set(result['cases']) == {'linear', 'network'}, 'Missing retained case')
    for row in result['cases'].values():
        require(row == {'status': 'asset_unavailable', 'acceptance_met': False, 'charge_generation_executed': False, 'dynamics_executed': False, 'downstream': 'not_executed'}, 'Unexecuted RESP dynamics cannot be claimed')
    require(result['native_am1bcc'] == {'status': 'blocked', 'classification': 'native_charge_generation_unavailable', 'observed_component_total_e': '-0.003', 'acceptance_met': False}, 'Native AM1-BCC diagnostic changed')
    require(set(result['external_controls']) == {'linear', 'network'}, 'Missing external controls')
    for kind, row in result['external_controls'].items():
        path = Path(f'docs/evidence/phase_4j35/references/{kind}.json')
        ref = json.loads(path.read_text())
        require(row == {'charge_origin': 'external_oplsaa', 'charge_method': 'oplsaa_provided', 'reference_identity': ref['identity'], 'reference_sha256': file_hash(path), 'historical_method': ref['method'], 'resp_psi4': False}, 'External OPLS-AA is not RESP/Psi4')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()), args.receipt.parent)
    except Exception as error:  # noqa: BLE001 -- preserve a nonzero diagnostic
        print('FAIL:', error)
        return 1
    print('VALID BLOCKED DIAGNOSTIC; RESP/Psi4 acceptance is not established')
    return 0 if args.diagnostic_only else 1


if __name__ == '__main__':
    raise SystemExit(main())
