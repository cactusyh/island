"""Strict C3 raw, workflow and receipt validation; diagnostic is not acceptance."""

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from island.analysis.npt import FLAGS
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from scripts.phase_4j38_c3_input import build_input
from scripts.phase_4j38_c3_reference import audit, sha
from scripts.phase_4j38_contract import require, validate_identity, validate_row
from scripts.phase_4j38_fixture import build
from scripts.validate_phase_4j38_c2_receipt import historical_hashes
from scripts.validate_phase_4j38_c2_receipt import validate_bundle as validate_c2
from scripts.validate_phase_4j38_c3 import BASE, KINDS, QUANTUM, matrix
from scripts.validate_phase_4j38_receipt import validate_case_bundle

SOURCES = {str(p) for p in Path('scripts').glob('*4j38_c3*.py')} | {
    'scripts/phase_4j38_c3_anchors.json', 'tests/test_phase_4j38_c3.py',
}


def read(path):
    return json.loads(Path(path).read_text())


def validate_raw(root):
    """Offline trust rooted in newly reviewed C3 raw anchors in tested source."""
    anchors = read('scripts/phase_4j38_c3_anchors.json')
    require(read(root / 'anchors.json') == anchors, 'Changed reviewed raw-source anchors')
    runtime = read(QUANTUM)
    require(read(root / 'runtime.json') == runtime, 'Changed executable environment')
    reviews = {}
    for kind in KINDS:
        require(read(root / ('quantum-' + kind) / 'input.json') == build_input(kind), 'Wrong exact benchmark graph/coordinates/input')
        data = read(root / (kind + '-charges.json'))
        review = audit(root, kind, data, runtime, anchors)
        require(review == read(root / (kind + '-review.json')), 'Changed raw-source independent audit')
        charges = FiniteRESPCharges.load(root / (kind + '-charges.json'), build(kind)[0].packed.graph,
                                         expected_identity=review['charge_identity'])
        require(charges.payload['policy']['accepted'], 'Charge conservation policy failed')
        reviews[kind] = review
    return reviews


def validate_bundle(root):
    root = Path(root)
    validate_c2(root / 'c2-recheck')
    reviews = validate_raw(root / 'raw')
    outcomes = read(root / 'outcomes.json')
    validate_identity(outcomes)
    require(outcomes['schema'] == 'island.j38-c3-outcomes.v1' and set(outcomes['cases']) == set(KINDS), 'Wrong C3 outcome inventory')
    for kind in KINDS:
        row = outcomes['cases'][kind]
        validate_row('gaff2/' + kind, row)
        require(all(row[k] == v for k, v in FLAGS.items()), 'Changed readiness')
        require(row['software_acceptance_met'] is True, 'Missing required corrected-cell acceptance')
        require(row['provided_charge_identity'] == reviews[kind]['charge_identity'] and row['independent_review_identity'] == reviews[kind]['identity'], 'Wrong source-to-workflow binding')
        case = root / ('gaff2-' + kind)
        require(row == read(case / 'result.json'), 'Outcome/result disagreement')
        validate_case_bundle(case, row)
        require(read(case / 'resp-charges.json') == read(root / 'raw' / (kind + '-charges.json')), 'Changed workflow charge vector')
        require(read(case / 'raw-review.json') == reviews[kind], 'Changed workflow provenance')
        backend = read(case / 'uninterrupted/backend.json')
        vector = {r['atom_id']: r['charge_e'] for r in read(case / 'resp-charges.json')['rows']}
        error = max(abs(r['charge'] - vector[r['id']]) for r in backend['numerical_data']['sites'])
        require(error <= 1e-14 and read(case / 'charge-roundtrip.json') == {'authoritative_vector_unchanged': True, 'max_native_prmtop_roundtrip_error_e': error, 'existing_serialization_bound_e': 1e-14}, 'Changed parameter charge roundtrip')
        for name in ('resume',):
            require(read(case / (name + '.command.json'))['exit_code'] == 0, 'Child reconstruction failed')
    expected = matrix(read(root / 'c2-recheck/matrix.json'), outcomes['cases'])
    require(read(root / 'matrix.json') == expected, 'Changed matrix source, graph, stage or acceptance')
    require(expected['unmet_charge_cells'] == ['oplsaa/aromatic'], 'Native aromatic source must remain blocked')
    return expected


def validate_receipt(path):
    path = Path(path); root = path.parent
    receipt, manifest = read(path), read(root / 'manifest.json')
    for value in (receipt, manifest):
        validate_identity(value)
        require(value['base_commit'] == BASE, 'Wrong base')
    require(receipt['schema'] == 'island.j38-c3-receipt.v1' and manifest['schema'] == 'island.j38-c3-manifest.v1', 'Wrong receipt schema')
    require(receipt['manifest_identity'] == manifest['identity'], 'Receipt/manifest mismatch')
    require(receipt['acceptance_met'] is False and receipt['corrected_cells_software_acceptance_met'] is True and receipt['native_am1bcc_success'] is False, 'False acceptance/origin claim')
    require(all(receipt[k] == v for k, v in FLAGS.items()), 'Readiness changed')
    inventory = {str(p.relative_to(root)): sha(p) for p in root.rglob('*') if p.is_file() and p.name not in ('receipt.json', 'manifest.json')}
    require(manifest['files'] == inventory, 'Changed or unbound artifact')
    require(read(root / 'historical-hashes.json') == historical_hashes(BASE), 'Changed historical evidence')
    tested = receipt['tested_implementation_commit']
    require(tested == manifest['tested_implementation_commit'], 'Wrong tested SHA')
    require(set(manifest['implementation_files']) == SOURCES, 'Incomplete implementation inventory')
    for name, expected in manifest['implementation_files'].items():
        blob = subprocess.check_output(['git', 'show', f'{tested}:{name}'])
        require(hashlib.sha256(blob).hexdigest() == expected == sha(name), 'Untested code/source anchors')
    with tempfile.TemporaryDirectory(prefix='j38-c3-receipt-') as tmp:
        with tarfile.open(root / 'artifacts.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        value = validate_bundle(Path(tmp) / 'run')
    require(receipt['matrix_identity'] == value['identity'], 'Wrong matrix binding')
    require(set(receipt['tests']) == {'focused', 'j22-j38-c2', 'full_supported'}, 'Missing test groups')
    for name, counts in receipt['tests'].items():
        suites = ET.parse(root / 'validation' / (name + '.xml')).getroot().findall('testsuite')
        actual = {k: sum(int(s.attrib.get(k, 0)) for s in suites) for k in ('tests', 'failures', 'errors', 'skipped')}
        actual['passed'] = actual['tests'] - actual['failures'] - actual['errors'] - actual['skipped']
        require(counts == actual and actual['passed'] > 0 and actual['failures'] == actual['errors'] == 0, 'Missing/failed test evidence')
        require(read(root / 'validation' / (name + '-command.json'))['implementation_commit'] == tested, 'Tests ran against another commit')
    required = {'ruff', 'compileall', 'pip_supported', 'pip_foyer', 'pip_analysis', 'pip_amber', 'pip_quantum', 'diff_check', 'historical_hashes', 'historical_receipts', 'c2_diagnostic', 'c2_acceptance', 'raw_driver', 'relocated_bundle', 'independent_raw_audit'}
    require(set(receipt['command_exit_codes']) == required, 'Missing required validation command')
    for name, code in receipt['command_exit_codes'].items():
        require(code == (1 if name in ('c2_acceptance', 'raw_driver') else 0), 'Invalid declared command outcome')
        result = read(root / 'validation' / (name + '-command.json'))
        require(result['exit_code'] == code and result['implementation_commit'] == tested, 'Unexpected validation status or untested code: ' + name)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt', type=Path, nargs='?', default=Path('docs/evidence/phase_4j38_c3/receipt.json'))
    p.add_argument('--bundle', type=Path)
    p.add_argument('--diagnostic-only', action='store_true')
    args = p.parse_args()
    try:
        if args.bundle:
            validate_bundle(args.bundle)
        else:
            validate_receipt(args.receipt)
    except (ValueError, KeyError, OSError, TypeError, AssertionError) as error:
        print('INVALID:', error)
        return 2
    print('C3 corrected GAFF2 software gates pass; native aromatic OPLS remains source_incompatible; equilibrium/scientific validation not established')
    return 0 if args.diagnostic_only else 1


if __name__ == '__main__':
    raise SystemExit(main())
