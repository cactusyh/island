"""Validate additive J33-C1 evidence, including every nested continuation bundle.

A valid diagnostic receipt may honestly record unavailable GAFF2 assets. It is
not an all-family acceptance success; --require-all makes that distinction an
exit-status gate.
"""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

from island.analysis.npt import FLAGS, _body, _digest
from island.analysis.npt_continuation import checksum, compare_continuations, require
from scripts.validate_phase_4j33_c1 import CASES, validate_case


def signed(record):
    require(record.get('identity') == _digest(_body(record)), 'Evidence identity mismatch')
    return record


def validate(receipt, directory):
    directory = Path(directory)
    signed(receipt)
    require(receipt.get('schema') == 'island.j33-c1-receipt.v1', 'Wrong corrective receipt schema')
    require(receipt.get('base_commit') == '6f34bf370746ff01fc0491490c0b3dda4a696ec0', 'Wrong reviewed base')
    for key in ('implementation_commit',):
        value = receipt.get(key, '')
        require(len(value) == 40 and all(c in '0123456789abcdef' for c in value), 'Missing exact tested implementation commit')
    require(all(receipt.get(k) == v and type(receipt.get(k)) is type(v) for k, v in FLAGS.items()), 'Readiness constraints changed')
    require(receipt.get('historical_evidence_preserved') is True, 'Historical preservation missing')
    for name, digest in receipt['artifacts'].items():
        path = directory / name
        require(path.resolve().is_relative_to(directory.resolve()), 'Unsafe evidence path')
        require(checksum(path) == digest, 'Evidence artifact changed: ' + name)
    archive = directory / receipt['continuation_archive']
    require(receipt['continuation_archive'] in receipt['artifacts'], 'Unbound continuation archive')
    with tempfile.TemporaryDirectory(prefix='j33-c1-readonly-') as tmp:
        with tarfile.open(archive) as tar:
            members = tar.getmembers()
            require(all(m.isfile() or m.isdir() for m in members), 'Unsupported archive entry')
            require(all(not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in members), 'Unsafe archive path')
            tar.extractall(tmp, filter='data')
        root = Path(tmp)
        inventory = {str(p.relative_to(root)): checksum(p) for p in root.rglob('*') if p.is_file()}
        require(inventory == receipt['bundle_files'], 'Bundle content inventory mismatch')
        matrix = signed(json.loads((root / 'matrix.json').read_text()))
        require(set(matrix['cases']) == set(CASES) == set(receipt['cases']), 'Wrong acceptance case inventory')
        for name in CASES:
            row = signed(json.loads((root / name.replace('/', '-') / 'result.json').read_text()))
            require(_body(row) == matrix['cases'][name] == receipt['cases'][name], 'Matrix/result/receipt disagreement')
            validate_case(name, row)
            if name.startswith('oplsaa/'):
                case = root / name.replace('/', '-')
                comparison = signed(json.loads((case / 'comparison.json').read_text()))
                require(compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated') == comparison, 'Continuation comparison mismatch')
                require(comparison['passed'] is True and comparison['identity'] == row['comparison_identity'], 'Continuation did not pass')
                require(json.loads((case / 'child-command.json').read_text())['exit_code'] == 0, 'Fresh process did not succeed')
        expected = all(row['acceptance_met'] for row in matrix['cases'].values())
        require(receipt['acceptance_met'] is matrix['acceptance_met'] is expected, 'False all-family acceptance claim')
        require(matrix['opls_pcff_contract_met'] is True and receipt['opls_pcff_contract_met'] is True, 'Mandatory continuation cases unmet')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--require-all', action='store_true')
    args = parser.parse_args()
    try:
        record = validate(json.loads(args.receipt.read_text()), args.receipt.parent)
        require(not args.require_all or record['acceptance_met'], 'All-family acceptance unmet: see asset-unavailable diagnostics')
    except Exception as error:  # noqa: BLE001 -- diagnostic CLI
        print('FAIL:', error)
        return 1
    print('PASS: corrective evidence validated; all-family acceptance=' + str(record['acceptance_met']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
