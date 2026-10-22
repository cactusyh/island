"""Validate J34-C1 blocked diagnostic integrity separately from acceptance.

Default exit remains nonzero when fresh AM1-BCC acceptance is unmet.
--diagnostic-only checks the diagnostic record, without an acceptance claim.
"""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

from island.analysis.npt import FLAGS
from island.forcefields.ambertools.charge_audit import require
from island.forcefields.ambertools.runtime import (
    checksum,
    declared_environment,
    digest,
    validate_runtime,
)
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j34_c1 import collect_audits

BASE = '74ff640ccad8d5f0c1dd88e172f061b800489d84'


def signed(record):
    require(record.get('identity') == digest({k: v for k, v in record.items() if k != 'identity'}), 'Changed record identity')
    return record


def validate_outcome(record, output, runtime):
    signed(record)
    require(record['schema'] == 'island.j34-c1-fresh-outcome.v1', 'Wrong fresh outcome schema')
    require(record['case'] == 'gaff2/linear' and record['preparation_mode'] == 'fresh_am1bcc' and record['frozen_input_used'] is False, 'Frozen or wrong fresh case')
    require(record['runtime_identity'] == runtime['identity'], 'Wrong runtime binding')
    require(all(record.get(k) == v and type(record.get(k)) is type(v) for k, v in FLAGS.items()), 'Readiness changed')
    row = signed(json.loads((output / 'result.json').read_text()))
    require({k: v for k, v in row.items() if k != 'identity'} == record['j34_result'], 'Fresh outcome/result mismatch')
    fixture = build_fixture('linear', 'GAFF2').fixture
    require(json.loads((output / 'fixture.json').read_text()) == json.loads(json.dumps(fixture.evidence())), 'Stale final graph or fixture')
    require(row['final_graph_identity'] == fixture.packed.graph.identity and row['runtime_identity'] == runtime['identity'], 'Wrong result graph/runtime identity')
    require(all(row.get(k) == v and type(row.get(k)) is type(v) for k, v in FLAGS.items()), 'Result readiness changed')
    audits = collect_audits(output, runtime)
    require(record['charge_audits'] == audits and audits, 'Native charge audits changed or missing')
    failures = [a for a in audits.values() if a['failure'] is not None]
    require(failures and all(a['classification'] == 'native_sqm_printed_charge_precision_loss' for a in failures), 'Wrong native failure classification')
    require(all(a['published_success_marker'] is False for a in failures), 'Failed charges were published')
    require(record['status'] == 'blocked' and record['acceptance_met'] is False, 'Diagnostic cannot claim acceptance')
    require(row['status'] == 'failed' and row['acceptance_met'] is False and row['npt_executed'] is False, 'Contradictory execution claims')
    require(row['diagnostic']['stage'] == 'preparation' and row['diagnostic']['exception_type'] == 'AmberToolsStageError', 'Wrong failure stage')
    require((output / 'failure.log').is_file() and not any((output / d).exists() for d in ('uninterrupted', 'split', 'resumed')), 'Missing failure diagnostic or dynamics after failed charge generation')
    return record


def validate(receipt, directory, *, check_runtime=True):
    directory = Path(directory)
    signed(receipt)
    require(receipt['schema'] == 'island.j34-c1-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong schema/base')
    require(len(receipt['implementation_commit']) == 40 and all(c in '0123456789abcdef' for c in receipt['implementation_commit']), 'Missing tested implementation commit')
    require(receipt['status'] == 'blocked' and receipt['acceptance_met'] is False, 'Unestablished fresh acceptance')
    require(all(receipt.get(k) == v and type(receipt.get(k)) is type(v) for k, v in FLAGS.items()), 'Readiness changed')
    for name, sha in receipt['artifacts'].items():
        path = directory / name
        require(path.resolve().is_relative_to(directory.resolve()) and checksum(path) == sha, 'Changed artifact: ' + name)
    for key in ('runtime_manifest', 'fresh_archive', 'historical_hash_manifest'):
        require(receipt[key] in receipt['artifacts'], 'Unbound artifact: ' + key)
    runtime = signed(json.loads((directory / receipt['runtime_manifest']).read_text()))
    require(runtime['identity'] == receipt['runtime_identity'], 'Wrong declared runtime')
    if check_runtime:
        validate_runtime(runtime, expected_identity=receipt['runtime_identity'], environment=declared_environment(runtime))
    hashes = json.loads((directory / receipt['historical_hash_manifest']).read_text())
    require(hashes.get('docs/evidence/phase_4j34/receipt.json') == receipt['frozen_j34_receipt_sha256'], 'Frozen control was not separately bound')
    for path, sha in hashes.items():
        require(checksum(path) == sha, 'Historical evidence changed: ' + path)
    with tempfile.TemporaryDirectory(prefix='j34-c1-receipt-') as tmp:
        with tarfile.open(directory / receipt['fresh_archive']) as tar:
            require(all((m.isfile() or m.isdir()) and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in tar.getmembers()), 'Unsafe native archive')
            tar.extractall(tmp, filter='data')
        output = Path(tmp)
        record = signed(json.loads((output / 'fresh-outcome.json').read_text()))
        validate_outcome(record, output, runtime)
        require(receipt['fresh_outcome_identity'] == record['identity'], 'Wrong fresh outcome identity')
    require(receipt['downstream_gates'] and all(v == 'not_executed_fresh_charge_generation_failed' for v in receipt['downstream_gates'].values()), 'False downstream gate claim')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    try:
        receipt = validate(json.loads(args.receipt.read_text()), args.receipt.parent, check_runtime=not args.offline)
    except Exception as error:  # noqa: BLE001 -- preserve exact validator diagnostic
        print('INVALID:', error)
        return 1
    if args.diagnostic_only:
        print('VALID BLOCKED DIAGNOSTIC: fresh GAFF2 acceptance remains unmet')
        return 0
    print('BLOCKED: native fresh AM1-BCC charge conservation; no downstream acceptance')
    return int(not receipt['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
