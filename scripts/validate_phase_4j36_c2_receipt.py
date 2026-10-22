"""Validate finite precision policy plus real RESP-only GAFF2 continuation."""

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from island.analysis.npt import _digest
from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.provided_preparation import validate_backend_record
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges, create
from island.forcefields.ambertools.resp_precision import audit_in_memory
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.phase_4j36_c2_reference import check
from scripts.phase_4j36_input import build_input
from scripts.validate_phase_4j34_c2_receipt import junit_counts
from scripts.validate_phase_4j35 import validate_case
from scripts.validate_phase_4j36_c1 import audit_case
from scripts.validate_phase_4j36_c2 import BASE


def verify_outcome_identity(matrix, expected):
    require(matrix['identity'] == expected == _digest({k: v for k, v in matrix.items() if k != 'identity'}), 'Wrong outcome binding')


def validate(receipt, directory):
    directory = Path(directory).resolve()
    require(receipt['identity'] == reference_digest({k: v for k, v in receipt.items() if k != 'identity'}), 'Changed C2 receipt identity')
    require(receipt['schema'] == 'island.j36-c2-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong C2 schema/base')
    require(receipt['acceptance_met'] is True and receipt['charge_origin'] == 'RESP/Psi4' and receipt['native_am1bcc_success'] is False, 'Wrong acceptance/origin')
    require(receipt['production_validated'] is False and receipt['simulation_readiness'] == 'not_established', 'Changed readiness')
    for name, expected in receipt['artifacts'].items():
        path = (directory / name).resolve()
        require(path.is_relative_to(directory) and file_hash(path) == expected, 'Changed corrective artifact: ' + name)
    require({'artifacts/workflows.tar.gz', 'historical-hashes.json', 'outcomes.json'} <= set(receipt['artifacts']), 'Missing C2 evidence')
    matrix = json.loads((directory / 'outcomes.json').read_text())
    verify_outcome_identity(matrix, receipt['outcomes_identity'])
    require(matrix['acceptance_met'] is True and set(matrix['cases']) == {'linear', 'network'}, 'Missing required workflow')
    runtime = json.loads(Path('docs/evidence/phase_4j36/environment.json').read_text())
    parent = json.loads(Path('docs/evidence/phase_4j36/receipt.json').read_text())
    with tempfile.TemporaryDirectory(prefix='j36-c2-receipt-') as tmp:
        root = Path(tmp)
        with tarfile.open(directory / 'artifacts/workflows.tar.gz') as archive:
            archive.extractall(root, filter='data')
        run = root / 'runs'
        for label in ('linear', 'network', 'memory-linear', 'memory-network'):
            kind = label.removeprefix('memory-');graph = build_fixture(kind, 'GAFF2').fixture.packed.graph
            audit = json.loads((run / (label + '-audit.json')).read_text())
            require(audit == json.loads(Path('docs/evidence/phase_4j36_c1', label + '.json').read_text()), 'Unreviewed native charge source')
            native = run / 'quantum-evidence' / (('native-esp-' if label.startswith('memory-') else 'quantum-') + kind)
            original = run / 'quantum-evidence' / ('quantum-' + kind)
            if label.startswith('memory-'):
                rebuilt = audit_in_memory(native, graph, native_directory=original, in_memory_worker_sha256=file_hash('scripts/phase_4j36_c1_native_esp.py'),
                    expected_identity=parent['cases'][kind]['execution_identity'], expected_input=build_input(kind), runtime_identity=runtime['identity'], worker_sha256=file_hash('scripts/phase_4j36_quantum_reference.py'))
            else:
                rebuilt = audit_case(original, kind, runtime, parent)
            require(rebuilt == audit, 'Source geometry/intermediate validation failed')
            independent = check(native, audit)
            require(independent == json.loads((run / (label + '-independent.json')).read_text()), 'Independent raw parser mismatch')
            provided = create(audit, graph, runtime, independent)
            loaded = FiniteRESPCharges.load(run / (label + '-charges.json'), graph, expected_identity=matrix['charge_identities'][label])
            require(loaded.payload == provided.payload, 'Changed source/policy/vector')
        for kind in ('linear', 'network'):
            row = matrix['cases'][kind];validate_case('gaff2/' + kind, row)
            require(row['charge_origin'] == 'resp_psi4' and row['charge_method'] == 'RESP/Psi4', 'External charges substituted')
            case = run / kind
            comparison = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
            require(comparison['passed'] and comparison == json.loads((case / 'comparison.json').read_text()) and row['comparison_identity'] == comparison['identity'], 'Continuation comparison mismatch')
            backend = json.loads((case / 'uninterrupted/backend.json').read_text());validate_backend_record(backend)
            native = backend['parameter_assignment']['native_payload']['record']
            provided = json.loads((case / 'resp-charges.json').read_text())
            require(provided == json.loads(native['provided_charge_raw']) and provided['identity'] == matrix['charge_identities']['memory-' + kind] == row['provided_charge_identity'], 'Stale charge/policy preparation binding')
            require(native['parmchk2_estimation_lines'] == [], 'Estimated force field')
            for name, expected in native['artifact_sha256'].items():
                require(file_hash(case / 'native-components' / name) == expected, 'Changed native parameter asset')
            workflow, _ = load_continuation(case / 'uninterrupted')
            require(workflow['manifest']['config']['family'] == 'GAFF2', 'Wrong family')
    hashes = json.loads((directory / 'historical-hashes.json').read_text())
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json'], text=True).splitlines()
    require(set(hashes) == set(names), 'Incomplete historical inventory')
    for name, expected in hashes.items():
        require(file_hash(name) == expected == hashlib.sha256(subprocess.check_output(['git', 'show', BASE + ':' + name])).hexdigest(), 'Historical evidence mutation')
    for name, expected in receipt['implementation_files'].items():
        require(file_hash(name) == expected == hashlib.sha256(subprocess.check_output(['git', 'show', receipt['implementation_commit'] + ':' + name])).hexdigest(), 'Wrong tested implementation')
    require(set(receipt['tests']) == {'focused', 'j22-j36', 'full_supported'}, 'Missing validation gate')
    for gate in receipt['tests'].values():
        require({gate['junit'], gate['command']} <= set(receipt['artifacts']), 'Unbound test evidence')
        counts = junit_counts(directory / gate['junit']);command = json.loads((directory / gate['command']).read_text())
        require(counts == gate['counts'] and counts['failures'] == counts['errors'] == 0 and counts['passed'] > 0, 'Failed/changed test counts')
        require(command['exit_code'] == 0 and command['implementation_commit'] == receipt['implementation_commit'], 'Wrong tested implementation')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__);parser.add_argument('receipt', type=Path);args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()), args.receipt.parent)
    except Exception as error:  # noqa: BLE001 -- retain nonzero integrity diagnostics
        print('FAIL:', error);return 1
    print('PASS: finite-precision RESP/Psi4 GAFF2 acceptance; no AM1-BCC or production-readiness claim')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
