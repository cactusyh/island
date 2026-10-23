"""Validate J37 acceptance, independent raw provenance, workflows and evidence."""

import argparse
import hashlib
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from island.analysis.npt import _digest
from island.analysis.npt_continuation import compare_continuations, load_continuation
from island.forcefields.ambertools.provided_preparation import validate_backend_record
from scripts.phase_4j37_reference import digest, read, require, sha
from scripts.validate_phase_4j34_c2_receipt import junit_counts
from scripts.validate_phase_4j35 import validate_case
from scripts.validate_phase_4j37 import BASE, validate_bundle


def validate_structure(receipt, matrix):
    require(receipt['schema'] == 'island.j37-receipt.v1' and receipt['base_commit'] == BASE, 'Wrong schema/base')
    require(receipt['identity'] == digest({k: v for k, v in receipt.items() if k != 'identity'}), 'Changed receipt identity')
    require(matrix['identity'] == receipt['outcomes_identity'] == digest({k: v for k, v in matrix.items() if k != 'identity'}), 'Changed outcomes binding')
    require(matrix['schema'] == 'island.j37-outcomes.v1' and matrix['base_commit'] == BASE, 'Wrong outcomes scope')
    for record in (receipt, matrix):
        require(record['acceptance_met'] is True and record['charge_origin'] == 'RESP/Psi4' and record['native_am1bcc_success'] is False, 'Wrong acceptance/charge origin')
        require(record['production_validated'] is False and record['simulation_readiness'] == 'not_established', 'Changed readiness')
    require(matrix['provenance_relocation'] is matrix['provenance_fresh_process'] is True, 'Missing provenance reconstruction')
    require(set(matrix['cases']) == set(matrix['reviews']) == {'linear', 'network'}, 'Incomplete required matrix')
    for kind, row in matrix['cases'].items():
        validate_case('gaff2/' + kind, row)
        require(row['case'] == 'gaff2/' + kind and row['charge_origin'] == 'resp_psi4' and row['charge_method'] == 'RESP/Psi4' and row['native_am1bcc_success'] is False, 'Mismatched fixture/origin')
        require(row['worker_exit_code'] == 0 and row['periodic_reference_passed'] is True and row['independent_review_identity'] == matrix['reviews'][kind], 'Unmet independent verification')
    require(set(receipt['tests']) == {'focused', 'j22-j36', 'full_supported'}, 'Missing required test suites')
    require(len(receipt['implementation_commit']) == 40 and receipt['implementation_commit'] != BASE, 'Wrong implementation commit')


def validate(receipt, directory):
    directory = Path(directory).resolve()
    for name, expected in receipt['artifacts'].items():
        path = (directory / name).resolve()
        require(path.is_relative_to(directory) and sha(path) == expected, 'Changed J37 artifact: ' + name)
    require({'artifacts/workflows.tar.gz', 'outcomes.json', 'historical-hashes.json'} <= set(receipt['artifacts']), 'Missing evidence')
    matrix = read(directory / 'outcomes.json');validate_structure(receipt, matrix)
    require(receipt['cases'] == matrix['cases'] and receipt['reviews'] == matrix['reviews'], 'Receipt/matrix scientific outcome mismatch')
    require(receipt['charge_vector_preservation']['normalization'] is receipt['charge_vector_preservation']['redistribution'] is False and receipt['charge_vector_preservation']['raw_npy_to_imported_json'] == 'bit_exact' and receipt['charge_vector_preservation']['native_prmtop_serialization_bound_e'] == 1e-14, 'Unproven charge preservation claim')
    required_sources = {'scripts/phase_4j37_anchors.json', 'scripts/phase_4j37_reference.py', 'scripts/phase_4j37_refit.py', 'scripts/phase_4j37_physics.py', 'scripts/validate_phase_4j37.py', 'scripts/validate_phase_4j37_receipt.py', 'tests/test_phase_4j37.py'}
    require(required_sources <= set(receipt['implementation_files']), 'Incomplete tested source inventory')
    required_checks = {'acceptance', 'ruff', 'compileall', 'pip_supported', 'pip_foyer', 'pip_analysis', 'pip_amber', 'pip_quantum', 'diff_check', 'j36_c2_receipt'}
    require(required_checks <= set(receipt['commands']), 'Missing required verification commands')
    historical_nonzero = {'j34_c1_acceptance', 'j34_c2_acceptance', 'j35_c1_acceptance', 'j36_acceptance', 'j36_c1_acceptance'}
    for name, command in receipt['commands'].items():
        artifact = 'artifacts/' + name + '-command.json'
        require(artifact in receipt['artifacts'] and read(directory / artifact) == command, 'Receipt/command evidence mismatch')
        require(command['exit_code'] == (1 if name in historical_nonzero else 0), 'Unmet verification gate: ' + name)
    manifest = read(directory / 'manifest.json')
    require(manifest == {'schema': 'island.j37-manifest.v1', 'receipt_identity': receipt['identity'], 'artifacts': receipt['artifacts'], 'production_validated': False, 'simulation_readiness': 'not_established'}, 'Receipt/manifest mismatch')
    with tempfile.TemporaryDirectory(prefix='j37-receipt-') as tmp:
        with tarfile.open(directory / 'artifacts/workflows.tar.gz') as archive:
            archive.extractall(tmp, filter='data')
        run = Path(tmp) / 'runs'
        require(read(run / 'result.json') == matrix, 'Archived matrix mismatch')
        for kind in ('linear', 'network'):
            review = validate_bundle(run, kind)
            require(review['identity'] == matrix['reviews'][kind], 'Wrong independent review')
            require(validate_bundle(run / 'relocated-provenance', kind) == review == read(run / (kind + '-child-review.json')), 'Missing relocation/child proof')
            for stage in ('refit', 'reconstruction', 'workflow', 'physics'):
                require(read(run / (kind + '-' + stage + '-command.json'))['exit_code'] == 0, 'Failed child stage')
            refit = read(run / (kind + '-refit.json'))
            require(refit['passed'] is True and refit['native_resp_refit_executed'] is True and refit['quantum_scf_reexecuted'] is False and refit['psi4_version'] == '1.9.1' and refit['resp_module_version'] == '0.8', 'Missing native refit proof')
            require(len(refit['components']) == len(review['mapping']), 'Incomplete refit components')
            for index, (component, mapping) in enumerate(zip(refit['components'], review['mapping'], strict=True)):
                require(component['bit_exact'] is True and component['atom_ids'] == mapping['atom_ids'] and component['charges_sha256'] == sha(run / 'raw' / ('native-esp-' + kind) / f'component-{index}/charges-native.npy') and component['inverse_distance_max_error_bohr_inverse'] <= 1e-13, 'Refit/source association changed')
            case = run / kind
            worker = read(case / 'result.json')
            require(worker['identity'] == _digest({k: v for k, v in worker.items() if k != 'identity'}), 'Changed worker diagnostic identity')
            require({k: v for k, v in matrix['cases'][kind].items() if k not in ('worker_exit_code', 'periodic_reference_passed', 'independent_review_identity')} == worker, 'Matrix/worker mismatch')
            comparison = compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated')
            require(comparison['passed'] is True and comparison == read(case / 'comparison.json') and comparison['identity'] == matrix['cases'][kind]['comparison_identity'], 'Continuation mismatch')
            workflow, _ = load_continuation(case / 'uninterrupted')
            backend = read(case / 'uninterrupted/backend.json');validate_backend_record(backend)
            native = backend['parameter_assignment']['native_payload']['record']
            require(read(case / 'resp-charges.json') == json.loads(native['provided_charge_raw']) == read(run / (kind + '-charges.json')), 'Wrong imported charge vector')
            require(native['parmchk2_estimation_lines'] == [], 'Estimated force-field terms')
            for name, expected in native['artifact_sha256'].items():
                path = (case / 'native-components' / name).resolve()
                require(path.is_relative_to((case / 'native-components').resolve()) and sha(path) == expected, 'Changed native assignment asset')
            require(workflow['manifest']['config']['family'] == 'GAFF2', 'Wrong workflow family')
            physics = read(run / (kind + '-physics.json'))
            require(physics['passed'] is True and physics['lammps_kspace_executed'] is False and physics['independent_engine'] is False, 'Unproven reference claim')
            require([s['frame_index'] for s in physics['snapshots']] == [0, len(workflow['frames']) - 1], 'Wrong reference snapshots')
            for snapshot in physics['snapshots']:
                require(snapshot['energy_error_kj_mol'] <= 1e-6 and snapshot['force_error_kj_mol_angstrom'] <= 1e-5 and set(snapshot['component_errors_kj_mol']) == {'bond', 'angle', 'torsion', 'nonbonded'} and max(snapshot['component_errors_kj_mol'].values()) <= 1e-6, 'Reference check failed')
    historical = read(directory / 'historical-hashes.json')
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'docs/evidence', 'docs/paper/reproducibility_manifest.json', '.gitattributes'], text=True).splitlines()
    require(set(historical) == set(names), 'Incomplete historical inventory')
    for name, expected in historical.items():
        require(sha(name) == expected == hashlib.sha256(subprocess.check_output(['git', 'show', BASE + ':' + name])).hexdigest(), 'Historical bytes changed')
    for name, expected in receipt['implementation_files'].items():
        require(sha(name) == expected == hashlib.sha256(subprocess.check_output(['git', 'show', receipt['implementation_commit'] + ':' + name])).hexdigest(), 'Wrong implementation source')
    for test in receipt['tests'].values():
        require({test['junit'], test['command']} <= set(receipt['artifacts']), 'Unbound JUnit/command')
        counts = junit_counts(directory / test['junit']);command = read(directory / test['command'])
        require(counts == test['counts'] and counts['failures'] == counts['errors'] == 0 and counts['passed'] > 0, 'Failed/mutated tests')
        require(command['exit_code'] == 0 and command['implementation_commit'] == receipt['implementation_commit'], 'Wrong tested commit')
    require(_digest(matrix) != matrix['identity'], 'Outcome identity domain confused')
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('receipt', type=Path);args = p.parse_args()
    try:
        validate(read(args.receipt), args.receipt.parent)
    except Exception as error:  # noqa: BLE001 -- fail closed with diagnostic
        print('FAIL:', error);return 1
    print('PASS: independent raw RESP/Psi4 provenance and finite-run GAFF2 reproducibility; readiness not established')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
