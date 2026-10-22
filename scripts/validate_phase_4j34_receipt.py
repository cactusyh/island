"""Strict J34 acceptance receipt and nested workflow/runtime binding validation."""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

from island.analysis.npt import FLAGS
from island.analysis.npt_continuation import (
    checksum,
    compare_continuations,
    load_continuation,
    require,
)
from island.forcefields.ambertools.multimolecule import validate_backend_record
from island.forcefields.ambertools.runtime import (
    declared_environment,
    digest,
    validate_runtime,
)
from scripts.validate_phase_4j33_c1_receipt import signed
from scripts.validate_phase_4j34 import CASES, SUCCESS, validate_case


def validate(receipt, directory, *, check_runtime=True):
    directory = Path(directory)
    signed(receipt)
    require(receipt.get('schema') == 'island.j34-receipt.v1', 'Unsupported J34 receipt schema')
    require(receipt.get('base_commit') == '3866e5366092816df2238e96ef35a779ed2dcf5b', 'Wrong J34 reviewed base')
    require(len(receipt.get('implementation_commit', '')) == 40, 'Missing tested implementation commit')
    require(all(receipt.get(k) == v and type(receipt.get(k)) is type(v) for k, v in FLAGS.items()), 'Readiness changed')
    for name, expected in receipt['artifacts'].items():
        path = directory / name
        require(path.resolve().is_relative_to(directory.resolve()) and checksum(path) == expected, 'Artifact changed: ' + name)
    require(receipt.get('workflow_archive') in receipt['artifacts'], 'Workflow archive is not bound')
    runtime = json.loads((directory / 'runtime.json').read_text())
    require('runtime.json' in receipt['artifacts'], 'Runtime manifest is not bound')
    require(runtime['identity'] == receipt['runtime_identity'] == digest({k: v for k, v in runtime.items() if k != 'identity'}), 'Runtime manifest changed')
    if check_runtime:
        validate_runtime(runtime, expected_identity=receipt['runtime_identity'], environment=declared_environment(runtime))
    require(set(receipt['cases']) == set(CASES), 'Wrong J34 required inventory')
    with tempfile.TemporaryDirectory(prefix='island-j34-receipt-') as tmp:
        with tarfile.open(directory / receipt['workflow_archive']) as tar:
            require(all((m.isfile() or m.isdir()) and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in tar.getmembers()), 'Unsafe archive member')
            tar.extractall(tmp, filter='data')
        root = Path(tmp)
        require({str(p.relative_to(root)): checksum(p) for p in root.rglob('*') if p.is_file()} == receipt['bundle_files'], 'Bundle inventory changed')
        matrix = signed(json.loads((root / 'matrix.json').read_text()))
        require(set(matrix['cases']) == set(CASES) and matrix['cases'] == receipt['cases'], 'Matrix/receipt mismatch')
        for name, row in matrix['cases'].items():
            case = root / name.replace('/', '-')
            saved = signed(json.loads((case / 'result.json').read_text()))
            require({k: v for k, v in saved.items() if k != 'identity'} == row, 'Matrix/result mismatch')
            validate_case(name, row)
            if name not in SUCCESS:
                continue
            workflow, _ = load_continuation(case / 'uninterrupted')
            comparison = signed(json.loads((case / 'comparison.json').read_text()))
            require(compare_continuations(case / 'uninterrupted', case / 'resumed', case / 'relocated') == comparison and comparison['passed'] is True, 'Continuation comparison failed')
            require(row['comparison_identity'] == comparison['identity'] and row['density_frame_identity'] == comparison['density_frame_identity'], 'Stale comparison/density identity')
            binding = signed(json.loads((case / 'runtime-binding.json').read_text()))
            require(row['runtime_binding_identity'] == binding['identity'], 'Runtime/workflow binding changed')
            backend = json.loads((case / 'uninterrupted/backend.json').read_text())
            for key in ('final_graph_identity', 'assignment_identity', 'source_identity', 'npt_config_identity'):
                expected = workflow['manifest']['config']['identity'] if key == 'npt_config_identity' else comparison[key]
                require(binding[key] == expected, 'Stale workflow binding: ' + key)
            require(binding['backend_identity'] == backend['identity'] and binding['family'] == workflow['manifest']['config']['family'], 'Wrong backend/family binding')
            require(json.loads((case / 'child-command.json').read_text())['exit_code'] == 0, 'Child continuation failed')
            independent = json.loads((case / 'independent-assignment.json').read_text())
            if name.startswith('gaff2/'):
                require(binding['runtime_identity'] == row['runtime_identity'] == runtime['identity'], 'Wrong GAFF2 runtime identity')
                record = backend['parameter_assignment']['native_payload']['record']
                require(record['force_field_data']['sha256'] == runtime['sources']['dat/leap/parm/gaff2.dat'] and record['leaprc']['sha256'] == runtime['sources']['dat/leap/cmd/leaprc.gaff2'], 'Wrong GAFF2 source')
                require(record['tool_package']['version'] == runtime['package']['version'] and record['tool_package']['build'] == runtime['package']['build'], 'Wrong AmberTools package version')
                require(all(record['executable_sha256'][tool] == runtime['executables'][tool]['sha256'] for tool in record['executables']), 'Wrong preparation executable')
                require(independent['parmchk2_estimation_lines'] == [] and independent['max_force_error_kj_mol_angstrom'] <= 1e-5, 'Independent GAFF2 verification unmet')
                if name == 'gaff2/linear':
                    require(record['runtime_identity'] == runtime['identity'], 'Stale component runtime identity')
                    validate_backend_record(backend)
        complete = all(r['acceptance_met'] for r in matrix['cases'].values())
        require(receipt.get('acceptance_met') is complete is matrix['acceptance_met'] is True, 'Incomplete all-family NPT acceptance')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--offline', action='store_true', help='Validate archived identities only; does not attest present runtime availability')
    args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()), args.receipt.parent, check_runtime=not args.offline)
    except Exception as error:  # noqa: BLE001 -- precise nonzero diagnostic
        print('FAIL:', error)
        return 1
    print('PASS: J34 declared all-family NPT matrix and receipt validated')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
