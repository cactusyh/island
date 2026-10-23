"""Gate fresh GAFF2 NPT on independent reviewed raw RESP/Psi4 provenance."""

import argparse
import json
import shutil
import subprocess
import sys
import tarfile
import traceback
from pathlib import Path

from island.forcefields.ambertools.psi4_runtime import environment, validate_runtime
from island.forcefields.ambertools.runtime import declared_environment
from scripts.phase_4j37_reference import ANCHORS, audit, digest, read, require, sha

BASE = '4e0457e0dafcd2c6e1338f64e8143ebcf7a83236'
FLAGS = {'production_validated': False, 'simulation_readiness': 'not_established'}


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False);stream.write('\n')


def command(argv, output, env=None, cwd=None):
    child = subprocess.run(argv, env=env, cwd=cwd, capture_output=True, text=True, check=False)
    write(str(output) + '-command.json', {'argv': argv, 'cwd': str(cwd) if cwd else None, 'exit_code': child.returncode})
    Path(str(output) + '.log').write_text(child.stdout + child.stderr)
    return child.returncode


def validate_bundle(root, kind):
    """Rebuild the path-independent offline review; self-signing is insufficient."""
    root = Path(root)
    saved = read(root / (kind + '-review.json'))
    charges = read(root / (kind + '-charges.json'))
    rebuilt = audit(root / 'raw', kind, charges, read(root / 'runtime.json'))
    require(saved == rebuilt, 'Independent bundle/provenance changed')
    return rebuilt


def run(output, audit_only=False):
    output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'island.j37-outcomes.v1', 'base_commit': BASE, 'acceptance_met': False,
              'charge_origin': 'RESP/Psi4', 'native_am1bcc_success': False, 'cases': {}, 'reviews': {}, **FLAGS}
    try:
        raw = output / 'raw';raw.mkdir()
        for name, expected in read(ANCHORS)['archives'].items():
            require(sha(name) == expected, 'Changed reviewed raw archive')
            with tarfile.open(name) as archive:
                archive.extractall(raw, filter='data')
        runtime = read('docs/evidence/phase_4j36/environment.json')
        write(output / 'runtime.json', runtime)
        for kind in ('linear', 'network'):
            charges = read(f'docs/evidence/phase_4j36_c2/memory-{kind}-charges.json')
            write(output / (kind + '-charges.json'), charges)
            review = audit(raw, kind, charges, runtime)
            write(output / (kind + '-review.json'), review)
            result['reviews'][kind] = review['identity']
            require(validate_bundle(output, kind) == review, 'Offline reconstruction failed')
        relocated = output / 'relocated-provenance';relocated.mkdir()
        shutil.copytree(raw, relocated / 'raw')
        shutil.copyfile(output / 'runtime.json', relocated / 'runtime.json')
        for kind in ('linear', 'network'):
            for suffix in ('-charges.json', '-review.json'):
                shutil.copyfile(output / (kind + suffix), relocated / (kind + suffix))
            require(validate_bundle(relocated, kind) == read(output / (kind + '-review.json')), 'Relocated raw provenance differs')
            argv = [sys.executable, '-m', 'scripts.phase_4j37_reference', '--root', str(relocated / 'raw'), '--fixture', kind, '--charges', str(relocated / (kind + '-charges.json')), '--runtime', str(relocated / 'runtime.json'), '--output', str(output / (kind + '-child-review.json'))]
            require(command(argv, output / (kind + '-reconstruction')) == 0 and read(output / (kind + '-child-review.json')) == read(output / (kind + '-review.json')), 'Fresh-process provenance reconstruction failed')
        result['provenance_relocation'] = result['provenance_fresh_process'] = True
        if audit_only:
            result['scope'] = 'provenance_inventory_only_no_workflow_acceptance'
        else:
            # Independent actual external fit replay supplements hashes and log
            # checks; no importer is used in the external process.
            validate_runtime(runtime, expected_identity=runtime['identity'])
            for kind in ('linear', 'network'):
                argv = [str(Path(runtime['prefix']) / 'bin/python'), str(Path('scripts/phase_4j37_refit.py').resolve()), '--directory', str((raw / ('native-esp-' + kind)).resolve()), '--output', str((output / (kind + '-refit.json')).resolve())]
                require(command(argv, output / (kind + '-refit'), environment(runtime['prefix']), cwd=output.resolve()) == 0, 'Independent native RESP refit failed')
                require(read(output / (kind + '-refit.json'))['passed'], 'Missing raw fitting evidence')
            # Both independent gates complete before either molecular workflow.
            amber = read('docs/evidence/phase_4j34/runtime.json')
            for kind in ('linear', 'network'):
                review = validate_bundle(output, kind)
                argv = [sys.executable, '-m', 'scripts.validate_phase_4j36_c2', '--output', str(output / kind), '--worker', kind,
                        '--charges', str(output / (kind + '-charges.json')), '--expected', review['charge_identity']]
                code = command(argv, output / (kind + '-workflow'), declared_environment(amber))
                path = output / kind / 'result.json'
                row = read(path) if path.is_file() else {'acceptance_met': False, 'diagnostic': 'No worker result'}
                if code == 0 and row.get('acceptance_met'):
                    physics = [sys.executable, '-m', 'scripts.phase_4j37_physics', '--case', str(output / kind), '--output', str(output / (kind + '-physics.json'))]
                    physics_code = command(physics, output / (kind + '-physics'))
                    row['periodic_reference_passed'] = physics_code == 0 and read(output / (kind + '-physics.json'))['passed']
                    row['acceptance_met'] = row['acceptance_met'] and row['periodic_reference_passed']
                row['independent_review_identity'] = review['identity']
                row['worker_exit_code'] = code
                result['cases'][kind] = row
            result['acceptance_met'] = set(result['cases']) == {'linear', 'network'} and all(r['acceptance_met'] and r['worker_exit_code'] == 0 for r in result['cases'].values())
    except Exception as error:  # noqa: BLE001 -- preserve partial outcomes, fail closed
        (output / 'failure.log').write_text(traceback.format_exc())
        result['diagnostic'] = str(error)
    result['identity'] = digest(result)
    write(output / 'result.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--audit-only', action='store_true')
    args = p.parse_args()
    return int(not run(args.output, args.audit_only)['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
