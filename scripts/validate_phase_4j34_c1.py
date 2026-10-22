"""Fresh-only J34-C1 acceptance; frozen preparation is never an input option."""

import argparse
import json
import subprocess
import sys
import traceback
from pathlib import Path

from island.analysis.npt import FLAGS, _signed, _write_json
from island.forcefields.ambertools.charge_audit import audit_component
from island.forcefields.ambertools.multimolecule import unwrapped_component
from island.forcefields.ambertools.runtime import declared_environment
from island.periodic import periodic_chemical_view
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j34 import validate_case, worker


def collect_audits(output, runtime):
    """Associate every completed native output with its original final component."""
    fixture = build_fixture('linear', 'GAFF2').fixture
    system = periodic_chemical_view(fixture.packed.system)
    components = sorted(system.topology.connected_components(), key=min)
    audits = {}
    for index, ids in enumerate(components):
        component, _ = unwrapped_component(system, sorted(ids), fixture.packed.system.box.lengths)
        root = Path(output) / 'native-components' / f'component-{index}'
        for native in sorted(root.glob('island-*')):
            if (native / 'result.prmtop').is_file():
                audits[str(native.relative_to(output))] = audit_component(component, native, runtime_identity=runtime['identity'])
    return audits


def run(output, runtime, environment):
    # The existing driver already prevents dynamics after failed preparation.
    # No retained/frozen preparation argument is supplied on this path.
    row = worker(output, 'gaff2/linear', runtime, runtime['identity'], environment)
    record = {'schema': 'island.j34-c1-fresh-outcome.v1', 'preparation_mode': 'fresh_am1bcc',
              'frozen_input_used': False, 'runtime_identity': runtime['identity'],
              'case': 'gaff2/linear', 'j34_result': row, 'acceptance_met': False, **FLAGS}
    try:
        record['charge_audits'] = collect_audits(output, runtime)
        if row['acceptance_met']:
            validate_case('gaff2/linear', row)
            if len(record['charge_audits']) != 2 or any(a['failure'] is not None for a in record['charge_audits'].values()):
                raise ValueError('Fresh charge-path audit is incomplete or inconsistent')
            record['acceptance_met'] = True
        record['status'] = 'complete' if record['acceptance_met'] else 'blocked'
    except Exception as error:  # noqa: BLE001 -- preserve diagnostic failures, never accept them
        (Path(output) / 'charge-audit-failure.log').write_text(traceback.format_exc())
        record.update(status='failed', diagnostic={'exception': type(error).__name__, 'message': str(error)})
    _write_json(Path(output) / 'fresh-outcome.json', _signed(record))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--runtime-identity', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--activated-child', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    runtime = json.loads(args.runtime.read_text())
    if runtime.get('identity') != args.runtime_identity:
        parser.error('Declared runtime identity differs')
    environment = declared_environment(runtime)
    if not args.activated_child:
        return subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j34_c1', *sys.argv[1:], '--activated-child'], env=environment, check=False).returncode
    return int(not run(args.output, runtime, environment)['acceptance_met'])


if __name__ == '__main__':
    raise SystemExit(main())
