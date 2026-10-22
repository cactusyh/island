"""Read-only native precision investigation. Completion is not acceptance.

Default exit is nonzero: no native high-precision charge vector was found.
--diagnostic-only validates the investigation without claiming fresh generation.
No frozen preparation, charge vector or success override is accepted.
"""

import argparse
import json
import tarfile
import tempfile
import traceback
from pathlib import Path

from island.forcefields.ambertools.charge_audit import require
from island.forcefields.ambertools.multimolecule import unwrapped_component
from island.forcefields.ambertools.native_precision import (
    inspect_component,
    validate_study,
)
from island.forcefields.ambertools.runtime import (
    checksum,
    declared_environment,
    signed,
    validate_runtime,
)
from island.periodic import periodic_chemical_view
from scripts.phase_4j34_real_fixture import build_fixture

STUDY = Path('docs/evidence/phase_4j34_c2/study.json')
RUNTIME = Path('docs/evidence/phase_4j34/runtime.json')


def reconstruct(study, runtime, archive):
    require(checksum(archive) == study['native_archive_sha256'], 'Changed retained native archive')
    fixture = build_fixture('linear', 'GAFF2').fixture
    system = periodic_chemical_view(fixture.packed.system)
    components = sorted(system.topology.connected_components(), key=min)
    result = {}
    with tempfile.TemporaryDirectory(prefix='j34-c2-native-replay-') as tmp:
        with tarfile.open(archive) as tar:
            require(all((m.isfile() or m.isdir()) and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in tar.getmembers()), 'Unsafe native archive')
            tar.extractall(tmp, filter='data')
        for case in sorted(study['audits']):
            index = int(Path(case).parts[1].split('-')[1])
            component, _ = unwrapped_component(system, sorted(components[index]), fixture.packed.system.box.lengths)
            result[case] = inspect_component(component, Path(tmp) / case,
                                            final_graph_identity=fixture.packed.graph.identity,
                                            runtime=runtime, study=study, case=case)
    return signed({'schema': 'island.j34-c2-native-outcome.v1', 'status': 'blocked',
                   'acceptance_met': False, 'study_identity': study['identity'],
                   'runtime_identity': runtime['identity'], 'components': result,
                   'authoritative_charge_vector': None,
                   'evidence_mode': study['evidence_mode'],
                   'fresh_preparation_executed': False, 'npt_executed': False,
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--diagnostic-only', action='store_true')
    args = parser.parse_args()
    # New, empty directory only; never overwrite retained evidence.
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        study = validate_study(json.loads(STUDY.read_text()), STUDY)
        runtime = json.loads(RUNTIME.read_text())
        validate_runtime(runtime, expected_identity=study['runtime_identity'], environment=declared_environment(runtime))
        record = reconstruct(study, runtime, Path(study['native_archive']))
    except Exception as error:  # noqa: BLE001 -- preserve failed diagnostic, no success fallback
        (args.output / 'failure.log').write_text(traceback.format_exc())
        record = signed({'schema': 'island.j34-c2-native-outcome.v1', 'status': 'failed',
                         'acceptance_met': False, 'diagnostic': str(error),
                         'fresh_preparation_executed': False, 'npt_executed': False,
                         'production_validated': False, 'simulation_readiness': 'not_established'})
    (args.output / 'outcome.json').write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(record['status'] + ': no authoritative high-precision charge vector; fresh GAFF2 acceptance unmet')
    return int(not (args.diagnostic_only and record['status'] == 'blocked'))


if __name__ == '__main__':
    raise SystemExit(main())
