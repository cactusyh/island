"""Real retained native failure plus software negative controls, never MD evidence."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.exceptions import AmberToolsInputError, AmberToolsStageError
from island.forcefields.ambertools.charge_audit import audit_component, validate_audit
from island.forcefields.ambertools.multimolecule import unwrapped_component
from island.forcefields.ambertools.runtime import checksum, signed
from island.periodic import periodic_chemical_view
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j34_c1 import collect_audits
from scripts.validate_phase_4j34_c1_receipt import validate_outcome

RUNTIME = Path('docs/evidence/phase_4j34/runtime.json')


@pytest.fixture(scope='module')
def retained(tmp_path_factory):
    root = tmp_path_factory.mktemp('j34-c1-native-control')
    archive = Path('docs/evidence/phase_4j34/artifacts/initial-validation.tar.gz')
    with tarfile.open(archive) as tar:
        tar.extractall(root, filter='data')
    fixture = build_fixture('linear', 'GAFF2').fixture
    system = periodic_chemical_view(fixture.packed.system)
    ids = sorted(system.topology.connected_components(), key=min)[1]
    component, _ = unwrapped_component(system, sorted(ids), fixture.packed.system.box.lengths)
    native = next((root / 'native-from-scratch/gaff2-linear/native-components/component-1').glob('island-*'))
    runtime = json.loads(RUNTIME.read_text())
    record = audit_component(component, native, runtime_identity=runtime['identity'])
    return component, native, runtime, record


def test_exact_native_failure_without_charge_repair(retained):
    component, native, runtime, record = retained
    assert record['classification'] == 'native_sqm_printed_charge_precision_loss'
    assert record['failure']['stage'] == 'import/component_total_validation'
    assert record['totals']['expected_e'] == 0
    assert record['totals']['sqm_reported_total_e'] == '0.000'
    assert record['totals']['sqm_printed_atom_sum_e'] == '-0.003'
    assert record['totals']['mol2_e'] == '-0.003000'
    assert abs(record['totals']['prmtop_e'] + .003) < 3e-9
    assert record['tolerance_e'] == 1e-4
    assert record['published_success_marker'] is False
    assert record['native_pre_bcc_changed_atom_count'] > 0
    assert validate_audit(record, component, native, runtime_identity=runtime['identity']) == record
    assert all(checksum(native / name) == sha for name, sha in record['files'].items())


@pytest.mark.parametrize('mutation', ['input', 'sqm_formal_charge', 'sqm_input_name', 'sqm_input_coordinate', 'sqm_element', 'sqm_order', 'sqm_missing', 'ac_name', 'mol2_name', 'missing_prmtop'])
def test_malformed_charge_path_never_classified_as_supported_limitation(retained, tmp_path, mutation):
    component, native, runtime, _ = retained
    target = tmp_path / 'native'
    shutil.copytree(native, target)
    if mutation == 'missing_prmtop':
        (target / 'result.prmtop').unlink()
    else:
        filename, old, new = {
            'input': ('input.mol2', 'C000', 'X000'),
            'sqm_formal_charge': ('sqm.in', 'qmcharge=0', 'qmcharge=1'),
            'sqm_input_name': ('sqm.in', 'C000', 'X000'),
            'sqm_input_coordinate': ('sqm.in', '3.6260', '9.6260'),
            'sqm_element': ('sqm.out', '1      C', '1      O'),
            'sqm_order': ('sqm.out', '1      C', '9      C'),
            'sqm_missing': ('sqm.out', 'Mulliken Charge', 'unavailable charge'),
            'ac_name': ('ANTECHAMBER_AM1BCC_PRE.AC', 'C000', 'X000'),
            'mol2_name': ('typed.mol2', 'C000', 'X000'),
        }[mutation]
        path = target / filename
        text = path.read_text()
        assert old in text
        path.write_text(text.replace(old, new))
    with pytest.raises((ValueError, AmberToolsInputError, AmberToolsStageError, OSError)):
        audit_component(component, target, runtime_identity=runtime['identity'])


@pytest.mark.parametrize('field', ['totals', 'runtime_identity', 'failure', 'tolerance_e', 'files', 'published_success_marker'])
def test_resigned_diagnostic_mutation_rejected(retained, field):
    component, native, runtime, record = retained
    changed = deepcopy(record)
    changed[field] = {'software_negative_control': True}
    changed = signed({k: v for k, v in changed.items() if k != 'identity'})
    with pytest.raises(ValueError, match='no longer reconstructs'):
        validate_audit(changed, component, native, runtime_identity=runtime['identity'])


def test_relocation_and_fresh_process_audit_reconstruction(retained, tmp_path):
    component, native, runtime, record = retained
    output = tmp_path / 'different-location'
    destination = output / 'native-components/component-1' / native.name
    shutil.copytree(native, destination)
    assert audit_component(component, destination, runtime_identity=runtime['identity']) == record
    code = '''import json,sys
from pathlib import Path
from scripts.validate_phase_4j34_c1 import collect_audits
r=collect_audits(Path(sys.argv[1]),json.loads(Path(sys.argv[2]).read_text()))
print(json.dumps(r,sort_keys=True))
'''
    child = subprocess.run([sys.executable, '-c', code, str(output), str(RUNTIME)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src'})
    assert child.returncode == 0, child.stderr
    assert list(json.loads(child.stdout).values()) == [record]


def test_missing_runtime_cli_nonzero_with_diagnostic(tmp_path):
    runtime = json.loads(RUNTIME.read_text())
    runtime['amberhome'] = str(tmp_path / 'absent-ambertools')
    runtime = signed({k: v for k, v in runtime.items() if k != 'identity'})
    manifest = tmp_path / 'runtime.json'
    manifest.write_text(json.dumps(runtime))
    output = tmp_path / 'failed'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j34_c1', '--runtime', str(manifest), '--runtime-identity', runtime['identity'], '--output', str(output)], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src'})
    assert child.returncode != 0
    record = json.loads((output / 'fresh-outcome.json').read_text())
    assert record['acceptance_met'] is False
    assert record['j34_result']['status'] == 'asset_unavailable'
    assert record['j34_result']['npt_executed'] is False
    assert (output / 'failure.log').is_file()


def test_frozen_preparation_is_not_a_fresh_cli_option(tmp_path):
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j34_c1', '--runtime', str(RUNTIME), '--runtime-identity', '0' * 64, '--output', str(tmp_path / 'run'), '--linear-preparation', 'historical.tar.gz'], capture_output=True, text=True, check=False, env={**os.environ, 'PYTHONPATH': '.:src'})
    assert child.returncode != 0
    assert 'unrecognized arguments' in child.stderr


@pytest.fixture(scope='module')
def diagnostic_control(retained):
    # Historical failure replay is a SOFTWARE schema control, never a new
    # fresh-generation success or evidence that dynamics executed.
    _, native, runtime, _ = retained
    output = native.parents[2]
    row = json.loads((output / 'result.json').read_text())
    record = signed({'schema': 'island.j34-c1-fresh-outcome.v1',
                     'case': 'gaff2/linear', 'preparation_mode': 'fresh_am1bcc',
                     'frozen_input_used': False, 'runtime_identity': runtime['identity'],
                     'status': 'blocked', 'acceptance_met': False,
                     'production_validated': False, 'simulation_readiness': 'not_established',
                     'j34_result': {k: v for k, v in row.items() if k != 'identity'},
                     'charge_audits': collect_audits(output, runtime)})
    return output, runtime, record


def test_blocked_diagnostic_readable_but_not_acceptance(diagnostic_control):
    output, runtime, record = diagnostic_control
    assert validate_outcome(record, output, runtime)['acceptance_met'] is False


@pytest.mark.parametrize('field,value', [('acceptance_met', True), ('status', 'complete'), ('frozen_input_used', True), ('runtime_identity', '0' * 64), ('charge_audits', {}), ('simulation_readiness', 'established')])
def test_resigned_outcome_mutations(diagnostic_control, field, value):
    output, runtime, record = diagnostic_control
    changed = deepcopy(record)
    changed[field] = value
    changed = signed({k: v for k, v in changed.items() if k != 'identity'})
    with pytest.raises(ValueError):
        validate_outcome(changed, output, runtime)
