"""Software negative controls; no synthetic charge vector is scientific evidence."""

import json
import os
import shutil
import subprocess
import sys
import tarfile
from copy import deepcopy
from pathlib import Path

import pytest

from island.forcefields.ambertools.multimolecule import unwrapped_component
from island.forcefields.ambertools.native_precision import (
    inspect_component,
    require_authoritative_vector,
    validate_diagnostic,
    validate_study,
)
from island.forcefields.ambertools.runtime import checksum, signed
from island.periodic import periodic_chemical_view
from scripts.phase_4j34_real_fixture import build_fixture
from scripts.validate_phase_4j34_c2 import RUNTIME, STUDY, reconstruct


def resign(record):
    return signed({k: v for k, v in record.items() if k != 'identity'})


@pytest.fixture(scope='module')
def retained(tmp_path_factory):
    study = validate_study(json.loads(STUDY.read_text()), STUDY)
    runtime = json.loads(RUNTIME.read_text())
    root = tmp_path_factory.mktemp('native-precision-replay')
    with tarfile.open(study['native_archive']) as tar:
        tar.extractall(root, filter='data')
    fixture = build_fixture('linear', 'GAFF2').fixture
    system = periodic_chemical_view(fixture.packed.system)
    ids = sorted(system.topology.connected_components(), key=min)[1]
    component, _ = unwrapped_component(system, sorted(ids), fixture.packed.system.box.lengths)
    case = next(c for c in study['audits'] if 'component-1/' in c)
    bindings = {'final_graph_identity': fixture.packed.graph.identity, 'runtime': runtime, 'study': study, 'case': case}
    return component, root / case, bindings


def test_actual_native_failure_never_promoted_to_authoritative_vector(retained):
    component, native, bindings = retained
    before = {p.name: checksum(p) for p in native.iterdir() if p.is_file()}
    record = inspect_component(component, native, **bindings)
    assert record['observed_component_totals_e']['sqm_printed_atom_sum_e'] == '-0.003'
    assert record['observed_component_totals_e']['expected_e'] == 0
    assert record['original_import_failure']['stage'] == 'import/component_total_validation'
    assert record['authoritative_charge_vector'] is None
    assert record['charge_tolerance_e'] == 1e-4
    assert record['charge_record_published'] is False
    assert record['npt_executed'] is False
    assert validate_diagnostic(record, component, native, **bindings) == record
    with pytest.raises(ValueError, match='native_charge_precision'):
        require_authoritative_vector(record)
    assert before == {p.name: checksum(p) for p in native.iterdir() if p.is_file()}


@pytest.mark.parametrize('mutation', ['reorder', 'truncate', 'zero_pad', 'rounded', 'stale_input', 'mutated_output', 'unexpected_file'])
def test_intermediate_mutations_rejected(retained, tmp_path, mutation):
    component, native, bindings = retained
    target = tmp_path / 'mutated'
    shutil.copytree(native, target)
    if mutation == 'unexpected_file':
        (target / 'unsupported-high-precision.qout').write_text('0.000000000000\n')
    else:
        filename = 'input.mol2' if mutation == 'stale_input' else 'sqm.out'
        path = target / filename
        text = path.read_text()
        if mutation == 'reorder':
            lines = text.splitlines(keepends=True)
            index = next(i for i, line in enumerate(lines) if 'Mulliken Charge' in line and 'Total' not in line)
            lines[index + 1], lines[index + 2] = lines[index + 2], lines[index + 1]
            changed = ''.join(lines)
        elif mutation == 'truncate':
            changed = text[:len(text) // 2]
        elif mutation == 'zero_pad':
            import re
            changed = re.sub(r'(\s[-+]?\d+\.\d{3})(\s*\n)', r'\g<1>000000\2', text)
        elif mutation == 'rounded':
            import re
            changed = re.sub(r'([-+]?\d+\.\d{2})\d', r'\1', text)
        else:
            changed = text + '\nSOFTWARE NEGATIVE CONTROL\n'
        assert changed != text
        path.write_text(changed)
    with pytest.raises(ValueError, match='Changed native|inventory changed'):
        inspect_component(component, target, **bindings)


@pytest.mark.parametrize('mutation', ['graph', 'case', 'wrong_executable', 'source', 'method', 'study_atom_mapping'])
def test_stale_bindings_rejected_even_when_resigned(retained, mutation):
    component, native, bindings = retained
    changed = deepcopy(bindings)
    if mutation == 'graph':
        changed['final_graph_identity'] = '0' * 64
    elif mutation == 'case':
        changed['case'] = 'native-components/component-9/stale'
    elif mutation in ('wrong_executable', 'source'):
        runtime = changed['runtime']
        if mutation == 'wrong_executable':
            runtime['executables']['sqm']['sha256'] = '0' * 64
        else:
            runtime['sources']['dat/leap/parm/gaff2.dat'] = '0' * 64
        changed['runtime'] = resign(runtime)
    else:
        if mutation == 'method':
            changed['study']['method'] = 'invented_precision_repair'
        else:
            changed['study']['audits'][changed['case']]['atoms'].reverse()
        changed['study'] = resign(changed['study'])
    with pytest.raises(ValueError):
        inspect_component(component, native, **changed)


@pytest.mark.parametrize('field,value', [('authoritative_charge_vector', [0.0] * 30), ('acceptance_met', True), ('status', 'complete'), ('charge_tolerance_e', .01), ('charge_record_published', True), ('npt_executed', True), ('atom_mapping', []), ('observed_component_totals_e', {'expected_e': 0}), ('fresh_preparation_executed', True)])
def test_resigned_diagnostic_claims_rejected(retained, field, value):
    component, native, bindings = retained
    record = inspect_component(component, native, **bindings)
    record[field] = value
    with pytest.raises(ValueError, match='no longer reconstructs'):
        validate_diagnostic(resign(record), component, native, **bindings)


def test_native_probes_do_not_export_more_precision():
    study = validate_study(json.loads(STUDY.read_text()), STUDY)
    assert study['supported_native_high_precision_export_found'] is False
    probes = json.loads((STUDY.parent / 'artifacts/native-probes.json').read_text())
    assert len(probes['results']) == 7
    for probe in probes['results'][:6]:
        assert probe['exit_code'] == 0
        assert probe['charge_tables']
        for table, _total in probe['charge_tables']:
            assert all(len(line.split()[-1].split('.')[-1]) == 3 for line in table.splitlines() if line.strip())
    assert probes['results'][-1]['total_e'] == '-0.003000'


def test_relocation_and_fresh_process_diagnostic_reconstruction(retained, tmp_path):
    component, native, bindings = retained
    target = tmp_path / 'different-directory'
    shutil.copytree(native, target)
    assert inspect_component(component, target, **bindings) == inspect_component(component, native, **bindings)
    expected = reconstruct(bindings['study'], bindings['runtime'], Path(bindings['study']['native_archive']))
    code = '''import json
from scripts.validate_phase_4j34_c2 import reconstruct,STUDY,RUNTIME
from pathlib import Path
s=json.loads(STUDY.read_text());r=json.loads(RUNTIME.read_text())
print(json.dumps(reconstruct(s,r,Path(s['native_archive'])),sort_keys=True))
'''
    child = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == expected


@pytest.mark.parametrize('mode,exit_code', [([], 1), (['--diagnostic-only'], 0)])
def test_real_cli_exit_and_preserved_diagnostic(tmp_path, mode, exit_code):
    output = tmp_path / 'cli-diagnostic'
    child = subprocess.run([sys.executable, '-m', 'scripts.validate_phase_4j34_c2', '--output', str(output), *mode], capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': '.:src'}, check=False)
    assert child.returncode == exit_code, child.stderr
    outcome = json.loads((output / 'outcome.json').read_text())
    assert outcome['status'] == 'blocked', outcome
    assert outcome['acceptance_met'] is False
    assert outcome['authoritative_charge_vector'] is None
    assert outcome['fresh_preparation_executed'] is False


@pytest.mark.parametrize('field,value', [('acceptance_met', True), ('status', 'complete'), ('authoritative_charge_vector', [0.0] * 60), ('fresh_preparation_executed', True), ('npt_executed', True), ('components', {}), ('runtime_identity', '0' * 64)])
def test_resigned_full_outcome_mutations(retained, field, value):
    from scripts.validate_phase_4j34_c2_receipt import validate_outcome

    _, _, bindings = retained
    study, runtime = bindings['study'], bindings['runtime']
    record = reconstruct(study, runtime, Path(study['native_archive']))
    record[field] = value
    with pytest.raises(ValueError, match='does not reconstruct'):
        validate_outcome(resign(record), study, runtime)


def test_official_source_and_patch_provenance_match_native_observations():
    report = json.loads((STUDY.parent / 'artifacts/source-investigation.json').read_text())
    assert report['source_download']['sha256'] == '52fb4fb3370a89b7ce738a2dc3e513c2fc1943fde4b4381846d9e75cc48d840f'
    assert len(report['updates']) == 8
    assert all(not row['relevant_source_changed'] for row in report['updates'])
    import hashlib

    with tarfile.open(STUDY.parent / 'artifacts/native-source-review.tar.gz') as tar:
        for name, sha in report['relevant_sources'].items():
            content = tar.extractfile('AmberTools/src/' + name).read()
            assert hashlib.sha256(content).hexdigest() == sha
        printer = tar.extractfile('AmberTools/src/sqm/qm2_print_charges.F90').read().decode()
        reader = tar.extractfile('AmberTools/src/antechamber/charge.c').read().decode()
    assert 'F14.3' in printer and 'F12.3' in printer
    assert 'atom[number].charge = tmpfloat1' in reader and '%d%s%lf' in reader
    assert '%10.6lf' in reader


def test_modified_study_bytes_rejected(tmp_path):
    study = json.loads(STUDY.read_text())
    study['supported_native_high_precision_export_found'] = True
    changed = tmp_path / 'study.json'
    changed.write_text(json.dumps(resign(study)))
    with pytest.raises(ValueError, match='Unreviewed'):
        validate_study(resign(study), changed)
