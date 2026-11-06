"""Archive actual J91 fenced checks; synthetic controls are not scientific data."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j91_intake import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo = Path.cwd()
root = repo / ROOT
v = root / 'validation'
v.mkdir(exist_ok=False)
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
history = historical(repo)
write(v / 'tested-implementation.json', {'commit': commit, 'files': {n: pin(repo/n) for n in IMPLEMENTATION}})
write(v / 'repository-state.json', {
    'base_commit': BASE, 'implementation_commit': commit, 'branch': BRANCH,
    'verified_remote_base_branch': 'codex/opls-provider-registry-provisional-polymers',
    'verified_remote_base_head': BASE,
    'verified_remote_main_before_work': '63ec95951e5f54f8a8ecdc43c5777b9b4bad3222',
    'local_main': subprocess.check_output(['git','rev-parse','main'],text=True).strip(),
    'main_modified': False, 'historical_files_preserved': len(history),
    'source_bundle_search': 'None. ISLAND_OPLS_SOURCE_DIR and all three control pins were absent. No fallback assets inspected.',
    'unrelated_untracked_files': 'preserved',
})
shutil.copyfile('/tmp/j91-terminal-collection.json', v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j91-terminal.xml').read_bytes()))
active = Path(os.environ['ISLAND_C1_POLICY']).parent
runs = []
for directory in sorted(Path('/tmp').glob('island_j91_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory == active:
        continue
    assert (directory/'result.json').is_file(), directory
    dest = v/'fences'/directory.name
    shutil.copytree(directory, dest)
    events = dest/'events.jsonl'
    (dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes()))
    events.unlink()
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json', runs)
logs = v/'logs'
logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j91_*.log')):
    if p.stem != active.name:
        (logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
names = ['terminal_tests', 'ruff', 'compileall', 'pip_supported', 'diff', 'offline_replay', 'relocated_replay', 'C1_receipt']
names += ['J'+n+'_receipt' for n in ('76','78','79','80','81','82','83','84','85','86','87','88','89','90')]
names += ['J'+n+'_replay' for n in ('81','82','83','84','85','86','87','88','89','90')]
checks = []
for name in names:
    attempt = 'island_j91_final_' + ('terminal' if name == 'terminal_tests' else name)
    result = read(Path('/tmp')/attempt/'result.json')
    assert result['returncode'] == 0, attempt
    checks.append({'name': name, 'returncode': 0, 'log': 'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json', checks)
reasons = {'island_j91_focused1': 'Unsupported tacticity raised IslandError outside the structured blocker path. Fixed by catching package topology exceptions; no scientific execution.'}
failures = []
for n in runs:
    r = read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons, n
        failures.append({'run': n, 'returncode': r['returncode'], 'reason': reasons[Path(n).name], 'scientific_execution': False})
write(v/'development-outcomes.json', {
    'failed_fenced_software_attempts': failures,
    'preliminary_unfenced_lint': {'status': 'failed_then_corrected', 'diagnostic': 'Ruff C401 set comprehension; observed in development tool output before implementation freeze.', 'scientific_execution': False},
    'all_pre_archive_fenced_attempts_with_logs': True,
    'maintenance_fences': 'Archive, seal and post-seal checks use separate fences. Their logs stay in /tmp to avoid a circular seal.',
})
cases = list(ET.parse('/tmp/j91-terminal.xml').iter('testcase'))
groups = {}
skips = []
for case in cases:
    group = 'J91_focused' if case.attrib['classname'] == 'tests.test_phase_4j91' else 'regressions'
    counts = groups.setdefault(group, {k: 0 for k in ('passed','skipped','failure','error')})
    outcome = next((k for k in ('failure','error','skipped') if case.find(k) is not None), 'passed')
    counts[outcome] += 1
    if outcome == 'skipped':
        skips.append({'test': case.attrib['classname']+'::'+case.attrib['name'], 'reason': case.find('skipped').attrib.get('message','')})
required = ['test_complete_source_contract_software_only', 'test_missing_source_does_not_construct_fake_graph_or_inspect', 'test_missing_interaction_no_zero_or_derived_exclusion', 'test_fragment_exports_roundtrip_no_engine_execution', 'test_lammps_conversion_requires_independent_audit', 'test_deterministic_relocated_external_replay_and_missing_originals', 'test_strict_receipt_reconstructs_semantics_even_after_resigning']
for name in required:
    found = [c for c in cases if c.attrib['classname'] == 'tests.test_phase_4j91' and c.attrib['name'].split('[')[0] == name]
    assert found and all(all(c.find(k) is None for k in ('failure','error','skipped')) for c in found), name
write(v/'opls-control-outcomes.json', {
    'tests': required, 'passed': True, 'source_kind': 'manufactured_software_controls_only',
    'external_manifest_census_and_pins': True, 'explicit_exclusion_coverage_required': True,
    'complete_audit_required_before_parameter_inspection': True, 'scientific_validation': False,
    'three_engine_serialization_round_trips': True, 'LAMMPS_conversion_requires_independent_audit': True,
    'OpenMM_artifact': 'data-only XML, loadableForceField=false; never loaded in an engine',
    'authorized_real_sources_supplied': 0, 'scientific_graphs_qualified': 0,
})
for name in ('archive','relocate','checks'):
    shutil.copyfile('/tmp/island_j91_'+name+'.py', v/(name+'.py'))
real_case = next(c for c in cases if c.attrib['classname'] == 'tests.test_phase_4j91' and c.attrib['name'] == 'test_real_external_intake_only_when_explicitly_supplied')
write(v/'real-source-integration-status.json', {
    'status': 'not_run_missing_user_inputs' if real_case.find('skipped') is not None else 'executed',
    'scientific_qualification_claimed': False,
    'reason': real_case.find('skipped').attrib.get('message','') if real_case.find('skipped') is not None else None,
    'synthetic_controls_are_not_real_sources': True,
})
write(v/'scope-skips.json', {
    'explicit_skips': skips, 'scientific_checks_passed': 0,
    'scientific_execution': {k: 0 for k in ('QM','MD','converter','physical_fitting','parameter_generation','energy_force','minimization','OpenMM_Context','reserved_validation','production')},
    'full_scientific_suite': 'Not run; outside authorized scope. Inherited full portable regression inventory executed.',
})
write(v/'summary.json', {
    'tested_implementation_commit': commit, 'groups': groups, 'unique_terminal_tests': len(cases),
    'historical_files_preserved': len(history), 'archived_fenced_attempts': len(runs),
    'failed_fenced_software_attempts': len(failures), 'preliminary_lint_failures': 1,
    'required_checks': len(checks), 'J81_J90_unchanged': True,
    'authorized_real_sources_supplied': 0, 'scientific_graphs_qualified': 0,
    'runtime_registration': False, 'accepted_coefficients': 0,
    'scientific_validation': 'not_established',
})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
