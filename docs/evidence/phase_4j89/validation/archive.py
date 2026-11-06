"""Archive actual fenced J89 software verification without scientific execution."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j89_opls import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo=Path.cwd(); root=repo/ROOT; v=root/'validation'; v.mkdir(exist_ok=False)
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
history=historical(repo)
write(v/'tested-implementation.json', {'commit':commit,'files':{n:pin(repo/n) for n in IMPLEMENTATION}})
write(v/'repository-state.json', {'base_commit':BASE,'implementation_commit':commit,'branch':BRANCH,'verified_remote_base_branch':'codex/phase-4j88-authorized-opls-bundle-qualification','verified_remote_base_head':'4ed7377a66c02362d69da0840b7701ba2b6a797f','fetched_origin_main':subprocess.check_output(['git','rev-parse','origin/main'],text=True).strip(),'local_main':subprocess.check_output(['git','rev-parse','main'],text=True).strip(),'main_modified':False,'historical_files_preserved':len(history),'source_bundle_search':'None. No authorized scientific source was supplied; untracked source assets were not inspected or imported.','unrelated_untracked_files':'preserved'})
shutil.copyfile('/tmp/j89-terminal-collection.json',v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j89-terminal.xml').read_bytes()))
active=Path(os.environ['ISLAND_C1_POLICY']).parent
runs=[]
for directory in sorted(Path('/tmp').glob('island_j89_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory==active: continue
    assert (directory/'result.json').is_file(), directory
    dest=v/'fences'/directory.name; shutil.copytree(directory,dest)
    events=dest/'events.jsonl'; (dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes())); events.unlink()
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json',runs)
logs=v/'logs'; logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j89_*.log')):
    if p.stem!=active.name: (logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
checks=[]
for name in ('terminal_tests','ruff','compileall','pip_supported','J81_receipt','J82_receipt','J83_receipt','J84_receipt','J85_receipt','J86_receipt','J87_receipt','J88_receipt','J81_replay','J82_replay','J83_replay','J84_replay','J85_replay','J86_replay','J87_replay','J88_replay','J80_receipt','J79_receipt','J78_receipt','J76_receipt','C1_receipt','offline_replay','relocated_replay','diff'):
    attempt='island_j89_final_'+('terminal' if name=='terminal_tests' else name)
    result=read(Path('/tmp')/attempt/'result.json'); assert result['returncode']==0, attempt
    checks.append({'name':name,'returncode':0,'log':'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json',checks)
reasons={}

failures=[]
for n in runs:
    r=read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons,n
        failures.append({'run':n,'returncode':r['returncode'],'reason':reasons[Path(n).name],'scientific_execution':False})
write(v/'development-outcomes.json',{'failed_software_attempts':failures,'all_pre_archive_fenced_attempts_with_logs':True,'scientific_execution':False,'maintenance_fences':'Archive, seal and post-seal verification use separate inherited-fence runs. Their /tmp logs remain outside this non-circular evidence seal.'})
cases=list(ET.parse('/tmp/j89-terminal.xml').iter('testcase')); groups={}; skips=[]
for case in cases:
    group='J89_focused' if case.attrib['classname']=='tests.test_phase_4j89' else 'regressions'
    counts=groups.setdefault(group,{k:0 for k in ('passed','skipped','failure','error')})
    outcome=next((n for n in ('failure','error','skipped') if case.find(n) is not None),'passed'); counts[outcome]+=1
    if outcome=='skipped': skips.append({'test':case.attrib['classname']+'::'+case.attrib['name'],'reason':case.find('skipped').attrib.get('message','')})
required=['test_complete_authorized_input_contract_all_eight_graphs_software_only', 'test_relocated_external_replay_requires_originals_not_metadata', 'test_no_source_is_deterministic_and_creates_no_scientific_coverage', 'test_persistence_contains_no_raw_tables_or_coefficient_arrays', 'test_external_strict_receipt_requires_sources_and_pins', 'test_invalid_pin_does_not_open_any_external_input']
for name in required:
    case=next(c for c in cases if c.attrib['classname']=='tests.test_phase_4j89' and c.attrib['name']==name)
    assert all(case.find(n) is None for n in ('failure','error','skipped'))
write(v/'opls-control-outcomes.json',{'tests':required,'passed':True,'source_kind':'manufactured_software_controls_only','eight_graph_software_contract_qualified':True,'source_bound_profile_matching':True,'other_families_rejected':True,'relocated_source_and_fixture_bytes_rehashed':True,'qualification_identity_preserved':True,'scientific_validation':False,'authorized_scientific_bundles_supplied':0,'qualified_scientific_profiles':0})
for name in ('archive','relocate','checks'): shutil.copyfile('/tmp/island_j89_'+name+'.py',v/(name+'.py'))
real_case=next(c for c in cases if c.attrib['classname']=='tests.test_phase_4j89' and c.attrib['name']=='test_real_user_supplied_authorized_bundle')
write(v/'real-source-integration-status.json',{'status':'not_run_missing_user_inputs' if real_case.find('skipped') is not None else 'executed','scientific_qualification_claimed':False,'reason':real_case.find('skipped').attrib.get('message','') if real_case.find('skipped') is not None else None,'synthetic_controls_are_not_real_sources':True})
write(v/'scope-skips.json',{'explicit_skips':skips,'scientific_checks_passed':0,'scientific_execution':{'QM':0,'MD':0,'converter':0,'physical_fitting':0,'parameter_generation':0,'energy_force':0,'reserved_validation':0,'production':0},'synthetic_controls':'Manufactured source tables/audits and software values only; no scientific values accepted. Fake prohibited child operations are rejected before execution.','full_scientific_suite':'Not run; scientific execution is outside authorized scope.'})
write(v/'summary.json',{'tested_implementation_commit':commit,'groups':groups,'unique_terminal_tests':len(cases),'historical_files_preserved':len(history),'archived_fenced_attempts':len(runs),'failed_software_attempts':len(failures),'required_checks':len(checks),'J81_J88_unchanged':True,'authorized_scientific_bundles_supplied':0,'qualified_scientific_profiles':0,'runtime_registration':False,'accepted_coefficients':0,'scientific_validation':'not_established'})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
