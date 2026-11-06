"""Archive exact J84 fenced software checks; no source or scientific execution."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j84_bundles import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo=Path.cwd();root=repo/ROOT;v=root/'validation';v.mkdir(exist_ok=False)
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
write(v/'tested-implementation.json',{'commit':commit,'files':{n:pin(repo/n) for n in IMPLEMENTATION}})
write(v/'repository-state.json',{'base_commit':BASE,'implementation_commit':commit,'branch':BRANCH,'fetched_origin_main':subprocess.check_output(['git','rev-parse','origin/main'],text=True).strip(),'local_main':subprocess.check_output(['git','rev-parse','main'],text=True).strip(),'main_modified':False,'historical_files_preserved':len(historical(repo)),'source_bundle_search':'None; no authorized scientific source bundle supplied; untracked source assets were not inspected or imported.','unrelated_untracked_files':'preserved'})
shutil.copyfile('/tmp/j84-terminal-collection.json',v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j84-terminal.xml').read_bytes()))
active=Path(os.environ['ISLAND_C1_POLICY']).parent
runs=[]
for directory in sorted(Path('/tmp').glob('island_j84_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory==active:continue
    assert (directory/'result.json').is_file(),directory
    dest=v/'fences'/directory.name;shutil.copytree(directory,dest)
    events=dest/'events.jsonl';(dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes()));events.unlink()
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json',runs)
logs=v/'logs';logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j84_*.log')):
    if p.stem != active.name:(logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
checks=[]
for name in ('terminal_tests','ruff','compileall','pip_supported','J81_receipt','J82_receipt','J83_receipt','J81_replay','J82_replay','J83_replay','J80_receipt','J79_receipt','J78_receipt','J76_receipt','C1_receipt','offline_replay','relocated_replay','diff'):
    label='terminal' if name=='terminal_tests' else name
    attempt='island_j84_final_'+label
    result=read(Path('/tmp')/attempt/'result.json');assert result['returncode']==0,attempt
    checks.append({'name':name,'returncode':0,'log':'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json',checks)
reasons={'island_j84_lint_01':'Initial new-file formatting, dictionary style and Python 3.11 UTC parsing style checks; corrected in J84 files only.','island_j84_focused_01':'One mutation control selected the opaque attachment (raw_text=None) instead of the parsed record source; corrected the control to select a parsed record explicitly.'}
failures=[]
for n in runs:
    r=read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons,n
        failures.append({'run':n,'returncode':r['returncode'],'reason':reasons[Path(n).name],'scientific_execution':False})
write(v/'development-outcomes.json',{'failed_software_attempts':failures,'scientific_execution':False,'all_pre_archive_fenced_attempts_with_logs':True,'maintenance_fences':'Archive, seal and post-seal validation are separate inherited-fence maintenance runs; their /tmp logs remain outside this non-circular seal.'})
cases=list(ET.parse('/tmp/j84-terminal.xml').iter('testcase'));groups={};skips=[]
for case in cases:
    group='J84_focused' if case.attrib['classname']=='tests.test_phase_4j84' else 'regressions'
    counts=groups.setdefault(group,{k:0 for k in ('passed','skipped','failure','error')})
    outcome=next((n for n in ('failure','error','skipped') if case.find(n) is not None),'passed');counts[outcome]+=1
    if outcome=='skipped':skips.append({'test':case.attrib['classname']+'::'+case.attrib['name'],'reason':case.find('skipped').attrib.get('message','')})
relocation=next(c for c in cases if c.attrib['classname']=='tests.test_phase_4j84' and c.attrib['name']=='test_bundle_reingestion_after_directory_relocation')
assert all(relocation.find(n) is None for n in ('failure','error','skipped'))
write(v/'bundle-reingestion-relocation.json',{'test':'tests/test_phase_4j84.py::test_bundle_reingestion_after_directory_relocation','passed':True,'source_kind':'manufactured_software_control_only','source_files':'Canonical synthetic records and a synthetic original attachment; complete independent audit fixture.','bundle_identity_preserved':True,'profile_identity_preserved':True,'coverage_replay_identity_preserved':True,'different_capture_timestamp_does_not_change_bundle_or_profile_identity':True,'scientific_validation':False})
for name in ('relocate','archive'):shutil.copyfile('/tmp/island_j84_'+name+'.py',v/(name+'.py'))
write(v/'scope-skips.json',{'explicit_skips':skips,'scientific_checks_passed':0,'scientific_execution':{'QM':0,'MD':0,'converter':0,'physical_fitting':0,'parameter_generation':0,'energy_force':0,'reserved_validation':0,'production':0},'synthetic_controls':'Manufactured source tables/audits and software values only. No scientific force-field values accepted. Fake executable child operations are rejected before execution.','full_scientific_suite':'Not run: no scientific execution authorized.'})
write(v/'summary.json',{'tested_implementation_commit':commit,'groups':groups,'unique_terminal_tests':len(cases),'historical_files_preserved':len(historical(repo)),'archived_fenced_attempts':len(runs),'failed_software_attempts':len(failures),'required_checks':len(checks),'J81_J82_J83_unchanged':True,'authorized_scientific_bundles_supplied':0,'runtime_registration':False,'scientific_validation':'not_established'})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
