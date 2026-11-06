"""Archive exact J85 fenced software checks; no source or scientific execution."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j85_qualification import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo=Path.cwd();root=repo/ROOT;v=root/'validation';v.mkdir(exist_ok=False)
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
write(v/'tested-implementation.json',{'commit':commit,'files':{n:pin(repo/n) for n in IMPLEMENTATION}})
write(v/'repository-state.json',{'base_commit':BASE,'implementation_commit':commit,'branch':BRANCH,'fetched_origin_main':subprocess.check_output(['git','rev-parse','origin/main'],text=True).strip(),'local_main':subprocess.check_output(['git','rev-parse','main'],text=True).strip(),'main_modified':False,'historical_files_preserved':len(historical(repo)),'source_bundle_search':'None; no authorized scientific source bundle supplied; untracked source assets were not inspected or imported.','unrelated_untracked_files':'preserved'})
shutil.copyfile('/tmp/j85-terminal-verified-collection.json',v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j85-terminal-verified.xml').read_bytes()))
active=Path(os.environ['ISLAND_C1_POLICY']).parent
runs=[]
for directory in sorted(Path('/tmp').glob('island_j85_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory==active:continue
    assert (directory/'result.json').is_file(),directory
    dest=v/'fences'/directory.name;shutil.copytree(directory,dest)
    events=dest/'events.jsonl';(dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes()));events.unlink()
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json',runs)
logs=v/'logs';logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j85_*.log')):
    if p.stem != active.name:(logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
checks=[]
for name in ('terminal_tests','ruff','compileall','pip_supported','J81_receipt','J82_receipt','J83_receipt','J84_receipt','J81_replay','J82_replay','J83_replay','J84_replay','J80_receipt','J79_receipt','J78_receipt','J76_receipt','C1_receipt','offline_replay','relocated_replay','diff'):
    label={'terminal_tests':'terminal_verified','ruff':'ruff_verified','compileall':'compileall_verified','diff':'diff_verified','offline_replay':'offline_replay_verified','relocated_replay':'relocated_replay_verified','J81_replay':'J81_replay_corrected'}.get(name,name)
    attempt='island_j85_final_'+label
    result=read(Path('/tmp')/attempt/'result.json');assert result['returncode']==0,attempt
    checks.append({'name':name,'returncode':0,'log':'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json',checks)
reasons={'island_j85_format_01':'Initial new-file Ruff checks found set-comprehension and loop-variable binding issues; corrected only new J85 files.','island_j85_focused_01':'Two acceptance-branch controls correctly blocked inherited H-capped ether/ester graphs because terminal alcohol/acid environments are outside the unchanged J83 motif catalog. J85 now declares methyl caps for its benchmark, with additional controls preserving H-capped blocking. Original run: 62 passed, 2 failed.','island_j85_final_J81_replay':'A replay command used the nonexistent module scripts.phase_4j81_source_ingestion. Corrected to the historical scripts.phase_4j81_ingest; original failure retained.'}
failures=[]
for n in runs:
    r=read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons,n
        failures.append({'run':n,'returncode':r['returncode'],'reason':reasons[Path(n).name],'scientific_execution':False})
old=Path('/tmp/island_j85_prequalification_evidence');superseded=v/'superseded-benchmark';superseded.mkdir()
for source in sorted(old.iterdir()):
    assert source.is_file(),source
    (superseded/(source.name+'.xz')).write_bytes(lzma.compress(source.read_bytes()))
write(v/'development-outcomes.json',{'failed_software_attempts':failures,'scientific_execution':False,'all_pre_archive_fenced_attempts_with_logs':True,'superseded_passing_terminal_run':{'run':'validation/fences/island_j85_final_terminal','implementation_commit':'21fc0f0ba309cd76a74d5b2336605ecc2c0ec648','passed':1246,'skipped':58,'reason':'Follow-up source-authorization status refinement; final terminal verification is a complete new run at the recorded tested implementation. Tests are not double-counted.','superseded_unsealed_benchmark':{p.name:pin(p) for p in old.iterdir()}},'maintenance_fences':'Archive, seal and post-seal validation are separate inherited-fence maintenance runs; their /tmp logs remain outside this non-circular seal.'})
write(v/'superseded-terminal-runs.json',{'runs':[{'run':'validation/fences/island_j85_final_terminal','implementation_commit':'21fc0f0ba309cd76a74d5b2336605ecc2c0ec648','passed':1246,'skipped':58,'reason':'Source authorization status refinement.'},{'run':'validation/fences/island_j85_final_terminal_corrected','implementation_commit':'2b9be188121f50b161992d2237e22d8f624c5a08','passed':1246,'skipped':58,'reason':'Final graph/ranking iterator materialization and empty-input controls.'}],'terminal_counts_include_only':'final_terminal_verified; superseded runs are not double-counted'})
cases=list(ET.parse('/tmp/j85-terminal-verified.xml').iter('testcase'));groups={};skips=[]
for case in cases:
    group='J85_focused' if case.attrib['classname']=='tests.test_phase_4j85' else 'regressions'
    counts=groups.setdefault(group,{k:0 for k in ('passed','skipped','failure','error')})
    outcome=next((n for n in ('failure','error','skipped') if case.find(n) is not None),'passed');counts[outcome]+=1
    if outcome=='skipped':skips.append({'test':case.attrib['classname']+'::'+case.attrib['name'],'reason':case.find('skipped').attrib.get('message','')})
relocation=next(c for c in cases if c.attrib['classname']=='tests.test_phase_4j85' and c.attrib['name']=='test_relocated_reingestion_preserves_qualification')
assert all(relocation.find(n) is None for n in ('failure','error','skipped'))
write(v/'bundle-reingestion-relocation.json',{'test':'tests/test_phase_4j85.py::test_relocated_reingestion_preserves_qualification','passed':True,'source_kind':'manufactured_software_control_only','source_files':'Canonical synthetic records and a synthetic original attachment; complete manufactured audit fixture exercising the trusted-input software branch.','bundle_identity_preserved':True,'qualification_identity_preserved':True,'fresh_source_bytes_rehashed':True,'scientific_validation':False})
for name in ('relocate','archive'):shutil.copyfile('/tmp/island_j85_'+name+'.py',v/(name+'.py'))
write(v/'scope-skips.json',{'explicit_skips':skips,'scientific_checks_passed':0,'scientific_execution':{'QM':0,'MD':0,'converter':0,'physical_fitting':0,'parameter_generation':0,'energy_force':0,'reserved_validation':0,'production':0},'synthetic_controls':'Manufactured source tables/audits and software values only. No scientific force-field values accepted. Fake executable child operations are rejected before execution.','full_scientific_suite':'Not run: no scientific execution authorized.'})
write(v/'summary.json',{'tested_implementation_commit':commit,'groups':groups,'unique_terminal_tests':len(cases),'historical_files_preserved':len(historical(repo)),'archived_fenced_attempts':len(runs),'failed_software_attempts':len(failures),'required_checks':len(checks),'J81_J82_J83_J84_unchanged':True,'authorized_scientific_bundles_supplied':0,'runtime_registration':False,'scientific_validation':'not_established'})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
