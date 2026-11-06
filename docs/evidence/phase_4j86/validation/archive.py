"""Archive actual fenced J86 software verification without scientific execution."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j86_qualification import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo=Path.cwd(); root=repo/ROOT; v=root/'validation'; v.mkdir(exist_ok=False)
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
history=historical(repo)
write(v/'tested-implementation.json', {'commit':commit,'files':{n:pin(repo/n) for n in IMPLEMENTATION}})
write(v/'repository-state.json', {'base_commit':BASE,'implementation_commit':commit,'branch':BRANCH,'verified_remote_base_branch':'codex/phase-4j85-authorized-bundle-qualification','verified_remote_base_head':'19fa2e76cd822987618cf03e42fca9dd9d884975','fetched_origin_main':subprocess.check_output(['git','rev-parse','origin/main'],text=True).strip(),'local_main':subprocess.check_output(['git','rev-parse','main'],text=True).strip(),'main_modified':False,'historical_files_preserved':len(history),'source_bundle_search':'None. No authorized scientific source was supplied; untracked source assets were not inspected or imported.','unrelated_untracked_files':'preserved'})
shutil.copyfile('/tmp/j86-terminal-collection.json',v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j86-terminal.xml').read_bytes()))
active=Path(os.environ['ISLAND_C1_POLICY']).parent
runs=[]
for directory in sorted(Path('/tmp').glob('island_j86_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory==active: continue
    assert (directory/'result.json').is_file(), directory
    dest=v/'fences'/directory.name; shutil.copytree(directory,dest)
    events=dest/'events.jsonl'; (dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes())); events.unlink()
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json',runs)
logs=v/'logs'; logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j86_*.log')):
    if p.stem!=active.name: (logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
checks=[]
for name in ('terminal_tests','ruff','compileall','pip_supported','J81_receipt','J82_receipt','J83_receipt','J84_receipt','J85_receipt','J81_replay','J82_replay','J83_replay','J84_replay','J85_replay','J80_receipt','J79_receipt','J78_receipt','J76_receipt','C1_receipt','offline_replay','relocated_replay','diff'):
    attempt='island_j86_final_'+('terminal' if name=='terminal_tests' else name)
    result=read(Path('/tmp')/attempt/'result.json'); assert result['returncode']==0, attempt
    checks.append({'name':name,'returncode':0,'log':'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json',checks)
reasons={
 'island_j86_focused_01':'A test helper consumed a one-shot result iterator while constructing declarations; aggregate correctly rejected the omitted results. Fixed only the new test helper. 44 passed, 1 failed.',
 'island_j86_focused_02':'After missing sources were excluded from candidate ranking, a mutation control copied one empty list over another and made no mutation. Replaced the control with a fabricated cross-family ranking entry. 49 passed, 1 failed.',
}
failures=[]
for n in runs:
    r=read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons,n
        failures.append({'run':n,'returncode':r['returncode'],'reason':reasons[Path(n).name],'scientific_execution':False})
write(v/'development-outcomes.json',{'failed_software_attempts':failures,'all_pre_archive_fenced_attempts_with_logs':True,'scientific_execution':False,'maintenance_fences':'Archive, seal and post-seal verification use separate inherited-fence runs. Their /tmp logs remain outside this non-circular evidence seal.'})
cases=list(ET.parse('/tmp/j86-terminal.xml').iter('testcase')); groups={}; skips=[]
for case in cases:
    group='J86_focused' if case.attrib['classname']=='tests.test_phase_4j86' else 'regressions'
    counts=groups.setdefault(group,{k:0 for k in ('passed','skipped','failure','error')})
    outcome=next((n for n in ('failure','error','skipped') if case.find(n) is not None),'passed'); counts[outcome]+=1
    if outcome=='skipped': skips.append({'test':case.attrib['classname']+'::'+case.attrib['name'],'reason':case.find('skipped').attrib.get('message','')})
required=['test_authorized_opls_and_four_missing_profiles_are_partial','test_incomplete_pcff_cannot_hide_behind_complete_opls','test_opls_and_gaff2_rank_separately','test_partial_report_deterministic_and_relocated','test_aggregate_rechecks_complete_counts_and_authorized_identity']
for name in required:
    case=next(c for c in cases if c.attrib['classname']=='tests.test_phase_4j86' and c.attrib['name']==name)
    assert all(case.find(n) is None for n in ('failure','error','skipped'))
write(v/'scoped-control-outcomes.json',{'tests':required,'passed':True,'source_kind':'manufactured_software_controls_only','partial_bundle_decision':'partial_bundle_set','opls_with_four_missing_profiles_qualified_globally':False,'incomplete_pcff_preserved':True,'cross_family_ranking_rejected':True,'relocated_source_bytes_rehashed':True,'qualification_identity_preserved':True,'scientific_validation':False,'authorized_scientific_bundles_supplied':0})
for name in ('archive','relocate','checks'): shutil.copyfile('/tmp/island_j86_'+name+'.py',v/(name+'.py'))
write(v/'scope-skips.json',{'explicit_skips':skips,'scientific_checks_passed':0,'scientific_execution':{'QM':0,'MD':0,'converter':0,'physical_fitting':0,'parameter_generation':0,'energy_force':0,'reserved_validation':0,'production':0},'synthetic_controls':'Manufactured source tables/audits and software values only; no scientific values accepted. Fake prohibited child operations are rejected before execution.','full_scientific_suite':'Not run; scientific execution is outside authorized scope.'})
write(v/'summary.json',{'tested_implementation_commit':commit,'groups':groups,'unique_terminal_tests':len(cases),'historical_files_preserved':len(history),'archived_fenced_attempts':len(runs),'failed_software_attempts':len(failures),'required_checks':len(checks),'J81_J85_unchanged':True,'authorized_scientific_bundles_supplied':0,'qualified_scientific_profiles':0,'runtime_registration':False,'accepted_coefficients':0,'scientific_validation':'not_established'})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
