"""Archive all J83 fenced software attempts and exact terminal evidence."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import pin, read, write
from scripts.phase_4j83_registry import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, scope

scope()
repo = Path.cwd()
root = repo / ROOT
v = root / 'validation'
v.mkdir(exist_ok=False)
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
write(v / 'tested-implementation.json', {'commit': commit, 'files': {n: pin(repo / n) for n in IMPLEMENTATION}})
write(v / 'repository-state.json', {
    'base_commit': BASE, 'implementation_commit': commit, 'branch': BRANCH,
    'fetched_origin_main': subprocess.check_output(['git','rev-parse','origin/main'],text=True).strip(),
    'local_main': subprocess.check_output(['git','rev-parse','main'],text=True).strip(),
    'main_modified': False, 'historical_files_preserved': len(historical(repo)),
    'unrelated_untracked_files': 'preserved; no literature/source assets read or imported',
})
shutil.copyfile('/tmp/j83-terminal-final-collection.json', v/'collection.json')
(v/'terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j83-terminal-final.xml').read_bytes()))
(v/'superseded-terminal.xml.xz').write_bytes(lzma.compress(Path('/tmp/j83-terminal.xml').read_bytes()))
active = Path(os.environ['ISLAND_C1_POLICY']).parent
runs = []
for directory in sorted(Path('/tmp').glob('island_j83_*')):
    if not directory.is_dir() or not (directory/'policy.json').is_file() or directory == active:
        continue
    assert (directory/'result.json').is_file(), directory
    dest = v/'fences'/directory.name
    shutil.copytree(directory, dest)
    events = dest/'events.jsonl'
    (dest/'events.jsonl.xz').write_bytes(lzma.compress(events.read_bytes()))
    events.unlink()  # Archive copy only; original events remain intact in /tmp.
    runs.append(dest.relative_to(root).as_posix())
write(v/'fenced-runs.json', runs)
logs=v/'logs';logs.mkdir()
for p in sorted(Path('/tmp').glob('island_j83_*.log')):
    if p.stem != active.name:
        (logs/(p.name+'.xz')).write_bytes(lzma.compress(p.read_bytes()))
checks=[]
for name in ('terminal_tests','ruff','compileall','pip_supported','J81_receipt','J82_receipt','J81_replay','J82_replay','J80_receipt','J79_receipt','J78_receipt','J76_receipt','C1_receipt','offline_replay','relocated_replay','diff'):
    label='terminal' if name=='terminal_tests' else name
    revision='_2' if label in ('terminal','ruff','compileall','J81_replay') else ''
    attempt='island_j83_final_'+label+revision
    result=read(Path('/tmp')/attempt/'result.json')
    assert result['returncode']==0, attempt
    checks.append({'name':name,'returncode':0,'log':'validation/logs/'+attempt+'.log.xz'})
write(v/'checks.json',checks)
reasons={
    'island_j83_focused_01':'Four software-control failures: uncapped single-repeat ether/ester were alcohol/acid, ethane lacked a five-site torsion-torsion path, and the class mutation initially wrote the existing II value. Corrected only J83 fixtures/tests; no scientific tolerance changed.',
    'island_j83_lint_01':'Initial new-file import formatting, __all__ ordering, dataclass default factory and dictionary style corrections.',
    'island_j83_final_J81_replay':'Incorrect module name in a software command; retried with the unchanged existing scripts.phase_4j81_ingest entry point.',
}
failures=[]
for n in runs:
    r=read(root/n/'result.json')
    if r['returncode']:
        assert Path(n).name in reasons, n
        failures.append({'run':n,'returncode':r['returncode'],'reason':reasons[Path(n).name],'scientific_execution':False})
write(v/'development-outcomes.json',{'failed_software_attempts':failures,'scientific_execution':False,'superseded_terminal':'The first terminal suite started at 4888f53 before report readiness and duplicate-cell schema hardening. It is preserved and superseded by the complete terminal_2 run at the tested implementation commit.','all_pre_archive_fenced_attempts_with_logs':True,'maintenance_fences':'Archive, seal and post-seal validation are separate inherited-fence maintenance runs; their /tmp logs remain outside this non-circular seal.'})
cases=list(ET.parse('/tmp/j83-terminal-final.xml').iter('testcase'))
groups={};skips=[]
for case in cases:
    group='J83_focused' if case.attrib['classname']=='tests.test_phase_4j83' else 'regressions'
    counts=groups.setdefault(group,{k:0 for k in ('passed','skipped','failure','error')})
    outcome=next((n for n in ('failure','error','skipped') if case.find(n) is not None),'passed')
    counts[outcome]+=1
    if outcome=='skipped':
        skips.append({'test':case.attrib['classname']+'::'+case.attrib['name'],'reason':case.find('skipped').attrib.get('message','')})
for name in ('relocate','archive'):
    shutil.copyfile('/tmp/island_j83_'+name+'.py',v/(name+'.py'))
write(v/'scope-skips.json',{'explicit_skips':skips,'scientific_checks_passed':0,'scientific_execution':{'QM':0,'MD':0,'converter':0,'physical_fitting':0,'parameter_generation':0,'energy_force':0,'reserved_validation':0,'production':0},'synthetic_controls':'Manufactured graph/rule/parameter/audit/receipt controls only; fake executables cannot execute. No scientific parameter data accepted.','full_scientific_suite':'Not run: no scientific execution authorized.'})
write(v/'summary.json',{'tested_implementation_commit':commit,'groups':groups,'unique_terminal_tests':len(cases),'historical_files_preserved':len(historical(repo)),'archived_fenced_attempts':len(runs),'failed_software_attempts':len(failures),'required_checks':len(checks),'native_pcff_unchanged':True,'runtime_registration':False,'scientific_validation':'not_established'})
print(json.dumps(read(v/'summary.json'),sort_keys=True))
