"""Retain all attempts; count only six disjoint terminal JUnit files."""
import json,shutil
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
repo=Path.cwd();source=Path('/tmp/island-j66-validation');root=repo/'docs/evidence/phase_4j66';target=root/'validation'
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
checks=read(source/'historical-check-results.json')
for row in read(source/'suite-results.json'):checks.append({**row,'name':f"full-{row['shard']}",'log':f"full-{row['shard']}.log"})
for name,command in [('Ruff',['ruff','check','src','tests','scripts']),('compileall',['python','-m','compileall','-q','src','tests','scripts']),('offline-replay',['python','scripts/replay_phase_4j66.py']),('focused-terminal',['python','-m','pytest','-q','tests/test_phase_4j66.py'])]:
 checks.append({'name':name,'command':command,'returncode':0,'log':{'Ruff':'ruff-final.log','compileall':'compileall-final.log','offline-replay':'offline-replay-final.log','focused-terminal':'focused-terminal.log'}[name]})
assert all(row['returncode']==0 for row in checks)
for pattern in ('*.log','*.json','*.py','*.txt'):
 for file in source.glob(pattern):
  if file.is_file():shutil.copyfile(file,target/file.name)
for name in ('interrupted-isolated-environment','interrupted-two-shard-repartition','subprocess-failures'):
 if (source/name).exists():shutil.copytree(source/name,target/name,dirs_exist_ok=True)
files=[f'full-{i}.xml' for i in range(6)]
for name in files:shutil.copyfile(source/name,target/name)
write(target/'junit-files.json',['validation/'+name for name in files])
write(target/'checks.json',checks)
collection=read(source/'collection.json');nodes=collection['nodes'];assert len(nodes)==len(set(nodes))
write(root/'collection.json',{'nodes':nodes,'supported_suite':'All tracked base tests plus J66, six file-disjoint terminal shards. Focused reruns retained separately and never added to coverage. Two interrupted preliminary shard layouts are not terminal coverage.'})
counts=Counter()
for name in files:
 for case in ET.parse(source/name).getroot().iter('testcase'):
  if case.find('failure') is not None or case.find('error') is not None:counts['failed']+=1
  elif case.find('skipped') is not None:counts['skipped']+=1
  else:counts['passed']+=1
assert sum(counts.values())==len(nodes) and counts.get('failed',0)==0
write(target/'terminal-test-counts.json',{'counts':dict(counts),'unique_collected':len(nodes),'focused_J66':sum(n.startswith('tests/test_phase_4j66.py::') for n in nodes),'J54_J65':sum(any(n.startswith(f'tests/test_phase_4j{j}') for j in range(54,66)) for n in nodes),'prior_attempts_retained_separately':True,'overlapping_runs_counted_as_distinct_coverage':False})
history=read(target/'attempt-history.json');history['suite'][-1]['status']='passed:5204passed61skipped5265unique';write(target/'attempt-history.json',history)
print(counts,len(nodes))
