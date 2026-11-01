import json, shutil, xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

repo=Path.cwd(); source=Path('/tmp/island-j65-validation'); root=repo/'docs/evidence/phase_4j65'; target=root/'validation'; target.mkdir(exist_ok=True)
def read(p): return json.loads(p.read_text())
def write(p,r): p.write_text(json.dumps(r,indent=2)+'\n')
checks=read(source/'historical-check-results.json')
for r in read(source/'suite-results.json'):
 checks.append({**r,'name':f"full-{r['shard']}",'log':f"full-{r['shard']}.log"})
checks.extend(read(source/'final-check-results.json'))
assert all(r['returncode']==0 for r in checks), checks
files=['full-0.xml','full-1.xml','focused-final.xml']
for name in files: shutil.copyfile(source/name,target/name)
write(target/'junit-files.json',['validation/'+n for n in files])
for pattern in ['*.log','*.json','*.py','*.txt','focused-*.xml']:
 for p in source.glob(pattern):
  if p.name not in files: shutil.copyfile(p,target/p.name)
if (source/'subprocess-failures').exists():
 shutil.copytree(source/'subprocess-failures',target/'subprocess-failures',dirs_exist_ok=True)
controls_log=Path('/tmp/island-j65-controls-attempt1.log')
if controls_log.exists(): shutil.copyfile(controls_log,target/'controls-coordinator-attempt1.log')
write(target/'checks.json',checks)
nodes=read(source/'collection.json')['nodes']
focused=read(source/'final-focused-collection.json')['nodes']
nodes += focused
assert len(nodes)==len(set(nodes))
# Existing receipt JUnit parser checks exact terminal test identity coverage.
write(root/'collection.json',{'nodes':nodes,'supported_suite':'All tracked base test files plus final J65 tests, two file-disjoint base shards and one J65 shard. Native QM is not called by tests.'})
counts=Counter()
for name in files:
 for case in ET.parse(source/name).getroot().iter('testcase'):
  if case.find('failure') is not None or case.find('error') is not None: counts['failed']+=1
  elif case.find('skipped') is not None: counts['skipped']+=1
  else: counts['passed']+=1
write(target/'terminal-test-counts.json',{'counts':dict(counts),'unique_collected':len(nodes),'J54_J64_regressions':sum(any(n.startswith(f'tests/test_phase_4j{j}') for j in range(54,65)) for n in nodes),'final_J65_cases':len(focused),'terminal_results_only':True,'prior_attempts_retained_separately':True})
print(json.dumps({'counts':dict(counts),'unique_collected':len(nodes)}))
