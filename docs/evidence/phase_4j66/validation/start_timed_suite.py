"""Bounded six-worker terminal suite, balanced with retained J65 file timings."""
import collections,json,os,subprocess,sys
from pathlib import Path
from xml.etree import ElementTree as ET
root=Path.cwd();out=Path('/tmp/island-j66-validation')
files=sorted(set(subprocess.check_output(['git','ls-files','tests/test*.py'],text=True).splitlines()+['tests/test_phase_4j66.py']))
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':str(root.parent/'island-validation/phase4j66-solvers-only')+':/tmp/island-j66-validation:.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','J62_DIAGNOSTICS':str(out/'subprocess-failures')}
r=subprocess.run([python,'-m','pytest','--collect-only','-q',*files],env=env,capture_output=True,text=True,check=False)
(out/'base-collection.txt').write_text(r.stdout+r.stderr)
if r.returncode:sys.exit(r.returncode)
nodes=[x for x in r.stdout.splitlines() if x.startswith('tests/') and '::' in x]
times=collections.Counter()
for name in ('full-0.xml','full-1.xml','focused-final.xml'):
 for row in ET.parse(root/'docs/evidence/phase_4j65/validation'/name).iter('testcase'):
  file=row.attrib['classname'].replace('.','/')+'.py'
  times[file]+=float(row.attrib.get('time',0))
times['tests/test_phase_4j66.py']=30
counts=collections.Counter(x.split('::')[0] for x in nodes)
groups=[[] for _ in range(6)];loads=[0.]*6
for name in sorted(counts,key=lambda f:(-times.get(f,1.),f)):
 i=min(range(6),key=lambda i:loads[i]);groups[i].append(name);loads[i]+=times.get(name,1.)
(out/'collection.json').write_text(json.dumps({'nodes':nodes,'file_shards':groups,'estimated_seconds_from_J65':loads,'workers':6,'reason':'J65 measured110CPUminutes; six single-thread workers fit8CPU cgroup with unrelated jobs preserved.'},indent=2)+'\n')
processes=[]
for i,group in enumerate(groups):
 log=(out/f'full-{i}.log').open('w')
 cmd=[python,'-m','pytest','-p','failure_capture','-q',*group,'--basetemp',str(out/f'pytest-terminal-{i}'),'--junitxml',str(out/f'full-{i}.xml')]
 processes.append((subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT),log,cmd))
print('Started6 shards',loads,flush=True);results=[]
for i,(p,log,cmd) in enumerate(processes):
 code=p.wait();log.close();results.append({'shard':i,'returncode':code,'command':cmd});(out/'suite-results.json').write_text(json.dumps(results,indent=2)+'\n');print(i,code,flush=True)
sys.exit(any(x['returncode'] for x in results))
