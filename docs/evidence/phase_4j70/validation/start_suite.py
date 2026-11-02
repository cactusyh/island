import collections,json,os,subprocess,sys
from pathlib import Path
from xml.etree import ElementTree as ET
root=Path.cwd();out=Path('/tmp/island-j70-validation');old=Path('/tmp/island-j69-validation');fixtures=Path('/tmp/island-j68-validation')
files=sorted(set(subprocess.check_output(['git','ls-files','tests/test*.py'],text=True).splitlines()+['tests/test_phase_4j70.py']))
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':str(root.parent/'island-validation/phase4j66-solvers-only')+':.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
reuse={'ISLAND_J32_TEST_RUN':'pytest-terminal-1/j32-bundles0','ISLAND_J33_C1_RUN':'pytest-terminal-1/j33-c10/runs','ISLAND_J34_RUN':'pytest-terminal-5/j340/runs','ISLAND_J35_RUN':'pytest-terminal-1/j350/runs','ISLAND_J39_RUN':'pytest-terminal-3/j39-retained0/run','ISLAND_J40_RUN':'pytest-terminal-2/j400/run','ISLAND_J40_C2_RUN':'pytest-terminal-0/j40-c20/run'}
for k,v in reuse.items():
 p=fixtures/v
 if p.is_dir():env[k]=str(p)
resources={n:Path('/sys/fs/cgroup/'+n).read_text() for n in ('memory.max','memory.current','cpu.max')}
(out/'suite-policy.json').write_text(json.dumps({'plugin':'scripts.phase_4j69_test_policy','reused_fixture_paths':{k:env[k] for k in reuse if k in env},'workers':2,'resources':resources,'research_QM':0,'MD_steps':0,'scope_skips_not_passes':True},indent=2)+'\n')
r=subprocess.run([python,'-m','pytest','-p','scripts.phase_4j69_test_policy','--collect-only','-q',*files],env=env,capture_output=True,text=True,check=False)
(out/'base-collection.txt').write_text(r.stdout+r.stderr)
if r.returncode:sys.exit(r.returncode)
nodes=[x for x in r.stdout.splitlines() if x.startswith('tests/') and '::' in x]
times=collections.Counter()
for i in range(2):
 for row in ET.parse(old/f'full-{i}.xml').iter('testcase'):
  times[row.attrib['classname'].replace('.','/')+'.py']+=float(row.attrib.get('time',0))
times['tests/test_phase_4j70.py']=90
groups=[[] for _ in range(2)];loads=[0.,0.]
for name in sorted(set(n.split('::')[0] for n in nodes),key=lambda f:(-times.get(f,1.),f)):
 i=min(range(2),key=lambda i:loads[i]);groups[i].append(name);loads[i]+=times.get(name,1.)
(out/'collection.json').write_text(json.dumps({'nodes':nodes,'file_shards':groups,'workers':2,'estimated_seconds_from_previous_run':loads,'focused_reruns_add_unique_coverage':False},indent=2)+'\n')
processes=[]
for i,group in enumerate(groups):
 log=(out/f'full-{i}.log').open('w');cmd=[python,'-m','pytest','-p','scripts.phase_4j69_test_policy','-q',*group,'--basetemp',str(out/f'pytest-terminal-{i}'),'--junitxml',str(out/f'full-{i}.xml')]
 processes.append((subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT),log,cmd))
print('Started two scope-protected full shards',loads,flush=True);results=[]
for i,(p,log,cmd) in enumerate(processes):
 code=p.wait();log.close();results.append({'shard':i,'returncode':code,'command':cmd});(out/'suite-results.json').write_text(json.dumps(results,indent=2)+'\n');print(i,code,flush=True)
sys.exit(any(x['returncode'] for x in results))
