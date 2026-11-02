import collections,json,os,subprocess,sys
from pathlib import Path
from xml.etree import ElementTree as ET
repo=Path.cwd();out=Path('/tmp/island-j72-validation/guarded-final');old=Path('/tmp/island-j71-validation');base=json.load(open('/tmp/island-j72-validation/collection.json'))
python=str(repo.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':str(repo.parent/'island-validation/phase4j66-solvers-only')+':.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','ISLAND_J72_FENCE_DIRECTORY':'/tmp/island-j72-validation/fence'}
reuse=json.load(open('/tmp/island-j72-validation/suite-policy.json'))['reused_fixture_paths'];env.update(reuse)
times={}
for i in (0,1):
 for c in ET.parse(old/f'full-{i}.xml').iter('testcase'):
  node=c.attrib['classname'].replace('.','/')+'.py::'+c.attrib['name'];times[node]=float(c.attrib.get('time',1))
groups=[[] for _ in range(8)];loads=[0.]*8
for node in sorted(base['nodes'],key=lambda n:(-times.get(n,1),n)):
 i=min(range(8),key=lambda i:loads[i]);groups[i].append(node);loads[i]+=times.get(node,1)
(out/'collection.json').write_text(json.dumps({'nodes':base['nodes'],'node_shards':groups,'workers':8,'estimated_seconds':loads,'reason':'8 single-thread software workers within8-CPU/32GiB cgroup; independent retained static fixtures; no scientific workers'},indent=2)+'\n')
(out/'suite-policy.json').write_text(json.dumps({'plugin':'scripts.phase_4j72_test_policy','child_fence':'scripts/phase_4j72_execution_fence.py','workers':8,'threads_per_worker':1,'resources':{n:Path('/sys/fs/cgroup/'+n).read_text() for n in ('memory.max','memory.current','cpu.max')},'reuse':reuse,'corrected_run_authorized_QM':0,'corrected_run_authorized_MD_steps':0,'first_run_scope_violation':'docs/evidence/phase_4j72/scope-deviation/execution-deviation.json'},indent=2)+'\n')
processes=[]
for i,group in enumerate(groups):
 log=(out/f'full-{i}.log').open('w');cmd=[python,'-m','pytest','-p','scripts.phase_4j72_test_policy','-q',*group,'--basetemp',str(out/f'pytest-{i}'),'--junitxml',str(out/f'full-{i}.xml')];processes.append((subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT),log,cmd))
print('Started8 guarded software shards',loads,flush=True);results=[]
for i,(p,log,cmd) in enumerate(processes):
 code=p.wait();log.close();results.append({'shard':i,'returncode':code,'command':cmd});(out/'suite-results.json').write_text(json.dumps(results,indent=2)+'\n');print(i,code,flush=True)
sys.exit(any(x['returncode'] for x in results))
