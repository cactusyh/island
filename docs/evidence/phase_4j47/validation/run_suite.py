import collections,json,os,subprocess,sys,time
from pathlib import Path
root=Path('/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island')
out=Path('/tmp/island-j47-validation')
nodes=[x for x in Path('/tmp/island-j47-validation/final-collection.txt').read_text().splitlines() if x.startswith('tests/') and '::' in x and not x.startswith('tests/test_phase_4j47.py::')]
files=collections.Counter(n.split('::')[0] for n in nodes)
groups=[[] for _ in range(8)]
loads=[0]*8
for name,count in sorted(files.items(),key=lambda x:(-x[1],x[0])):
 i=min(range(8),key=lambda i:loads[i]);groups[i].append(name);loads[i]+=count
(out/'collection.json').write_text(json.dumps({'nodes':nodes,'file_shards':groups,'loads':loads},indent=2)+'\n')
env={**os.environ,'PYTHONPATH':'.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
procs=[]
for i,group in enumerate(groups):
 log=(out/f'full-{i}.log').open('w')
 command=[str(root.parent/'island-validation/phase4g1-env/bin/python'),'-m','pytest','-q',*group,'--junitxml',str(out/f'full-{i}.xml')]
 procs.append((subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT),log,command))
print('Started 8 file shards:',loads,flush=True)
results=[]
for i,(p,log,cmd) in enumerate(procs):
 code=p.wait();log.close();results.append({'shard':i,'returncode':code,'command':cmd});(out/'suite-results.json').write_text(json.dumps(results,indent=2)+'\n')
 print(i,code,flush=True)
sys.exit(any(r['returncode'] for r in results))
