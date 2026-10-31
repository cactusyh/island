import collections,json,os,subprocess,sys
from pathlib import Path
root=Path.cwd(); out=Path('/tmp/island-j63-validation')
files=subprocess.check_output(['git','ls-files','tests/test*.py'],text=True).splitlines()
env={**os.environ,'PYTHONPATH':'/tmp/island-j63-validation:.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','J62_DIAGNOSTICS':str(out/'subprocess-failures')}
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
r=subprocess.run([python,'-m','pytest','--collect-only','-q',*files],env=env,capture_output=True,text=True)
(out/'base-collection.txt').write_text(r.stdout+r.stderr)
if r.returncode: sys.exit(r.returncode)
nodes=[x for x in r.stdout.splitlines() if x.startswith('tests/') and '::' in x]
counts=collections.Counter(x.split('::')[0] for x in nodes); groups=[[] for _ in range(8)];loads=[0]*8
for name,n in sorted(counts.items(),key=lambda x:(-x[1],x[0])):
 i=min(range(8),key=lambda i:loads[i]);groups[i].append(name);loads[i]+=n
(out/'collection.json').write_text(json.dumps({'nodes':nodes,'file_shards':groups,'loads':loads},indent=2)+'\n')
procs=[]
for i,group in enumerate(groups):
 log=(out/f'full-{i}.log').open('w');cmd=[python,'-m','pytest','-p','failure_capture','-q',*group,'--basetemp',str(out/f'pytest-shard-{i}'),'--junitxml',str(out/f'full-{i}.xml')]
 procs.append((subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT),log,cmd))
print('Started8 shards',loads,flush=True);results=[]
for i,(p,log,cmd) in enumerate(procs):
 code=p.wait();log.close();results.append({'shard':i,'returncode':code,'command':cmd});(out/'suite-results.json').write_text(json.dumps(results,indent=2)+'\n');print(i,code,flush=True)
sys.exit(any(x['returncode'] for x in results))
