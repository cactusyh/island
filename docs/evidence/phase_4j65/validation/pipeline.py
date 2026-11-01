import json, os, subprocess, sys, time
from pathlib import Path
repo=Path.cwd(); out=Path('/tmp/island-j65-validation'); work=repo.parent/'island-validation/phase4j65-qm'; evidence=repo/'docs/evidence/phase_4j65'
python=repo.parent/'island-validation/phase4g1-env/bin/python'; native=repo.parent/'island-validation/phase4j36-native-env/bin/python'
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
results=[]
def read(p):return json.loads(p.read_text())
def run(name,cmd):
 start=time.time()
 with (out/(name+'.log')).open('w') as log:r=subprocess.run([str(x) for x in cmd],env=env,stdout=log,stderr=subprocess.STDOUT)
 results.append({'name':name,'command':[str(x) for x in cmd],'returncode':r.returncode,'started_unix':start,'ended_unix':time.time(),'log':name+'.log'})
 (out/'pipeline-results.json').write_text(json.dumps(results,indent=2)+'\n');print(name,r.returncode,flush=True)
 if r.returncode:sys.exit(r.returncode)
def execute(stage):run(stage+'-attempt1',[python,'scripts/phase_4j65.py',stage,'--work',work,'--python',native])
def analyze(stage):run('analyze-'+stage+'-attempt1',[python,'scripts/analyze_phase_4j65.py',stage,'--work',work])
keys=[k for k,j in read(work/'jobs.json')['jobs'].items() if j['stage']=='controls']
start=time.time()
while True:
 paths=[work/'processes'/key/'result.json' for key in keys]
 records=[read(p) for p in paths if p.exists()]
 if any(r['returncode']!=0 for r in records):print('controls failure; no dependent jobs launched',flush=True);sys.exit(1)
 if len(records)==len(keys):break
 if time.time()-start>36000:print('coordinator controls wait cap reached',flush=True);sys.exit(124)
 time.sleep(20)
analyze('controls');execute('primary');analyze('primary')
if not read(evidence/'primary-numerical-gates.json')['joints_authorized']:
 print('primary numerical gate failed; no joints, fit or validation activated',flush=True);sys.exit(3)
execute('joints');execute('precision');analyze('trial')
if read(evidence/'trial-fit.json')['validation_authorized']:
 analyze('validation-preflight');execute('validation');analyze('validation')
else:print('training candidate ineligible; three reserved validation configurations untouched',flush=True)
