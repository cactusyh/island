import os,time,json,subprocess,sys
from pathlib import Path
repo=Path('/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island')
work=repo.parent/'island-validation/phase4j61-qm'
check=repo.parent/'island-validation/phase4g1-env/bin/python'
qm=repo.parent/'island-validation/phase4j36-native-env/bin/python'
env={**os.environ,'PYTHONPATH':str(repo)+':'+str(repo/'src'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
start=time.monotonic()
while not (work/'optimization/result.json').exists():
 if time.monotonic()-start>7200:sys.exit('optimization wait limit')
 time.sleep(5)
assert json.loads((work/'optimization/result.json').read_text())['status']=='completed'
common=[str(check),'scripts/phase_4j61.py']
subprocess.run(common+['prepare','--work',str(work)],cwd=repo,env=env,check=True)
for group in ('curvature','sensitivity','finite-difference','fresh'):
 subprocess.run(common+['execute','--work',str(work),'--python',str(qm),'--group',group],cwd=repo,env=env,check=True)
subprocess.run(common+['import','--work',str(work)],cwd=repo,env=env,check=True)
print('J61 counted native campaign completed',flush=True)
