import json,os,subprocess,sys,time
from pathlib import Path
repo=Path('/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island')
work=repo.parent/'island-validation/phase4j61-qm'
python=repo.parent/'island-validation/phase4g1-env/bin/python'
evidence=repo/'docs/evidence/phase_4j61'
env={**os.environ,'PYTHONPATH':str(repo)+':'+str(repo/'src'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
start=time.monotonic()
while not (evidence/'normalized-reference.json').exists():
    if time.monotonic()-start>29000: sys.exit('bounded native campaign wait expired')
    time.sleep(20)
with Path('/tmp/island-j61-validation/analysis-attempt1.log').open('w') as log:
    result=subprocess.run([str(python),'scripts/analyze_phase_4j61.py','fit','--work',str(work)],cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
if result.returncode:sys.exit(result.returncode)
with Path('/tmp/island-j61-validation/retain-attempt1.log').open('w') as log:
    result=subprocess.run([str(python),'scripts/replay_phase_4j61.py','retain','--work',str(work)],cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
print('J61 analysis and retention',result.returncode,flush=True)
sys.exit(result.returncode)
