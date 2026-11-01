import json,os,subprocess
from pathlib import Path
repo=Path.cwd();out=Path('/tmp/island-j63-validation');old=json.loads((repo/'docs/evidence/phase_4j62/validation/independent-checks.json').read_text());old.append({'name':'J62_receipt','command':[str(repo.parent/'island-validation/phase4g1-env/bin/python'),'scripts/validate_phase_4j62_receipt.py','docs/evidence/phase_4j62','--repo','.']})
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'};results=[]
for row in old:
 with (out/(row['name']+'.log')).open('w') as f:r=subprocess.run(row['command'],env=env,stdout=f,stderr=subprocess.STDOUT)
 results.append({'name':row['name'],'command':row['command'],'returncode':r.returncode,'log':row['name']+'.log'});(out/'independent-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(row['name'],r.returncode,flush=True)
