import concurrent.futures,json,os,subprocess
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j62-validation');prior=json.loads((root/'docs/evidence/phase_4j61/validation/checks.json').read_text());tasks=[r for r in prior if r['name'].startswith('pip-') or r['name'].endswith('_receipt')];tasks.append({'name':'J61_receipt','command':[str(root.parent/'island-validation/phase4g1-env/bin/python'),'scripts/validate_phase_4j61_receipt.py','docs/evidence/phase_4j61','--repo','.']})
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
def one(t):
 cmd=t['command'];cmd=cmd if isinstance(cmd,list) else [str(root.parent/'island-validation/phase4g1-env/bin/python'),*cmd.split()]
 with (out/(t['name']+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
 return {'name':t['name'],'command':cmd,'returncode':r.returncode,'log':t['name']+'.log'}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,tasks))
(out/'independent-checks.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows)
