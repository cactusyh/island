import concurrent.futures,json,os,subprocess,time
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j71-validation');old=json.load(open('/tmp/island-j70-validation/historical-and-environment-checks.json'))
env={**os.environ,'PYTHONPATH':str(root.parent/'island-validation/phase4j66-solvers-only')+':.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
def run(item):
 item={k:item[k] for k in ('name','command','cwd')};name=item['name'];e=env.copy()
 if name.startswith('pip-'): e.pop('PYTHONPATH',None)
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as f:
  p=subprocess.run(item['command'],cwd=item['cwd'],env=e,stdout=f,stderr=subprocess.STDOUT,check=False)
 return {**item,'returncode':p.returncode,'log':name+'.log','runtime_seconds':time.monotonic()-start}
old.append({'name':'J70_receipt','command':[str(root.parent/'island-validation/phase4g1-env/bin/python'),'scripts/validate_phase_4j70_receipt.py','docs/evidence/phase_4j70'],'cwd':str(root)})
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for result in pool.map(run,old):
  results.append(result);(out/'historical-and-environment-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(result['name'],result['returncode'],flush=True)
