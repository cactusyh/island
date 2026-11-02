import json,os,subprocess,time
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j70-validation')
rows=json.loads(Path('/tmp/island-j69-validation/historical-and-environment-checks.json').read_text())
rows.append({'name':'J69_receipt','command':[str(root.parent/'island-validation/phase4g1-env/bin/python'),'scripts/validate_phase_4j69_receipt.py','docs/evidence/phase_4j69'],'cwd':str(root),'log':'J69_receipt.log'})
results=[]
for row in rows:
 env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
 if row['name']=='pip-analysis-solver-bridge':env['PYTHONPATH']=str(root.parent/'island-validation/phase4j66-solvers-only')+':.:src'
 if row['name'].startswith('pip-') and row['name']!='pip-analysis-solver-bridge':env.pop('PYTHONPATH',None)
 start=time.monotonic()
 with (out/row['log']).open('w') as f:
  p=subprocess.run(row['command'],cwd=row['cwd'],env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
 result={**row,'returncode':p.returncode,'runtime_seconds':time.monotonic()-start}
 results.append(result);(out/'historical-and-environment-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(row['name'],p.returncode,flush=True)
raise SystemExit(any(r['returncode'] for r in results))
