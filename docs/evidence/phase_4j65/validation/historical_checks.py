import json,os,subprocess,sys
from pathlib import Path
out=Path('/tmp/island-j65-validation');repo=Path.cwd();hist=repo.parent/'island-validation/phase4j65-historical-base'
checks=json.loads((out/'historical-check-commands.json').read_text()); results=[]
for c in checks:
 cwd=repo if c['name'].startswith('pip-') else hist
 with (out/(c['name']+'.log')).open('w') as f:
  r=subprocess.run(c['command'],cwd=cwd,env={**os.environ,'PYTHONPATH':'.:src','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'},stdout=f,stderr=subprocess.STDOUT)
 results.append({'name':c['name'],'command':c['command'],'cwd':str(cwd),'returncode':r.returncode})
 (out/'historical-check-results.json').write_text(json.dumps(results,indent=2)+'\n'); print(c['name'],r.returncode,flush=True)
sys.exit(any(c['returncode'] for c in results))
