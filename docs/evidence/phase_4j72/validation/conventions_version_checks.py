import json,os,subprocess
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j72-validation');py=str(root.parent/'island-validation/phase4g1-env/bin/python')
files=['src/island/pcff_qm_pilot/source_adjudication.py','scripts/adjudicate_phase_4j72.py','scripts/check_phase_4j72_sources.py','scripts/validate_phase_4j72_receipt.py','tests/test_phase_4j72.py']
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
results=[]
for name,cmd in [('ruff-conventions-final-version',['ruff','check','src','scripts','tests']),('format-conventions-final-version',['ruff','format','--check',*files]),('compileall-conventions-final-version',[py,'-m','compileall','-q','src','scripts','tests']),('diff-conventions-final-version',['git','diff','--cached','--check'])]:
 with (out/(name+'.log')).open('w') as log:
  p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
 results.append({'name':name,'command':cmd,'cwd':str(root),'returncode':p.returncode,'log':name+'.log'})
 (out/'conventions-final-version-checks.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,p.returncode,flush=True)
