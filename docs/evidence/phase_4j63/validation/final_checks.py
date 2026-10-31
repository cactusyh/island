"""Run the terminal J63 checks; retain each attempt and its process status."""
import json,os,subprocess,sys
from pathlib import Path
repo=Path.cwd();out=Path('/tmp/island-j63-validation');python=repo.parent/'island-validation/phase4g1-env/bin/python'
env={**os.environ,'PYTHONPATH':str(out)+':.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1','J62_DIAGNOSTICS':str(out/'subprocess-failures')}
commands=[('focused-final',[str(python),'-m','pytest','-p','failure_capture','-q','tests/test_phase_4j63.py','--basetemp',str(out/'pytest-focused-final'),'--junitxml',str(out/'focused-final.xml')]),('ruff-final',['ruff','check','src','tests','scripts']),('compileall-final',[str(python),'-m','compileall','-q','src','scripts','tests']),('diff-final',['git','diff','--check'])]
results=[]
for name,command in commands:
 with (out/(name+'.log')).open('w') as log: p=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
 results.append({'name':name,'command':command,'returncode':p.returncode,'log':name+'.log'});(out/'final-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(name,p.returncode,flush=True)
 if p.returncode:sys.exit(p.returncode)
command=[str(python),'-m','pytest','--collect-only','-q']
p=subprocess.run(command,env=env,capture_output=True,text=True);(out/'source-collection-final.txt').write_text(p.stdout+p.stderr)
if p.returncode:sys.exit(p.returncode)
nodes=[s for s in p.stdout.splitlines() if s.startswith('tests/') and '::' in s]
(repo/'docs/evidence/phase_4j63/collection.json').write_text(json.dumps({'nodes':nodes},indent=2)+'\n');print('source collection',len(nodes),flush=True)
