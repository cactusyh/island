import json,os,subprocess,sys
from pathlib import Path
repo=Path.cwd();out=Path('/tmp/island-j65-validation');py=repo.parent/'island-validation/phase4g1-env/bin/python';native=repo.parent/'island-validation/phase4j36-native-env/bin/python';work=repo.parent/'island-validation/phase4j65-qm'
env={**os.environ,'PYTHONPATH':'/tmp/island-j65-validation:.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','J62_DIAGNOSTICS':str(out/'subprocess-failures')};records=[]
def run(name,cmd,logname=None):
 logname=logname or name+'.log'
 with (out/logname).open('w') as log:r=subprocess.run([str(x) for x in cmd],env=env,stdout=log,stderr=subprocess.STDOUT)
 records.append({'name':name,'command':[str(x) for x in cmd],'cwd':str(repo),'returncode':r.returncode,'log':logname});(out/'final-check-results.json').write_text(json.dumps(records,indent=2)+'\n');print(name,r.returncode,flush=True)
 if r.returncode:sys.exit(r.returncode)
run('native-relocation-replay',[py,'scripts/replay_phase_4j65.py','replay'])
run('focused-final',[py,'-m','pytest','-p','failure_capture','-q','tests/test_phase_4j65.py','--basetemp',out/'pytest-j65-terminal','--junitxml',out/'focused-final.xml'])
run('final-focused-collection',[py,'-m','pytest','--collect-only','-q','tests/test_phase_4j65.py'],'focused-collection.txt')
nodes=[x for x in (out/'focused-collection.txt').read_text().splitlines() if x.startswith('tests/') and '::' in x];(out/'final-focused-collection.json').write_text(json.dumps({'nodes':nodes},indent=2)+'\n')
run('ruff-final',['ruff','check','src','tests','scripts'])
run('compileall-final',[py,'-m','compileall','-q','src','tests','scripts'])
before=(work/'attempts.json').read_bytes()
run('verified-resume-controls',[py,'scripts/phase_4j65.py','controls','--work',work,'--python',native])
assert (work/'attempts.json').read_bytes()==before,'Verified resume changed attempt ledger'
run('diff-final',['git','diff','--check'])
