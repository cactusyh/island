import json, os, subprocess, time
from pathlib import Path
repo=Path.cwd(); out=Path('/tmp/island-j72-validation/scope-corrected-checks'); out.mkdir(exist_ok=True)
py=str(repo.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
files=['src/island/pcff_qm_pilot/source_adjudication.py','scripts/adjudicate_phase_4j72.py','scripts/check_phase_4j72_sources.py','scripts/phase_4j72_execution_fence.py','scripts/phase_4j72_test_policy.py','scripts/validate_phase_4j72_receipt.py','tests/test_phase_4j72.py']
commands=[
 ('ruff-scope-corrected',['ruff','check','src','scripts','tests']),
 ('format-scope-corrected',['ruff','format','--check',*files]),
 ('compileall-scope-corrected',[py,'-m','compileall','-q','src','scripts','tests']),
 ('offline-deep-scope-corrected',[py,'scripts/adjudicate_phase_4j72.py','--replay','--deep-parent']),
 ('independent-PDF-scope-corrected',['/tmp/island-j71-validation/pdf-env/bin/python','scripts/check_phase_4j72_sources.py']),
 ('original-sources-scope-corrected',[py,'scripts/read_phase_4j71_sources.py']),
 ('pip-analysis-scope-corrected',[py,'-m','pip','check']),
 ('pip-PDF-scope-corrected',['/tmp/island-j71-validation/pdf-env/bin/python','-m','pip','check']),
 ('historical-bytes-scope-corrected',[py,'-c','from pathlib import Path; from scripts.validate_phase_4j47_receipt import history,pin; from scripts.adjudicate_phase_4j72 import BASE; r=Path.cwd(); h=history(r,BASE); assert all(pin(r/n)==v for n,v in h.items()); print({"files_preserved":len(h),"base":BASE})']),
 ('diff-unstaged-scope-corrected',['git','diff','--check']),
 ('diff-staged-scope-corrected',['git','diff','--cached','--check']),
]
results=[]
for name,cmd in commands:
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as log:
  p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
 results.append({'name':name,'command':cmd,'cwd':str(repo),'returncode':p.returncode,'log':name+'.log','runtime_seconds':time.monotonic()-start})
 (out/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,p.returncode,flush=True)
assert all(r['returncode']==0 for r in results), results
