import json,os,subprocess,time
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j72-validation/terminal-checks')
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':str(root.parent/'island-validation/phase4j66-solvers-only')+':.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
files=['src/island/pcff_qm_pilot/source_adjudication.py','scripts/adjudicate_phase_4j72.py','scripts/check_phase_4j72_sources.py','scripts/validate_phase_4j72_receipt.py','tests/test_phase_4j72.py']
commands=[
('ruff-full-final',['ruff','check','src','scripts','tests']),
('ruff-format-final',['ruff','format','--check',*files]),
('compileall-final',[python,'-m','compileall','-q','src','scripts','tests']),
('diff-check-final',['git','diff','--cached','--check']),
('offline-deep-final',[python,'scripts/adjudicate_phase_4j72.py','--replay','--deep-parent']),
('independent-PDF-tables',['/tmp/island-j71-validation/pdf-env/bin/python','scripts/check_phase_4j72_sources.py']),
('original-sources-final',[python,'scripts/read_phase_4j71_sources.py']),
('pip-isolated-pdf-final',['/tmp/island-j71-validation/pdf-env/bin/python','-m','pip','check']),
('historical-bytes-final',[python,'-c','from pathlib import Path; from scripts.validate_phase_4j47_receipt import history,pin; from scripts.adjudicate_phase_4j72 import BASE; r=Path.cwd(); h=history(r,BASE); assert all(pin(r/n)==v for n,v in h.items()); print({"files_preserved":len(h),"base":BASE})']),
]
results=[]
for name,cmd in commands:
 start=time.monotonic()
 e=env.copy()
 if name=='independent-PDF-tables':e.pop('PYTHONPATH',None)
 if name.startswith('pip-'):e.pop('PYTHONPATH',None)
 with (out/(name+'.log')).open('w') as log:
  p=subprocess.run(cmd,cwd=root,env=e,stdout=log,stderr=subprocess.STDOUT,check=False)
 results.append({'name':name,'command':cmd,'cwd':str(root),'returncode':p.returncode,'log':name+'.log','runtime_seconds':time.monotonic()-start})
 (out/'final-checks.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,p.returncode,flush=True)
