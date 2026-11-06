"""Run independent J86 offline checks under separate inherited execution fences."""
import os
import subprocess
import sys
from pathlib import Path

repo=Path.cwd()
py=str(repo.parent/'island-validation/phase4g1-env/bin/python')
checks={
 'ruff':['-c','from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["ruff","check","src/island/forcefields/bundle_adapters/scoped_qualification.py","scripts/phase_4j86_qualification.py","scripts/validate_phase_4j86_receipt.py","tests/test_phase_4j86.py"],check=True); subprocess.run(["ruff","format","--check","src/island/forcefields/bundle_adapters/scoped_qualification.py","scripts/phase_4j86_qualification.py","scripts/validate_phase_4j86_receipt.py","tests/test_phase_4j86.py"],check=True)'],
 'compileall':['-m','compileall','-q','src/island','scripts','tests'],
 'pip_supported':['-c','from scripts.phase_4j76_source_resolution import scope; scope(); import runpy; import sys; sys.argv=["pip","check"]; runpy.run_module("pip",run_name="__main__")'],
 'diff':['-c','from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["git","diff","--check"],check=True)'],
}
for phase in ('81','82','83','84','85','80','79','78','76','72_c1'):
 label='C1_receipt' if phase=='72_c1' else 'J'+phase+'_receipt'
 checks[label]=['-m','scripts.validate_phase_4j'+phase+'_receipt']
for phase,module in [('81','ingest'),('82','coverage'),('83','registry'),('84','bundles'),('85','qualification')]:
 checks['J'+phase+'_replay']=['-m','scripts.phase_4j'+phase+'_'+module,'replay']
checks['offline_replay']=['-m','scripts.phase_4j86_qualification','replay']
checks['relocated_replay']=['/tmp/island_j86_relocate.py']
name=sys.argv[1]
with open('/tmp/island_j86_final_'+name+'.log','w') as log:
 run=subprocess.run([py,'-m','scripts.run_phase_4j72_c1','--output','/tmp/island_j86_final_'+name,'--',*checks[name]],env={**os.environ,'PYTHONPATH':'.:src'},stdout=log,stderr=subprocess.STDOUT,check=False)
print(name,run.returncode)
raise SystemExit(run.returncode)
