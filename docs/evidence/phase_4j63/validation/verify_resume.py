import json,os,subprocess
from pathlib import Path
from island.pcff_completion import pin,signed,write
repo=Path.cwd();work=repo.parent/'island-validation/phase4j63-qm';logs=Path('/tmp/island-j63-validation');py=repo.parent/'island-validation/phase4g1-env/bin/python';engine=repo.parent/'island-validation/phase4j36-native-env/bin/python'
before=pin(work/'ledger.json');rows=[];env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
for stage in ('optimize','deform'):
 command=[str(py),'scripts/phase_4j63.py',stage,'--work',str(work),'--python',str(engine)]
 with (logs/f'resume-{stage}.log').open('w') as log: result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
 rows.append({'name':'resume-'+stage,'command':command,'returncode':result.returncode,'log':f'resume-{stage}.log'})
 if result.returncode:raise RuntimeError(rows[-1])
after=pin(work/'ledger.json');assert before==after
write(repo/'docs/evidence/phase_4j63/validation/exact-resume.json',signed({'checks':rows,'ledger_before_sha256':before,'ledger_after_sha256':after,'new_gradient_calls':0,'all_exact_completed_outputs_verified':True}));print('Exact optimize/deform resume verified; zero native calls.')
