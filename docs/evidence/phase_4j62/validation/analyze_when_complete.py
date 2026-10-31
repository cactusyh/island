"""Wait for the one bounded campaign, then reconstruct/analyze without new QM."""
import json,os,subprocess,sys,time
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j62-validation');e=root/'docs/evidence/phase_4j62';work=root.parent/'island-validation/phase4j62-qm';python=str(root.parent/'island-validation/phase4g1-env/bin/python');start=time.monotonic()
while not (e/'execution-outcome.json').exists():
 if time.monotonic()-start>14400:raise TimeoutError('No campaign outcome within four-hour analysis wait; inspect retained ledger, do not relaunch QM')
 time.sleep(15)
status=json.loads((e/'execution-outcome.json').read_text());require_success=not status['failures']
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
commands=[('retain',[python,'scripts/replay_phase_4j62.py','retain','--work',str(work)]),('analysis',[python,'scripts/analyze_phase_4j62.py','--work',str(work)]),('finalize',[python,'scripts/finalize_phase_4j62.py','--work',str(work)]),('replay',[python,'scripts/replay_phase_4j62.py','replay']),('focused-final',[python,'-m','pytest','-q','tests/test_phase_4j62.py','--junitxml',str(out/'focused-final.xml')])]
results=[]
for name,cmd in commands:
 if name!='retain' and not require_success:break
 with (out/(name+'-attempt1.log')).open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 results.append({'name':name,'command':cmd,'returncode':r.returncode});(out/'analysis-results.json').write_text(json.dumps(results,indent=2)+'\n');print(name,r.returncode,flush=True)
 if r.returncode:sys.exit(r.returncode)
sys.exit(not require_success)
