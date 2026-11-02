import json,os,shutil,subprocess,time
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j70-validation');python=str(root.parent/'island-validation/phase4g1-env/bin/python')
checks=[('ruff',[shutil.which('ruff'),'check','src','tests','scripts']),('compileall',[python,'-m','compileall','-q','src','tests','scripts']),('offline_replay',[python,'scripts/replay_phase_4j70.py']),('offline_deep_parent',[python,'scripts/replay_phase_4j70.py','--deep-parent']),('historical_bytes',[python,'-c',"from pathlib import Path;from scripts.validate_phase_4j47_receipt import history,pin;from scripts.intake_phase_4j70 import BASE;r=Path.cwd();h=history(r,BASE);changed=[n for n,v in h.items() if pin(r/n)!=v];assert not changed,changed;print({'base':BASE,'files_preserved':len(h),'changed':changed})"])]
results=[]
for name,command in checks:
 env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as f:p=subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
 results.append({'name':name,'command':command,'cwd':str(root),'returncode':p.returncode,'runtime_seconds':time.monotonic()-start,'log':name+'.log'})
 (out/'current-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(name,p.returncode,flush=True)
raise SystemExit(any(r['returncode'] for r in results))
