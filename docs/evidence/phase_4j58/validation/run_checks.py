import concurrent.futures,json,os,subprocess,sys
from pathlib import Path
root=Path.cwd();out=Path('/tmp/island-j58-validation');py=str(root.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':'.:src','OPENMM_CPU_THREADS':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
commands=[('ruff',['ruff','check','src','tests','scripts']),('compileall',[py,'-m','compileall','-q','src','tests','scripts'])]
for name,path in [('phase4e2','island-validation/phase4e2-env'),('phase4g1','island-validation/phase4g1-env'),('phase4j32','island-validation/phase4j32-analysis-env'),('phase4j36','island-validation/phase4j36-native-env'),('md_rg','miniforge3/envs/md_rg')]:
 commands.append(('pip-'+name,[str(root.parent/path/'bin/python'),'-m','pip','check']))
for phase in range(48,54):
 commands.append((f'J{phase}_receipt',[py,f'scripts/validate_phase_4j{phase}_receipt.py',f'docs/evidence/phase_4j{phase}','--repo','.']))
commands.append(('J55_receipt',[py,'scripts/validate_phase_4j55_receipt.py','docs/evidence/phase_4j55','--repo','.']))
commands.append(('J54_receipt',[py,'scripts/validate_phase_4j54_receipt.py','docs/evidence/phase_4j54','--repo','.', '--distribution',str(root.parent/'island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files'), '--archive','island_ffref/pcff/PDMS_benchmark-main.zip','--references','island_ffref']))
commands.append(('J56_receipt',[py,'scripts/validate_phase_4j56_receipt.py','docs/evidence/phase_4j56','--repo','.', '--distribution',str(root.parent/'island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files')]))
commands.append(('J57_receipt',[py,'scripts/validate_phase_4j57_receipt.py','docs/evidence/phase_4j57','--repo','.', '--distribution',str(root.parent/'island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files')]))
commands.append(('diff',['git','diff','--check']))
commands.append(('historical_hashes',[py,'-c',"from scripts.validate_phase_4j47_receipt import history,pin; PREREQUISITE='f471e06894275b5b6485461d2bc1dcd594f7dfa7'; from pathlib import Path; h=history(Path.cwd(),PREREQUISITE); assert all(pin(Path(n))==v for n,v in h.items()); print(len(h),'historical files byte-identical')"]))
def run(item):
 name,cmd=item
 with (out/(name+'.log')).open('w') as f:r=subprocess.run(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
 return {'name':name,'command':cmd,'returncode':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(run,commands))
(out/'checks-preliminary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results));sys.exit(any(r['returncode'] for r in results))
