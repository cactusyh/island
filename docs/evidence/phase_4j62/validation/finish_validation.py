"""Collect unique suite coverage and final independent checks; no QM invoked."""
import json,os,shutil,subprocess,sys
from pathlib import Path
from collections import Counter
from xml.etree import ElementTree as ET
root=Path.cwd();out=Path('/tmp/island-j62-validation');e=root/'docs/evidence/phase_4j62';v=e/'validation';v.mkdir(exist_ok=True)
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
checks=json.loads((out/'independent-checks.json').read_text());assert all(c['returncode']==0 for c in checks)
commands=[('Ruff',['ruff','check','src','tests','scripts']),('compileall',[python,'-m','compileall','-q','src','scripts','tests']),('source-collection',[python,'-m','pytest','--collect-only','-q']),('diff-check',['git','diff','--check']),('cached-diff-check',['git','diff','--cached','--check'])]
for name,cmd in commands:
 log=out/(name+'-final.log')
 with log.open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 checks.append({'name':name,'command':cmd,'returncode':r.returncode,'log':log.name})
 if r.returncode:raise RuntimeError(f'{name} failed; preserve original log')
suite=json.loads((out/'suite-results.json').read_text());assert len(suite)==8 and all(s['returncode']==0 for s in suite)
analysis=json.loads((out/'analysis-results.json').read_text());assert len(analysis)==5 and all(a['returncode']==0 for a in analysis)
checks.extend(analysis)
nodes=[];counts=Counter();junit=[]
for name in [*[f'full-{i}.xml' for i in range(8)],'focused-final.xml']:
 p=out/name;cases=list(ET.parse(p).iter('testcase'));junit.append('validation/'+name)
 for c in cases:
  nodes.append(c.attrib['classname'].replace('.','/')+'.py::'+c.attrib['name']);counts['tests']+=1
  for kind in ('failure','error','skipped'):counts[kind]+=int(c.find(kind)is not None)
 assert not any(c.find('failure')is not None or c.find('error')is not None for c in cases)
collected=[x for x in (out/'source-collection-final.log').read_text().splitlines() if x.startswith('tests/') and '::'in x]
assert Counter(nodes)==Counter(collected) and len(nodes)==len(set(nodes))
(e/'collection.json').write_text(json.dumps({'nodes':collected},indent=2)+'\n')
(v/'junit-files.json').write_text(json.dumps(junit,indent=2)+'\n')
(v/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
(v/'suite-accounting.json').write_text(json.dumps({'software_counts':dict(counts),'original_base_nodes':5001,'source_collected_nodes':len(collected),'unique_node_coverage':True,'original_full_suite_failures':0,'full_suite_retries':[],'early_focused_and_inventory_failures':'early-attempts.json','negative_control_subprocess_nonzero_results_are_not_pytest_failures':True},indent=2)+'\n')
for p in out.iterdir():
 if p.is_file() and (p.suffix in ('.log','.xml') or p.name in ('suite-results.json','independent-checks.json','analysis-results.json','collection.json','software-priority.json','start_suite.py','failure_capture.py','run_checks.py','analyze_when_complete.py','finish_validation.py')):
  shutil.copyfile(p,v/p.name)
if (out/'subprocess-failures').exists():shutil.copytree(out/'subprocess-failures',v/'subprocess-failures',dirs_exist_ok=True)
print(json.dumps({'counts':dict(counts),'nodes':len(nodes),'checks':len(checks)},sort_keys=True))
