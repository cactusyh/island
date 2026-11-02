import json,os,shutil,subprocess,time
from pathlib import Path
from hashlib import sha256
from xml.etree import ElementTree as ET
root=Path.cwd();out=Path('/tmp/island-j71-validation');e=root/'docs/evidence/phase_4j71';v=e/'validation'
deadline=time.monotonic()+1800
while time.monotonic()<deadline:
 p=out/'suite-results.json'
 if p.exists():
  results=json.loads(p.read_text())
  if len(results)==2:break
 time.sleep(5)
else:raise RuntimeError('full suite did not finish in bounded wait')
assert all(r['returncode']==0 for r in results),results
for name in ('full-0.xml','full-1.xml','full-0.log','full-1.log','suite-results.json'):
 shutil.copyfile(out/name,v/name)
nodes=[];skips=0
for name in ('full-0.xml','full-1.xml','j71-final.xml'):
 for c in ET.parse(v/name).iter('testcase'):
  assert c.find('failure') is None and c.find('error') is None
  skips+=c.find('skipped') is not None
  nodes.append(c.attrib['classname'].replace('.','/')+'.py::'+c.attrib['name'])
expected=json.load(open(e/'collection.json'))['nodes']
assert len(nodes)==len(set(nodes)) and sorted(nodes)==sorted(expected)
counts={'unique_tests':len(nodes),'passed':len(nodes)-skips,'skipped':skips,'failure':0,'error':0,'historical_base_nodes':5567,'J71_nodes':93,'focused_overlap_coverage_added':0}
(v/'terminal-test-counts.json').write_text(json.dumps(counts,indent=2)+'\n')
doc=root/'docs/phase_4j71.md';s=doc.read_text()
marker='The supported full suite, J54–J70 receipt regressions,'
s=s.replace(marker,f"The unique terminal total is **{len(nodes):,} tests: {len(nodes)-skips:,} passed, {skips:,} skipped, zero failures/errors** (5,567 historical-suite nodes plus93 J71 nodes). Historical byte checks preserve8,722 files. Focused reruns add zero unique coverage.\n\n"+marker)
doc.write_text(s)
files=['src/island/pcff_qm_pilot/literature_intake.py','scripts/extract_phase_4j71.py','scripts/read_phase_4j71_sources.py','scripts/replay_phase_4j71.py','scripts/validate_phase_4j71_receipt.py','tests/test_phase_4j71.py','docs/phase_4j71.md']
(v/'implementation-files.json').write_text(json.dumps({n:sha256((root/n).read_bytes()).hexdigest() for n in files},indent=2)+'\n')
shutil.copyfile(__file__,v/'finish_evidence.py')
subprocess.run(['git','add',*files,'docs/evidence/phase_4j71'],check=True)
with (v/'diff-seal-final.log').open('w') as f:
 r=subprocess.run(['git','diff','--cached','--check'],stdout=f,stderr=subprocess.STDOUT,check=False)
assert r.returncode==0
checks=json.load(open(v/'checks.json'))
checks.append({'name':'diff-seal-final','command':['git','diff','--cached','--check'],'cwd':str(root),'returncode':r.returncode,'log':'diff-seal-final.log'})
(v/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
env={**os.environ,'PYTHONPATH':'../island-validation/phase4j66-solvers-only:.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
python=str(root.parent/'island-validation/phase4g1-env/bin/python')
command=[python,'scripts/validate_phase_4j71_receipt.py','docs/evidence/phase_4j71']
with (out/'seal.log').open('w') as f:
 r=subprocess.run(command+['--seal'],env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
assert r.returncode==0,(out/'seal.log').read_text()
with (out/'strict-receipt.log').open('w') as f:
 r=subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
assert r.returncode==0,(out/'strict-receipt.log').read_text()
print(counts,flush=True)
print((out/'strict-receipt.log').read_text(),flush=True)
