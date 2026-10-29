import json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
from island.forcefields.selection.inventory import digest
from scripts.validate_phase_4j49_receipt import validate
from scripts.validate_phase_4j47_receipt import pin
root=Path('docs/evidence/phase_4j49').resolve();repo=Path.cwd();expected=pin(root/'receipt.json')
checks={}
def rejected(name,call):
 try:call()
 except ValueError as error:checks[name]={'rejected':True,'reason':str(error)}
 else:raise AssertionError('Accepted tamper: '+name)
def resign(path,payload):
 payload.pop('identity',None);payload['identity']=digest(payload)
 path.write_text(json.dumps(payload,sort_keys=True,indent=2)+'\n')
rejected('independent_receipt_hash',lambda:validate(root,repo=repo,expected_sha256='0'*64))
with tempfile.TemporaryDirectory(prefix='island-j49-receipt-relocated-') as tmp:
 moved=Path(tmp)/'evidence';shutil.copytree(root,moved)
 proc=subprocess.run([sys.executable,'scripts/validate_phase_4j49_receipt.py',str(moved),'--repo',str(repo),'--sha256',expected],capture_output=True,text=True,check=True,env={**os.environ,'PYTHONHASHSEED':'31231'})
 checks['fresh_process_receipt_relocation']={'passed':True,'process_hash_seed':31231,'result':json.loads(proc.stdout)}
 receipt=json.loads((moved/'receipt.json').read_text());receipt['author_profile_MD_executed']=True
 resign(moved/'receipt.json',receipt)
 rejected('resigned_author_MD_promotion',lambda:validate(moved,repo=repo))
 shutil.copyfile(root/'receipt.json',moved/'receipt.json')
 name='audit/components-dp0.json';p=json.loads((moved/name).read_text())
 p['exact-erfc-ewald']['finite_difference']['raw_records'][0]['components']['coulomb']['positive_energy_kj_mol']+=1
 (moved/name).write_text(json.dumps(p,sort_keys=True,indent=2)+'\n')
 rejected('raw_FD_artifact_mutation',lambda:validate(moved,repo=repo))
 manifest=json.loads((moved/'manifest.json').read_text());manifest['artifacts'][name]=pin(moved/name);resign(moved/'manifest.json',manifest)
 receipt=json.loads((moved/'receipt.json').read_text());receipt['manifest_sha256']=pin(moved/'manifest.json');resign(moved/'receipt.json',receipt)
 rejected('resigned_raw_FD_forgery',lambda:validate(moved,repo=repo))
record={'tested_receipt_sha256':expected,'checks':checks,'all_passed':True}
Path('/tmp/island-j49-validation/receipt-controls.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps(record,sort_keys=True))
