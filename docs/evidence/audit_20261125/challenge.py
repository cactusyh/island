import json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.assignment_core import normalize, PolymerSpec
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.emc_profiles import install
CASES={
'PE':'*CC*','PS':'*CC(c1ccccc1)*','PEG':'*CCO*','PA6':'*NCCCCCC(=O)*',
'PET':'*OCCOC(=O)c1ccc(C(=O)*)cc1',
'PU':'*OCCOC(=O)NCCCCCCNC(=O)*',
'PAN':'*CC(C#N)*','PVC':'*CC(Cl)*','PTFE':'*C(F)(F)C(F)(F)*',
'PDMS':'*O[Si](C)(C)*','PPS':'*Sc1ccc(*)cc1',
'PEEK':'*Oc1ccc(Oc2ccc(C(=O)c3ccc(*)cc3)cc2)cc1',
'phosphate':'*OCCOP(=O)(O)*','sulfonated_PS':'*CC(c1ccc(S(=O)(=O)O)cc1)*'}
ROOT=Path('/tmp/island-audit-20261125')
def probe(case):
 name,smiles,field,mode=case;t=time.monotonic();r=dict(name=name,psmiles=smiles,field=field,mode=mode,dp=3)
 try:
  structure=smiles if mode=='repeat' else molecule(normalize(PolymerSpec((('A',smiles),),3)))
  options={'dp':3} if mode=='repeat' else {}
  a=assign_polymer_parameters(structure,field,field_dir=ROOT/(field+'-field'),**options)
  r.update(assigned=True,hashes=a.field_sha256)
 except Exception as e:r.update(assigned=False,error_type=type(e).__name__,diagnostic=str(e))
 r['seconds']=time.monotonic()-t
 return r
if __name__=='__main__':
 for f,p in [('pcff','pcff-carboxamide-v1'),('opls','opls-polymer-v1')]:
  if not (ROOT/(f+'-field')).exists():install(p,ROOT/(f+'-field'))
 jobs=[(n,s,f,m) for n,s in CASES.items() for f in ('pcff','opls') for m in ('repeat','graph')]
 with ProcessPoolExecutor(3) as pool,(ROOT/'challenge.jsonl').open('w') as out:
  for r in pool.map(probe,jobs):
   out.write(json.dumps(r)+'\n');out.flush();print(r['name'],r['field'],r['mode'],r['assigned'],flush=True)
