import json,tempfile
from pathlib import Path
from rdkit import Chem
from island.workflows.recipe_config import PolymerRecipe
from island.workflows.recipe import _build_chain
from island.forcefields.assignment import assign_polymer_parameters
cases=[('PP-explicit-H-caps','*C[C@H](C)*',['*[H]','*[H]'],None),('PP-isotactic-methyl-caps','*C[C@H](C)*',['*C','*C'],'isotactic'),('PDMS-TMS-caps','*O[Si](C)(C)*',['*[Si](C)(C)C','*C'],None)]
rows=[]
for name,smi,ends,tacticity in cases:
 r=dict(name=name,psmiles=smi,end_groups=ends,tacticity=tacticity)
 p={'repeat_units':{'A':smi},'dp':3,'end_groups':ends}
 if tacticity:p['tacticity']=tacticity
 try:
  recipe=PolymerRecipe.from_dict(dict(architecture='homopolymer',polymer=p,forcefield={'name':'gaff2','charge':'gasteiger'}))
  with tempfile.TemporaryDirectory() as d:m=_build_chain(recipe,Path(d))
  r['graph']='passed';r['smiles_final']=Chem.MolToSmiles(Chem.RemoveHs(m))
  if name.startswith('PDMS'):
   for field in ('pcff','opls'):
    try:a=assign_polymer_parameters(m,field,field_dir=Path('/tmp/island-audit-20261125')/(field+'-field'));r[field]={'assigned':True}
    except Exception as e:r[field]={'assigned':False,'diagnostic':str(e)}
 except Exception as e:r.update(graph='failed',error=type(e).__name__,diagnostic=str(e))
 rows.append(r)
Path('/tmp/island-audit-20261125/caps-probes.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
