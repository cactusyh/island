import json,sys,platform,tempfile
from pathlib import Path
import rdkit
from rdkit import Chem
from island.workflows.recipe_config import PolymerRecipe
from island.workflows.recipe import _build_chain
from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.assignment import assign_polymer_parameters
cases=[('PP-unspecified','*CC(C)*',None),('PP-explicit','*C[C@H](C)*',None),('PP-isotactic','*CC(C)*','isotactic'),('PP-syndiotactic','*CC(C)*','syndiotactic'),('stereo-internal','*C[C@H](F)O*',None),('PI-EZ','*C/C=C\\C*',None),('PA6-default-ends','*NCCCCCC(=O)*',None)]
rows=[]
for name,smiles,tacticity in cases:
 r={'name':name,'psmiles':smiles,'tacticity':tacticity}
 config={'architecture':'homopolymer','polymer':{'repeat_units':{'A':smiles},'dp':3},'forcefield':{'name':'gaff2','charge':'gasteiger'}}
 if tacticity:config['polymer']['tacticity']=tacticity
 try:
  recipe=PolymerRecipe.from_dict(config);g=normalize(recipe.polymer_spec());m=molecule(g)
  r.update(blockers=list(g.blockers),radicals=[a.GetIdx() for a in m.GetAtoms() if a.GetNumRadicalElectrons()],smiles_final=Chem.MolToSmiles(Chem.RemoveHs(m)),aldehyde_sites=len(m.GetSubstructMatches(Chem.MolFromSmarts('[CX3H1](=O)'))))
  with tempfile.TemporaryDirectory() as out:_build_chain(recipe,Path(out))
  r['recipe_graph']='passed'
 except Exception as e:r.update(recipe_graph='failed',error_type=type(e).__name__,diagnostic=str(e))
 if not tacticity:
  try:
   a=assign_polymer_parameters(smiles,'gaff2',dp=3,charge='gasteiger',atom_order='radonpy')
   r['direct_assignment']='passed';r['direct_smiles']=Chem.MolToSmiles(Chem.RemoveHs(a.mol))
  except Exception as e:r['direct_assignment']='failed';r['direct_diagnostic']=str(e)
 rows.append(r)
result={'python':platform.python_version(),'rdkit':rdkit.__version__,'results':rows}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
