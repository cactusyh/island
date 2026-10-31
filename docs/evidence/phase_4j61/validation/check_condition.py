import sys,json
import numpy as np
from pathlib import Path
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.stretch_baseline import Configuration,source_baseline,stretching_projection
from island.pcff_qm_pilot.stretch_fit import Problem
m=load_model(Path.cwd());old=Path('docs/evidence/phase_4j60')
x=json.loads((old/'plan/reference-plan.json').read_text())['geometries']['center']
p=np.r_[source_baseline(m),np.linspace(-2.4,3.1,20)]
r=Configuration(m,x).evaluate(p)
s={'coordinates':x,'energy':r['energy'],'forces':r['forces'].tolist(),'graph_identity':m['graph_identity'],'atom_order':m['atom_order']}
q=stretching_projection(x,[[1,2],[2,3],[1,4],[1,5],[1,6],[3,7],[3,8],[3,9]])['Q']
b=Problem(m,[s]*4,s,q,np.eye(8));b.curvature=b.mm_curvature(p)[0]
try:
 print(b.fit(p)['scaled_condition'])
except Exception as e: print(str(e))
