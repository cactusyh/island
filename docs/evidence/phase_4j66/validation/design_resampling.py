"""Geometry-only stress test of the bounded proposal; no reference values."""
from pathlib import Path
import numpy as np
from island.pcff_completion import signed,write
from island.pcff_qm_pilot.hmds_baseline import extended_parts
from island.pcff_qm_pilot.hmds_constrained import physical_problem
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
r=Path('docs/evidence/phase_4j66');proposal=read(r/'geometry-only-next-design.json');plan=read(Path('docs/evidence/phase_4j65/plan/scientific-plan.json'));model=load_model(Path.cwd())
with np.load(r/'physical-arrays.npz') as stored:data={k:stored[k] for k in stored.files}
data['names']=data['names'].tolist();u=read(r/'block-uncertainty.json')
parts=[extended_parts(model,x['coordinates']) for x in proposal['selected']]
rows=[]
for record in u['attempts']:
 p=physical_problem(data,record['indices']);n=p['n'];extra=[]
 for item in parts:extra.extend([np.r_[item['energy_jacobian'],1]/(.1*np.sqrt(n)), *np.c_[item['force_jacobian'].reshape(81,58),np.zeros(81)]/(.5*np.sqrt(81*n))])
 after=metrics(np.vstack([p['w'],extra]));rows.append({'replicate':record['replicate'],'old_scaled_condition':record['actual_weighted_design']['scaled_condition'],'after_geometry_only_condition':after['scaled_condition'],'after_rank':after['rank']})
write(r/'design-block-stress-check.json',signed({'rows':rows,'maximum_old_condition':max(x['old_scaled_condition'] for x in rows),'maximum_after_condition':max(x['after_geometry_only_condition'] for x in rows),'new_QM':0,'scope':'Append all4 new geometry feature blocks as guaranteed training information to each retained old draw. Not a new bootstrap distribution; not evidence of predictive accuracy, physical feasibility or qualified uncertainty. Future draws must also resample actual new observations.','protected_existing_errors':'Any future fitting protocol must keep all original22 errors/limits visible as a separate protected subset; adding new configurations cannot dilute their failed errors.','new_independently_sampled_conformational_basins':False}))
