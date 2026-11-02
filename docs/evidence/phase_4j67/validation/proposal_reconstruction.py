import json,time
from pathlib import Path
import numpy as np
from island.pcff_completion import require,signed,write
from island.pcff_qm_pilot.hmds_baseline import extended_parts
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
repo=Path.cwd();root=repo/'docs/evidence/phase_4j67';parent=repo/'docs/evidence/phase_4j66'
p=read(parent/'geometry-only-next-design.json');audit=read(parent/'block-identifiability-diagnosis.json')
directions=np.array([r['physical_direction'] for r in audit['nominal_weak']]+[r['weak'][-1]['physical_direction'] for r in audit['weak_draws']])
require(np.array_equal(directions,p['weak_directions']),'J67_proposed_weak_direction_mapping')
m=load_model(repo);rows=[]
for row in p['selected']:
 f=extended_parts(m,row['coordinates'])
 e=np.r_[f['energy_jacobian'],1]/(.1*np.sqrt(22))
 forces=np.c_[f['force_jacobian'].reshape(81,58),np.zeros(81)]/(.5*np.sqrt(81*22))
 measured=np.linalg.norm(np.vstack([e,forces])@directions.T,axis=0)
 error=float(np.max(abs(measured-row['target_weak_direction_norms'])))
 require(error<1e-12,'J67_proposal_sensitivity_reconstruction')
 rows.append({'id':row['id'],'geometry_identity':row['geometry_identity'],'maximum_absolute_sensitivity_error':error,'weak_direction_norms':measured.tolist()})
write(root/'proposal-sensitivity-reconstruction.json',signed({'rows':rows,'directions':directions.tolist(),'reconstructed_unix':time.time(),
 'timing':'Additional independent geometry-only sensitivity replay after the first native gradient completed. The exact selected coordinates, parent identities, graph leakage, requests, fit rules and64draws were frozen before any new gradient; no selections or rules changed.',
 'new_QM_used_for_sensitivity_reconstruction':False}))
print(rows)
