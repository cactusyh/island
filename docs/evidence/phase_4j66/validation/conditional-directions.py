from pathlib import Path
import numpy as np
from island.pcff_completion import signed,write
from island.pcff_qm_pilot.hmds_constrained import physical_problem
from island.pcff_qm_pilot.hmds_precision import read
r=Path('docs/evidence/phase_4j66')
with np.load(r/'physical-arrays.npz') as f:data={k:f[k] for k in f.files}
names=data['names'].tolist()+['energy_offset'];data['names']=data['names'].tolist();p=physical_problem(data);a=p['w'];norms=np.linalg.norm(a,axis=0)
audit=read(r/'block-identifiability-diagnosis.json');rows=[]
for name,_ in sorted(audit['J66_bound_parameter_recurrence'].items(),key=lambda item:-item[1])[:5]:
 j=names.index(name);others=[k for k in range(59) if k!=j];v=np.zeros(59);v[j]=1;v[others]=np.linalg.lstsq(a[:,others]/norms[others],-a[:,j],rcond=1e-10)[0]/norms[others]
 scaled=v*norms;normalization=np.linalg.norm(scaled);v/=normalization;scaled/=normalization
 effects=[]
 for block in sorted({x['block'] for x in read(Path('docs/evidence/phase_4j65/fitting-dataset.json'))['observations']}):
  ids=[i for i,x in enumerate(read(Path('docs/evidence/phase_4j65/fitting-dataset.json'))['observations']) if x['block']==block]
  e=(p['e']@v)[ids];f=(p['f']@v).reshape(22,27,3)[ids];effects.append(dict(block=block,energy_RMS=float(np.sqrt(np.mean(e**2))),force_RMS=float(np.sqrt(np.mean(f**2))),maximum_local_vector=float(np.linalg.norm(f,axis=2).max())))
 top=np.argsort(abs(scaled))[-10:][::-1]
 rows.append(dict(fixed_direction_parameter=name,scaled_direction=scaled.tolist(),physical_direction=v.tolist(),conditional_weighted_residual_norm=float(np.linalg.norm(a@v)),top=[dict(name=names[k],scaled=float(scaled[k])) for k in top],observations=effects))
write(r/'bound-combination-sensitivities.json',signed(dict(directions=rows,normalization='For each recurring bound coefficient, minimize the actual weighted prediction change with its physical derivative fixed1 and other coefficients free, then normalize in column-scaled coordinates. Descriptive local linear dependence; not new fitting or parameter constraints.',new_QM=0)))
