import json
from pathlib import Path
import numpy as np
from island.pcff_completion import signed,write,digest,pin
from island.pcff_qm.geometry import _adjust
from island.pcff_qm_pilot.campaign import selected
from island.pcff_qm_pilot.stretch_fit import leakage
root=Path.cwd();e=root/'docs/evidence/phase_4j62';read=lambda p:json.loads(p.read_text());plan=read(e/'plan/scientific-plan.json');x=np.array(plan['control']['center']['coordinates'])
system,_,_=selected(root,root.parent/'island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files',x)
rng=np.random.default_rng(4062);rows=[]
coords=[('bond',[1,4],.03),('bond',[1,2],.03),('bond',[2,3],.03),('bond',[3,7],.03),('angle',[4,1,2],.08),('angle',[4,1,5],.08),('angle',[1,2,3],.08),('angle',[2,3,7],.08)]
for k,t in enumerate([(1.31,.61),(-.55,1.45),(1.02,-1.17)]):
 y=x.copy();offsets=[]
 for kind,sites,extent in coords:
  z=float(rng.uniform(-extent,extent));offsets.append(z);_adjust(system,y,list(range(1,10)),{'kind':kind,'sites':sites},z)
 for sites,z in zip(([4,1,2,3],[1,2,3,7]),t):_adjust(system,y,list(range(1,10)),{'kind':'torsion','sites':sites},z)
 rows.append({'key':f'j62-fresh/validation-{k}','coordinates':y.tolist(),'geometry_identity':digest(y.tolist()),'local_offsets':offsets,'torsion_offsets':t})
history=[]
d=read(root/'docs/evidence/phase_4j60/joint-target-design.json');history.extend(d['selected_training']+d['held_out']);history.extend(read(root/'docs/evidence/phase_4j61/plan/fresh-configurations.json')['selections'])
report=signed({'schema':'island.j62-conditional-fresh-selection.v1','seed':4062,'selection_rules':coords,'selections':rows,'separation':leakage(history,rows),'plan_sha256':pin(e/'plan/scientific-plan.json'),'status':'reserved_only_no_QM_or_errors_inspected','scope':'New deformations, not independently sampled conformational basins; only a training/uncertainty-qualified frozen model activates execution. No model tuning from fresh values.'})
write(e/'plan/conditional-fresh-validation.json',report)
print(report['separation'])
