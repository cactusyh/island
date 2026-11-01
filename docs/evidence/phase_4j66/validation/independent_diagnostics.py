from pathlib import Path
import numpy as np
from island.pcff_completion import signed,write
from island.pcff_qm_pilot.hmds_precision import read
r=Path('docs/evidence/phase_4j66');u=read(r/'block-uncertainty.json');names=read(r/'reconstruction.json')['parameters'];values=np.array([row['values'] for row in u['attempts'] if row['values'] is not None]);bounds=read(r/'plan/scientific-plan.json')['parameter_specifications'];trial=read(r/'constrained-epigraph.json');v=np.array(trial['solution']['values']);w=read(r/'nominal-rho.json')
intervals=[]
for i,name in enumerate(names+['energy_offset']):
 q=np.quantile(values[:,i],[.025,.975]);intervals.append({'name':name,'minimum':float(values[:,i].min()),'maximum':float(values[:,i].max()),'lower':float(q[0]),'upper':float(q[1]),'absolute_sd':float(values[:,i].std(ddof=1))})
locals=trial['original_unit_verification']['local_vector_force_errors'];active=[]
for i,row in enumerate(locals):
 for j,error in enumerate(row):
  if 2-error<1e-5:active.append({'configuration':read(r/'plan/scientific-plan.json')['configuration_order'][i],'stable_atom_id':j+1,'local_error':error,'slack':2-error})
write(r/'all-trial-sensitivity-and-active-constraints.json',signed({'all_feasible_conditional_trials':len(values),'interval_scope':'Descriptive failed-and-surviving trial ranges, not qualified confidence intervals. Infeasible draw vectors are not silently turned into fitted coefficients. Exact64 outcomes separately retained.','absolute_intervals':intervals,'active_local_constraints':active,'force_RMS_slack':.5-trial['original_unit_verification']['force_RMS'],'parameter_bound_distances':[{'name':b['name'],'lower_distance':float(v[i]-b['lower']),'upper_distance':float(b['upper']-v[i]),'units':b['units']} for i,b in enumerate(bounds)],'no_active_bound_maximum_eligible_fraction_ignoring_status':sum(all(c=='solver_status' for c in row['failure_categories']) for row in u['attempts'])/64,'global_problem_infeasible':False,'new_QM':0}))
