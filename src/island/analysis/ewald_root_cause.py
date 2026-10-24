"""Executed root-cause isolation records for polymer PME/KSPACE comparisons."""
from island.exceptions import AnalysisError

CONTROL_ORDER=('point_charges','counterion_system','pairwise_control','small_fragment','linear_polymer','aromatic_polymer','packed_network')
VARIABLES=('charge_units','atom_mapping','exclusions','one_four_scaling','box_convention','cutoff','reciprocal_tolerance','precision','periodic_images','network_topology')

def isolate(cases):
 missing=[x for x in CONTROL_ORDER if x not in cases]
 first=next((x for x in CONTROL_ORDER if x in cases and cases[x].get('status')!='passed'),None)
 return {'schema':'island.ewald-root-cause.v1','control_order':CONTROL_ORDER,'cases':cases,'missing_controls':missing,'first_failing_system':first,'candidate_variables':VARIABLES,'status':'blocked' if missing or first else 'closed','root_cause':None if first else 'none demonstrated','ewald_unvalidated':bool(missing or first)}

def one_variable_experiment(*, baseline, variable, changed):
 if variable not in VARIABLES: raise AnalysisError('Unknown Ewald diagnosis variable')
 fields=('coordinates','box','atom_order','charges','charge_origin','exclusions','one_four_scaling','cutoff','switching','reciprocal_tolerance','units','precision')
 return {'variable':variable,'baseline_identity':baseline.get('identity'),'changed_identity':changed.get('identity'),'metadata_changed':[x for x in fields if baseline.get(x)!=changed.get(x)],'status':'blocked' if baseline.get('status')!='passed' or changed.get('status')!='passed' else 'passed','repair_justified':False}
