"""Construct and validate an explicit safe recovery record; no dynamics run."""
import argparse
import json
from pathlib import Path

from island.analysis.periodic_preflight import PeriodicPreflightPolicy, preflight_box
from island.exceptions import AnalysisError


def recovery(*, original_graph_identity, original_parameter_identity, charge_origin, original_box_nm, recovered_box_nm, cutoff_nm, target_density_g_cm3, actual_density_g_cm3):
    d=preflight_box(initial_box_nm=recovered_box_nm, predicted_min_box_nm=recovered_box_nm, cutoff_nm=cutoff_nm, policy=PeriodicPreflightPolicy())
    if d['status']!='ok': raise AnalysisError('Recovery box fails periodic preflight')
    return {'schema':'island.j40-c4-safe-recovery.v1','original_run_preserved':True,'transformation':'recorded_box_repack_only','original_graph_identity':original_graph_identity,'recovered_graph_identity':original_graph_identity,'parameter_identity':original_parameter_identity,'charge_origin':charge_origin,'original_box_nm':list(original_box_nm),'recovered_box_nm':list(recovered_box_nm),'target_density_g_cm3':target_density_g_cm3,'actual_initial_density_g_cm3':actual_density_g_cm3,'cutoff_nm':cutoff_nm,'switching_distance_nm':cutoff_nm,'neighbor_skin_nm':0.0,'preflight':d,'dynamics_executed':False,'scientific_acceptance':False,'blocked_reason':'pressure_controls_and_full_recovery_not_completed'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(recovery(original_graph_identity='declared-source-graph',original_parameter_identity='declared-parameter-source',charge_origin='historical-control',original_box_nm=(1,1,1),recovered_box_nm=(3,3,3),cutoff_nm=.5,target_density_g_cm3=1.0,actual_density_g_cm3=1.0),indent=2)+'\n');return 1
if __name__=='__main__': raise SystemExit(main())
