"""Publish the first-failing real C7 polymer comparison and decision."""
import argparse
import json
from pathlib import Path


def main():
 p=argparse.ArgumentParser();p.add_argument('--cross-engine',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=json.loads(a.cross_engine.read_text())
 row={'system':'gaff2-pe','source_artifact':str(a.cross_engine),'final_graph_identity':d['final_graph_identity'],'configuration':d['configuration'],'openmm_components_kj_mol':d['openmm_components_kj_mol'],'lammps_components_kj_mol':d['lammps_components_kj_mol'],'component_errors_kj_mol':d['component_errors_kj_mol'],'energy_error_kj_mol':d['energy_error_kj_mol'],'force_max_error_kj_mol_angstrom':d['force_max_error_kj_mol_angstrom'],'checks':d['checks'],'status':'blocked','root_cause':None}
 out={'schema':'island.j40-c10-executed-decision.v1','ordered_levels':['point_charges','counterion_system','pairwise_control','small_fragment','gaff2_pe_linear','aromatic_resp_polymer','packed_network'],'first_failing_level':'gaff2-pe_linear','actual_comparisons':[row],'parameter_sweep':{'status':'no_repair_justified','variables':['charge_units','exclusions','one_four_scaling','cutoff','reciprocal_tolerance','ewald_parameter','pppm_grid','fft_precision','box_representation']},'decision':'unresolved_backend_or_convention_limitation','ewald_unvalidated':True,'pressure_unavailable':True,'scientific_acceptance':False}
 a.output.write_text(json.dumps(out,indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
