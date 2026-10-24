"""Bounded numerical characterization of the retained GAFF2-PE discrepancy."""
import argparse
import json
import math
from pathlib import Path


def main():
 p=argparse.ArgumentParser();p.add_argument('--comparison',type=Path,required=True);p.add_argument('--forces',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=json.loads(a.comparison.read_text());cfg=d['configuration'];charges=[];forces=[]
 if a.forces and a.forces.exists():
  for line in a.forces.read_text().splitlines():
   cols=line.split()
   if len(cols)>=9 and cols[0].isdigit(): charges.append(float(cols[2]));forces.append(math.sqrt(sum(float(x)**2 for x in cols[6:9])))
 atom_count=len(charges) or None; total=d['openmm_energy_kj_mol']; err=d['energy_error_kj_mol']
 out={'schema':'island.j40-c11-characterization.v1','system':'GAFF2 linear PE','graph_identity':d['final_graph_identity'],'atom_count':atom_count,'box_lengths_angstrom':cfg['box_lengths'],'total_energy_kj_mol':total,'coulomb_energy_kj_mol':d['openmm_components_kj_mol']['coulomb'],'absolute_energy_error_kj_mol':err,'relative_energy_error':abs(err/total),'force_rms_error_kj_mol_nm':None,'force_max_error_kj_mol_nm':d['force_max_error_kj_mol_angstrom']/10,'force_error_percentiles':None,'maximum_lammps_force_kj_mol_nm':max(forces) if forces else None,'charge_range':([min(charges),max(charges)] if charges else None),'cutoff_angstrom':cfg['nonbonded_cutoff'],'reciprocal_tolerance':cfg['pme_tolerance'],'sweep':{'status':'bounded_difference_not_converged','variables':['reciprocal_tolerance','cutoff','mesh','interpolation_order','ewald_parameter','fft_precision','box_size'],'repair_applied':False},'decision':'bounded_numerical_difference_candidate_but_cross_engine_scientific_gate_unvalidated'}
 a.output.write_text(json.dumps(out,indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
