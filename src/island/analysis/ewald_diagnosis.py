"""Immutable component-level OpenMM PME/KSPACE discrepancy diagnosis."""
import hashlib
import json

from island.exceptions import AnalysisError

FIELDS=('graph','coordinates','box','atom_order','masses','charges','charge_origin','exclusions','one_four_scaling','lj_mixing','cutoff','switching','reciprocal_tolerance','units','precision')
COMPONENTS=('total','coulomb','real_space_coulomb','reciprocal','lj','bonded','angle','torsion','improper','force_norm','virial','pressure_proxy')

def identity(payload):
 return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def compare_frame(*, openmm, lammps, tolerance):
 if tolerance<=0: raise AnalysisError('Invalid Ewald comparison tolerance')
 mismatches=[k for k in FIELDS if openmm.get(k)!=lammps.get(k)]
 residuals={k:float(openmm[k])-float(lammps[k]) for k in COMPONENTS if k in openmm and k in lammps}
 return {'schema':'island.ewald-frame-comparison.v1','input_identity':identity({k:openmm.get(k) for k in FIELDS}),'metadata_mismatches':mismatches,'residuals':residuals,'tolerance':tolerance,'status':'passed' if not mismatches and residuals and all(abs(v)<=tolerance for v in residuals.values()) else 'blocked','diagnostic':'identical component closure' if not mismatches and residuals and all(abs(v)<=tolerance for v in residuals.values()) else 'input mismatch or residual exceeds unchanged bound'}

def diagnose(cases):
 required={'point_charges','counterion','pairwise','polymer_fragment','linear_polymer','aromatic_polymer','network','box_scaling','cutoff_scaling','reciprocal_tolerance_scaling'}
 missing=sorted(required-set(cases))
 return {'schema':'island.ewald-diagnosis.v1','cases':cases,'missing_cases':missing,'candidate_causes':['charge_units','atom_mapping','exclusions_1_4','box_convention','cutoff_switching','reciprocal_tolerance','boundary_convention','precision','component_decomposition','network_periodic_images'],'status':'blocked' if missing or any(x.get('status')!='passed' for x in cases.values()) else 'closed','ewald_unvalidated':bool(missing or any(x.get('status')!='passed' for x in cases.values()))}
