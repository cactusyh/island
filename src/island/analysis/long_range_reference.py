"""Provenance and fail-closed closure records for PME/Ewald comparisons."""
import hashlib
import subprocess
from math import isfinite
from pathlib import Path

from island.exceptions import AnalysisError


def executable_manifest(path: str, *, openmm_version: str, precision: str, platform: str) -> dict:
 p=Path(path)
 if not p.is_file(): return {'status':'asset_unavailable','path':str(p),'reason':'executable missing'}
 digest=hashlib.sha256(p.read_bytes()).hexdigest()
 try: version=subprocess.run([str(p),'-h'],capture_output=True,text=True,check=False).stdout.splitlines()[:2]
 except OSError as e: version=[str(e)]
 return {'status':'available','executable':p.name,'path':str(p),'sha256':digest,'version_lines':version,'kspace_capability':'Ewald/PPPM capability must be confirmed by executed control','openmm_version':openmm_version,'precision':precision,'platform':platform}

def pairwise_coulomb_control(q1, q2, distance_nm, dielectric=1.0):
 if distance_nm<=0 or dielectric<=0 or not all(isfinite(float(x)) for x in (q1,q2,distance_nm,dielectric)): raise AnalysisError('Invalid Coulomb control')
 factor=138.93545764438198
 energy=factor*q1*q2/(dielectric*distance_nm)
 return {'control':'pairwise_coulomb','energy_kj_mol':energy,'force_norm_kj_mol_nm':abs(energy/distance_nm),'dielectric':dielectric,'conversion_factor':factor}

def closure_record(*, controls, polymer_cells, executable_manifest, tolerance):
 if tolerance<=0: raise AnalysisError('Invalid closure tolerance')
 polymer_ok=all(c.get('status')=='passed' for c in polymer_cells)
 controls_ok=all(c.get('status')=='passed' for c in controls)
 return {'schema':'island.long-range-closure.v1','controls':controls,'polymer_cells':polymer_cells,'executable':executable_manifest,'tolerance':tolerance,'toy_controls_passed':controls_ok,'polymer_closure_passed':polymer_ok,'status':'validated' if polymer_ok and controls_ok else 'ewald_unvalidated','scientific_acceptance':False if not polymer_ok else controls_ok}
