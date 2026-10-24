"""Inventory and audit the executed KSPACE reference without fabricating closure."""
import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path


def inventory(path):
    p=Path(path)
    if not p.exists(): return {'status':'asset_unavailable','path':str(p)}
    raw=p.read_bytes(); proc=subprocess.run([str(p),'-h'],capture_output=True,text=True,check=False)
    lines=proc.stdout.splitlines()
    try:
        import openmm
        openmm_version=openmm.version.version
    except (ImportError, AttributeError):
        openmm_version='unavailable'
    return {'status':'available','executable':p.name,'path':str(p),'sha256':hashlib.sha256(raw).hexdigest(),'version':next((x.strip() for x in lines if 'Large-scale' in x),None),'build_info':[x.strip() for x in lines[:3]],'kspace_capability':'executed Ewald reference; PPPM capability requires separate input','precision':'double (declared reference)','platform':platform.platform(),'machine':platform.machine(),'openmm_version':openmm_version,'python_version':platform.python_version()}
def main():
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--executable',default='../island-validation/phase4j22/lmp-kspace');x=a.parse_args()
    manifest=inventory(x.executable)
    controls=['neutral_point_charge','counterion_system','pairwise_reference','small_polymer_fragment','linear_polymer','aromatic_resp_polymer','packed_network','box_scaling','cutoff_scaling','reciprocal_tolerance_scaling']
    result={'schema':'island.j40-c7-kspace-audit.v1','reference':manifest,'controls':{k:{'status':'blocked','reason':'existing component-level Ewald discrepancy or unavailable raw paired PME record'} for k in controls},'polymer_pressure':{'status':'pressure_unavailable','reason':'no independently closed force-field virial and finite-volume pair'},'scientific_acceptance':False}
    x.output.write_text(json.dumps(result,indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
