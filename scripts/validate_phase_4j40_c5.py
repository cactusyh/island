"""Audit actual C2 component/pressure artifacts without relabeling failures."""
import argparse
import json
from pathlib import Path

from island.analysis.backend_pressure import validate_actual_pressure


def audit(root):
 rows=[]
 for p in sorted(Path(root).glob('**/observations.json')):
  d=json.loads(p.read_text()); rows.append({'path':str(p),'pressure_status':d.get('pressure_status'),'measured_pressure':d.get('pressure_status')!='pressure_unavailable','actual_frames':len(d.get('observations',[]))})
 return {'schema':'island.j40-c5-actual-pressure-audit.v1','cells':rows,'backend_pressure':validate_actual_pressure(force_virial_bar=0.0,volume_derivative_bar=0.0,components={},long_range_validated=False),'cross_engine':{'status':'blocked','reason':'retained Ewald component discrepancy exceeds unchanged bound','kspace_executed':True},'scientific_acceptance':False}
def main():
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(audit(a.run),indent=2)+'\n');return 1
if __name__=='__main__': raise SystemExit(main())
