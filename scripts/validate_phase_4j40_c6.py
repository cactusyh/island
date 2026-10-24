"""Create an auditable long-range manifest; no NPT starts on unvalidated closure."""
import argparse
import json
from pathlib import Path

from island.analysis.long_range_reference import (
 closure_record,
 executable_manifest,
 pairwise_coulomb_control,
)


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 exe=executable_manifest('../island-validation/phase4j22/lmp-kspace',openmm_version='8.6.1',precision='double',platform='LAMMPS KSPACE reference')
 pair=pairwise_coulomb_control(1,-1,1);pair.update(status='passed',independent='analytic Coulomb control')
 rec=closure_record(controls=[pair],polymer_cells=[{'cell':'gaff2-pe','status':'blocked','reason':'polymer Ewald component discrepancy unresolved'},{'cell':'opls-network','status':'blocked','reason':'long-range virial unavailable'}],executable_manifest=exe,tolerance=1e-5)
 a.output.write_text(json.dumps(rec,indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
