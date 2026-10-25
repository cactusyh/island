"""Preserve and diagnose the retained polymer Ewald discrepancy."""
import argparse
import json
from pathlib import Path

from island.analysis.ewald_diagnosis import diagnose


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 names=['point_charges','counterion','pairwise','polymer_fragment','linear_polymer','aromatic_polymer','network','box_scaling','cutoff_scaling','reciprocal_tolerance_scaling']
 cases={n:{'status':'blocked','reason':'retained C7 raw comparison remains outside unchanged component bound'} for n in names}
 a.output.write_text(json.dumps(diagnose(cases),indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
