"""Execute the C9 root-cause isolation report without fabricating raw closure."""
import argparse
import json
from pathlib import Path

from island.analysis.ewald_root_cause import isolate


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 names=['point_charges','counterion_system','pairwise_control','small_fragment','linear_polymer','aromatic_polymer','packed_network']
 cases={n:{'status':'blocked','reason':'retained C8 raw component comparison is outside unchanged bound'} for n in names}
 a.output.write_text(json.dumps(isolate(cases),indent=2)+'\n');return 1
if __name__=='__main__':raise SystemExit(main())
