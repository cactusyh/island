import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c9-receipt.v1': raise AnalysisError('Wrong C9 schema')
 if r.get('scientific_acceptance') is not False or r.get('root_cause') is not None: raise AnalysisError('C9 unresolved result relabeled')
 if r.get('ewald_status')!='ewald_unvalidated' or r.get('pressure_status')!='pressure_unavailable': raise AnalysisError('Unsafe C9 gate')
 if r.get('repair_applied') or r.get('cutoff_recovery')!='blocked': raise AnalysisError('Unjustified repair or recovery')
 if r['readiness']!={'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C9 receipt integrity passed')
