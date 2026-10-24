import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c11-receipt.v1': raise AnalysisError('Wrong C11 schema')
 if r.get('scientific_acceptance') is not False or r.get('ewald_status')!='ewald_unvalidated': raise AnalysisError('C11 blocked result relabeled')
 if r.get('pressure_status')!='pressure_unavailable' or r.get('repair_applied') or r.get('cutoff_recovery')!='blocked': raise AnalysisError('Unsafe C11 gate')
 if r['readiness']!={'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C11 receipt integrity passed')
