import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c4-receipt.v1': raise AnalysisError('Wrong C4 schema')
 if r.get('scientific_acceptance') is not False: raise AnalysisError('C4 blocked result relabeled')
 if r.get('pressure_status')!='pressure_unavailable' or r['pressure_controls']['backend_virial_validated']: raise AnalysisError('Unvalidated backend pressure')
 if r['recovery']['dynamics_executed'] or not r['recovery']['original_failures_preserved']: raise AnalysisError('Unsafe recovery record')
 if r['readiness'] != {'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C4 receipt integrity passed')
