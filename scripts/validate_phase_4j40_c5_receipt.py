import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c5-receipt.v1': raise AnalysisError('Wrong C5 schema')
 if r.get('scientific_acceptance') is not False: raise AnalysisError('C5 blocked result relabeled')
 if r.get('backend_pressure')!='pressure_unavailable': raise AnalysisError('Backend pressure incorrectly validated')
 if r['cross_engine'].get('ewald_closure') is not False: raise AnalysisError('Cross-engine closure incorrectly passed')
 if r['readiness'] != {'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C5 receipt integrity passed')
