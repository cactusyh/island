import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c7-receipt.v1': raise AnalysisError('Wrong C7 schema')
 if r.get('scientific_acceptance') is not False: raise AnalysisError('C7 blocked result relabeled')
 if r.get('polymer_ewald')!='ewald_unvalidated': raise AnalysisError('Unclosed polymer Ewald')
 if r.get('backend_pressure')!='pressure_unavailable': raise AnalysisError('Unvalidated backend pressure')
 if r['readiness']!={'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C7 receipt integrity passed')
