import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 if r.get('schema')!='island.j40-c6-receipt.v1': raise AnalysisError('Wrong C6 schema')
 if r.get('scientific_acceptance') is not False: raise AnalysisError('C6 relabeled blocked result')
 if r['kspace'].get('polymer_closure')!='ewald_unvalidated': raise AnalysisError('Polymer closure incorrectly passed')
 if r.get('backend_pressure')!='pressure_unavailable' or r.get('cutoff_recovery')!='blocked': raise AnalysisError('Unsafe C6 gate')
 if r['readiness']!={'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('receipt');a=p.parse_args();validate(a.receipt);print('C6 receipt integrity passed')
