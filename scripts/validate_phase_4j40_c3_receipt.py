"""Validate the additive C3 receipt without rewriting C2 evidence."""
import argparse
import json
from pathlib import Path

from island.exceptions import AnalysisError


def validate(path):
 r=json.loads(Path(path).read_text())
 required=['schema','base_commit','scientific_acceptance','preflight_policy','pressure','blocked_gates','readiness']
 if any(k not in r for k in required): raise AnalysisError('Incomplete C3 receipt')
 if r['schema']!='island.j40-c3-receipt.v1': raise AnalysisError('Wrong C3 receipt schema')
 if r['scientific_acceptance'] is not False: raise AnalysisError('Blocked C3 receipt relabeled')
 p=r['preflight_policy']
 if not (p['safety_margin_nm']>=0 and 'half smallest predicted box' in p['condition']): raise AnalysisError('Changed preflight contract')
 if r['pressure']['status']!='pressure_unavailable' or not r['pressure']['target_is_not_measurement']: raise AnalysisError('Pressure target copied as measurement')
 if r['readiness']!={'production_validated':False,'simulation_readiness':'not_established'}: raise AnalysisError('Changed readiness')
 return r
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('receipt');x=a.parse_args();validate(x.receipt);print('C3 receipt integrity passed')
