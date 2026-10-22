"""Diagnose RESP/Psi4 acceptance without substituting external OPLS-AA charges.

The present declared environments have no quantum/RESP assets. This diagnostic
entry point deliberately cannot claim workflow acceptance. A validated loader is
provided separately in island.forcefields.ambertools.resp for reviewed results.
"""

import argparse
import json
import subprocess
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest
from island.forcefields.ambertools.resp import file_hash

BASE = '7f1a4201c3e8bb82628acf5f01860005874827dd'
PYTHONS = ('../island-validation/phase4e2-env/bin/python',
           '../island-validation/phase4g1-env/bin/python',
           '../island-validation/phase4j32-analysis-env/bin/python',
           '../miniforge3/envs/md_rg/bin/python')
CODE = 'import sys,importlib.util,shutil,json; print(json.dumps({"python":sys.executable,"version":sys.version,"psi4_module":str(importlib.util.find_spec("psi4")),"resp_module":str(importlib.util.find_spec("resp")),"psi4_executable":shutil.which("psi4")}))'


def signed(body):
    return {**body, 'identity': reference_digest(body)}


def diagnose(pythons=PYTHONS):
    audit = []
    for python in pythons:
        try:
            run = subprocess.run([python, '-c', CODE], capture_output=True, text=True, check=False)
            audit.append({'argv': [python, '-c', CODE], 'exit_code': run.returncode,
                          'stdout': run.stdout, 'stderr': run.stderr})
        except OSError as error:
            audit.append({'argv': [python, '-c', CODE], 'exit_code': None, 'error': str(error)})
    controls = {}
    for kind in ('linear', 'network'):
        path = Path(f'docs/evidence/phase_4j35/references/{kind}.json')
        ref = json.loads(path.read_text())
        controls[kind] = {'charge_origin': 'external_oplsaa', 'charge_method': 'oplsaa_provided',
                          'reference_identity': ref['identity'], 'reference_sha256': file_hash(path),
                          'historical_method': ref['method'], 'resp_psi4': False}
    return signed({'schema': 'island.j35-c1-diagnostic.v1', 'base_commit': BASE,
                   'classification': 'asset_unavailable', 'acceptance_met': False,
                   'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
                   'diagnostic': 'No pinned Psi4/RESP runtime and independently validated Psi4 ESP/RESP fit results are declared for either retained graph. Required: Psi4 executable/library, version and environment fingerprints; a pinned RESP fitter; quantum input/output, ESP grid, fit input/output and full-precision mapped charges with independent review.',
                   'environment_audit': audit, 'external_controls': controls,
                   'cases': {kind: {'status': 'asset_unavailable', 'acceptance_met': False,
                                    'charge_generation_executed': False, 'dynamics_executed': False,
                                    'downstream': 'not_executed'} for kind in ('linear', 'network')},
                   'native_am1bcc': {'status': 'blocked', 'classification': 'native_charge_generation_unavailable',
                                     'observed_component_total_e': '-0.003', 'acceptance_met': False},
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = diagnose()
    (args.output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print('BLOCKED: ' + result['diagnostic'])
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
