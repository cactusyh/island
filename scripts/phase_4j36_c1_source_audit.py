"""Reproduce pinned source/format inspection; no charge generation or repair."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from island.forcefields.ambertools.provided import require
from island.forcefields.ambertools.psi4_runtime import validate_runtime
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_precision import signed


def collect():
    runtime = json.loads(Path('docs/evidence/phase_4j36/environment.json').read_text())
    validate_runtime(runtime, expected_identity=runtime['identity'])
    names = ['lib/python3.11/site-packages/resp/driver.py', 'lib/python3.11/site-packages/resp/espfit.py',
             'lib/python3.11/site-packages/resp/vdw_surface.py', 'lib/python3.11/site-packages/psi4/tests/python/mints14/input.py']
    sources, excerpts = {}, {}
    for name in names:
        path = Path(runtime['prefix']) / name
        sources[name] = file_hash(path)
        require(runtime['files'][name] == sources[name], 'Pinned source changed: ' + name)
        excerpts[name] = [{'line': i, 'text': line} for i, line in enumerate(path.read_text().splitlines(), 1)
                          if any(token in line for token in ('np.linalg.solve', 'np.zeros', 'savetxt', 'loadtxt', 'return qf', '%12.8f', 'compute_esp_over_grid_in_memory'))]
    amber = json.loads(Path('docs/evidence/phase_4j35_c1/artifacts/standalone-resp-audit.json').read_text())[0]
    require(file_hash(amber['path']) == amber['sha256'], 'Changed pinned Amber RESP executable')
    with tempfile.TemporaryDirectory(prefix='j36-c1-amber-help-') as tmp:
        help_run = subprocess.run([amber['path'], '-h'], cwd=tmp, capture_output=True, text=True, check=False)
    strings = subprocess.run(['strings', amber['path']], capture_output=True, text=True, check=True)
    return signed({'schema': 'island.resp-precision-source-audit.v1', 'runtime_identity': runtime['identity'],
                   'source_files': sources, 'source_excerpts': excerpts, 'psi4_version': '1.9.1',
                   'resp_package_version': '1.0', 'resp_module_version': '0.8',
                   'native_solver_dtype': 'float64', 'native_solver_significand_bits': 53,
                   'higher_precision_native_export_identified': False,
                   'scope': 'No charge output exceeding the 53-bit native ndarray was identified in the inspected pinned fitter. Native Psi4 in-memory ESP improves the intermediate transfer and is used separately; it does not add charge-solver bits.',
                   'esp_text_quantization': 'Original CDS driver loads grid_esp.dat text; new worker transfers native ESPPropCalc Vector directly to espfit.fit.',
                   'charge_import_precision_loss': False, 'native_constraint_algorithm': 'NumPy float64 constrained linear solves with iterative hyperbolic restraints',
                   'amber_resp': {'path': amber['path'], 'sha256': amber['sha256'], 'package': amber['package'],
                                  'help_command': [amber['path'], '-h'], 'help_exit_code': help_run.returncode,
                                  'help_stdout': help_run.stdout, 'help_stderr': help_run.stderr,
                                  'binary_format_strings': [line for line in strings.stdout.splitlines() if '(8f10.6)' in line],
                                  'interpretation': 'Standalone RESP qout uses six decimals; not a higher-precision replacement for native float64 charges. No Amber fit is substituted or claimed.'}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.output.open('x') as stream:
        json.dump(collect(), stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
