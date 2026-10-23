"""Pinned executable Psi4/CDS RESP runtime; no system PATH fallback."""

import json
import os
import subprocess
import tempfile
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest
from island.forcefields.ambertools.resp import file_hash

SCHEMA = 'island.psi4-resp-runtime.v1'
PSI4_VERSION = '1.9.1'
RESP_PACKAGE_VERSION = '1.0'
RESP_MODULE_VERSION = '0.8'
PROBE = '''import sys,json,psi4,resp,numpy,importlib.metadata as m
print(json.dumps({'python_version':sys.version,'python_executable':sys.executable,'psi4_version':psi4.__version__,'psi4_module':psi4.__file__,'resp_module_version':resp.__version__,'resp_module':resp.__file__,'numpy_version':numpy.__version__,'packages':sorted((d.metadata['Name'],d.version) for d in m.distributions())}))'''


class Psi4RuntimeError(ValueError):
    """An explicit runtime gate, distinct from a charge or scientific result."""

    def __init__(self, message, classification='runtime_integrity_failure'):
        super().__init__(message)
        self.classification = classification


def require(ok, message, classification='runtime_integrity_failure'):
    if not ok:
        raise Psi4RuntimeError(message, classification)


def environment(prefix):
    prefix = Path(prefix).resolve()
    env = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'PYTHONHOME', 'LD_LIBRARY_PATH', 'AMBERHOME')}
    env.update(PATH=str(prefix / 'bin') + os.pathsep + env.get('PATH', ''), OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', ISLAND_J36_QM_THREADS='4')
    return env


def probe(prefix):
    python = Path(prefix).resolve() / 'bin/python'
    require(python.is_file(), 'Missing declared Psi4 Python: ' + str(python), 'asset_unavailable')
    require((python.parent / 'psi4').is_file(), 'Missing declared Psi4 executable', 'asset_unavailable')
    with tempfile.TemporaryDirectory(prefix='island-psi4-probe-') as scratch:
        run = subprocess.run([str(python), '-c', PROBE], env=environment(prefix), cwd=scratch, capture_output=True, text=True, check=False, timeout=120)
    require(run.returncode == 0, 'Psi4/RESP import or dependency execution failed:\n' + run.stdout + run.stderr, 'runtime_execution_failure')
    try:
        result = json.loads(run.stdout)
    except ValueError as error:
        raise Psi4RuntimeError('Malformed runtime probe output', 'runtime_execution_failure') from error
    require(result['psi4_version'] == PSI4_VERSION and result['resp_module_version'] == RESP_MODULE_VERSION, 'Wrong Psi4/RESP version')
    return result


def capture_runtime(prefix):
    prefix = Path(prefix).resolve()
    observed = probe(prefix)
    metadata = {}
    for path in sorted((prefix / 'conda-meta').glob('*.json')):
        data = json.loads(path.read_text())
        metadata[data['name']] = {k: data.get(k) for k in ('name', 'version', 'build', 'url', 'sha256', 'md5')}
    require(metadata['psi4']['version'] == PSI4_VERSION and metadata['resp']['version'] == RESP_PACKAGE_VERSION, 'Wrong native package versions')
    # Pin the compatible LibXC/Libint builds, not merely Psi4's broad solver bounds.
    require(metadata['libxc-c']['version'] == '7.0.0' and metadata['libint']['version'] == '2.9.0', 'Unsupported quantum integral/functional libraries')
    sites = Path(observed['psi4_module']).parent.parent
    files = {prefix / 'bin/python', prefix / 'bin/psi4'}
    files.update((prefix / 'conda-meta').glob('*.json'))
    for package in ('psi4', 'resp'):
        files.update((sites / package).rglob('*.py'))
        files.update((sites / package).rglob('*.so'))
    files.update((prefix / 'share/psi4/basis').glob('*.gbs'))
    for pattern in ('libxc.so*', 'libint2.so*', 'libmkl_rt.so*', 'libiomp5.so*', 'libomp.so*', 'libstdc++.so*', 'libgcc_s.so*', 'libpython3.11.so*'):
        files.update((prefix / 'lib').glob(pattern))
    hashes = {str(p.relative_to(prefix)): file_hash(p) for p in sorted(files) if p.is_file()}
    body = {'schema': SCHEMA, 'prefix': str(prefix), 'probe': observed, 'packages': metadata,
            'files': hashes, 'psi4_version': PSI4_VERSION, 'resp_package_version': RESP_PACKAGE_VERSION,
            'resp_module_version': RESP_MODULE_VERSION, 'resp_implementation': 'cdsgroup/resp',
            'version_note': 'Conda package 1.0 contains module __version__ 0.8; both pinned without relabeling.',
            'probe_command': [str(prefix / 'bin/python'), '-c', PROBE]}
    return {**body, 'identity': reference_digest(body)}


def validate_runtime(record, *, expected_identity, execute_probe=True):
    require(record.get('schema') == SCHEMA and record.get('identity') == expected_identity == reference_digest({k: v for k, v in record.items() if k != 'identity'}), 'Changed Psi4 runtime identity')
    require(record.get('psi4_version') == PSI4_VERSION and record.get('resp_package_version') == RESP_PACKAGE_VERSION and record.get('resp_module_version') == RESP_MODULE_VERSION, 'Wrong frozen Psi4/RESP version')
    prefix = Path(record['prefix']).resolve()
    require({'bin/python', 'bin/psi4'} <= set(record['files']), 'Incomplete quantum executable inventory')
    for relative, expected in record['files'].items():
        path = prefix / relative
        require(not Path(relative).is_absolute() and '..' not in Path(relative).parts and path.is_file(), 'Missing/unsafe quantum asset: ' + relative, 'asset_unavailable')
        require(file_hash(path) == expected, 'Changed quantum source/environment asset: ' + relative)
    if execute_probe:
        require(probe(prefix) == record['probe'], 'Quantum runtime probe/version/dependencies changed')
    return record
