"""Explicit, content-bound GAFF2 runtime inventory; no system-install fallback."""

import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
from hashlib import sha256
from pathlib import Path

from island.exceptions import AmberToolsUnavailableError

SCHEMA = 'island.gaff2-runtime.v1'
TOOLS = ('antechamber', 'parmchk2', 'tleap', 'sqm')
SOURCES = ('dat/leap/cmd/leaprc.gaff2', 'dat/leap/parm/gaff2.dat')


def digest(value):
    return sha256((json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()).hexdigest()


def signed(value):
    return {**value, 'identity': digest(value)}


def require(ok, message):
    if not ok:
        raise AmberToolsUnavailableError(message)


def checksum(path):
    with Path(path).open('rb') as stream:
        return sha256(stream.read()).hexdigest()


def declared_environment(record, inherited=None):
    """Explicit activation, returned to the caller; never alters process globals."""
    env = dict(os.environ if inherited is None else inherited)
    home = Path(record['amberhome'])
    env['AMBERHOME'] = str(home)
    env['PATH'] = str(home / 'bin') + os.pathsep + env.get('PATH', '')
    return env


def _run(argv, env, *, cwd=None):
    result = subprocess.run(list(map(str, argv)), env=env, cwd=cwd, capture_output=True, text=True, check=False, timeout=120)
    return result


def inventory(home, *, environment):
    """Audit actual binaries, source files, package records and resolved libraries."""
    home = Path(home).resolve()
    require(platform.system() == 'Linux', 'Declared GAFF2 runtime requires Linux')
    require(environment.get('AMBERHOME') == str(home), 'AMBERHOME differs from declared runtime')
    require((home / 'conda-meta').is_dir(), 'Missing declared AmberTools environment: ' + str(home))
    for name in TOOLS:
        path = home / 'bin' / name
        require(path.is_file() and os.access(path, os.X_OK), 'Missing required executable: ' + str(path))
        found = shutil.which(name, path=environment.get('PATH', ''))
        require(found is not None and Path(found).resolve() == path.resolve(), 'Wrong PATH executable resolution for ' + name + ': ' + str(found))
    metadata_paths = sorted((home / 'conda-meta').glob('*.json'))
    packages = [json.loads(p.read_text()) for p in (home / 'conda-meta').glob('ambertools-*.json')]
    require(len(packages) == 1 and packages[0].get('name') == 'ambertools', 'Missing or ambiguous AmberTools package record')
    package = packages[0]
    package_files = {}
    for name in package['files']:
        path = home / name
        # conda records compiled Python bytecode, which is not a runtime asset.
        if name.endswith('.pyc'):
            continue
        require(path.is_file() and path.resolve().is_relative_to(home), 'Incomplete AmberTools package: ' + name)
        package_files[name] = checksum(path)
    sources = {}
    for name in SOURCES:
        require((home / name).is_file(), 'Missing GAFF2 source: ' + name)
        sources[name] = checksum(home / name)
    require('gaff2.dat' in (home / SOURCES[0]).read_text(), 'GAFF2 leaprc does not load declared source')
    executables, libraries = {}, {}
    with tempfile.TemporaryDirectory(prefix='island-gaff2-probe-') as tmp:
        for name in TOOLS:
            path = home / 'bin' / name
            native = home / 'bin' / 'teLeap' if name == 'tleap' else home / 'bin' / 'wrapped_progs' / name
            if not native.is_file():
                native = path
            require(os.access(native, os.X_OK), 'Missing native AmberTools executable: ' + str(native))
            ldd = _run(['ldd', native], environment)
            require(ldd.returncode == 0 and 'not found' not in ldd.stdout + ldd.stderr, 'Required AmberTools library unavailable for ' + name + ': ' + ldd.stdout + ldd.stderr)
            deps = sorted(set(re.findall(r'(?:=>\s+|^\s*)(/\S+)\s+\(', ldd.stdout, re.MULTILINE)))
            require(deps, 'No dynamic-library inventory for ' + name)
            for dep in deps:
                libraries[dep] = checksum(dep)
            help_result = _run([path, '-h'], environment, cwd=tmp)
            # Package release and executable content are authoritative. Several
            # Amber programs do not expose the package version in their banner.
            banner = next((line.strip() for line in (help_result.stdout + help_result.stderr).splitlines() if re.search(r'(?:version|Version|antechamber|parmchk)', line)), None)
            executables[name] = {'path': str(path), 'sha256': checksum(path), 'native_path': str(native), 'native_sha256': checksum(native), 'package_version': package['version'], 'package_build': package['build'], 'observed_banner': banner, 'help_exit_code': help_result.returncode, 'libraries': deps}
    # Prevent a caller's PYTHONPATH from influencing the declared environment audit.
    probe_env = {k: v for k, v in environment.items() if k not in {'PYTHONPATH', 'PYTHONHOME'}}
    probe = "import importlib.metadata as m,json,hashlib; print(json.dumps(sorted([(d.metadata['Name'],d.version,hashlib.sha256((d.read_text('RECORD') or '').encode()).hexdigest()) for d in m.distributions()])))"
    py = home / 'bin' / 'python'
    require(py.is_file(), 'Missing AmberTools environment Python')
    installed = _run([py, '-c', probe], probe_env)
    require(installed.returncode == 0, 'Incomplete AmberTools Python environment: ' + installed.stderr)
    pip = _run([py, '-m', 'pip', 'check'], probe_env)
    require(pip.returncode == 0, 'Incomplete AmberTools dependencies: ' + pip.stdout + pip.stderr)
    return {
        'schema': SCHEMA, 'amberhome': str(home),
        'platform': {'system': platform.system(), 'machine': platform.machine()},
        'package': {k: package.get(k) for k in ('name', 'version', 'build', 'channel', 'url', 'sha256')},
        'conda_records': {p.name: checksum(p) for p in metadata_paths},
        'package_files': package_files, 'executables': executables, 'sources': sources,
        'libraries': libraries, 'python_sha256': checksum(py),
        'python_distributions': json.loads(installed.stdout), 'pip_check': 'passed',
    }


def capture_runtime(home):
    env = declared_environment({'amberhome': str(Path(home).resolve())})
    return signed(inventory(home, environment=env))


def validate_runtime(record, *, expected_identity, environment=None):
    """Never accept a re-signed replacement for an independently pinned manifest."""
    try:
        require(record.get('schema') == SCHEMA, 'Wrong GAFF2 runtime schema')
        body = {k: v for k, v in record.items() if k != 'identity'}
        require(record.get('identity') == digest(body) == expected_identity, 'GAFF2 runtime manifest identity mismatch')
        actual = inventory(record['amberhome'], environment=dict(os.environ if environment is None else environment))
        for key, value in body.items():
            require(value == actual.get(key), 'Changed GAFF2 runtime ' + key + ' fingerprint')
    except AmberToolsUnavailableError:
        raise
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        raise AmberToolsUnavailableError('Cannot validate declared GAFF2 runtime: ' + str(error)) from error
    return record
