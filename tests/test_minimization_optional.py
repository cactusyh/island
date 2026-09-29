"""Core and evaluation remain importable without SciPy or chemistry extras."""

import subprocess
import sys


def test_minimization_dependencies_are_lazy():
    code = """
import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'scipy','openmm','parmed','rdkit'}:
            raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0,Block())
import island
import island.evaluation
import island.minimization
from island.minimization.engine import _scipy
from island.exceptions import MinimizationUnavailableError
assert not {'scipy','openmm','parmed','rdkit'} & sys.modules.keys()
try:
    _scipy()
except MinimizationUnavailableError as error:
    assert 'island[minimization]' in str(error)
else:
    raise AssertionError('SciPy unavailability was not reported')
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
