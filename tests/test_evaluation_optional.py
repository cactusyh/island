"""Optional backend import contract, runs without any scientific extras."""

import subprocess
import sys


def test_evaluation_import_without_scientific_extras():
    code = """
import importlib.abc
import sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'openmm', 'rdkit', 'parmed'}:
            raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0, Block())
import island
from island.evaluation import PotentialEvaluator, EvaluationResult, OpenMMSinglePointEvaluator
from island.evaluation.openmm import _openmm
from island.exceptions import EvaluationUnavailableError
assert not {'openmm', 'rdkit', 'parmed'} & sys.modules.keys()
try:
    _openmm()
except EvaluationUnavailableError as error:
    assert 'island[evaluation]' in str(error)
else:
    raise AssertionError('Missing OpenMM was not diagnosed')
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
