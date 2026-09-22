import os
import subprocess
import sys


def test_core_import_does_not_eagerly_import_rdkit() -> None:
    code = """
import importlib.abc
import sys

class BlockRDKit(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'rdkit' or fullname.startswith('rdkit.'):
            raise RuntimeError('island.core attempted to import RDKit')
        return None

sys.meta_path.insert(0, BlockRDKit())
import island.core
from island.core import MolecularSystem
assert 'rdkit' not in sys.modules
assert MolecularSystem.__name__ == 'MolecularSystem'
"""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=False,
        env=environment,
        text=True,
    )
    assert result.returncode == 0, result.stderr
