"""The analytical NVE layer needs no optional scientific backend."""

import subprocess
import sys


def test_optional_import_isolation():
    code = """
import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'scipy', 'openmm', 'parmed', 'rdkit'}:
            raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0, Block())
import island
import island.evaluation
import island.minimization
import island.dynamics
from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import EvaluationResult
from island.core.coordinate_provenance import coordinate_hash
graph = Topology()
graph.add_site(AtomSite(17, 'C', 12., element='C', atomic_number=6))
system = MolecularSystem(graph, Coordinates({17: (0.,0.,0.)}))
class Free:
    def evaluate(self, coordinates=None, *, coordinate_unit='angstrom'):
        return EvaluationResult(0., {'free':0.}, {17:(0.,0.,0.)}, coordinate_hash(coordinates), 'params', 'model', 'frame', 'free', '1', 'python', {})
result = run_nve(system, Free(), {17:(1.,0.,0.)}, DynamicsOptions(1.,2,4,3))
assert result.completed
result.to_system(system).validate()
from island.dynamics import LangevinOptions, initialize_velocities, run_langevin
initial = initialize_velocities(system, temperature_kelvin=300, seed=5)
thermal = run_langevin(system, Free(), initial.velocities, LangevinOptions(1.,300.,2.,17,2,4,3))
assert thermal.completed
thermal.to_system(system).validate()
from island.dynamics import run_dynamics_segment, create_dynamics_checkpoint, resume_dynamics, DynamicsSegmentOptions
segment = run_dynamics_segment(system, Free(), initial.velocities, LangevinOptions(1.,300.,2.,17,2,4,3))
checkpoint = create_dynamics_checkpoint(segment)
continued = resume_dynamics(checkpoint, system, Free(), DynamicsSegmentOptions(2,4,3))
assert continued.completed
assert not {'scipy', 'openmm', 'parmed', 'rdkit'} & sys.modules.keys()
"""
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
