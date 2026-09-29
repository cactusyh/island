"""Short isolated NVE trajectory of synthetic Amber bond records; no AmberTools."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import import_amber_prmtop

source = pmd.Structure()
graph = Topology()
for index, site in enumerate((17, 91)):
    kind = pmd.AtomType(f"C{index}", index + 1, 12.01, 6)
    kind.set_lj_params(0.1, 1.9)
    atom = pmd.Atom(name=f"C{index}", type=f"C{index}", atomic_number=6, mass=12.01)
    atom.atom_type = kind
    source.add_atom(atom, "SYN", 1)
    graph.add_site(AtomSite(site, atom.name, 12.01, element="C", atomic_number=6))
source.bond_types.append(pmd.BondType(100, 1.5))
source.bonds.append(pmd.Bond(*source.atoms, type=source.bond_types[0]))
graph.add_bond(17, 91)
system = MolecularSystem(graph, Coordinates({17: (0, 0, 0), 91: (1.9, 0.1, 0)}))
with TemporaryDirectory() as directory:
    path = Path(directory) / "synthetic.prmtop"
    pmd.amber.AmberParm.from_structure(source).save(str(path))
    imported = import_amber_prmtop(
        system, path, {0: 17, 1: 91}, source="synthetic software test"
    )
    evaluator = OpenMMSinglePointEvaluator(system, imported, platform="Reference")
options = DynamicsOptions(
    timestep_fs=0.1,
    steps=20,
    max_evaluations=22,
    max_frames=5,
    recording_interval=5,
    max_energy_deviation=0.01,
)
velocities = {17: (0.2, 0.1, 0.0), 91: (-0.1, 0.0, 0.1)}  # angstrom/ps
result = run_nve(system, evaluator, velocities, options)
for frame in result.frames:
    print(
        "step/time(ps)/U/K/E(kJ/mol):",
        frame.step,
        frame.time_ps,
        frame.potential_energy,
        frame.kinetic_energy,
        frame.total_energy,
    )
print("completed/requested:", result.completed_steps, result.requested_steps)
print("termination/calls:", result.termination_reason, result.evaluations)
print("max absolute energy deviation (kJ/mol):", result.max_abs_energy_deviation)
print(
    "synchronized final velocities (angstrom/ps):", dict(result.final_state.velocities)
)
print("dynamics fingerprint:", result.dynamics_fingerprint)
print("simulation readiness:", result.simulation_readiness)
if result.completed:
    moved = result.to_system(system)
    print("current coordinate source:", moved.metadata["coordinate_source"])
    print("Keep result.final_state.velocities: MolecularSystem alone is not a restart.")
