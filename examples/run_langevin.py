"""Thermostatted isolated-molecule trajectory of synthetic Amber bond records; no AmberTools."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics import LangevinOptions, initialize_velocities, run_langevin
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import import_amber_prmtop
from island.minimization import minimize_geometry

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
with evaluator.open_session() as session:
    minimum = minimize_geometry(system, session)
    if not minimum.converged:
        raise RuntimeError(minimum.termination_reason)
    starting = minimum.to_system(system)
    initialization = initialize_velocities(starting, temperature_kelvin=300, seed=78123)
    options = LangevinOptions(
        timestep_fs=0.1,
        temperature_kelvin=300,
        friction_per_ps=5,
        thermostat_seed=99181,
        steps=200,
        max_evaluations=202,
        max_frames=5,
        recording_interval=50,
    )
    result = run_langevin(starting, session, initialization.velocities, options)
for frame in result.frames:
    print(
        "step/time(ps)/U/K/E(kJ/mol)/T(K):",
        frame.step,
        frame.time_ps,
        frame.potential_energy,
        frame.kinetic_energy,
        frame.total_energy,
        frame.instantaneous_temperature_kelvin,
    )
print(
    "termination/completed/requested:",
    result.termination_reason,
    result.completed_steps,
    result.requested_steps,
)
print(
    "calls/final independently verified:",
    result.evaluations,
    result.final_evaluation_verified,
)
print(
    "thermostat seed/normal draws:", result.options.thermostat_seed, result.normal_draws
)
print("dynamics fingerprint:", result.dynamics_fingerprint)
print("simulation readiness:", result.simulation_readiness)
if result.completed:
    moved = result.to_system(starting)
    print("current coordinate source:", moved.metadata["coordinate_source"])
print("Velocities and RNG identity remain separate; MolecularSystem is not a restart.")
