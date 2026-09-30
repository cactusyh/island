"""Explicitly owned reusable session for single points, minimization, and NVE."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics import DynamicsOptions, run_nve
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
options = DynamicsOptions(
    timestep_fs=0.1,
    steps=20,
    max_evaluations=22,
    max_frames=5,
    recording_interval=5,
    max_energy_deviation=0.01,
)
velocities = {17: (0.2, 0.1, 0.0), 91: (-0.1, 0.0, 0.1)}  # angstrom/ps
with evaluator.open_session() as session:
    first = session.evaluate()
    minimum = minimize_geometry(system, session)
    assert minimum.converged, minimum.optimizer_message
    minimized = minimum.to_system(system)
    result = run_nve(minimized, session, velocities, options)
    print(
        "initial potential / minimized potential (kJ/mol):",
        first.potential_energy,
        minimum.final_energy,
    )
    print("NVE completion / calls:", result.termination_reason, result.evaluations)
    print(
        "final energy / max absolute deviation:",
        result.final_state.total_energy,
        result.max_abs_energy_deviation,
    )
    print("model:", result.model_fingerprint)
    print("session Contexts: 1; separate final verification Contexts: 2")
print("closed:", session.closed)
print(
    "unchanged parent coordinate fingerprint:",
    evaluator.evaluate().coordinate_fingerprint,
)
print("simulation readiness:", result.simulation_readiness)
