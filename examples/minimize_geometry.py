"""Bounded local minimization of synthetic Amber records (software demonstration)."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import import_amber_prmtop
from island.minimization import MinimizationOptions, minimize_geometry

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
result = minimize_geometry(system, evaluator, MinimizationOptions(force_tolerance=0.01))
print(
    "energy before/after:",
    result.initial_energy,
    result.final_energy,
    result.energy_unit,
)
print("fmax before/after:", result.initial_fmax, result.final_fmax, result.force_unit)
print(
    "RMS atomic force before/after:", result.initial_rms_force, result.final_rms_force
)
print("termination:", result.termination_reason, result.optimizer_message)
print("iterations/evaluator calls:", result.iterations, result.evaluations)
print("model/coordinates:", result.model_fingerprint, result.coordinate_fingerprint)
print(
    "backend/optimizer:",
    result.final_evaluation.backend_version,
    result.optimizer_version,
)
if result.converged:
    optimized = result.to_system(system)
    print("coordinate provenance:", optimized.metadata["local_minimization"])
print("simulation readiness:", result.simulation_readiness)
