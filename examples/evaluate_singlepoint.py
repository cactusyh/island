"""Two synthetic charged sites: a nonperiodic single point, no integration steps."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import import_amber_prmtop

source = pmd.Structure()
graph = Topology()
for index, (site_id, charge) in enumerate(((17, 1), (91, -1))):
    atom_type = pmd.AtomType(f"C{index}", index + 1, 12.01, 6)
    atom_type.set_lj_params(0.1, 1.9)
    atom = pmd.Atom(
        name=f"C{index}", type=f"C{index}", atomic_number=6, mass=12.01, charge=charge
    )
    atom.atom_type = atom_type
    source.add_atom(atom, "SYN", 1)
    graph.add_site(
        AtomSite(
            site_id,
            atom.name,
            12.01,
            element="C",
            atomic_number=6,
            formal_charge=charge,
        )
    )
system = MolecularSystem(graph, Coordinates({17: (0, 0, 0), 91: (4.2, 0, 0)}))
with TemporaryDirectory() as directory:
    path = Path(directory) / "synthetic.prmtop"
    pmd.amber.AmberParm.from_structure(source).save(str(path))
    imported = import_amber_prmtop(
        system,
        path,
        {0: 17, 1: 91},
        source="synthetic software test, not a physical model",
    )
    evaluator = OpenMMSinglePointEvaluator(system, imported, platform="Reference")

# No source files or ParmEd are consulted after binding.
first = evaluator.evaluate()
second = evaluator.evaluate(
    {17: (0, 0, 0), 91: (4.3, 0.1, 0)}, coordinate_unit="angstrom"
)
print("potential energy:", first.potential_energy, first.energy_unit)
print("components:", dict(first.energy_components))
print("forces by stable site ID:", dict(first.forces), first.force_unit)
print("backend:", first.backend_name, first.backend_version, first.platform)
print("new coordinate fingerprint:", second.coordinate_fingerprint)
print("unchanged model:", first.model_fingerprint == second.model_fingerprint)
print("simulation readiness:", first.simulation_readiness)
