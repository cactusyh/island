"""Import a hand-checkable synthetic prmtop; NOT a GAFF/GAFF2 example."""

from pathlib import Path
from tempfile import TemporaryDirectory

import parmed as pmd

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.forcefields import ParameterizedSystem, import_amber_prmtop


def synthetic_prmtop(path: Path) -> None:
    """Create a tiny already-parameterized Amber-format chain for demonstration."""
    structure = pmd.Structure()
    atom_type = pmd.AtomType("CX", 1, 12.01, 6)
    atom_type.set_lj_params(0.1, 1.9)
    atoms = []
    for index in range(4):
        atom = pmd.Atom(
            name=f"C{index}", type="CX", atomic_number=6, mass=12.01, charge=0
        )
        atom.atom_type = atom_type
        structure.add_atom(atom, "SYN", 1)
        atoms.append(atom)
    bond_type = pmd.BondType(100, 1.5)
    angle_type = pmd.AngleType(50, 109.5)
    torsion_type = pmd.DihedralType(1, 2, 180, scee=1.2, scnb=2)
    structure.bond_types.append(bond_type)
    structure.angle_types.append(angle_type)
    structure.dihedral_types.append(torsion_type)
    for left, right in ((0, 1), (1, 2), (2, 3)):
        structure.bonds.append(pmd.Bond(atoms[left], atoms[right], type=bond_type))
    for left, center, right in ((0, 1, 2), (1, 2, 3)):
        structure.angles.append(
            pmd.Angle(atoms[left], atoms[center], atoms[right], type=angle_type)
        )
    structure.dihedrals.append(pmd.Dihedral(*atoms, type=torsion_type))
    structure.bond_types.claim()
    structure.angle_types.claim()
    structure.dihedral_types.claim()
    pmd.amber.AmberParm.from_structure(structure).save(str(path), overwrite=True)


topology = Topology()
for source_index in range(4):
    topology.add_site(
        AtomSite(10 + source_index * 10, f"C{source_index}", 12.01,
                 element="C", atomic_number=6)
    )
for left, right in ((10, 20), (20, 30), (30, 40)):
    topology.add_bond(left, right, order=1)
system = MolecularSystem(topology, Coordinates())  # Import needs no coordinates.
mapping = {source_index: 10 + source_index * 10 for source_index in range(4)}

with TemporaryDirectory() as directory:
    path = Path(directory) / "synthetic.prmtop"
    synthetic_prmtop(path)
    imported = import_amber_prmtop(
        system, path, mapping, source="SYNTHETIC SOFTWARE-TEST TOPOLOGY"
    )
    snapshot = ParameterizedSystem.from_amber_import(system, imported)

print("stable site IDs:", sorted(imported.site_assignments))
print("bond/angle/proper counts:", len(imported.bond_assignments),
      len(imported.angle_assignments), len(imported.proper_torsion_assignments))
print("source checksum:", imported.source_sha256)
print("1-4 pairs:", imported.source_14_pairs)
print("production validated:", snapshot.metadata["aggregate"]["production_validated"])
print("simulation readiness:", snapshot.metadata["aggregate"]["simulation_readiness"])
