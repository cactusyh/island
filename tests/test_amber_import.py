"""Hand-checkable synthetic Amber fixtures; not GAFF reference data."""

import builtins
import copy
import subprocess
import sys
from dataclasses import replace
from math import cos, pi

import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.exceptions import (
    AmberImportError,
    IncompleteChargeAssignmentError,
    InvalidAmberImportResultError,
    InvalidAmberMappingError,
    UnsupportedAmberFeatureError,
)
from island.forcefields.amber import import_amber_prmtop
from island.forcefields.parameterized import ParameterizedSystem

pmd = pytest.importorskip("parmed")


def synthetic_source(
    tmp_path, *, improper=False, reverse=False, multi=False, epsilon=0.1,
    zero_first=False,
):
    """Four carbon atoms, zero charges, manually chosen source parameters."""
    structure = pmd.Structure()
    carbon = pmd.AtomType("CX", 1, 12.01, 6)
    carbon.set_lj_params(epsilon, 1.9)
    zero_carbon = pmd.AtomType("CZ", 2, 12.01, 6)
    zero_carbon.set_lj_params(0, 0)
    atoms = [
        pmd.Atom(name=f"C{i}", type="CX", atomic_number=6, mass=12.01, charge=0)
        for i in range(4)
    ]
    source_order = tuple(reversed(range(4))) if reverse else tuple(range(4))
    for logical_index in source_order:
        atom = atoms[logical_index]
        atom.atom_type = zero_carbon if zero_first and logical_index == 0 else carbon
        if zero_first:
            atom.charge = 0.2 if logical_index == 0 else -0.2 if logical_index == 1 else 0
        structure.add_atom(atom, "SYN", 1)
    bond_type = pmd.BondType(100.0, 1.5)
    angle_type = pmd.AngleType(50.0, 109.5)
    torsion_type = pmd.DihedralType(1.5 if improper else 0.0, 2, 180, scee=1.2, scnb=2)
    structure.bond_types.append(bond_type)
    structure.angle_types.append(angle_type)
    structure.dihedral_types.append(torsion_type)
    if multi:
        other_type = pmd.DihedralType(0.75, 3, 0, scee=1.2, scnb=2)
        structure.dihedral_types.append(other_type)
    edges = ((0, 2), (1, 2), (2, 3)) if improper else ((0, 1), (1, 2), (2, 3))
    for a, b in edges:
        structure.bonds.append(pmd.Bond(atoms[a], atoms[b], type=bond_type))
    angle_indices = ((0, 2, 1), (0, 2, 3), (1, 2, 3)) if improper else ((0, 1, 2), (1, 2, 3))
    for a, b, c in angle_indices:
        structure.angles.append(pmd.Angle(atoms[a], atoms[b], atoms[c], type=angle_type))
    if improper:
        structure.dihedrals.append(
            pmd.Dihedral(*atoms, improper=True, ignore_end=True, type=torsion_type)
        )
    else:
        structure.dihedrals.append(pmd.Dihedral(*atoms, type=torsion_type))
        if multi:
            structure.dihedrals.append(
                pmd.Dihedral(*atoms, ignore_end=True, type=other_type)
            )
    structure.bond_types.claim()
    structure.angle_types.claim()
    structure.dihedral_types.claim()
    path = tmp_path / "synthetic.prmtop"
    pmd.amber.AmberParm.from_structure(structure).save(str(path), overwrite=True)
    mapping = {source_index: 10 + 10 * logical_index
               for source_index, logical_index in enumerate(source_order)}
    topology = Topology()
    for i in reversed(range(4)) if reverse else range(4):
        topology.add_site(AtomSite(10 + 10 * i, f"C{i}", 12.01, element="C", atomic_number=6))
    for a, b in reversed(edges) if reverse else edges:
        topology.add_bond(10 + 10 * a, 10 + 10 * b, order=1)
    system = MolecularSystem(topology, Coordinates({site_id: (float(site_id), 0, 0) for site_id in mapping.values()}))
    return system, path, mapping


@pytest.mark.parametrize("reverse", [False, True])
def test_resolved_import_mapping_units_snapshot(tmp_path, reverse):
    system, path, mapping = synthetic_source(tmp_path, reverse=reverse)
    original = copy.deepcopy(system.to_dict())
    result = import_amber_prmtop(system, path, mapping, source="synthetic test")
    assert result.is_compatible_with(system)
    assert system.to_dict() == original
    assert set(result.site_assignments) == set(mapping.values())
    bond = next(iter(result.bond_assignments.values())).parameter
    angle = next(iter(result.angle_assignments.values())).parameter
    torsion = next(iter(result.proper_torsion_assignments.values())).parameter
    lj = next(iter(result.site_assignments.values())).parameter
    assert bond.force_constant == pytest.approx(83680)
    assert bond.equilibrium_length == pytest.approx(0.15)
    assert angle.force_constant == pytest.approx(418.4)
    assert lj.epsilon == pytest.approx(0.4184)
    assert lj.sigma == pytest.approx(0.38 / 2 ** (1 / 6))
    source_lj = 0.1 * ((3.8 / 4.1) ** 12 - 2 * (3.8 / 4.1) ** 6) * 4.184
    converted_lj = 4 * lj.epsilon * ((lj.sigma / 0.41) ** 12 - (lj.sigma / 0.41) ** 6)
    assert source_lj == pytest.approx(converted_lj)
    assert torsion.terms[0].force_constant == 0
    assert result.nonbonded_policy.lj_scale_14 == pytest.approx(0.5)
    assert result.nonbonded_policy.coulomb_scale_14 == pytest.approx(1 / 1.2)
    assert len(result.source_14_pairs) == 1
    snapshot = result.to_parameterized_system(system)
    second_snapshot = ParameterizedSystem.from_amber_import(system, result)
    assert second_snapshot.aggregate_signature == snapshot.aggregate_signature
    assert snapshot.metadata["aggregate"]["simulation_readiness"] == "not_established"
    system.coordinates.set(mapping[0], (99, 99, 99))
    assert result.is_compatible_with(system)
    assert snapshot.system.coordinates.get(mapping[0])[0] != 99
    assert copy.deepcopy(result).is_compatible_with(system)
    assert replace(result).is_compatible_with(system)
    corrupted = replace(result, bond_assignments={})
    assert not corrupted.is_compatible_with(system)
    with pytest.raises(InvalidAmberImportResultError):
        corrupted.to_parameterized_system(system)
    selected = result.site_assignments[mapping[0]]
    changed_site = replace(
        selected, parameter=replace(selected.parameter, epsilon=selected.parameter.epsilon + 1)
    )
    changed = replace(result, site_assignments={**result.site_assignments, mapping[0]: changed_site})
    assert not changed.is_compatible_with(system)
    with pytest.raises(InvalidAmberImportResultError):
        ParameterizedSystem.from_amber_import(system, changed)


def test_improper_order_and_geometry_energy(tmp_path):
    system, path, mapping = synthetic_source(tmp_path, improper=True)
    result = import_amber_prmtop(system, path, mapping, source="synthetic improper")
    assert tuple(result.improper_assignments) == (tuple(mapping[i] for i in range(4)),)
    improper = next(iter(result.improper_assignments.values())).parameter
    assert improper.central_atom_position == 3
    assert improper.terms[0].force_constant == pytest.approx(1.5 * 4.184)
    # Independently evaluate non-equilibrium source versus converted formulas.
    bond = next(iter(result.bond_assignments.values())).parameter
    angle = next(iter(result.angle_assignments.values())).parameter
    source_bond = 100.0 * (1.7 - 1.5) ** 2 * 4.184
    converted_bond = 0.5 * bond.force_constant * (0.17 - bond.equilibrium_length) ** 2
    assert source_bond == pytest.approx(converted_bond)
    source_angle = 50.0 * ((120 - 109.5) * pi / 180) ** 2 * 4.184
    converted_angle = 0.5 * angle.force_constant * ((120 - angle.equilibrium_angle) * pi / 180) ** 2
    assert source_angle == pytest.approx(converted_angle)
    source_torsion = 1.5 * (1 + cos(2 * pi / 3 - pi)) * 4.184
    converted_torsion = improper.terms[0].force_constant * (1 + cos(2 * pi / 3 - pi))
    assert source_torsion == pytest.approx(converted_torsion)


def test_multiterm_proper_is_one_resolved_record(tmp_path):
    system, path, mapping = synthetic_source(tmp_path, multi=True)
    result = import_amber_prmtop(system, path, mapping, source="synthetic multi")
    parameter = next(iter(result.proper_torsion_assignments.values())).parameter
    assert len(parameter.terms) == 2
    assert [term.force_constant for term in parameter.terms] == pytest.approx([0, 0.75 * 4.184])
    assert len(result.source_14_pairs) == 1
    phi = 40 * pi / 180
    source_energy = (0.0 * (1 + cos(2 * phi - pi)) + 0.75 * (1 + cos(3 * phi))) * 4.184
    converted_energy = sum(
        term.force_constant * (1 + cos(term.periodicity * phi - term.phase * pi / 180))
        for term in parameter.terms
    )
    assert source_energy == pytest.approx(converted_energy)


def test_zero_lj_is_explicit_and_signed(tmp_path):
    system, path, mapping = synthetic_source(tmp_path, epsilon=0)
    result = import_amber_prmtop(system, path, mapping, source="zero LJ fixture")
    assert all(s.parameter.epsilon == 0 for s in result.site_assignments.values())
    assert all(s.parameter.sigma == 0 for s in result.site_assignments.values())
    assert result.is_compatible_with(system)
    changed = replace(result, source_sha256="another")
    assert not changed.is_compatible_with(system)


def test_on_disk_zero_lj_mixes_with_nonzero_type_and_keeps_charge(tmp_path):
    system, path, mapping = synthetic_source(tmp_path, zero_first=True)
    parsed = pmd.load_file(str(path))
    assert parsed.atoms[0].epsilon == parsed.atoms[0].rmin == 0
    original = system.to_dict()
    result = import_amber_prmtop(system, path, mapping, source="mixed zero LJ source")
    zero = result.site_assignments[mapping[0]].parameter
    regular = result.site_assignments[mapping[1]].parameter
    assert (zero.epsilon, zero.sigma) == (0, 0)
    assert regular.epsilon > 0 and regular.sigma > 0
    assert result.nonbonded_policy.mix_lj(zero, regular).epsilon == 0
    geometric = replace(result.nonbonded_policy, mixing_rule="geometric")
    assert geometric.mix_lj(zero, regular).sigma == 0
    assert geometric.mix_lj(zero, regular).epsilon == 0
    assert result.charge_result.assignments[mapping[0]].charge == pytest.approx(0.2)
    assert result.charge_result.assignments[mapping[1]].charge == pytest.approx(-0.2)
    selection = result.site_assignments[mapping[0]]
    changed_record = replace(selection.parameter, epsilon=0.1, sigma=0.1)
    changed_result = replace(
        result,
        site_assignments={
            **result.site_assignments,
            mapping[0]: replace(selection, parameter=changed_record),
        },
    )
    assert changed_result.content_signature() != result.result_signature
    assert not changed_result.is_compatible_with(system)
    assert system.to_dict() == original


def test_zero_lj_with_positive_epsilon_is_rejected(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path, epsilon=0)
    parm = pmd.load_file(str(path))
    parm.atoms[0].epsilon = 0.1
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    with pytest.raises(UnsupportedAmberFeatureError, match="Positive LJ epsilon"):
        import_amber_prmtop(system, path, mapping, source="invalid zero LJ")


def test_mapping_and_chemistry_rejections(tmp_path):
    system, path, mapping = synthetic_source(tmp_path)
    with pytest.raises(InvalidAmberMappingError, match="bijection"):
        import_amber_prmtop(system, path, {**mapping, 4: 10}, source="test")
    with pytest.raises(InvalidAmberMappingError, match="bijection"):
        import_amber_prmtop(system, path, {**mapping, 0: mapping[1]}, source="test")
    system.topology.sites[mapping[0]].atomic_number = 8
    with pytest.raises(InvalidAmberMappingError, match="element"):
        import_amber_prmtop(system, path, mapping, source="test")
    system.topology.sites[mapping[0]].atomic_number = 6
    system.topology.remove_bond(mapping[0], mapping[1])
    with pytest.raises(InvalidAmberMappingError, match="bonds"):
        import_amber_prmtop(system, path, mapping, source="test")


def test_unsupported_source_pair_override(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    parm = pmd.load_file(str(path))
    parm.parm_data["LENNARD_JONES_ACOEF"][0] *= 1.25
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    with pytest.raises(UnsupportedAmberFeatureError, match="pair override"):
        import_amber_prmtop(system, path, mapping, source="test")


def test_on_disk_unlike_pair_override_for_zero_lj_rejected(tmp_path):
    system, path, mapping = synthetic_source(tmp_path, zero_first=True)
    original = system.to_dict()
    parm = pmd.load_file(str(path))
    ntypes = parm.ptr("NTYPES")
    first, second = parm.atoms[0].nb_idx, parm.atoms[1].nb_idx
    index = parm.parm_data["NONBONDED_PARM_INDEX"][
        (first - 1) * ntypes + second - 1
    ] - 1
    parm.parm_data["LENNARD_JONES_ACOEF"][index] = 0.125
    # Write the explicit source coefficient array, bypassing ParmEd's
    # parameter regeneration, then verify the real on-disk parser sees it.
    pmd.amber.AmberFormat.write_parm(parm, str(path))
    assert pmd.load_file(str(path)).parm_data["LENNARD_JONES_ACOEF"][index] == pytest.approx(0.125)
    with pytest.raises(UnsupportedAmberFeatureError, match="unlike-pair override"):
        import_amber_prmtop(system, path, mapping, source="bad unlike pair")
    assert system.to_dict() == original


@pytest.mark.parametrize("coefficients,diagnostic", [
    ([0.5], "LENNARD_JONES_CCOEF"),
    ([float("nan")], "LENNARD_JONES_CCOEF"),
])
def test_on_disk_two_component_1264_rejected(tmp_path, coefficients, diagnostic):
    _, path, _ = synthetic_source(tmp_path)
    source = pmd.load_file(str(path))
    combined = pmd.amber.AmberParm.from_structure(source + source)
    combined.add_flag("LENNARD_JONES_CCOEF", "5E16.8", data=coefficients)
    combined.save(str(path), overwrite=True)
    reloaded = pmd.load_file(str(path))
    assert len(reloaded.atoms) == 8
    assert "LENNARD_JONES_CCOEF" in reloaded.parm_data
    graph = Topology()
    for index in range(8):
        graph.add_site(AtomSite(101 + 7 * index, f"C{index}", 12.01,
                                element="C", atomic_number=6))
    for offset in (0, 4):
        for left, right in ((0, 1), (1, 2), (2, 3)):
            graph.add_bond(101 + 7 * (offset + left),
                           101 + 7 * (offset + right), order=1)
    paired_system = MolecularSystem(graph, Coordinates())
    mapping = {i: 101 + 7 * i for i in range(8)}
    original_sites = tuple(sorted(paired_system.topology.sites))
    original_bonds = tuple(sorted(paired_system.topology.bonds))
    with pytest.raises(UnsupportedAmberFeatureError, match=diagnostic):
        import_amber_prmtop(paired_system, path, mapping, source="synthetic 12-6-4")
    assert tuple(sorted(paired_system.topology.sites)) == original_sites
    assert tuple(sorted(paired_system.topology.bonds)) == original_bonds
    assert len(paired_system.coordinates) == 0


def test_malformed_ccoef_file_has_source_flag_diagnostic(tmp_path):
    system, path, mapping = synthetic_source(tmp_path)
    source = pmd.load_file(str(path))
    source.add_flag("LENNARD_JONES_CCOEF", "5E16.8", data=[0.0, 0.0])
    source.save(str(path), overwrite=True)
    with pytest.raises(AmberImportError, match="LENNARD_JONES_CCOEF"):
        import_amber_prmtop(system, path, mapping, source="malformed CCOEF")


def test_all_zero_ccoef_is_inert_and_recorded(tmp_path):
    system, path, mapping = synthetic_source(tmp_path)
    source = pmd.load_file(str(path))
    source.add_flag("LENNARD_JONES_CCOEF", "5E16.8", data=[0.0])
    source.save(str(path), overwrite=True)
    assert pmd.load_file(str(path)).parm_data["LENNARD_JONES_CCOEF"] == [0.0]
    result = import_amber_prmtop(system, path, mapping, source="inert zero CCOEF")
    assert result.provenance["ignored_inert_zero_energy_flags"] == ["LENNARD_JONES_CCOEF"]


def test_unknown_lj_energy_extension_rejected(tmp_path):
    system, path, mapping = synthetic_source(tmp_path)
    source = pmd.load_file(str(path))
    source.add_flag("LENNARD_JONES_DCOEF", "5E16.8", data=[0.2])
    source.save(str(path), overwrite=True)
    with pytest.raises(UnsupportedAmberFeatureError, match="LENNARD_JONES_DCOEF"):
        import_amber_prmtop(system, path, mapping, source="unknown energy term")


def test_source_exclusions_are_checked(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    parm = pmd.load_file(str(path))
    parm.parm_data["EXCLUDED_ATOMS_LIST"][0] = 0
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    original = system.to_dict()
    with pytest.raises(UnsupportedAmberFeatureError, match="exclusions"):
        import_amber_prmtop(system, path, mapping, source="test")
    assert system.to_dict() == original


def test_source_charge_consistency_and_no_mutation(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    parm = pmd.load_file(str(path))
    parm.atoms[0].charge = 0.25
    parm.atoms[1].charge = -0.25
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    result = import_amber_prmtop(system, path, mapping, source="synthetic charges")
    assert result.charge_result.assignments[mapping[0]].charge == pytest.approx(0.25)
    original = system.to_dict()
    parm.atoms[1].charge = 0
    with pytest.raises(IncompleteChargeAssignmentError):
        import_amber_prmtop(system, path, mapping, source="bad charges")
    assert system.to_dict() == original


def test_unsupported_interaction_family(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    parm = pmd.load_file(str(path))
    parm.cmaps.append(object())
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    with pytest.raises(UnsupportedAmberFeatureError, match="cmaps"):
        import_amber_prmtop(system, path, mapping, source="unsupported CMAP")


def test_unsupported_source_scaling(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    parm = pmd.load_file(str(path))
    parm.dihedrals[0].type.scee = 0.5
    monkeypatch.setattr(pmd, "load_file", lambda _: parm)
    with pytest.raises(UnsupportedAmberFeatureError, match="scaling"):
        import_amber_prmtop(system, path, mapping, source="unsupported scaling")


def test_core_import_does_not_eagerly_import_parmed_or_rdkit():
    script = (
        "import sys; import island.forcefields.amber; "
        "assert 'parmed' not in sys.modules; assert 'rdkit' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", script], check=True)


def test_missing_optional_parser_is_actionable(tmp_path, monkeypatch):
    system, path, mapping = synthetic_source(tmp_path)
    original_import = builtins.__import__

    def missing_parmed(name, *args, **kwargs):
        if name == "parmed":
            raise ImportError("simulated optional dependency absence")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", missing_parmed)
    with pytest.raises(AmberImportError, match=r"island\[amber\]"):
        import_amber_prmtop(system, path, mapping, source="test")
