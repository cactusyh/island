"""Pinned external Amber phenol fixture; force-field identity is not established."""

import hashlib
from math import atan2, cos, pi
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.forcefields import ParameterizedSystem, import_amber_prmtop

FIXTURE = Path(__file__).parent / "fixtures/phenol_parmed_13239c2.prmtop"
SHA256 = "4722fe1f53d89e9576841b74c2ad3a18494a466be8fb0a2995cd8d956d3d0c28"
SOURCE = "ParmEd upstream test/files/phenol.prmtop at 13239c2516dd371317ffedeecb411864a5ccb365"
ELEMENTS = ("C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H", "H")
EDGES = (
    (0, 5), (0, 1), (1, 2), (2, 3), (3, 4), (3, 6), (4, 5),
    (0, 7), (1, 8), (2, 9), (4, 10), (5, 11), (6, 12),
)
IMPROPER_SOURCE_ORDER = (
    (2, 4, 3, 6), (7, 0, 5, 1), (0, 2, 1, 8),
    (1, 3, 2, 9), (3, 5, 4, 10), (0, 4, 5, 11),
)
MAP = {index: 101 + 7 * index for index in range(13)}


def phenol_system() -> MolecularSystem:
    """Reviewed phenol graph from chemical identity, independent of prmtop bonds."""
    graph = Topology()
    for index, element in enumerate(ELEMENTS):
        graph.add_site(AtomSite(
            MAP[index], f"{element}{index}", 12.01 if element == "C" else 15.999
            if element == "O" else 1.008,
            element=element, atomic_number={"C": 6, "O": 8, "H": 1}[element],
            formal_charge=0, metadata={"aromatic": index < 6},
        ))
    for left, right in EDGES:
        aromatic = left < 6 and right < 6
        graph.add_bond(MAP[left], MAP[right], order=1.5 if aromatic else 1,
                       aromatic=aromatic)
    return MolecularSystem(graph, Coordinates())


def ordered_dihedral(points: tuple[np.ndarray, ...]) -> float:
    """Signed angle from four ordered Cartesian positions in radians."""
    p0, p1, p2, p3 = points
    axis = p2 - p1
    axis /= np.linalg.norm(axis)
    first = p0 - p1
    second = p3 - p2
    first -= np.dot(first, axis) * axis
    second -= np.dot(second, axis) * axis
    return atan2(np.dot(np.cross(axis, first), second), np.dot(first, second))


def test_pinned_external_phenol_inventory_and_zero_lj():
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == SHA256
    source = pmd.load_file(str(FIXTURE))
    assert len(source.atoms) == 13
    assert len(source.bonds) == 13
    assert sum(term.improper for term in source.dihedrals) == 6
    assert source.atoms[12].type == "ho"
    assert source.atoms[12].epsilon == source.atoms[12].rmin == 0
    system = phenol_system()
    result = import_amber_prmtop(system, FIXTURE, MAP, source=SOURCE)
    assert result.is_compatible_with(system)
    assert dict(result.mapping) == MAP
    assert result.source_sha256 == SHA256
    assert len(result.site_assignments) == 13
    assert len(result.bond_assignments) == 13
    assert len(result.angle_assignments) == 19
    assert len(result.proper_torsion_assignments) == 26
    assert len(result.improper_assignments) == 6
    assert len(result.source_14_pairs) == 23
    assert len(result.source_exclusions) == 55  # distinct unordered pairs
    assert result.charge_result.total_assigned_charge == pytest.approx(0, abs=1e-8)
    assert result.charge_result.assignments[MAP[12]].charge == pytest.approx(0.4186)
    zero_site = result.site_assignments[MAP[12]].parameter
    carbon = result.site_assignments[MAP[0]].parameter
    assert (zero_site.epsilon, zero_site.sigma) == (0, 0)
    mixed = result.nonbonded_policy.mix_lj(zero_site, carbon)
    assert mixed.epsilon == 0
    assert mixed.sigma == pytest.approx(carbon.sigma / 2)
    assert result.nonbonded_policy.mix_lj(zero_site, zero_site).sigma == 0
    assert result.nonbonded_policy.coulomb_scale_14 == pytest.approx(1 / 1.2)
    assert result.nonbonded_policy.lj_scale_14 == pytest.approx(0.5)
    expected_impropers = tuple(tuple(MAP[i] for i in order)
                               for order in IMPROPER_SOURCE_ORDER)
    assert tuple(result.improper_assignments) == expected_impropers
    assert result.improper_assignments[expected_impropers[1]].parameter.central_atom_position == 2
    snapshot = ParameterizedSystem.from_amber_import(system, result)
    assert set(snapshot.interaction_assignments) == {
        "bond", "angle", "proper_torsion", "periodic_improper",
    }
    assert tuple(snapshot.interaction_assignments["periodic_improper"]) == expected_impropers
    assert snapshot.metadata["aggregate"]["production_validated"] is False


def test_external_source_to_converted_energies():
    source = pmd.load_file(str(FIXTURE))
    result = import_amber_prmtop(phenol_system(), FIXTURE, MAP, source=SOURCE)
    bond = source.bonds[0]
    bond_key = tuple(sorted((MAP[bond.atom1.idx], MAP[bond.atom2.idx])))
    converted_bond = result.bond_assignments[bond_key].parameter
    r_angstrom = 1.52
    expected = bond.type.k * (r_angstrom - bond.type.req) ** 2 * 4.184
    observed = 0.5 * converted_bond.force_constant * (
        r_angstrom / 10 - converted_bond.equilibrium_length
    ) ** 2
    assert observed == pytest.approx(expected)
    angle = source.angles[0]
    angle_key = tuple(MAP[a.idx] for a in (angle.atom1, angle.atom2, angle.atom3))
    angle_key = min(angle_key, angle_key[::-1])
    converted_angle = result.angle_assignments[angle_key].parameter
    theta = 113.0
    expected = angle.type.k * ((theta - angle.type.theteq) * pi / 180) ** 2 * 4.184
    observed = 0.5 * converted_angle.force_constant * (
        (theta - converted_angle.equilibrium_angle) * pi / 180
    ) ** 2
    assert observed == pytest.approx(expected)

    # Deliberately noncoplanar positions. Derive each torsion from its ordered
    # four Cartesian points, including the source improper with center at atom 2.
    geometry = (
        np.array((0.2, 0.1, -0.1)), np.array((1.1, 0.0, 0.3)),
        np.array((1.5, 1.0, 0.2)), np.array((2.4, 1.1, 1.2)),
    )
    for improper in (False, True):
        source_term = next(d for d in source.dihedrals if d.improper == improper
                           and (not improper or d.atom1.idx == 7))
        ordered = tuple(a.idx for a in (
            source_term.atom1, source_term.atom2, source_term.atom3, source_term.atom4,
        ))
        phi = ordered_dihedral(geometry)
        source_energy = source_term.type.phi_k * (
            1 + cos(source_term.type.per * phi - source_term.type.phase * pi / 180)
        ) * 4.184
        key = tuple(MAP[i] for i in ordered)
        if not improper:
            key = min(key, key[::-1])
            record = result.proper_torsion_assignments[key].parameter
        else:
            record = result.improper_assignments[key].parameter
            assert ordered == (7, 0, 5, 1)
            assert record.central_atom_position == 2
        converted_energy = sum(
            term.force_constant * (
                1 + cos(term.periodicity * phi - term.phase * pi / 180)
            ) for term in record.terms
        )
        assert converted_energy == pytest.approx(source_energy)

    # Non-excluded type pairs: source A/B and decoded charges are independent
    # of importer output. Distances are deliberately away from LJ minimum.
    first, second, distance_a = 0, 6, 4.3
    types = (source.atoms[first].nb_idx, source.atoms[second].nb_idx)
    ntypes = source.ptr("NTYPES")
    pair_index = source.parm_data["NONBONDED_PARM_INDEX"][
        (types[0] - 1) * ntypes + types[1] - 1
    ] - 1
    a = source.parm_data["LENNARD_JONES_ACOEF"][pair_index]
    b = source.parm_data["LENNARD_JONES_BCOEF"][pair_index]
    expected_lj = (a / distance_a ** 12 - b / distance_a ** 6) * 4.184
    mixed = result.nonbonded_policy.mix_lj(
        result.site_assignments[MAP[first]].parameter,
        result.site_assignments[MAP[second]].parameter,
    )
    observed_lj = 4 * mixed.epsilon * (
        (mixed.sigma / (distance_a / 10)) ** 12
        - (mixed.sigma / (distance_a / 10)) ** 6
    )
    assert observed_lj == pytest.approx(expected_lj, abs=1e-5)
    zero_type = source.atoms[12].nb_idx
    assert all(source.parm_data["LENNARD_JONES_ACOEF"][
        source.parm_data["NONBONDED_PARM_INDEX"][(zero_type - 1) * ntypes + i] - 1
    ] == 0 for i in range(ntypes))
    assert all(source.parm_data["LENNARD_JONES_BCOEF"][
        source.parm_data["NONBONDED_PARM_INDEX"][(zero_type - 1) * ntypes + i] - 1
    ] == 0 for i in range(ntypes))
    # Amber's Coulomb coefficient is kcal*angstrom/(mol*e^2).
    q1, q2 = source.atoms[6].charge, source.atoms[7].charge
    expected_coul = 332.06371 * q1 * q2 / 5.0 * 4.184
    observed_coul = 138.935455864 * (
        result.charge_result.assignments[MAP[6]].charge
        * result.charge_result.assignments[MAP[7]].charge
    ) / 0.5
    assert observed_coul == pytest.approx(expected_coul, rel=1e-6)


def test_family_keys_match_native_parameter_snapshot():
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer
    from island.forcefields import (
        ParameterAssignmentEngine,
        RDKitSmartsAtomTypingEngine,
        island_demo_parameters_v1,
        island_demo_v1_ruleset,
    )

    native_system = build_linear_polymer("[*]CC[*]", dp=2, generate_3d=False)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(native_system, ruleset)
    native = ParameterAssignmentEngine().assign(
        native_system, typing, ruleset, island_demo_parameters_v1()
    ).to_parameterized_system(native_system)
    imported_system = phenol_system()
    imported_result = import_amber_prmtop(imported_system, FIXTURE, MAP, source=SOURCE)
    direct = imported_result.to_parameterized_system(imported_system)
    wrapped = ParameterizedSystem.from_amber_import(imported_system, imported_result)
    for family in ("bond", "angle", "proper_torsion"):
        assert family in native.interaction_assignments
        assert family in direct.interaction_assignments
        assert family in wrapped.interaction_assignments
    assert "periodic_improper" in direct.interaction_assignments
    assert "improper_assignments" not in direct.interaction_assignments
    assert direct.interaction_assignments == wrapped.interaction_assignments
    # Both entry points own their data: subsequent caller mutation cannot
    # change either validated snapshot.
    imported_system.topology.sites[MAP[0]].name = "changed"
    assert direct.system.topology.sites[MAP[0]].name != "changed"
    assert wrapped.system.topology.sites[MAP[0]].name != "changed"
