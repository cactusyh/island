"""Whole-system reference and analytical acceptance, no AmberTools/network needed."""

import copy
from dataclasses import replace
from math import cos, pi

import numpy as np
import pytest

mm = pytest.importorskip("openmm")
pmd = pytest.importorskip("parmed")
from openmm import app, unit
from test_amber_external_reference import FIXTURE, MAP, ordered_dihedral, phenol_system
from test_amber_import import synthetic_source
from test_ambertools_integrity import _valid_result

from island import Coordinates
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    UnsupportedEvaluationError,
)
from island.forcefields import import_amber_prmtop

XYZ = np.array(((0.2, 0.1, -0.1), (1.6, 0.0, 0.3), (2.4, 1.3, 0.2), (3.6, 1.5, 1.2)))


def source_case(tmp_path, *, improper=False, zero_first=False, reverse=False):
    system, path, mapping = synthetic_source(
        tmp_path,
        improper=improper,
        multi=not improper,
        zero_first=zero_first,
        reverse=reverse,
    )
    source = pmd.load_file(str(path))
    # Independent original Amber data: nonzero charges, arbitrary phases,
    # multiple torsion terms; zero-LJ carbon still participates in Coulomb.
    for atom in source.atoms:
        atom.charge = {10: 0.2, 20: -0.1, 30: -0.3, 40: 0.2}[mapping[atom.idx]]
    for i, term in enumerate(source.dihedral_types):
        term.phi_k = 0.75 + i * 0.3
        term.phase = 37 + i * 41
    source.remake_parm()
    source.save(str(path), overwrite=True)
    system.coordinates = Coordinates(dict(zip((10, 20, 30, 40), XYZ, strict=True)))
    imported = import_amber_prmtop(
        system, path, mapping, source="analytical synthetic Amber"
    )
    return system, imported, path, mapping


def reference(path, mapping, coordinates):
    """Independent path B: OpenMM's reader, never converted ISLAND parameters."""
    model = app.AmberPrmtopFile(str(path)).createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=None,
        removeCMMotion=False,
        hydrogenMass=None,
    )
    groups = {
        "HarmonicBondForce": 0,
        "HarmonicAngleForce": 1,
        "PeriodicTorsionForce": 2,
        "NonbondedForce": 3,
    }
    for force in model.getForces():
        force.setForceGroup(groups[type(force).__name__])
        if isinstance(force, mm.NonbondedForce):
            force.setUseDispersionCorrection(False)
            force.setUseSwitchingFunction(False)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        positions = np.array([coordinates[mapping[i]] for i in range(len(mapping))])
        context.setPositions(positions * 0.1 * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = (
            state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer
            )
            * 0.1
        )
        components = [
            context.getState(getEnergy=True, groups=1 << i)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
            for i in range(4)
        ]
        return energy, {mapping[i]: forces[i] for i in range(len(mapping))}, components
    finally:
        del context, integrator


def assert_reference(evaluator, path, mapping, coords):
    result = evaluator.evaluate(coords)
    energy, forces, components = reference(path, mapping, coords)
    assert result.potential_energy == pytest.approx(energy, rel=2e-10, abs=2e-7)
    observed = result.energy_components
    assert [
        observed["bond"],
        observed["angle"],
        observed["proper_torsion"] + observed["periodic_improper"],
        observed["nonbonded"],
    ] == pytest.approx(components, rel=2e-10, abs=2e-7)
    for site in coords:
        np.testing.assert_allclose(
            result.forces[site], forces[site], rtol=2e-9, atol=2e-6
        )
    assert sum(result.energy_components.values()) == pytest.approx(
        result.potential_energy, abs=1e-8
    )
    return result


@pytest.mark.parametrize(
    "improper,zero_first,reverse",
    [
        (False, False, False),
        (False, True, False),
        (True, False, False),
        (False, True, True),
    ],
)
def test_whole_system_synthetic(tmp_path, improper, zero_first, reverse):
    system, imported, path, mapping = source_case(
        tmp_path,
        improper=improper,
        zero_first=zero_first,
        reverse=reverse,
    )
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    rng = np.random.default_rng(401)
    for _ in range(3):
        xyz = XYZ + rng.normal(0, 0.035, XYZ.shape)
        assert_reference(
            evaluator, path, mapping, dict(zip((10, 20, 30, 40), xyz, strict=True))
        )


def phenol_coordinates():
    # Source atom order follows the reviewed aromatic ring graph. This is a
    # deterministic geometry, not an optimized reference or a charge calculation.
    ring = np.array(
        [(1.4 * cos(i * pi / 3), 1.4 * np.sin(i * pi / 3), 0) for i in range(6)]
    )
    xyz = list(ring)
    xyz.append(ring[3] * 1.96)
    xyz.extend(ring[i] * (2.48 / 1.4) for i in (0, 1, 2, 4, 5))
    xyz.append(np.array((-3.1, 0.8, 0.2)))
    return np.array(xyz)


def test_whole_phenol_ring_charges_impropers_and_zero_lj():
    system = phenol_system()
    system.coordinates = Coordinates(
        {MAP[i]: x for i, x in enumerate(phenol_coordinates())}
    )
    imported = import_amber_prmtop(
        system,
        FIXTURE,
        MAP,
        source="upstream phenol: exact FF/charge provenance unknown",
    )
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    rng = np.random.default_rng(2026)
    for _ in range(3):
        xyz = phenol_coordinates() + rng.normal(0, 0.04, (13, 3))
        result = assert_reference(
            evaluator, FIXTURE, MAP, {MAP[i]: x for i, x in enumerate(xyz)}
        )
        assert result.energy_components["periodic_improper"] != 0


@pytest.mark.parametrize(
    "improper,zero_first", [(False, False), (False, True), (True, False)]
)
def test_analytical_component_formulas(tmp_path, improper, zero_first):
    system, imported, path, mapping = source_case(
        tmp_path, improper=improper, zero_first=zero_first
    )
    source = pmd.load_file(str(path))
    xyz = {site: system.coordinates.get(site) for site in mapping.values()}
    pos = lambda atom: xyz[mapping[atom.idx]]
    expected_bond = sum(
        b.type.k
        * (np.linalg.norm(pos(b.atom1) - pos(b.atom2)) - b.type.req) ** 2
        * 4.184
        for b in source.bonds
    )
    expected_angle = 0
    for a in source.angles:
        u, v = pos(a.atom1) - pos(a.atom2), pos(a.atom3) - pos(a.atom2)
        theta = np.arccos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
        expected_angle += a.type.k * (theta - a.type.theteq * pi / 180) ** 2 * 4.184
    expected_torsion = sum(
        d.type.phi_k
        * (
            1
            + cos(
                d.type.per
                * ordered_dihedral(
                    tuple(pos(a) for a in (d.atom1, d.atom2, d.atom3, d.atom4))
                )
                - d.type.phase * pi / 180
            )
        )
        * 4.184
        for d in source.dihedrals
    )
    result = OpenMMSinglePointEvaluator(system, imported).evaluate()
    assert result.energy_components["bond"] == pytest.approx(expected_bond, rel=1e-10)
    assert result.energy_components["angle"] == pytest.approx(expected_angle, rel=1e-10)
    assert result.energy_components[
        "periodic_improper" if improper else "proper_torsion"
    ] == pytest.approx(expected_torsion, abs=1e-9)
    # Four-site chain has exactly one 1-4 pair, independent of torsion term count.
    if improper:
        assert result.energy_components["nonbonded"] == 0
    else:
        left, right = sorted(source.atoms, key=lambda a: mapping[a.idx])[::3]
        distance = np.linalg.norm(pos(left) - pos(right))
        expected_lj = (
            np.sqrt(left.epsilon * right.epsilon)
            * (
                ((left.rmin + right.rmin) / distance) ** 12
                - 2 * ((left.rmin + right.rmin) / distance) ** 6
            )
            * 4.184
            / 2
        )
        expected_coulomb = (
            138.93545764438198 * left.charge * right.charge / (distance * 0.1) / 1.2
        )
        assert result.energy_components["nonbonded"] == pytest.approx(
            expected_lj + expected_coulomb, abs=2e-7
        )
        if zero_first:
            assert expected_lj == 0 and expected_coulomb != 0


def test_force_finite_difference_rigid_motion_and_owned_frames(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    coords = {i: system.coordinates.get(i) for i in (10, 20, 30, 40)}
    original = copy.deepcopy(system.to_dict())
    result = evaluator.evaluate(coords)
    for site, axis in ((10, 0), (20, 1), (40, 2)):
        estimates = []
        for step in (1e-4, 5e-5):
            left, right = copy.deepcopy(coords), copy.deepcopy(coords)
            left[site][axis] -= step
            right[site][axis] += step
            derivative = -(
                evaluator.evaluate(right).potential_energy
                - evaluator.evaluate(left).potential_energy
            ) / (2 * step)
            estimates.append(derivative)
            assert derivative == pytest.approx(
                result.forces[site][axis], rel=2e-6, abs=2e-5
            )
        assert estimates[1] == pytest.approx(estimates[0], rel=2e-6, abs=2e-5)
    rotation = np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
    moved = {i: rotation @ coords[i] + (3, -4, 5) for i in reversed(coords)}
    transformed = evaluator.evaluate(moved)
    assert transformed.potential_energy == pytest.approx(
        result.potential_energy, abs=1e-8
    )
    for site in coords:
        np.testing.assert_allclose(
            transformed.forces[site], rotation @ result.forces[site], atol=1e-8
        )
    assert transformed.coordinate_fingerprint != result.coordinate_fingerprint
    assert transformed.model_fingerprint == result.model_fingerprint
    assert result.force_unit == "kJ/(mol*angstrom)"
    np.testing.assert_allclose(
        np.sum(list(result.forces.values()), axis=0), 0, atol=1e-8
    )
    coords[10][0] += 0.1
    assert (
        evaluator.evaluate(coords).evaluation_fingerprint
        != result.evaluation_fingerprint
    )
    assert system.to_dict() == original
    with pytest.raises(TypeError):
        result.forces[10] = (0, 0, 0)
    assert copy.deepcopy(result) == result
    system.coordinates.set(10, (90, 80, 70))
    assert evaluator.evaluate().potential_energy == result.potential_energy


def test_parameter_content_and_mutable_snapshot_validation(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    snapshot = imported.to_parameterized_system(system)
    assert (
        OpenMMSinglePointEvaluator.from_parameterized_system(snapshot, imported)
        .evaluate()
        .potential_energy
    )
    key = next(iter(snapshot.interaction_assignments["bond"]))
    del snapshot.interaction_assignments["bond"][key]
    with pytest.raises(EvaluationInputError, match="Snapshot"):
        OpenMMSinglePointEvaluator.from_parameterized_system(snapshot, imported)
    selected = imported.bond_assignments[key]
    changed = replace(
        imported,
        bond_assignments={
            **imported.bond_assignments,
            key: replace(
                selected, parameter=replace(selected.parameter, force_constant=123)
            ),
        },
    )
    with pytest.raises(EvaluationInputError):
        OpenMMSinglePointEvaluator(system, changed)
    # Malformed numerics with a freshly recomputed hash must still fail binding.
    bad = copy.deepcopy(imported)
    object.__setattr__(bad.bond_assignments[key].parameter, "force_constant", -1)
    bad = replace(bad, result_signature=bad.content_signature())
    with pytest.raises(EvaluationInputError):
        OpenMMSinglePointEvaluator(system, bad)
    geometric = replace(
        imported,
        nonbonded_policy=replace(imported.nonbonded_policy, mixing_rule="geometric"),
    )
    geometric = replace(geometric, result_signature=geometric.content_signature())
    with pytest.raises(UnsupportedEvaluationError, match="mixing"):
        OpenMMSinglePointEvaluator(system, geometric)
    with pytest.raises(UnsupportedEvaluationError, match="nonperiodic"):
        OpenMMSinglePointEvaluator(system, imported, periodic=True)


@pytest.mark.parametrize(
    "bad",
    [
        {10: (0, 0, 0)},
        {10: (0, 0, 0), 20: (1, 2, 3), 30: (2, 3, 4), 40: (float("nan"), 0, 0)},
        {10: (0, 0), 20: (1, 2, 3), 30: (2, 3, 4), 40: (3, 4, 5)},
        {10: ("0", "0", "0"), 20: (1, 2, 3), 30: (2, 3, 4), 40: (3, 4, 5)},
        np.zeros((4, 3)),
    ],
)
def test_bad_coordinates(tmp_path, bad):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    with pytest.raises(EvaluationInputError):
        evaluator.evaluate(bad)


def test_singular_and_changed_graph(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    with pytest.raises(EvaluationInputError, match="coincident"):
        evaluator.evaluate({site: (0, 0, 0) for site in (10, 20, 30, 40)})
    with pytest.raises(EvaluationInputError, match="collinear"):
        evaluator.evaluate({site: (site, 0, 0) for site in (10, 20, 30, 40)})
    with pytest.raises(EvaluationInputError, match="angstrom"):
        evaluator.evaluate(coordinate_unit="nm")
    system.topology.add_bond(10, 40)
    with pytest.raises(EvaluationInputError):
        evaluator.evaluate_system(system)
    # Extreme finite coordinates must not become a successful NaN result.
    with pytest.raises(EvaluationError):
        evaluator.evaluate(
            {site: (site * 1e300, site * 1e300, 1) for site in (10, 20, 30, 40)}
        )


def test_preparation_history_separate_from_frame(tmp_path):
    system, preparation = _valid_result(tmp_path)
    original = copy.deepcopy(dict(preparation.record))
    evaluator = OpenMMSinglePointEvaluator.from_preparation(system, preparation)
    before = evaluator.evaluate()
    system.coordinates.set(101, (0.1, 0.2, 0.3))
    after = evaluator.evaluate_system(system)
    assert before.coordinate_fingerprint != after.coordinate_fingerprint
    assert before.model_fingerprint == after.model_fingerprint
    assert dict(preparation.record) == original
    from island.exceptions import InvalidAmberImportResultError

    with pytest.raises(InvalidAmberImportResultError, match="coordinates changed"):
        OpenMMSinglePointEvaluator.from_preparation(system, preparation)


@pytest.mark.parametrize("zero_lj", [False, True])
def test_isolated_pair_separate_lj_coulomb_and_analytic_force(tmp_path, zero_lj):
    from island import AtomSite, MolecularSystem, Topology

    structure = pmd.Structure()
    graph = Topology()
    for i, charge in enumerate((1.0, -1.0)):
        atom_type = pmd.AtomType(f"C{i}", i + 1, 12.01, 6)
        atom_type.set_lj_params(
            0 if zero_lj and i == 0 else 0.1 * (i + 1),
            0 if zero_lj and i == 0 else 1.8 + 0.3 * i,
        )
        atom = pmd.Atom(
            name=f"C{i}", type=f"C{i}", atomic_number=6, mass=12.01, charge=charge
        )
        atom.atom_type = atom_type
        structure.add_atom(atom, "PAIR", 1)
        graph.add_site(
            AtomSite(
                91 - i * 74,
                f"C{i}",
                12.01,
                element="C",
                atomic_number=6,
                formal_charge=int(charge),
            )
        )
    source = pmd.amber.AmberParm.from_structure(structure)
    path = tmp_path / "pair.prmtop"
    source.save(str(path))
    system = MolecularSystem(graph, Coordinates({91: (0, 0, 0), 17: (4.2, 0, 0)}))
    imported = import_amber_prmtop(
        system, path, {0: 91, 1: 17}, source="independent two-site example"
    )
    result = OpenMMSinglePointEvaluator(system, imported).evaluate()
    assert_reference(
        OpenMMSinglePointEvaluator(system, imported),
        path,
        {0: 91, 1: 17},
        {91: (0, 0, 0), 17: (4.2, 0, 0)},
    )
    r = 4.2
    epsilon = 0 if zero_lj else np.sqrt(0.1 * 0.2) * 4.184
    radius = 3.9
    lj = epsilon * ((radius / r) ** 12 - 2 * (radius / r) ** 6)
    coulomb = 138.93545764438198 * (-1.0) / (r * 0.1)
    assert result.potential_energy == pytest.approx(lj + coulomb, abs=2e-7)
    derivative_lj = epsilon * (-12 * (radius / r) ** 12 + 12 * (radius / r) ** 6) / r
    derivative_coulomb = -coulomb / r
    assert result.forces[91][0] == pytest.approx(
        derivative_lj + derivative_coulomb, abs=2e-7
    )
    assert result.forces[17][0] == pytest.approx(
        -derivative_lj - derivative_coulomb, abs=2e-7
    )
    # A second independently parsed source has no electrostatics: isolate LJ
    # rather than allowing errors in LJ and Coulomb to cancel in their sum.
    for atom in source.atoms:
        atom.charge = 0
    source.remake_parm()
    source.save(str(path), overwrite=True)
    for site in system.topology.sites.values():
        site.formal_charge = 0
    uncharged = import_amber_prmtop(
        system, path, {0: 91, 1: 17}, source="uncharged analytical pair"
    )
    lj_result = OpenMMSinglePointEvaluator(system, uncharged).evaluate()
    assert lj_result.potential_energy == pytest.approx(lj, abs=1e-9)
    assert result.potential_energy - lj_result.potential_energy == pytest.approx(
        coulomb, abs=2e-7
    )


def test_ordered_improper_with_center_second_and_nonzero_phase(tmp_path):
    system, _, path, mapping = source_case(tmp_path, improper=True)
    source = pmd.load_file(str(path))
    term = source.dihedrals[0]
    term.atom2, term.atom3 = term.atom3, term.atom2
    source.remake_parm()
    source.save(str(path), overwrite=True)
    imported = import_amber_prmtop(
        system, path, mapping, source="center-second ordered improper"
    )
    selection = next(iter(imported.improper_assignments.values()))
    assert selection.parameter.central_atom_position == 2
    coords = {site: system.coordinates.get(site) for site in mapping.values()}
    result = assert_reference(
        OpenMMSinglePointEvaluator(system, imported), path, mapping, coords
    )
    phi = ordered_dihedral(
        tuple(
            coords[mapping[a.idx]]
            for a in (term.atom1, term.atom2, term.atom3, term.atom4)
        )
    )
    expected = (
        term.type.phi_k
        * (1 + cos(term.type.per * phi - term.type.phase * pi / 180))
        * 4.184
    )
    assert result.energy_components["periodic_improper"] == pytest.approx(
        expected, abs=1e-9
    )


def test_insertion_order_independence_and_model_ownership(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    expected = evaluator.evaluate()
    reordered = copy.deepcopy(system)
    reordered.topology.sites = dict(reversed(list(reordered.topology.sites.items())))
    reordered.topology.bonds = dict(reversed(list(reordered.topology.bonds.items())))
    actual = OpenMMSinglePointEvaluator(reordered, imported).evaluate()
    assert actual == expected
    parameter = next(iter(imported.bond_assignments.values())).parameter
    object.__setattr__(parameter, "force_constant", 123456)
    assert evaluator.evaluate() == expected


def test_unsupported_snapshot_family_and_pair_overrides(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    snapshot = imported.to_parameterized_system(system)
    snapshot.interaction_assignments["urey_bradley"] = {}
    with pytest.raises(EvaluationInputError, match="interaction_assignments"):
        OpenMMSinglePointEvaluator.from_parameterized_system(snapshot, imported)
    system.pair_overrides = {(10, 40): (0.1, 0.2)}
    with pytest.raises(UnsupportedEvaluationError, match="pair overrides"):
        OpenMMSinglePointEvaluator(system, imported)


@pytest.mark.parametrize(
    "failure", ["energy_nan", "energy_inf", "force_nan", "force_shape"]
)
def test_invalid_backend_outputs_are_domain_errors(tmp_path, monkeypatch, failure):
    from types import SimpleNamespace

    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    original = mm.Context

    class ContextWithBadOutput:
        def __init__(self, *args):
            self.context = original(*args)

        def setPositions(self, positions):
            self.context.setPositions(positions)

        def getPlatform(self):
            return self.context.getPlatform()

        def getState(self, **kwargs):
            state = self.context.getState(**kwargs)

            def energy():
                if failure.startswith("energy"):
                    return (
                        float("nan") if failure.endswith("nan") else float("inf")
                    ) * unit.kilojoule_per_mole
                return state.getPotentialEnergy()

            def forces(**kwargs):
                shape = (4, 2) if failure == "force_shape" else (4, 3)
                values = np.full(shape, float("nan") if failure == "force_nan" else 0)
                return values * unit.kilojoule_per_mole / unit.nanometer

            return SimpleNamespace(getPotentialEnergy=energy, getForces=forces)

    monkeypatch.setattr(mm, "Context", ContextWithBadOutput)
    with pytest.raises(EvaluationError, match="nonfinite"):
        evaluator.evaluate()


@pytest.mark.parametrize("label", ["cip_label", "chiral_tag"])
def test_changed_stereochemical_graph_requires_rebinding(tmp_path, label):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    system.topology.sites[20].metadata[label] = (
        "R" if label == "cip_label" else "CHI_TETRAHEDRAL_CW"
    )
    with pytest.raises(EvaluationInputError, match="graph changed"):
        evaluator.evaluate_system(system)


def test_binding_and_evaluation_do_not_need_parmed_or_rdkit(tmp_path, monkeypatch):
    import builtins

    system, imported, _, _ = source_case(tmp_path)
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"parmed", "rdkit"}:
            raise ModuleNotFoundError(name)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    result = OpenMMSinglePointEvaluator(system, imported).evaluate()
    assert np.isfinite(result.potential_energy)
