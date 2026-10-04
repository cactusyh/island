"""Synthetic software tests; retained LAMMPS isolated terms are a separate oracle."""

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_model import parameters, policy  # noqa: F401

from island.evaluation import PCFFSinglePointEvaluator
from island.evaluation.pcff import build_model
from island.exceptions import EvaluationError, EvaluationInputError
from island.forcefields.pcff import define_pcff_model
from island.forcefields.pcff.terms import Class2Term, finite_difference_forces

mm = pytest.importorskip("openmm")
from openmm import unit

FIXTURES = json.loads(
    (Path(__file__).parent / "fixtures/pcff_class2_terms.json").read_text()
)["fixtures"]


def compiled(data, xyz):
    model = build_model(data, {i: 12 for i in range(len(xyz))}, mm)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(np.array(xyz) * 0.1)
    state = context.getState(getEnergy=True, getForces=True)
    return state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    ), state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    ) * 0.1


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda f: f["family"])
def test_all_compiled_terms_against_executed_lammps(fixture):
    data = {
        "terms": [
            {k: fixture[k] for k in ("family", "sites", "coefficients", "equilibria")}
        ],
        "nonbonded": [
            {"site": s, "rmin": 1, "epsilon": 0, "charge": 0} for s in fixture["sites"]
        ],
        "special_pair_inventory": [],
    }
    energy, forces = compiled(data, fixture["xyz_angstrom"])
    assert energy == pytest.approx(fixture["energy_kcal_mol"] * 4.184, abs=2e-10)
    np.testing.assert_allclose(
        forces,
        np.array(fixture["forces_kcal_mol_angstrom"]) * 4.184,
        atol=2e-9,
        rtol=2e-10,
    )


def test_nonzero_phase_and_fd_convergence():
    xyz = np.array([[0.1, 1.3, 0.2], [0, 0, 0], [1.5, 0.2, 0.1], [2, 1.1, 1.4]])
    term = Class2Term("torsion_3", (0, 1, 2, 3), (1.2, 0.3, -0.7, -0.5, 0.4, 0.9))
    data = {
        "terms": [
            {
                "family": term.family,
                "sites": list(term.sites),
                "coefficients": list(term.coefficients),
                "equilibria": [],
            }
        ],
        "nonbonded": [
            {"site": s, "rmin": 1, "epsilon": 0, "charge": 0} for s in term.sites
        ],
        "special_pair_inventory": [],
    }
    energy, forces = compiled(data, xyz)
    coords = dict(enumerate(xyz))
    assert energy == pytest.approx(term.energy(coords), abs=1e-12)
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        fd = finite_difference_forces(term, coords, displacement=h)
        errors.append(np.max(np.abs(np.array(list(fd.values())) - forces)))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < 1e-8


def test_pair_weights_zero_lj_and_mixing():
    data = {
        "terms": [],
        "nonbonded": [
            {"site": 0, "rmin": 2, "epsilon": 0.3, "charge": 0.2},
            {"site": 1, "rmin": 4, "epsilon": 0.7, "charge": -0.2},
        ],
        "special_pair_inventory": [
            {"sites": [0, 1], "lj_weight": 0.3, "coulomb_weight": 0.7}
        ],
    }
    r = 3.1
    rmin = ((2**6 + 4**6) / 2) ** (1 / 6)
    eps = 2 * np.sqrt(0.3 * 0.7) * 2**3 * 4**3 / (2**6 + 4**6)
    expected = (
        0.3 * eps * (2 * (rmin / r) ** 9 - 3 * (rmin / r) ** 6)
        + 0.7 * 1389.35456264 * 0.2 * (-0.2) / r
    )
    e, f = compiled(data, [[0, 0, 0], [r, 0, 0]])
    assert e == pytest.approx(expected, abs=1e-11)
    data["nonbonded"][0]["epsilon"] = 0
    e, f = compiled(data, [[0, 0, 0], [r, 0, 0]])
    assert e == pytest.approx(0.7 * 1389.35456264 * 0.2 * (-0.2) / r)
    assert f[0, 0] > 0


@pytest.fixture
def bound(parameters):  # noqa: F811
    system = explicit([(0, 1)], [3, 3])
    xyz = [
        [0, 0, 0],
        [1.5, 0, 0],
        [-0.4, 1, 0],
        [-0.4, -0.5, 0.87],
        [-0.4, -0.5, -0.87],
        [1.9, -1, 0],
        [1.9, 0.5, 0.87],
        [1.9, 0.5, -0.87],
    ]
    for s, v in zip(sorted(system.topology.sites), xyz):
        system.coordinates.set(s, v)
    spec = define_pcff_model(parameters, special_pairs=policy())
    return system, spec, PCFFSinglePointEvaluator(system, spec)


def test_frames_binding_ownership_and_covariance(bound):
    system, _, evaluator = bound
    first = evaluator.evaluate()
    coords = {s: system.coordinates.get(s) for s in system.topology.sites}
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    moved = {s: v @ rotation + [3, -4, 7] for s, v in reversed(list(coords.items()))}
    second = evaluator.evaluate_fresh(moved)
    assert second.potential_energy == pytest.approx(first.potential_energy, abs=1e-10)
    for s in coords:
        np.testing.assert_allclose(
            second.forces[s], np.array(first.forces[s]) @ rotation, atol=1e-10
        )
    assert sum(first.energy_components.values()) == pytest.approx(
        first.potential_energy
    )
    for field, value in (("mass", 13), ("element", "F")):
        changed = deepcopy(system)
        setattr(changed.topology.sites[11], field, value)
        with pytest.raises(EvaluationInputError):
            evaluator.validate_system(changed)
    system.coordinates.translate([0.2, 0, 0])
    evaluator.validate_system(system)
    assert evaluator.evaluate() == first
    assert second.model_fingerprint == first.model_fingerprint
    assert not first.production_validated


@pytest.mark.parametrize(
    "bad", [True, np.bool_(False), float("nan"), float("inf"), 10**400, 1e308]
)
def test_invalid_inputs_before_context(bound, monkeypatch, bad):
    system, _, evaluator = bound
    coords = {s: list(system.coordinates.get(s)) for s in system.topology.sites}
    coords[11][0] = bad

    def forbidden(*args, **kwargs):
        pytest.fail("Context constructed for invalid input")

    monkeypatch.setattr(mm, "Context", forbidden)
    with pytest.raises(EvaluationInputError):
        evaluator.evaluate(coords)


def test_fresh_count_and_backend_failure(bound, monkeypatch):
    _, _, evaluator = bound
    original = mm.Context
    calls = []

    def create(*args):
        calls.append(1)
        return original(*args)

    monkeypatch.setattr(mm, "Context", create)
    evaluator.evaluate()
    evaluator.evaluate_fresh()
    assert len(calls) == 2

    def fail(*args):
        raise RuntimeError("injected backend failure")

    monkeypatch.setattr(mm, "Context", fail)
    with pytest.raises(EvaluationError, match="injected backend failure"):
        evaluator.evaluate()


def test_wrong_units_coverage_singular_and_model(bound):
    system, spec, evaluator = bound
    from island.charge_references.records import pack, unpack
    from island.forcefields.pcff import PCFFModelSpecification

    coords = {s: system.coordinates.get(s) for s in system.topology.sites}
    for kwargs in (
        {"coordinate_unit": "nanometer"},
        {"coordinates": {}},
        {"coordinates": {s: [0, 0, 0] for s in coords}},
        {"coordinates": {**coords, 11: [0, 0]}},
    ):
        with pytest.raises(EvaluationInputError):
            evaluator.evaluate(**kwargs)
    data = unpack(spec.json_text)
    data["terms"][0]["coefficients"][1] += 1
    with pytest.raises(EvaluationInputError):
        PCFFSinglePointEvaluator(
            system, PCFFModelSpecification(pack(data), spec.assignment)
        )
    altered = deepcopy(system)
    altered.topology.sites[11].metadata["cip_label"] = "R"
    with pytest.raises(EvaluationInputError):
        evaluator.validate_system(altered)
    altered = deepcopy(system)
    altered.topology.bonds.pop(next(iter(altered.topology.bonds)))
    with pytest.raises(EvaluationInputError):
        evaluator.validate_system(altered)


def test_binding_remapping_and_whole_force_fd(parameters):  # noqa: F811
    from island import AtomSite, Coordinates, MolecularSystem, Topology
    from island.forcefields.pcff import (
        assign_automatic_pcff_charges,
        assign_pcff_parameters,
        type_pcff_atoms,
    )

    original = explicit([(0, 1)], [3, 3])
    xyz = [
        [0, 0, 0],
        [1.5, 0.1, 0.2],
        [-0.4, 1, 0],
        [-0.4, -0.5, 0.87],
        [-0.4, -0.5, -0.87],
        [1.9, -1, 0],
        [1.9, 0.5, 0.87],
        [1.9, 0.5, -0.87],
    ]
    ids = sorted(original.topology.sites)
    for s, v in zip(ids, xyz):
        original.coordinates.set(s, v)
    spec = define_pcff_model(parameters, special_pairs=policy())
    evaluator = PCFFSinglePointEvaluator(original, spec)
    result = evaluator.evaluate()
    remap = {s: 1000 - 7 * i for i, s in enumerate(ids)}
    topology = Topology()
    for s in reversed(ids):
        a = original.topology.sites[s]
        topology.add_site(
            AtomSite(
                remap[s],
                a.name,
                a.mass,
                element=a.element,
                atomic_number=a.atomic_number,
            )
        )
    for b in reversed(list(original.topology.bonds.values())):
        topology.add_bond(remap[b.site2], remap[b.site1], order=1)
    remapped = MolecularSystem(
        topology, Coordinates({remap[s]: original.coordinates.get(s) for s in ids})
    )
    typing = type_pcff_atoms(remapped, parameters.source)
    charges = assign_automatic_pcff_charges(remapped, typing)
    assigned = assign_pcff_parameters(remapped, typing, charges)
    other = PCFFSinglePointEvaluator(
        remapped, define_pcff_model(assigned, special_pairs=policy())
    ).evaluate()
    assert other.potential_energy == pytest.approx(result.potential_energy, abs=1e-9)
    for s in ids:
        np.testing.assert_allclose(other.forces[remap[s]], result.forces[s], atol=1e-9)
    frame = {s: original.coordinates.get(s) for s in ids}
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        plus = deepcopy(frame)
        minus = deepcopy(frame)
        plus[ids[2]][1] += h
        minus[ids[2]][1] -= h
        fd = -(
            evaluator.evaluate(plus).potential_energy
            - evaluator.evaluate(minus).potential_energy
        ) / (2 * h)
        errors.append(abs(fd - result.forces[ids[2]][1]))
    assert errors[-1] < errors[1] < errors[0]
    assert errors[-1] < 1e-5


def test_cleanup_preserves_failure():
    from island.evaluation.pcff import _cleanup

    class Broken:
        def close(self):
            raise RuntimeError("cleanup")

    original = ValueError("first")
    _cleanup(Broken(), original)
    assert original.__notes__ == ["PCFF cleanup also failed: cleanup"]
    with pytest.raises(EvaluationError, match="cleanup"):
        _cleanup(Broken())


def test_import_isolation():
    import subprocess
    import sys

    code = """
import sys
for name in ('openmm','rdkit','scipy','parmed','foyer'):
    sys.modules[name] = None
from island.evaluation import PCFFSinglePointEvaluator
assert PCFFSinglePointEvaluator
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_integrator_never_stepped(bound, monkeypatch):
    def forbidden(*args):
        pytest.fail("single-point evaluation stepped integrator")

    monkeypatch.setattr(mm.VerletIntegrator, "step", forbidden)
    bound[2].evaluate()


def test_invalid_backend_output_and_resource_disposal(bound, monkeypatch):
    from island.evaluation import pcff

    original = pcff._OpenMMResources
    closed = []

    class State:
        def getPotentialEnergy(self):
            return float("nan") * unit.kilojoule_per_mole

        def getForces(self, asNumpy):
            return np.zeros((8, 3)) * unit.kilojoule_per_mole / unit.nanometer

    class Context:
        def setPositions(self, positions):
            pass

        def getState(self, **kwargs):
            return State()

    class Resources:
        def __init__(self, *args):
            self.context = Context()

        def close(self):
            closed.append(1)

    monkeypatch.setattr(pcff, "_OpenMMResources", Resources)
    with pytest.raises(EvaluationError, match="Nonfinite"):
        bound[2].evaluate()
    assert closed == [1]
    monkeypatch.setattr(pcff, "_OpenMMResources", original)
    assert np.isfinite(bound[2].evaluate_fresh().potential_energy)


def test_independent_converter_equilibrium_roles():
    import importlib.util
    import sys

    scripts = str(Path(__file__).resolve().parents[1] / "scripts")
    sys.path.insert(0, scripts)
    try:
        spec = importlib.util.spec_from_file_location(
            "pcff_reference", Path(scripts) / "validate_pcff_singlepoint.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        sections = {
            "Angles": [
                ["1", "1", "1", "2", "3"],
                ["2", "2", "1", "2", "4"],
                ["3", "3", "3", "2", "4"],
            ],
            "Angle Coeffs": [["1", "110"], ["2", "109"], ["3", "107"]],
            "AngleAngle Coeffs": [["1", ".2", "-.3", ".4", "110", "107", "109"]],
            "Impropers": [["1", "1", "1", "2", "3", "4"]],
        }
        commands = module.aa_reference_coefficients(sections)
        values = commands[0].split()
        assert list(map(float, values[3:])) == [0.2, -0.3, 0.4, 110, 109, 107]
        sections["AngleAngle Coeffs"][0][-1] = "120"
        with pytest.raises(ValueError, match="dependency"):
            module.aa_reference_coefficients(sections)
    finally:
        sys.path.remove(scripts)
