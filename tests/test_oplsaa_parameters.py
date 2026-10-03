"""Analytical software fixtures; no source charges are synthesized as real OPLS."""

from copy import deepcopy
from math import cos, pi, sqrt
from xml.etree import ElementTree as ET

import numpy as np
import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack
from island.evaluation.oplsaa import OPLSSinglePointEvaluator, build_model
from island.exceptions import EvaluationInputError, OPLSAssignmentError
from island.forcefields.oplsaa import (
    FoyerOPLSSource,
    OPLSTypingResult,
    assign_native_charges,
    source,
)
from island.forcefields.oplsaa.models import chemical_graph, typing_data
from island.forcefields.oplsaa.parameters import (
    VERIFICATION,
    OPLSParameterizationResult,
    inventory,
    resolved_data,
    selected_record,
)


def numeric(data, xyz):
    mm = pytest.importorskip("openmm")
    from openmm import unit

    model = build_model(data, mm)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(np.array(xyz) * 0.1)
        state = context.getState(getEnergy=True, getForces=True)
        return state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ), state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ) * 0.1
    finally:
        del context, integrator, model


def simple(n):
    return {
        "sites": {
            i: {
                "mass_dalton": 12.0,
                "charge_e": 0.0,
                "sigma_angstrom": 1.0,
                "epsilon_kj_mol": 0.0,
            }
            for i in range(n)
        },
        "bonds": [],
        "angles": [],
        "rb": [],
        "pairs": [],
    }


def test_harmonic_half_convention_and_units():
    d = simple(2)
    d["bonds"] = [
        {
            "sites": (0, 1),
            "converted": {"length_angstrom": 1.0, "k_kj_mol_angstrom2": 200.0},
        }
    ]
    e, f = numeric(d, [[0, 0, 0], [1.2, 0, 0]])
    assert e == pytest.approx(0.5 * 200 * 0.2**2)
    assert f[1, 0] == pytest.approx(-40)
    d = simple(3)
    d["angles"] = [
        {
            "sites": (0, 1, 2),
            "converted": {"angle_degrees": 60.0, "k_kj_mol_radian2": 50.0},
        }
    ]
    e, _ = numeric(d, [[1, 0, 0], [0, 0, 0], [0, 1, 0]])
    assert e == pytest.approx(0.5 * 50 * (pi / 6) ** 2)


@pytest.mark.parametrize("theta", [0.2, 0.8, 1.7, 2.8])
def test_rb_scan_constant_sign_and_phase(theta):
    d = simple(4)
    coeff = [1.0, -2.0, 3.0, -0.7, 0.0, 0.3]
    d["rb"] = [{"sites": (0, 1, 2, 3), "converted": {"coefficients_kj_mol": coeff}}]
    xyz = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, cos(theta), np.sin(theta)]]
    e, f = numeric(d, xyz)
    assert e == pytest.approx(
        sum(c * cos(theta - pi) ** i for i, c in enumerate(coeff)), abs=1e-12
    )
    for h in (1e-4, 5e-5):
        minus = np.array(xyz)
        plus = minus.copy()
        minus[3, 2] -= h
        plus[3, 2] += h
        fd = -(numeric(d, plus)[0] - numeric(d, minus)[0]) / (2 * h)
        assert fd == pytest.approx(f[3, 2], abs=1e-6)


@pytest.mark.parametrize("scale", [1.0, 0.5, 0.0])
def test_geometric_lj_coulomb_and_pairs(scale):
    d = simple(2)
    d["sites"][0].update(sigma_angstrom=2.0, epsilon_kj_mol=0.3, charge_e=0.4)
    d["sites"][1].update(sigma_angstrom=4.0, epsilon_kj_mol=0.7, charge_e=-0.4)
    if scale != 1:
        d["pairs"] = [{"sites": (0, 1), "lj_scale": scale, "coulomb_scale": scale}]
    r = 5.0
    e, _ = numeric(d, [[0, 0, 0], [r, 0, 0]])
    lj = 4 * sqrt(0.3 * 0.7) * ((sqrt(8) / r) ** 12 - (sqrt(8) / r) ** 6)
    coulomb = 1389.35457644382 * 0.4 * (-0.4) / r
    assert e == pytest.approx(scale * (lj + coulomb), abs=1e-8)
    d["sites"][0]["epsilon_kj_mol"] = 0.0
    e, _ = numeric(d, [[0, 0, 0], [r, 0, 0]])
    assert e == pytest.approx(scale * coulomb, abs=1e-8)


def test_ring_shortest_paths_no_multiple_torsion_counting():
    graph = Topology()
    for i in range(6):
        graph.add_site(AtomSite(i, str(i), 12.0, element="C", atomic_number=6))
    for i in range(6):
        graph.add_bond(i, (i + 1) % 6, order=1)
    system = MolecularSystem(graph, Coordinates())
    paths, pairs = inventory(system)
    assert len(paths["rb"]) == 6 and len(pairs) == 15
    assert sum(p["lj_scale"] == 0.5 for p in pairs) == 3
    graph.add_bond(0, 3, order=1)
    _, pairs = inventory(system)
    assert next(p for p in pairs if p["sites"] == (0, 3))["lj_scale"] == 0
    assert len({p["sites"] for p in pairs}) == len(pairs)


def test_upstream_matching_precedence():
    types = {"A": {"class": "X"}, "B": {"class": "Y"}}
    xml = "<ForceField><RBTorsionForce>"
    for c1, c4, c0 in [("", "", 1), ("", "Y", 2), ("X", "Y", 3)]:
        xml += f'<Proper class1="{c1}" class2="X" class3="Y" class4="{c4}" c0="{c0}" c1="0" c2="0" c3="0" c4="0" c5="0"/>'
    root = ET.fromstring(xml + "</RBTorsionForce></ForceField>")
    assert selected_record(root, types, ("A", "A", "B", "B"), "rb")[0] == 2
    root.find("RBTorsionForce").remove(root.find("RBTorsionForce")[-1])
    assert selected_record(root, types, ("A", "A", "B", "B"), "rb")[0] == 0
    assert selected_record(root, types, ("B", "B", "A", "A"), "rb")[0] == 0


@pytest.fixture
def synthetic(monkeypatch):
    # Independent toy hydrogen molecule, not native OPLS scientific evidence.
    xml = b"""<ForceField name="OPLS-AA" version="0.1.0" combining_rule="geometric">
    <AtomTypes><Type name="H" class="H" element="H" def="H"/></AtomTypes>
    <HarmonicBondForce><Bond class1="H" class2="H" length="0.1" k="10000"/></HarmonicBondForce>
    <HarmonicAngleForce/><RBTorsionForce/>
    <NonbondedForce coulomb14scale="0.5" lj14scale="0.5"><Atom type="H" charge="0" sigma="0.1" epsilon="0"/></NonbondedForce></ForceField>"""
    monkeypatch.setitem(source.PIN, "xml_sha256", source.digest(xml))
    src = FoyerOPLSSource(xml)
    g = Topology()
    for i in (13, 39):
        g.add_site(AtomSite(i, str(i), 1.008, element="H", atomic_number=1))
    g.add_bond(13, 39, order=1)
    s = MolecularSystem(g, Coordinates({13: (0, 0, 0), 39: (1.2, 0, 0)}))
    p = {
        "schema": "island_foyer_typing_v1",
        "source": src.identity,
        "graph": chemical_graph(s),
        "index_to_site_id": {0: 13, 1: 39},
        "matches": {
            i: {"atomtype": "H", "whitelist": ["H"], "blacklist": []} for i in (13, 39)
        },
        "environment": {"foyer": source.PIN["foyer_version"]},
        "data": {},
    }
    p["data"] = typing_data(p, s, src)
    charges = assign_native_charges(s, OPLSTypingResult(pack(p)), src)
    result = OPLSParameterizationResult(
        pack(
            {
                "schema": "island_opls_parameters_v1",
                "charges": charges.payload,
                "resolved": resolved_data(s, src, charges),
                "upstream_verification": VERIFICATION,
            }
        )
    )
    return s, src, result


def test_owned_snapshot_and_coordinate_frames(synthetic):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    before = deepcopy(s.to_dict())
    p.validate_integrity(s, src)
    s.topology.angles[(13, 13, 39)] = "stale cache"
    p.validate_integrity(s, src)  # Only authoritative bonds matter.
    snap = p.to_parameterized_system(s, src)
    e = OPLSSinglePointEvaluator.from_parameterized_system(snap, src)
    result = e.evaluate()
    assert result.potential_energy == pytest.approx(2.0)
    assert (
        e.evaluate({39: (1.4, 0, 0), 13: (0, 0, 0)}).coordinate_fingerprint
        != result.coordinate_fingerprint
    )
    assert result.potential_energy == pytest.approx(2.0)
    snap.site_assignments[13]["charge_e"] = 1
    with pytest.raises(EvaluationInputError):
        OPLSSinglePointEvaluator.from_parameterized_system(snap, src)
    assert s.to_dict() == before


@pytest.mark.parametrize("defect", ["missing", "value", "pair", "rb", "source", "flag"])
def test_rechecksummed_parameters_rejected(synthetic, defect):
    s, src, result = synthetic
    p = result.payload
    if defect == "missing":
        p["resolved"]["bonds"] = []
    if defect == "value":
        p["resolved"]["bonds"][0]["converted"]["k_kj_mol_angstrom2"] = 1.0
    if defect == "pair":
        p["resolved"]["pairs"][0]["lj_scale"] = 0.5
    if defect == "rb":
        p["resolved"]["rb"] = [{"sites": (13, 39, 13, 39)}]
    if defect == "source":
        p["resolved"]["source"]["xml_sha256"] = "0" * 64
    if defect == "flag":
        p["resolved"]["production_validated"] = True
    corrupt = OPLSParameterizationResult(pack(p))
    with pytest.raises(OPLSAssignmentError):
        corrupt.to_parameterized_system(s, src)


@pytest.mark.parametrize(
    "coords",
    [
        {13: (True, 0, 0), 39: (1, 0, 0)},
        {13: (1e308, 0, 0), 39: (-1e308, 0, 0)},
        {13: (0, 0, 0), 39: (0, 0, 0)},
        {13: (0, 0, 0)},
    ],
)
def test_invalid_evaluation_coordinates(synthetic, coords):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    with pytest.raises(EvaluationInputError):
        e.evaluate(coords)
    assert np.isfinite(e.evaluate().potential_energy)


def test_numeric_coordinate_types(synthetic):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    with pytest.raises(EvaluationInputError):
        e.evaluate({13: ("0", 0, 0), 39: (1.2, 0, 0)})
    s.coordinates = Coordinates()
    without_coordinates = OPLSSinglePointEvaluator(s, p, src)
    with pytest.raises(EvaluationInputError, match="No coordinates"):
        without_coordinates.evaluate()


def test_real_source_noncontiguous_permutation():
    """Opt-in actual source test, independently remapped ethanol graph."""
    import os
    from dataclasses import replace

    path = os.environ.get("ISLAND_OPLSAA_XML")
    if not path:
        pytest.skip("Set ISLAND_OPLSAA_XML in the pinned Foyer environment")
    from island.chemistry import from_smiles
    from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa

    s = from_smiles("CCO", generate_3d=True, random_seed=20261003)
    src = load_oplsaa_source(path)
    p = parameterize_oplsaa(s, src)
    result = OPLSSinglePointEvaluator(s, p, src).evaluate()
    ids = sorted(s.topology.sites)
    mapping = {i: 1000 - 7 * j for j, i in enumerate(ids)}
    topology = Topology()
    for i in reversed(ids):
        topology.add_site(replace(deepcopy(s.topology.sites[i]), id=mapping[i]))
    for (a, b), bond in reversed(list(s.topology.bonds.items())):
        topology.add_bond(
            mapping[a], mapping[b], order=bond.order, aromatic=bond.aromatic
        )
    changed = MolecularSystem(
        topology, Coordinates({mapping[i]: s.coordinates.get(i) for i in ids})
    )
    changed.topology.angles[(mapping[ids[0]],) * 3] = "stale"
    other = parameterize_oplsaa(changed, src)
    output = OPLSSinglePointEvaluator(changed, other, src).evaluate()
    assert output.potential_energy == pytest.approx(
        result.potential_energy, abs=1e-9, rel=1e-12
    )
    for i in ids:
        assert output.forces[mapping[i]] == pytest.approx(
            result.forces[i], abs=1e-9, rel=1e-12
        )
    moved = deepcopy(changed)
    moved.coordinates.set(mapping[ids[0]], (20, 30, 40))
    other.validate_integrity(moved, src)
