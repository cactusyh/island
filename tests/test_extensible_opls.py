"""Extensible source contracts and analytical fixtures, never PDMS parameters."""

import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields import preparation as native
from island.forcefields.extensible import (
    ExternalOPLSError,
    ExternalOPLSOptions,
    ExternalOPLSPreparation,
    ForceFieldRequest,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.external_pdms import digest
from island.forcefields.opls_external import (
    extension_graph,
    load_source,
    write_source_manifest,
)
from island.forcefields.opls_external.assignment import _select
from island.forcefields.opls_external.source import inspect_xml
from island.forcefields.opls_external.storage import load_preparation, save_preparation

ROOT = Path(__file__).resolve().parents[1]


def pin(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def manifest_for(directory, xml, *, elements=("C", "H", "F", "Si")):
    """All values below are analytical test data, not scientific source rows."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "source.xml").write_bytes(xml)
    (directory / "LICENSE").write_text(
        "Analytical software test data only; no physical parameter authority.\n"
    )
    output = directory / "source.json"
    p = write_source_manifest(
        directory / "source.xml",
        output,
        profile_id="island_external_opls_analytical_test_v1",
        provenance={
            "citation": "Analytical software fixture; NOT a physical PDMS model",
            "doi": None,
            "charge_method": "explicit analytical test charges",
            "license_status": "authorized_local_use",
            "license_evidence": {
                "file": "LICENSE",
                "sha256": pin(directory / "LICENSE"),
                "line": 1,
                "raw_text": (directory / "LICENSE").read_text().strip(),
            },
            "origin": "test_extensible_opls.py",
            "conversion": "identity test XML units",
        },
        model={
            "nonbonded_method": "NoCutoff",
            "periodic": False,
            "improper_policy": "source_defines_no_impropers_for_declared_domain",
            "chemical_state": "neutral_closed_shell",
            "elements": list(elements),
        },
    )
    return output, p


@pytest.fixture
def analytical_xml():
    r = ET.Element(
        "ForceField",
        name="Analytical software test",
        version="1",
        combining_rule="geometric",
    )
    types = ET.SubElement(r, "AtomTypes")
    # Source-driven rules demonstrate multiple new elements without editing code.
    for element, mass, valence in (
        ("C", "12.011", 4),
        ("H", "1.008", 1),
        ("F", "18.998", 1),
        ("Si", "28.086", 4),
    ):
        ET.SubElement(
            types,
            "Type",
            name=element,
            **{
                "class": element,
                "element": element,
                "mass": mass,
                "def": f"[{element};X{valence}]",
            },
        )
    b = ET.SubElement(r, "HarmonicBondForce")
    ET.SubElement(b, "Bond", class1="", class2="", length="0.15", k="10000")
    a = ET.SubElement(r, "HarmonicAngleForce")
    ET.SubElement(a, "Angle", class1="", class2="", class3="", angle="1.91", k="50")
    t = ET.SubElement(r, "RBTorsionForce")
    ET.SubElement(
        t,
        "Proper",
        class1="",
        class2="",
        class3="",
        class4="",
        c0="1",
        c1="-1",
        c2="0",
        c3="0",
        c4="0",
        c5="0",
    )
    n = ET.SubElement(r, "NonbondedForce", lj14scale="0.5", coulomb14scale="0.5")
    for name in ("C", "H", "F", "Si"):
        ET.SubElement(n, "Atom", type=name, charge="0", sigma="0.3", epsilon="0.1")
    ET.indent(r)
    return ET.tostring(r)


@pytest.fixture
def source_path(tmp_path, analytical_xml):
    return manifest_for(tmp_path / "analytical", analytical_xml)


@pytest.fixture
def typing_engine():
    from island.exceptions import OPLSDependencyError
    from island.forcefields.oplsaa.adapter import check_foyer_installation

    try:
        check_foyer_installation()
    except OPLSDependencyError as error:
        pytest.skip(str(error))


def prepare(s, source_path=None):
    return prepare_forcefield(
        s, ForceFieldRequest("external_opls", ExternalOPLSOptions(*(source_path or ())))
    )


def test_native_dispatch_unchanged():
    path = Path("explicit.xml")
    a = ForceFieldRequest("oplsaa", native.OPLSOptions(path))
    b = native.ForceFieldRequest("oplsaa", native.OPLSOptions(path))
    assert a == b and type(a) is native.ForceFieldRequest
    with pytest.raises(native.ForceFieldRequestError):
        ForceFieldRequest("external_opls", ExternalOPLSOptions(), final_graph="native")


@pytest.mark.parametrize(
    "smiles", ["[SiH3]O[SiH3]", "C[Si](C)(C)O[Si](C)(C)C", "CF", "CS", "N"]
)
def test_graph_validation_is_not_an_element_whitelist(smiles):
    s = from_smiles(smiles, random_seed=2026)
    g = extension_graph(s)
    assert g["atom_order"] == list(s.topology.sites)
    r = prepare(s)
    assert r.payload["status"] == "asset_unavailable"
    assert r.payload["profile_identity"] is None
    assert not r.payload["execution_permitted"]


@pytest.mark.parametrize(
    "mutation", ["radical", "charge", "valence", "element", "implicit_H"]
)
def test_invalid_graph_classified_before_engine(mutation, monkeypatch):
    s = from_smiles("[SiH3]O[SiH3]", random_seed=2026)
    i = next(i for i, a in s.topology.sites.items() if a.element == "Si")
    a = s.topology.sites[i]
    if mutation == "radical":
        a.metadata["radical_electrons"] = 1
    elif mutation == "charge":
        s.topology.sites[i] = replace(a, formal_charge=1)
    elif mutation == "element":
        s.topology.sites[i] = replace(a, element="P")
    elif mutation == "implicit_H":
        a.metadata["explicit_hydrogen_count"] = 1
    else:
        s.topology.bonds.pop(next(k for k in s.topology.bonds if i in k))
    r = prepare(s)
    assert r.payload["status"] == "typing_blocked"
    assert r.payload["graph"] is None
    monkeypatch.setattr(
        "island.forcefields.opls_external.evaluator._openmm",
        lambda: pytest.fail("engine created"),
    )
    with pytest.raises(ExternalOPLSError) as error:
        create_evaluator(s, r)
    assert error.value.status == "typing_blocked"
    assert error.value.result is r
    assert error.value.diagnostics == r.payload["diagnostics"]


@pytest.mark.parametrize(
    "defect", ["missing_pin", "wrong_pin", "bad_json", "missing_asset"]
)
def test_explicit_source_failures_have_complete_result(tmp_path, source_path, defect):
    path, expected = source_path
    if defect == "missing_pin":
        expected = None
    if defect == "wrong_pin":
        expected = "0" * 64
    if defect == "bad_json":
        path.write_text("{")
        expected = pin(path)
    if defect == "missing_asset":
        path = tmp_path / "absent.json"
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, (path, expected))
    assert r.payload["status"] in {"asset_unavailable", "source_incompatible"}
    with pytest.raises(ExternalOPLSError) as e:
        create_evaluator(s, r)
    assert e.value.result is r


@pytest.mark.parametrize(
    "defect",
    [
        "units",
        "mixing",
        "scaling",
        "license",
        "cross_terms",
        "charge",
        "row",
        "authority",
    ],
)
def test_source_mutations_rejected(source_path, defect):
    path, expected = source_path
    src = load_source(path, expected)
    p = json.loads(path.read_text())
    if defect == "units":
        p["units"] = "kcal/mol"
    if defect == "mixing":
        x = path.parent / "source.xml"
        x.write_text(
            x.read_text().replace(
                'combining_rule="geometric"', 'combining_rule="lorentz"'
            )
        )
    if defect == "scaling":
        x = path.parent / "source.xml"
        x.write_text(x.read_text().replace('lj14scale="0.5"', 'lj14scale="0.75"'))
    if defect == "license":
        (path.parent / "LICENSE").write_text("altered")
    if defect == "cross_terms":
        x = path.parent / "source.xml"
        x.write_text(
            x.read_text().replace("</ForceField>", "<Class2BondBond/></ForceField>")
        )
    if defect == "charge":
        x = path.parent / "source.xml"
        x.write_text(x.read_text().replace('charge="0"', 'charge="0.001"'))
    if defect == "row":
        x = path.parent / "source.xml"
        x.write_text(x.read_text().replace('k="10000"', 'k="0"'))
    if defect == "authority":
        p["execution_authorized"] = True
    if defect in {"units", "authority"}:
        path.write_text(json.dumps(p))
    with pytest.raises((ValueError, ExternalOPLSError)):
        src.validate_integrity()


@pytest.mark.parametrize("smiles", ["C[SiH3]", "CF"])
def test_multiple_new_elements_require_source_rules(source_path, typing_engine, smiles):
    r = prepare(from_smiles(smiles, random_seed=2026), source_path)
    assert r.payload["status"] == "prepared", r.payload["diagnostics"]
    assert not r.payload["native_opls_coverage"]
    assert all(r.payload["assignment"]["typing_matches"].values())


@pytest.mark.parametrize(
    "family", ["HarmonicBondForce", "HarmonicAngleForce", "RBTorsionForce"]
)
def test_all_missing_rows_preserved(tmp_path, analytical_xml, typing_engine, family):
    root = ET.fromstring(analytical_xml)
    root.find(family).clear()
    path = manifest_for(tmp_path / family, ET.tostring(root))
    s = from_smiles("CCCC", random_seed=2026)
    r = prepare(s, path)
    assert r.payload["status"] == "parameter_incomplete"
    kind = {
        "HarmonicBondForce": "bonds",
        "HarmonicAngleForce": "angles",
        "RBTorsionForce": "rb",
    }[family]
    from island.forcefields.oplsaa.parameters import inventory

    assert len(r.payload["diagnostics"]) == len(inventory(s)[0][kind])
    assert r.payload["profile_identity"] is None and r.payload["assignment"] is None
    with pytest.raises(ExternalOPLSError) as e:
        create_evaluator(s, r)
    assert e.value.status == "parameter_incomplete"


def test_explicit_source_zero_retained_not_filled(source_path, typing_engine):
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    assert r.payload["status"] == "prepared"
    row = r.payload["assignment"]["rb"][0]
    assert row["original_values"]["c2"] == "0"
    ev = r.payload["assignment"]["source_row_evidence"]
    assert any('c2="0"' in e["raw_text"] for e in ev)
    assert all(
        e["source_sha256"] == load_source(*source_path).document["xml_sha256"]
        for e in ev
    )


def test_conflicting_wildcards_and_specific_precedence(source_path):
    src = load_source(*source_path)
    root, types, _, _ = inspect_xml(src.xml, src.document)
    b = root.find("HarmonicBondForce")
    ET.SubElement(b, "Bond", class1="", class2="", length="0.16", k="10000")
    with pytest.raises(ValueError, match="ambiguous_conflicting"):
        _select(root, types, ("C", "H"), "HarmonicBondForce")
    ET.SubElement(b, "Bond", type1="H", type2="C", length="0.14", k="10000")
    row = _select(root, types, ("C", "H"), "HarmonicBondForce")
    assert not row["selected_wildcard"]
    assert row["converted"]["length_angstrom"] == 1.4


@pytest.mark.parametrize(
    "change", ["status", "charge", "units", "graph", "identity", "source"]
)
def test_resigned_preparation_mutations_rejected(source_path, typing_engine, change):
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    p = r.payload
    if change == "status":
        p["status"] = "asset_unavailable"
    if change == "charge":
        p["assignment"]["sites"][min(s.topology.sites)]["charge_e"] = 1
    if change == "units":
        p["assignment"]["bonds"][0]["converted"]["length_angstrom"] *= 10
    if change == "graph":
        p["graph"]["atom_order"].reverse()
    if change == "identity":
        p["profile_identity"] = "native_oplsaa"
    if change == "source":
        p["source_identity"] = "0" * 64
    p.pop("identity")
    p["identity"] = digest(p)
    with pytest.raises(ExternalOPLSError):
        ExternalOPLSPreparation(pack(p), s, r._source)


@pytest.mark.parametrize("change", ["mass", "order", "component", "bond"])
def test_final_graph_binding(source_path, typing_engine, change):
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    moved = deepcopy(s)
    if change == "mass":
        moved.topology.sites[min(s.topology.sites)] = replace(
            moved.topology.sites[min(s.topology.sites)], mass=1.0
        )
    if change == "order":
        moved.topology.sites = dict(reversed(list(moved.topology.sites.items())))
    if change == "component":
        moved.topology.sites[min(s.topology.sites)].metadata["chain_id"] = "changed"
    if change == "bond":
        moved.topology.bonds.pop(next(iter(moved.topology.bonds)))
    with pytest.raises(ExternalOPLSError):
        create_evaluator(moved, r)
    moved = deepcopy(s)
    moved.coordinates.set(min(s.topology.sites), (10, 20, 30))
    r.validate_integrity(moved)


def test_unknown_license_and_unsupported_model_never_execute(source_path, monkeypatch):
    path, _ = source_path
    p = json.loads(path.read_text())
    p["provenance"]["license_status"] = "not_established"
    p["provenance"]["license_evidence"] = None
    p["model"]["nonbonded_method"] = "lj/charmm/coul/long+PPPM"
    path.write_text(json.dumps(p))
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, (path, pin(path)))
    assert r.payload["status"] == "scientifically_unqualified"
    assert {d["stage"] for d in r.payload["diagnostics"]} == {
        "source_use_provenance_unestablished",
        "unsupported_nonbonded_method",
    }
    monkeypatch.setattr(
        "island.forcefields.opls_external.evaluator._openmm",
        lambda: pytest.fail("engine reached"),
    )
    with pytest.raises(ExternalOPLSError):
        create_evaluator(s, r)


def test_missing_dependency_classified(source_path, monkeypatch):
    from island.exceptions import OPLSDependencyError

    def missing():
        raise OPLSDependencyError("missing pinned engine")

    monkeypatch.setattr(
        "island.forcefields.oplsaa.adapter.check_foyer_installation", missing
    )
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    assert r.payload["status"] == "asset_unavailable"
    assert r.payload["diagnostics"][0]["stage"] == "pinned_typing_engine_unavailable"


def test_stale_source_before_backend_execution(source_path, typing_engine, monkeypatch):
    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    e = create_evaluator(s, r)
    (source_path[0].parent / "source.xml").write_text("changed")
    monkeypatch.setattr(
        "island.evaluation.oplsaa._openmm", lambda: pytest.fail("engine reached")
    )
    with pytest.raises(ExternalOPLSError):
        e.evaluate()


def test_relocation_reconstruction(source_path, typing_engine, tmp_path):
    s = from_smiles("CF", random_seed=2026)
    r = prepare(s, source_path)
    original = source_path[0].parent
    saved_pin = save_preparation(r, original / "prepared.json")
    moved = tmp_path / "relocated"
    shutil.copytree(original, moved)
    shutil.rmtree(original)
    options = ExternalOPLSOptions(moved / "source.json", source_path[1])
    replay = load_preparation(moved / "prepared.json", saved_pin, s, options)
    assert replay.identity == r.payload["identity"]
    code = """import sys
from island.chemistry import from_smiles
from island.forcefields.extensible import ExternalOPLSOptions
from island.forcefields.opls_external.storage import load_preparation
r=load_preparation(sys.argv[1],sys.argv[2],from_smiles('CF',random_seed=2026),ExternalOPLSOptions(sys.argv[3],sys.argv[4]))
print(r.identity)
"""
    out = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(moved / "prepared.json"),
            saved_pin,
            str(moved / "source.json"),
            source_path[1],
        ],
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert out.stdout.strip() == replay.identity


@pytest.mark.parametrize("smiles", ["CC", "CCCC"])
def test_real_source_numerical_control(tmp_path, typing_engine, smiles):
    """Real unchanged Foyer rows, not artificial Si scientific qualification."""
    xml = os.environ.get("ISLAND_OPLSAA_XML")
    license_path = os.environ.get("ISLAND_FOYER_LICENSE")
    if not xml or not license_path:
        pytest.skip("Explicit real Foyer XML + license paths required")
    directory = tmp_path / "real"
    directory.mkdir()
    shutil.copy(xml, directory / "source.xml")
    shutil.copy(license_path, directory / "LICENSE")
    p = write_source_manifest(
        directory / "source.xml",
        directory / "source.json",
        profile_id="island_external_opls_foyer_control_v1",
        provenance={
            "citation": "Pinned Foyer OPLS finite software control",
            "doi": None,
            "charge_method": "source NonbondedForce Atom charge",
            "license_status": "authorized_local_use",
            "license_evidence": {
                "file": "LICENSE",
                "sha256": pin(directory / "LICENSE"),
                "line": 5,
                "raw_text": (directory / "LICENSE").read_text().splitlines()[4],
            },
            "origin": "Foyer dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd",
            "conversion": "unchanged XML; explicit finite software control",
        },
        model={
            "nonbonded_method": "NoCutoff",
            "periodic": False,
            "improper_policy": "source_defines_no_impropers_for_declared_domain",
            "chemical_state": "neutral_closed_shell",
            "elements": ["C", "H"],
        },
    )
    s = from_smiles(smiles, random_seed=2026)
    r = prepare(s, (directory / "source.json", p))
    assert r.payload["status"] == "prepared", r.payload["diagnostics"]
    e = create_evaluator(s, r)
    result = e.evaluate()
    n = native.prepare_forcefield(
        s, native.ForceFieldRequest("oplsaa", native.OPLSOptions(Path(xml)))
    )
    reference = native.create_evaluator(s, n).evaluate()
    assert result.potential_energy == pytest.approx(
        reference.potential_energy, abs=1e-10
    )
    for key, value in reference.energy_components.items():
        assert result.energy_components[
            "proper_torsion" if key == "rb_proper" else key
        ] == pytest.approx(value, abs=1e-10)
    assert np.array(list(result.forces.values())) == pytest.approx(
        np.array(list(reference.forces.values())), abs=1e-10
    )
    for h in (1e-4, 5e-5):
        xyz = {i: s.coordinates.get(i) for i in s.topology.sites}
        i = min(xyz)
        minus = deepcopy(xyz)
        plus = deepcopy(xyz)
        minus[i] = (xyz[i][0] - h, *xyz[i][1:])
        plus[i] = (xyz[i][0] + h, *xyz[i][1:])
        fd = -(
            e.evaluate(plus).potential_energy - e.evaluate(minus).potential_energy
        ) / (2 * h)
        assert fd == pytest.approx(result.forces[i][0], abs=1e-4)
    from island.minimization import minimize_geometry

    minimum = minimize_geometry(s, e)
    assert minimum.converged and minimum.final_evaluation_verified
    assert minimum.final_evaluation.potential_energy < result.potential_energy


def test_langevin_checkpoint_continuation(source_path, typing_engine, tmp_path):
    """Finite analytical software control, no PDMS physical MD claim."""
    from island import Coordinates
    from island.dynamics import (
        DynamicsSegmentOptions,
        LangevinOptions,
        create_dynamics_checkpoint,
        initialize_velocities,
        load_dynamics_checkpoint,
        resume_dynamics,
        run_dynamics_segment,
        save_dynamics_checkpoint,
    )

    s = from_smiles("CC", random_seed=2026)
    r = prepare(s, source_path)
    e = create_evaluator(s, r)
    v = initialize_velocities(s, temperature_kelvin=10, seed=2026).velocities
    opts = LangevinOptions(
        timestep_fs=0.01,
        temperature_kelvin=10,
        friction_per_ps=1,
        thermostat_seed=123,
        steps=4,
        max_evaluations=10,
        max_frames=10,
    )
    whole = run_dynamics_segment(s, e, v, opts)
    first = run_dynamics_segment(s, e, v, replace(opts, steps=2))
    assert whole.completed and first.completed
    cp = create_dynamics_checkpoint(first)
    path = tmp_path / "checkpoint.json"
    save_dynamics_checkpoint(cp, path)
    restored = load_dynamics_checkpoint(path)
    assert restored.payload_json == cp.payload_json
    assert restored.content_checksum == cp.content_checksum
    moved = first.to_system(s)
    e.validate_system(moved)  # Coordinate history is not chemical identity.
    end = resume_dynamics(restored, moved, e, DynamicsSegmentOptions(2, 10, 10))
    assert end.completed
    assert end.final_state == whole.final_state
    assert end.payload["rng"] == whole.payload["rng"]
    assert prepare(moved, source_path).identity == r.identity
    moved.coordinates = Coordinates(dict(end.final_state.coordinates))
    e.validate_system(moved)


@pytest.mark.parametrize(
    "status",
    [
        "asset_unavailable",
        "typing_blocked",
        "parameter_incomplete",
        "source_incompatible",
        "scientifically_unqualified",
    ],
)
def test_fresh_process_blocked_reconstruction(
    tmp_path, analytical_xml, typing_engine, status
):
    path, expected = manifest_for(tmp_path / "inputs", analytical_xml)
    s = from_smiles("CC", random_seed=2026)
    if status == "typing_blocked":
        s.topology.sites[min(s.topology.sites)].metadata["radical_electrons"] = 1
    elif status == "parameter_incomplete":
        root = ET.fromstring(analytical_xml)
        root.find("HarmonicAngleForce").clear()
        (path.parent / "source.xml").write_bytes(ET.tostring(root))
        p = json.loads(path.read_text())
        p["xml_sha256"] = pin(path.parent / "source.xml")
        path.write_text(json.dumps(p))
        expected = pin(path)
    elif status == "source_incompatible":
        expected = "0" * 64
    elif status == "scientifically_unqualified":
        p = json.loads(path.read_text())
        p["provenance"]["license_status"] = "not_established"
        path.write_text(json.dumps(p))
        expected = pin(path)
    source_spec = None if status == "asset_unavailable" else (path, expected)
    r = prepare(s, source_spec)
    assert r.payload["status"] == status
    saved = save_preparation(r, tmp_path / "blocked.json")
    code = """import sys
from island.chemistry import from_smiles
from island.forcefields.extensible import ExternalOPLSOptions,ExternalOPLSError,create_evaluator
from island.forcefields.opls_external.storage import load_preparation
s=from_smiles('CC',random_seed=2026)
if sys.argv[1]=='typing_blocked':s.topology.sites[min(s.topology.sites)].metadata['radical_electrons']=1
options=ExternalOPLSOptions() if sys.argv[1]=='asset_unavailable' else ExternalOPLSOptions(sys.argv[4],sys.argv[5])
r=load_preparation(sys.argv[2],sys.argv[3],s,options)
try:create_evaluator(s,r)
except ExternalOPLSError as e:
 assert e.result is r
 assert r.payload['profile_identity'] is None
 print(e.status)
else:raise AssertionError('blocked evaluator executed')
"""
    out = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            status,
            str(tmp_path / "blocked.json"),
            saved,
            str(path),
            expected,
        ],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert out.stdout.strip() == status


@pytest.mark.parametrize("theta", [0.2, 0.8, 1.7, 2.8])
def test_four_term_opls_equation(theta):
    """Independent analytic torsion scan; arbitrary coefficients are test data."""
    from math import cos, sin

    from test_oplsaa_parameters import simple

    from island.forcefields.opls_external.evaluator import build_model

    mm = pytest.importorskip("openmm")
    from openmm import unit

    d = simple(4)
    d["torsion_form"] = "OPLSTorsionForce"
    coeff = [1.2, -0.7, 0.3, 0.9]
    d["rb"] = [{"sites": (0, 1, 2, 3), "converted": {"coefficients_kj_mol": coeff}}]
    model = build_model(d, mm)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(
            np.array([[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, cos(theta), sin(theta)]])
            * 0.1
        )
        energy = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        reference = 0.5 * (
            coeff[0] * (1 + cos(theta))
            + coeff[1] * (1 - cos(2 * theta))
            + coeff[2] * (1 + cos(3 * theta))
            + coeff[3] * (1 - cos(4 * theta))
        )
        assert energy == pytest.approx(reference, abs=1e-12)
    finally:
        del context, integrator, model


def test_source_charge_residual_never_rounds_to_neutral(
    tmp_path, analytical_xml, typing_engine
):
    # A deliberately tiny analytical residual used to reproduce ambient Decimal
    # precision loss. It is neither a source charge nor a neutrality tolerance.
    root = ET.fromstring(analytical_xml)
    charges = {"C": "-0.3", "H": "0.1", "F": "1e-40", "Si": "0"}
    for row in root.find("NonbondedForce"):
        row.set("charge", charges[row.get("type")])
    path = manifest_for(tmp_path / "nonneutral", ET.tostring(root))
    r = prepare(from_smiles("CF", random_seed=2026), path)
    assert r.payload["status"] == "parameter_incomplete"
    assert r.payload["assignment"] is None
    assert (
        "source_component_charge_not_exactly_neutral"
        in r.payload["diagnostics"][0]["reason"]
    )


def test_unit_conversion_independent_of_global_decimal_context(
    tmp_path, analytical_xml, typing_engine
):
    from decimal import localcontext

    from island.forcefields.opls_external.assignment import _resolved

    root = ET.fromstring(analytical_xml)
    bond = root.find("HarmonicBondForce/Bond")
    bond.set("length", "0.1543")
    bond.set("k", "12345")
    source = manifest_for(tmp_path / "units", ET.tostring(root))
    s = from_smiles("CF", random_seed=2026)
    reference = prepare(s, source)
    assert reference.payload["status"] == "prepared"
    _resolved.cache_clear()
    with localcontext() as ctx:
        ctx.prec = 2
        other = prepare(s, source)
        assert other.identity == reference.identity
    row = other.payload["assignment"]["bonds"][0]["converted"]
    assert row == {"length_angstrom": 1.543, "k_kj_mol_angstrom2": 123.45}
