"""Independent energy derivatives, mutation rejection and retained regressions."""

import copy
import json
import tarfile
from pathlib import Path

import numpy as np
import pytest

from scripts.phase_4j40_c11_numeric import signed
from scripts.phase_4j40_c12 import completion, corrected


@pytest.fixture(scope="module")
def retained(tmp_path_factory):
    root = tmp_path_factory.mktemp("c12-native-source")
    with tarfile.open("docs/evidence/phase_4j40_c11/artifacts.tar.xz") as archive:
        archive.extractall(root, filter="data")
    root /= "study"
    backend = json.loads((root / "inputs/backend.json").read_text())

    def case(name):
        frame = json.loads((root / name / "frame.json").read_text())
        native = json.loads((root / name / "lammps/point.json").read_text())
        om = json.loads((root / name / "openmm/point.json").read_text())
        return backend, frame, native["setting"], native, om

    return case


def test_exact_strained_energy_and_force_repair(retained):
    b, f, s, n, om = retained("box-1.5")
    c = completion(b, f, s, n)
    result = corrected(n, c)
    assert len(c["pairs"]) == 18
    assert c["energy_lj_kj_mol"] == pytest.approx(-0.401939047921807, abs=1e-12)
    assert c["energy_coulomb_kj_mol"] == pytest.approx(-1.352337038039448, abs=1e-7)
    assert abs(om["energy_kj_mol"] - result["energy_kj_mol"]) < 1e-5
    assert (
        np.max(np.abs(np.array(om["forces_kj_mol_nm"]) - result["forces_kj_mol_nm"]))
        < 1e-4
    )
    assert n["engine"] == "LAMMPS"
    assert "composite" in result["engine"]


@pytest.mark.parametrize("name", ["baseline", "box-1.25", "translation"])
def test_unstrained_no_correction(retained, name):
    b, f, s, n, _ = retained(name)
    c = completion(b, f, s, n)
    assert not c["pairs"]
    assert corrected(n, c)["energy_kj_mol"] == n["energy_kj_mol"]


@pytest.mark.parametrize("field", ["scaling", "cutoff", "box", "charge", "ordering"])
def test_mutated_inputs_rejected_even_with_renewed_frame_digest(retained, field):
    b, f, s, n, _ = copy.deepcopy(retained("box-1.5"))
    if field == "scaling":
        b["numerical_data"]["pairs"][-1]["coulomb"] += 0.01
        f["numerical_data"] = b["numerical_data"]
    elif field == "cutoff":
        s["cutoff_nm"] = 0.4
        # Change only requested setting; native record stays independently bound.
        n["setting"] = {**s, "cutoff_nm": 0.5}
    elif field == "box":
        f["box_lengths_nm"][0] += 0.1
    elif field == "charge":
        b["numerical_data"]["sites"][0]["charge"] += 0.1
        f["numerical_data"] = b["numerical_data"]
    else:
        f["atom_ids"].reverse()
    f = signed(f)
    with pytest.raises(ValueError):
        completion(b, f, s, n)


def test_central_energy_derivative_force_and_volume_virial(retained):
    b, f, s, n, _ = retained("box-1.5")
    c = completion(b, f, s, n)

    def energy(frame):
        frame = signed(frame)
        native = dict(n, frame_identity=frame["identity"])
        r = completion(b, frame, s, native)
        return r["energy_lj_kj_mol"] + r["energy_coulomb_kj_mol"]

    errors = []
    for h in (1e-4, 1e-5, 1e-6):
        values = []
        for sign in (-1, 1):
            frame = copy.deepcopy(f)
            frame["positions_nm"][0][0] += sign * h
            values.append(energy(frame))
        errors.append(
            abs(-(values[1] - values[0]) / (2 * h) - c["forces_kj_mol_nm"][0][0])
        )
    assert errors[-1] < 1e-7
    assert errors[-1] < errors[0]
    h = 1e-6
    values = []
    for sign in (-1, 1):
        scale = (1 + sign * h) ** (1 / 3)
        frame = copy.deepcopy(f)
        frame["positions_nm"] = (np.array(f["positions_nm"]) * scale).tolist()
        frame["box_lengths_nm"] = (np.array(f["box_lengths_nm"]) * scale).tolist()
        values.append(energy(frame))
    derivative = -(values[1] - values[0]) / (2 * h)
    assert derivative == pytest.approx(np.trace(c["virial_kj_mol"]) / 3, abs=1e-8)
    assert np.max(np.abs(np.sum(c["forces_kj_mol_nm"], axis=0))) < 1e-14


def test_tampered_completion_and_native_rejected(retained):
    b, f, s, n, _ = retained("box-1.5")
    c = completion(b, f, s, n)
    c["energy_lj_kj_mol"] += 1
    with pytest.raises(ValueError):
        corrected(n, c)
    c = completion(b, f, s, n)
    n["energy_kj_mol"] += 1
    with pytest.raises(ValueError):
        corrected(n, c)


def test_independent_openmm_exception_energy_force_control(retained):
    mm = pytest.importorskip("openmm")
    from openmm import unit

    b, frame, setting, native, _ = retained("box-1.5")
    c = completion(b, frame, setting, native)
    model = mm.System()
    for _ in frame["atom_ids"]:
        model.addParticle(12.0)
    model.setDefaultPeriodicBoxVectors(
        *[mm.Vec3(*v) for v in np.diag(frame["box_lengths_nm"])]
    )
    force = mm.CustomBondForce("4*eps*((sig/r)^12-(sig/r)^6)+q*(f-erf(alpha*r))/r")
    for name in ("eps", "sig", "q", "f", "alpha"):
        force.addPerBondParameter(name)
    sites = {s["id"]: s for s in b["numerical_data"]["sites"]}
    for p in c["pairs"]:
        a, d = p["atom_ids"]
        x, y = sites[a], sites[d]
        force.addBond(
            a - 1,
            d - 1,
            [
                np.sqrt(x["epsilon"] * y["epsilon"]) * p["lj_scaling"],
                (x["sigma"] + y["sigma"]) / 20,
                138.935456264 * x["charge"] * y["charge"],
                p["coulomb_scaling"],
                c["alpha_nm_from_8_digit_log"],
            ],
        )
    force.setUsesPeriodicBoundaryConditions(True)
    model.addForce(force)
    integrator = mm.VerletIntegrator(0.0005)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(frame["positions_nm"])
    state = context.getState(getEnergy=True, getForces=True)
    assert state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    ) == pytest.approx(c["energy_lj_kj_mol"] + c["energy_coulomb_kj_mol"], abs=1e-13)
    assert np.array(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
    ) == pytest.approx(np.array(c["forces_kj_mol_nm"]), abs=1e-12)
    del context, integrator


def test_excluded_12_13_pairs_are_completed_without_inventing_lj(retained):
    b, f, s, n, _ = copy.deepcopy(retained("box-1.5"))
    f["positions_nm"][0][0] += 0.6
    f = signed(f)
    n["frame_identity"] = f["identity"]
    c = completion(b, f, s, n)
    excluded = [p for p in c["pairs"] if p["path_length"] < 3]
    assert {p["path_length"] for p in excluded} == {1, 2}
    assert all(p["lj_kj_mol"] == 0 and p["coulomb_scaling"] == 0 for p in excluded)


@pytest.fixture(scope="module")
def actual_study(tmp_path_factory):
    import os

    configured = os.environ.get("ISLAND_J40_C12_RUN")
    if configured:
        return Path(configured)
    archive_path = Path("docs/evidence/phase_4j40_c12/artifacts.tar.xz")
    if not archive_path.exists():
        pytest.skip("Execute C12 study or provide its native evidence archive")
    root = tmp_path_factory.mktemp("c12-actual-study")
    with tarfile.open(archive_path) as archive:
        archive.extractall(root, filter="data")
    return root / "study"


def test_actual_complete_matrix_and_raw_reconstruction(actual_study):
    from scripts.validate_phase_4j40_c12_receipt import validate_study

    r = validate_study(actual_study)
    assert r["additional_rows"] == 42
    assert r["scientific_acceptance"] is False


@pytest.mark.parametrize(
    "field", ["completion", "native_force", "cutoff", "box", "readiness"]
)
def test_actual_study_mutation_rejected(actual_study, tmp_path, field):
    from scripts.validate_phase_4j40_c12_receipt import validate_study

    for p in actual_study.iterdir():
        if p.name not in ("c11-completed.json", "study.json"):
            (tmp_path / p.name).symlink_to(p)
    r = json.loads((actual_study / "study.json").read_text())
    rows = json.loads((actual_study / "c11-completed.json").read_text())
    if field == "readiness":
        r["scientific_acceptance"] = True
    elif field == "completion":
        rows[-1]["completion"]["energy_lj_kj_mol"] += 1
        rows[-1]["completion"] = signed(rows[-1]["completion"])
    elif field == "native_force":
        rows[-1]["completion"]["forces_kj_mol_nm"][0][0] += 1
        rows[-1]["completion"] = signed(rows[-1]["completion"])
    elif field == "cutoff":
        rows[-1]["completion"]["setting"]["cutoff_nm"] += 0.01
        rows[-1]["completion"] = signed(rows[-1]["completion"])
    else:
        rows[-1]["completion"]["frame_identity"] = "mutated-box"
        rows[-1]["completion"] = signed(rows[-1]["completion"])
    (tmp_path / "study.json").write_text(json.dumps(signed(r)))
    (tmp_path / "c11-completed.json").write_text(json.dumps(rows))
    with pytest.raises(ValueError):
        validate_study(tmp_path)


def test_truncated_completion_cannot_broadcast_to_all_atoms(retained):
    b, f, s, n, _ = retained("box-1.5")
    c = completion(b, f, s, n)
    c["forces_kj_mol_nm"] = c["forces_kj_mol_nm"][:1]
    with pytest.raises(ValueError, match="force-vector mapping"):
        corrected(n, signed(c))
