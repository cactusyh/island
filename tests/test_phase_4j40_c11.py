"""Native units, strict raw-output mapping and actual convergence records."""

import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest

from scripts.phase_4j40_c11_numeric import compare, force_metrics, parse_dump, signed
from scripts.validate_phase_4j40_c11 import characterize
from scripts.validate_phase_4j40_c11_receipt import validate_study


@pytest.fixture(scope="module")
def study(tmp_path_factory):
    configured = os.environ.get("ISLAND_J40_C11_RUN")
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp("c11-evidence")
    with tarfile.open("docs/evidence/phase_4j40_c11/artifacts.tar.xz") as archive:
        archive.extractall(root, filter="data")
    return root / "study"


def dump(tmp_path, *, rows=None):
    path = tmp_path / "forces.dump"
    rows = rows or ["1 1 -1 0 0 0 1 2 3", "2 1 1 1 1 1 -1 -2 -3"]
    path.write_text(
        "\n".join(
            [
                "ITEM: TIMESTEP",
                "0",
                "ITEM: NUMBER OF ATOMS",
                "2",
                "ITEM: BOX BOUNDS pp pp pp",
                "0 10",
                "0 10",
                "0 10",
                "ITEM: ATOMS id type q x y z fx fy fz",
                *rows,
            ]
        )
        + "\n"
    )
    return path


def test_raw_force_unit_conversion_and_full_mapping(tmp_path):
    r = parse_dump(dump(tmp_path))
    assert r["forces_kj_mol_nm"][0] == [41.84, 83.68, 125.52]
    assert r["box_lengths_nm"] == [1.0, 1.0, 1.0]
    assert r["positions_nm"][1] == [0.1, 0.1, 0.1]
    assert r["charges_e"] == [-1, 1]


@pytest.mark.parametrize(
    "rows",
    [
        ["1 1 -1 0 0 0 1 2 3"],
        ["2 1 -1 0 0 0 1 2 3", "1 1 1 1 1 1 -1 -2 -3"],
        ["1 1 -1 0 0 0 nan 2 3", "2 1 1 1 1 1 -1 -2 -3"],
    ],
)
def test_truncation_reordering_nonfinite_rejected(tmp_path, rows):
    with pytest.raises(ValueError):
        parse_dump(dump(tmp_path, rows=rows))


def test_force_normalization_and_percentiles():
    r = force_metrics([[1, 0, 0], [2, 0, 0]], [[0, 0, 0], [2, 0, 0]])
    assert r["rms_atom_norm_error"] == pytest.approx(np.sqrt(0.5))
    assert r["maximum_atom_norm_error"] == 1
    assert r["atom_norm_error_percentiles"]["50"] == 0.5
    assert r["relative_rms_component_error"] == pytest.approx(1 / np.sqrt(5))


def test_native_baseline_reproduced_and_characterized(study):
    c = characterize(study)
    assert c["atom_count"] == 52
    assert c["baseline_electrostatics"].startswith("OpenMM Ewald")
    assert c["absolute_energy_error_kj_mol"] == pytest.approx(
        0.04294867896004462, abs=1e-12
    )
    assert c["force_errors"]["maximum_component_error"] == pytest.approx(
        0.419290025488408, abs=1e-10
    )
    assert not c["bounded_numerical_difference_established"]
    assert c["pressure_status"] == "pressure_unavailable"


def test_native_sweep_not_metadata_contract(study):
    r = json.loads((study / "study.json").read_text())
    assert len(r["rows"]) == 38
    rows = {x["name"]: x for x in r["rows"]}
    assert (
        rows["tolerance-1e-09"]["comparison"]["absolute_energy_error_kj_mol"]
        < rows["baseline"]["comparison"]["absolute_energy_error_kj_mol"]
    )
    # Tighter requests are diagnostic controls, not a replacement acceptance policy.
    assert not rows["tolerance-1e-09"]["comparison"]["passed"]
    assert rows["box-1.5"]["comparison"]["absolute_energy_error_kj_mol"] > 1
    for row in r["rows"]:
        root = study / row["name"]
        assert "run 0" in (root / "lammps/log.lammps").read_text()
        assert (root / "openmm/system.xml").exists()
        assert len(row["openmm"]["forces_kj_mol_nm"]) == 52


def test_actual_polymer_pressure_derivative_step_sweep(study):
    r = json.loads((study / "study.json").read_text())
    for method in ("ewald", "pme"):
        for engine in ("openmm", "lammps"):
            rows = r["pressure_strain"][method][engine]
            assert [v["relative_volume_step"] for v in rows] == [1e-4, 1e-5, 1e-6]
            assert all(np.isfinite(x["total_bar"]) for x in rows)
    assert len(r["retained_frames"]) == 8


def test_offline_actual_raw_reconstruction(study):
    assert validate_study(study)["scientific_acceptance"] is False


@pytest.mark.parametrize(
    "mutation", ("forces", "energy", "tolerance", "frame", "readiness")
)
def test_renewed_outer_identity_does_not_hide_mutation(study, tmp_path, mutation):
    # Link source directories read-only, mutate only the new study record.
    for p in study.iterdir():
        if p.name != "study.json":
            (tmp_path / p.name).symlink_to(p)
    r = json.loads((study / "study.json").read_text())
    if mutation == "forces":
        r["rows"][0]["lammps"]["forces_kj_mol_nm"][0][0] += 1
    elif mutation == "energy":
        r["rows"][0]["openmm"]["energy_kj_mol"] += 1
    elif mutation == "tolerance":
        r["rows"][0]["comparison"]["energy_bound_kj_mol"] = 1
    elif mutation == "frame":
        r["rows"][0]["openmm"]["frame_identity"] = "stale"
    else:
        r["scientific_acceptance"] = True
    (tmp_path / "study.json").write_text(json.dumps(signed(r)))
    with pytest.raises(ValueError):
        validate_study(tmp_path)


def test_statistic_comparison_rejects_mismatched_frames():
    with pytest.raises(ValueError, match="Different physical frames"):
        compare({"frame_identity": "a"}, {"frame_identity": "b"})


def test_relocated_fresh_process_statistics(study, tmp_path):
    output = tmp_path / "characterization.json"
    command = [
        sys.executable,
        "-m",
        "scripts.validate_phase_4j40_c11",
        "--study",
        str(study),
        "--output",
        str(output),
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    assert result.returncode == 1, result.stdout + result.stderr
    assert json.loads(output.read_text()) == characterize(study)


def test_comparison_rejects_nonfinite_force_values():
    with pytest.raises(ValueError, match="Nonfinite"):
        force_metrics([[float("nan"), 0, 0]], [[1, 0, 0]])
