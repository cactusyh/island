"""J33 finite-run NPT records and density diagnostics."""

from copy import deepcopy

import pytest

from island.analysis.npt import (
    NPTConfiguration,
    analyze_npt_frames,
    load_npt_record,
    relocate_npt_record,
    validate_npt_analysis,
    write_npt_record,
)
from island.exceptions import AnalysisError, ValidationError


def config(**changes):
    values = {
        "family": "GAFF2",
        "source_identity": "a" * 64,
        "final_graph_identity": "b" * 64,
        "packing_identity": "c" * 64,
        "box_lengths_nm": (3.0, 3.0, 3.0),
    }
    values.update(changes)
    return NPTConfiguration(**values)


def frames(count=16):
    return [
        {
            "stage": "npt",
            "step": i,
            "box_vectors_nm": [[3 + i / 100, 0, 0], [0, 3, 0], [0, 0, 3]],
            "measured_temperature_kelvin": 300 + i / 10,
            "measured_pressure_bar": None,
        }
        for i in range(count)
    ]


def test_configuration_is_immutable_and_identity_bound():
    value = config()
    assert len(value.identity) == 64
    assert NPTConfiguration.from_payload(value.payload) == value
    with pytest.raises(ValidationError, match="pressure"):
        config(pressure_bar=0)
    with pytest.raises(ValidationError, match="platform"):
        config(platform="CPU")
    with pytest.raises(ValidationError, match="barostat"):
        config(barostat_type="none")


def test_density_uses_actual_box_and_four_blocks():
    result = analyze_npt_frames(
        frames(), total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
    )
    assert result["status"] == "ok"
    assert result["density"]["blocks"] == 4
    assert result["density"]["mean"] != pytest.approx(100 * 1.66053906660 / 27)
    assert result["equilibrium_status"] == "not_established"
    assert result["pressure_measurement_status"] == "unavailable"


def test_fixed_box_does_not_demonstrate_density_relaxation():
    fixed = frames()
    for frame in fixed:
        frame["box_vectors_nm"] = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
    result = analyze_npt_frames(
        fixed, total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
    )
    assert result["status"] == "blocked"
    assert result["density_relaxation_status"] == "not_demonstrated"


def test_insufficient_sampling_and_nonperiodic_rejected():
    short = analyze_npt_frames(
        frames(4), total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
    )
    assert short["status"] == "insufficient_sampling"
    broken = frames()
    broken[0]["box_vectors_nm"][0][1] = 0.1
    with pytest.raises(AnalysisError, match="orthorhombic"):
        analyze_npt_frames(
            broken, total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
        )


def test_measured_pressure_requires_backend_origin():
    broken = frames()
    broken[0]["measured_pressure_bar"] = 1.0
    with pytest.raises(AnalysisError, match="backend observation"):
        analyze_npt_frames(
            broken, total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
        )


def test_record_identity_and_relocation(tmp_path):
    record = analyze_npt_frames(
        frames(), total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
    )
    record["total_mass_dalton"] = 100.0
    write_npt_record(tmp_path / "record", record)
    assert load_npt_record(tmp_path / "record")["frame_content_identity"] == record["frame_content_identity"]
    relocate_npt_record(tmp_path / "record", tmp_path / "moved")
    (tmp_path / "record/analysis.json").write_text(
        (tmp_path / "record/analysis.json").read_text().replace('"density_g_cm3"', '"density_g_cm3"', 1)
    )


def test_mutated_analysis_is_rejected():
    record = analyze_npt_frames(
        frames(), total_mass_dalton=100.0, target_temperature_kelvin=300, target_pressure_bar=1
    )
    changed = deepcopy(record)
    changed["equilibrium_status"] = "converged"
    with pytest.raises(AnalysisError, match="equilibrium"):
        validate_npt_analysis(changed)

