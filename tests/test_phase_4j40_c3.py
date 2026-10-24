import pytest

from island.analysis.periodic_preflight import (
    PeriodicPreflightPolicy,
    preflight_box,
    pressure_from_virial,
    pressure_unavailable,
)
from island.exceptions import AnalysisError


def test_preflight_accepts_safe_box():
    x=preflight_box(initial_box_nm=(3,3,3),predicted_min_box_nm=(2,2,2),cutoff_nm=.5,policy=PeriodicPreflightPolicy(.05,.1))
    assert x['status']=='ok' and x['safe_predicted']

def test_preflight_rejects_predicted_collapse():
    x=preflight_box(initial_box_nm=(3,3,3),predicted_min_box_nm=(1,1,1),cutoff_nm=.5)
    assert x['status']=='blocked' and not x['safe_predicted']

def test_preflight_rejects_bad_switching():
    with pytest.raises(AnalysisError): preflight_box(initial_box_nm=(3,3,3),cutoff_nm=.5,switching_distance_nm=.6)

def test_preflight_rejects_nonpositive_box():
    with pytest.raises(AnalysisError): preflight_box(initial_box_nm=(3,0,3),cutoff_nm=.5)

def test_virial_pressure_is_not_target():
    p=pressure_from_virial(kinetic_energy_kj_mol=10,virial_kj_mol=-2,volume_nm3=2,temperature_kelvin=300,particle_count=10)
    assert p['status']=='measured_virial' and p['independent_of_target']

def test_pressure_unavailable_is_explicit():
    p=pressure_unavailable(reason='backend exposes no validated virial')
    assert p['status']=='pressure_unavailable' and p['measured_pressure_bar'] is None

def test_invalid_virial_inputs():
    with pytest.raises(AnalysisError): pressure_from_virial(kinetic_energy_kj_mol=1,virial_kj_mol=1,volume_nm3=0,temperature_kelvin=300,particle_count=2)
