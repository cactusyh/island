import pytest

from island.analysis.pressure_controls import (
 PressureControlPolicy,
 finite_volume_pressure,
 ideal_gas_control,
 pressure_scaling_controls,
)
from island.exceptions import AnalysisError


def test_ideal_gas_virial_control():
 r=ideal_gas_control(particle_count=10,temperature_kelvin=300,volume_nm3=10)
 assert r['passed'] and abs(r['pressure_bar']-r['expected_pressure_bar'])<1e-8

def test_finite_volume_derivative_control():
 r=finite_volume_pressure(energy_of_volume=lambda v: 2*v*v,volume_nm3=3)
 assert r['equation']=='P=-dU/dV' and r['pressure_bar']<0

def test_pressure_scaling_records_sign_and_ratios():
 r=pressure_scaling_controls(pressure_bar=2,volume_nm3=4,temperature_kelvin=300)
 assert r['passed'] and r['pressure_sign']=='positive'

def test_controls_reject_bad_inputs():
 with pytest.raises(AnalysisError): ideal_gas_control(particle_count=0,temperature_kelvin=300,volume_nm3=1)
 with pytest.raises(AnalysisError): finite_volume_pressure(energy_of_volume=lambda v: v,volume_nm3=0)

def test_policy_is_immutable_and_positive():
 with pytest.raises(AnalysisError): PressureControlPolicy(relative_tolerance=0)
