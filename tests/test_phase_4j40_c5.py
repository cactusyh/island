import pytest

from island.analysis.backend_pressure import (
 component_pressure,
 require_safe_execution,
 validate_actual_pressure,
)
from island.exceptions import AnalysisError


def test_component_pressure_accounting():
 r=component_pressure(kinetic_kj_mol=3,bonded_virial_kj_mol=1,nonbonded_virial_kj_mol=-2,volume_nm3=2,long_range_virial_kj_mol=1)
 assert r['total_pressure_bar']>0 and r['kinetic_pressure_bar']>0

def test_long_range_unvalidated_is_blocked():
 r=validate_actual_pressure(force_virial_bar=1,volume_derivative_bar=1,components={},long_range_validated=False)
 assert r['status']=='pressure_unavailable' and not r['scientific_acceptance']

def test_validated_estimators_need_agreement():
 r=validate_actual_pressure(force_virial_bar=1,volume_derivative_bar=1.05,components={},tolerance_bar=.1,long_range_validated=True)
 assert r['status']=='validated'

def test_disagreement_blocks():
 r=validate_actual_pressure(force_virial_bar=1,volume_derivative_bar=2,components={},tolerance_bar=.1,long_range_validated=True)
 assert r['status']=='blocked'

def test_unsafe_execution_rejected():
 with pytest.raises(AnalysisError,match='cutoff_guard_blocked'): require_safe_execution({'status':'blocked'},{'status':'validated'})
 with pytest.raises(AnalysisError,match='pressure_unavailable'): require_safe_execution({'status':'ok'},{'status':'pressure_unavailable'})
