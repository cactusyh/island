"""Fail-closed pressure validation for actual periodic force-field systems."""
from math import isfinite

from island.exceptions import AnalysisError

BAR_PER_KJ_MOL_NM3 = 16.60539067

def component_pressure(*, kinetic_kj_mol, bonded_virial_kj_mol, nonbonded_virial_kj_mol, volume_nm3, long_range_virial_kj_mol=None):
    vals=(kinetic_kj_mol,bonded_virial_kj_mol,nonbonded_virial_kj_mol,volume_nm3)
    if any(not isfinite(float(v)) for v in vals) or volume_nm3<=0: raise AnalysisError('Invalid actual-system virial inputs')
    lr=0.0 if long_range_virial_kj_mol is None else float(long_range_virial_kj_mol)
    if not isfinite(lr): raise AnalysisError('Invalid long-range virial')
    total=(2*float(kinetic_kj_mol)+float(bonded_virial_kj_mol)+float(nonbonded_virial_kj_mol)+lr)/(3*float(volume_nm3))
    return {'kinetic_pressure_bar':2*float(kinetic_kj_mol)/(3*volume_nm3)*BAR_PER_KJ_MOL_NM3,'bonded_pressure_bar':float(bonded_virial_kj_mol)/(3*volume_nm3)*BAR_PER_KJ_MOL_NM3,'nonbonded_pressure_bar':float(nonbonded_virial_kj_mol)/(3*volume_nm3)*BAR_PER_KJ_MOL_NM3,'long_range_pressure_bar':lr/(3*volume_nm3)*BAR_PER_KJ_MOL_NM3,'total_pressure_bar':total*BAR_PER_KJ_MOL_NM3}

def validate_actual_pressure(*, force_virial_bar, volume_derivative_bar, components, tolerance_bar=0.1, long_range_validated=False, audit=None):
    if tolerance_bar<=0 or not isfinite(float(tolerance_bar)): raise AnalysisError('Invalid pressure tolerance')
    if not long_range_validated: return {'status':'pressure_unavailable','scientific_acceptance':False,'reason':'long-range electrostatic virial is not independently validated','force_virial_bar':force_virial_bar,'volume_derivative_bar':volume_derivative_bar,'components':components}
    residual=float(force_virial_bar)-float(volume_derivative_bar)
    passed=abs(residual)<=tolerance_bar
    return {'status':'validated' if passed else 'blocked','scientific_acceptance':passed,'force_virial_bar':force_virial_bar,'volume_derivative_bar':volume_derivative_bar,'residual_bar':residual,'tolerance_bar':tolerance_bar,'components':components,'audit':audit or {}}

def require_safe_execution(preflight, pressure):
    if preflight.get('status')!='ok': raise AnalysisError('cutoff_guard_blocked: unsafe periodic box')
    if pressure.get('status')!='validated': raise AnalysisError('pressure_unavailable: backend pressure is not validated')
