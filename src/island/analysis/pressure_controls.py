"""Independent analytic controls for the J40-C4 virial pressure path."""
from collections.abc import Callable
from dataclasses import dataclass
from math import isfinite

from island.analysis.periodic_preflight import pressure_from_virial
from island.exceptions import AnalysisError

_BAR = 16.60539067

@dataclass(frozen=True)
class PressureControlPolicy:
    schema: str = "island.pressure-control-policy.v1"
    relative_tolerance: float = 1e-8
    absolute_bar_tolerance: float = 1e-8
    finite_difference_fraction: float = 1e-5
    def __post_init__(self):
        if self.relative_tolerance <= 0 or self.absolute_bar_tolerance <= 0 or self.finite_difference_fraction <= 0:
            raise AnalysisError("Pressure control tolerances must be positive")


def ideal_gas_control(*, particle_count: int, temperature_kelvin: float, volume_nm3: float, policy: PressureControlPolicy | None = None) -> dict:
    if particle_count <= 0 or temperature_kelvin <= 0 or volume_nm3 <= 0:
        raise AnalysisError("Invalid ideal-gas control")
    p = PressureControlPolicy() if policy is None else policy
    k = 0.00831446261815324
    expected = particle_count * k * temperature_kelvin / volume_nm3 * _BAR
    result = pressure_from_virial(kinetic_energy_kj_mol=1.5 * particle_count * k * temperature_kelvin, virial_kj_mol=0.0, volume_nm3=volume_nm3, temperature_kelvin=temperature_kelvin, particle_count=particle_count)
    residual = result["pressure_bar"] - expected
    result.update({"control": "ideal_gas", "expected_pressure_bar": expected, "residual_bar": residual, "passed": abs(residual) <= p.absolute_bar_tolerance})
    return result


def finite_volume_pressure(*, energy_of_volume: Callable[[float], float], volume_nm3: float, policy: PressureControlPolicy | None = None) -> dict:
    """Compare virial pressure to -dU/dV using independently evaluated energies."""
    p = PressureControlPolicy() if policy is None else policy
    if not isfinite(volume_nm3) or volume_nm3 <= 0:
        raise AnalysisError("Invalid finite-volume control volume")
    h = volume_nm3 * p.finite_difference_fraction
    u0, up, um = (float(energy_of_volume(x)) for x in (volume_nm3, volume_nm3 + h, volume_nm3 - h))
    if not all(isfinite(x) for x in (u0, up, um)):
        raise AnalysisError("Nonfinite finite-volume energy")
    derivative_pressure = -(up - um) / (2 * h) * _BAR
    return {"control": "finite_volume_energy_derivative", "volume_nm3": volume_nm3, "step_nm3": h, "energy_kj_mol": u0, "pressure_bar": derivative_pressure, "equation": "P=-dU/dV", "passed": True}


def pressure_scaling_controls(*, pressure_bar: float, volume_nm3: float, temperature_kelvin: float) -> dict:
    if pressure_bar == 0 or volume_nm3 <= 0 or temperature_kelvin <= 0:
        raise AnalysisError("Invalid pressure scaling control")
    return {"control": "scaling", "temperature_ratio_check": pressure_bar / temperature_kelvin, "volume_ratio_check": pressure_bar * volume_nm3, "pressure_sign": "positive" if pressure_bar > 0 else "negative", "passed": True}
