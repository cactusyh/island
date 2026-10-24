"""Fail-closed periodic box/cutoff and pressure diagnostics for J40-C3.

The checks are pure and require the caller to provide measurements from the
actual engine. A barostat target is never treated as measured pressure.
"""
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from math import isfinite

from island.exceptions import AnalysisError


@dataclass(frozen=True)
class PeriodicPreflightPolicy:
    schema: str = "island.periodic-preflight-policy.v1"
    safety_margin_nm: float = 0.05
    neighbor_skin_nm: float = 0.0
    require_orthorhombic: bool = True

    def __post_init__(self):
        if not isfinite(self.safety_margin_nm) or self.safety_margin_nm < 0:
            raise AnalysisError("Invalid periodic safety margin")
        if not isfinite(self.neighbor_skin_nm) or self.neighbor_skin_nm < 0:
            raise AnalysisError("Invalid neighbor-list skin")

    @property
    def payload(self):
        return asdict(self)


def _lengths(lengths: Sequence[float]) -> tuple[float, ...]:
    if len(lengths) != 3 or any(not isfinite(float(x)) or float(x) <= 0 for x in lengths):
        raise AnalysisError("Periodic box lengths must be three positive finite values")
    return tuple(float(x) for x in lengths)


def preflight_box(*, initial_box_nm: Sequence[float], cutoff_nm: float,
                  switching_distance_nm: float | None = None,
                  policy: PeriodicPreflightPolicy | None = None,
                  predicted_min_box_nm: Sequence[float] | None = None) -> dict:
    """Validate initial and predicted NPT boxes before any dynamics starts."""
    policy = policy or PeriodicPreflightPolicy()
    box = _lengths(initial_box_nm)
    minimum = _lengths(predicted_min_box_nm or box)
    if not isfinite(cutoff_nm) or cutoff_nm <= 0:
        raise AnalysisError("Nonbonded cutoff must be positive and finite")
    switch = cutoff_nm if switching_distance_nm is None else float(switching_distance_nm)
    if not isfinite(switch) or switch <= 0 or switch > cutoff_nm:
        raise AnalysisError("Switching distance must be finite, positive and <= cutoff")
    required = float(cutoff_nm) + policy.safety_margin_nm + policy.neighbor_skin_nm
    half_initial = min(box) / 2
    half_predicted = min(minimum) / 2
    diagnostics = {
        "schema": "island.periodic-preflight.v1",
        "initial_box_nm": list(box), "predicted_min_box_nm": list(minimum),
        "minimum_box_length_nm": min(minimum), "half_initial_box_nm": half_initial,
        "half_predicted_box_nm": half_predicted, "cutoff_nm": float(cutoff_nm),
        "switching_distance_nm": switch, "neighbor_skin_nm": policy.neighbor_skin_nm,
        "safety_margin_nm": policy.safety_margin_nm, "required_half_box_nm": required,
        "safe_initial": required < half_initial, "safe_predicted": required < half_predicted,
        "status": "ok" if required < half_initial and required < half_predicted else "blocked",
        "diagnostic": "periodic cutoff preflight passed" if required < half_initial and required < half_predicted else "cutoff plus safety margin is not below half the predicted periodic box",
    }
    return diagnostics


def pressure_from_virial(*, kinetic_energy_kj_mol: float, virial_kj_mol: float,
                         volume_nm3: float, temperature_kelvin: float,
                         particle_count: int, boltzmann_kj_mol_kelvin: float = 0.00831446261815324) -> dict:
    """Return P=(2K+W)/(3V), with units kJ mol^-1 nm^-3 and bar conversion."""
    vals = (kinetic_energy_kj_mol, virial_kj_mol, volume_nm3, temperature_kelvin)
    if any(not isfinite(float(x)) for x in vals) or volume_nm3 <= 0 or particle_count <= 0:
        raise AnalysisError("Invalid virial pressure inputs")
    # Temperature is retained as an independent consistency diagnostic; it is
    # not used to manufacture pressure from the target.
    pressure_kj_mol_nm3 = (2.0 * kinetic_energy_kj_mol + virial_kj_mol) / (3.0 * volume_nm3)
    return {"status": "measured_virial", "equation": "P=(2K+W)/(3V)",
            "pressure_kj_mol_nm3": pressure_kj_mol_nm3,
            "pressure_bar": pressure_kj_mol_nm3 * 16.60539067,
            "temperature_kelvin": float(temperature_kelvin), "particle_count": int(particle_count),
            "volume_nm3": float(volume_nm3), "independent_of_target": True}


def pressure_unavailable(*, reason: str, target_pressure_bar: float | None = None) -> dict:
    if not reason.strip():
        raise AnalysisError("Pressure-unavailable diagnostic requires a reason")
    return {"status": "pressure_unavailable", "reason": reason,
            "target_pressure_bar": target_pressure_bar,
            "measured_pressure_bar": None, "scientific_acceptance": False}
