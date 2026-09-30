"""Versioned BAOAB contracts; synchronized endpoints, no NVE energy guard."""

from dataclasses import asdict, dataclass, field
from math import isfinite

from island.evaluation.models import fingerprint
from island.exceptions import DynamicsInputError
from island.minimization.integrity import number

from .models import DynamicsFrame, DynamicsResult
from .thermal import RNG_ALGORITHM, validate_temperature_seed

BAOAB = "baoab_full_step_v1"


@dataclass(frozen=True)
class LangevinOptions:
    timestep_fs: float
    temperature_kelvin: float
    friction_per_ps: float
    thermostat_seed: int
    steps: int
    max_evaluations: int
    max_frames: int
    recording_interval: int = 1

    def __post_init__(self):
        validate_temperature_seed(self.temperature_kelvin, self.thermostat_seed)
        try:
            for key in ("steps", "max_evaluations", "max_frames", "recording_interval"):
                if type(getattr(self, key)) is not int or getattr(self, key) <= 0:
                    raise DynamicsInputError(f"{key} must be a positive integer")
            if (
                not number(self.timestep_fs)
                or self.timestep_fs <= 0
                or not isfinite(self.timestep_fs * 0.001 * self.steps)
                or self.timestep_fs * 0.001 == 0
            ):
                raise DynamicsInputError(
                    "Require finite positive timestep and duration"
                )
            if not number(self.friction_per_ps, nonnegative=True):
                raise DynamicsInputError(
                    "friction_per_ps must be finite and nonnegative"
                )
            if not isfinite(self.friction_per_ps * (self.timestep_fs * 0.001)):
                raise DynamicsInputError("Friction times timestep must be finite")
            for key in ("timestep_fs", "temperature_kelvin", "friction_per_ps"):
                object.__setattr__(self, key, float(getattr(self, key)))
        except (TypeError, ValueError, OverflowError) as error:
            raise DynamicsInputError(str(error)) from error


@dataclass(frozen=True)
class LangevinFrame(DynamicsFrame):
    degrees_of_freedom: int = 0
    instantaneous_temperature_kelvin: float = 0.0
    temperature_unit: str = field(default="kelvin", init=False)


@dataclass(frozen=True)
class LangevinResult(DynamicsResult):
    options: LangevinOptions
    frames: tuple[LangevinFrame, ...]

    # The shared result envelope validates all existing bookkeeping. Thermal
    # validation adds RNG/DOF identity and omits only NVE's conservation guard.
    normal_draws: int = 0
    random_steps: int = 0
    numpy_version: str = ""
    integrator: str = field(default=BAOAB, init=False)
    rng_algorithm: str = field(default=RNG_ALGORITHM, init=False)

    def _coordinate_source(self):
        return "langevin_dynamics"


def langevin_identity(options, masses, initial, system_fingerprint, numpy_version):
    return fingerprint(
        {
            "integrator": BAOAB,
            "options": asdict(options),
            "masses_dalton": sorted(masses.items()),
            "system": system_fingerprint,
            "coordinates": initial.coordinate_fingerprint,
            "velocities": initial.velocity_fingerprint,
            "model": None
            if initial.evaluation is None
            else initial.evaluation.model_fingerprint,
            "parameters": None
            if initial.evaluation is None
            else initial.evaluation.parameter_fingerprint,
            "rng": RNG_ALGORITHM,
            "numpy": numpy_version,
            "dof_policy": "3N_no_removal_v1",
        }
    )
