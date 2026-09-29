"""Owned, full-step states for finite nonperiodic velocity Verlet."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from math import isfinite
from types import MappingProxyType

import numpy as np

from island.core import Coordinates, MolecularSystem
from island.core.coordinate_provenance import (
    coordinate_hash,
    previous_coordinate_source,
    updated_coordinate_metadata,
)
from island.evaluation import EvaluationResult
from island.evaluation.models import fingerprint
from island.exceptions import (
    DynamicsInputError,
    DynamicsUnavailableError,
    InvalidDynamicsResultError,
    IslandError,
)
from island.minimization.integrity import number
from island.minimization.models import force_metrics, system_identity

INTEGRATOR = "velocity_verlet_full_step_v1"


def owned_vectors(mapping):
    try:
        result = {}
        for site, vector in mapping.items():
            values = tuple(vector)
            if (
                type(site) is not int
                or len(values) != 3
                or not all(number(v) for v in values)
            ):
                raise InvalidDynamicsResultError(
                    "Require integer IDs and finite three-component vectors"
                )
            result[site] = tuple(float(v) for v in values)
        if not result:
            raise InvalidDynamicsResultError("State vectors cannot be empty")
        return MappingProxyType(result)
    except (AttributeError, TypeError, ValueError, OverflowError) as error:
        raise InvalidDynamicsResultError(f"Malformed state vectors: {error}") from error


def velocity_hash(velocities):
    return fingerprint(
        {
            "unit": "angstrom/ps",
            "sites": [(site, list(velocities[site])) for site in sorted(velocities)],
        }
    )


def kinetic_energy(masses, velocities):
    with np.errstate(over="raise", invalid="raise"):
        return float(
            0.005
            * sum(
                masses[site] * np.dot(velocities[site], velocities[site])
                for site in sorted(masses)
            )
        )


def dynamics_identity(options, masses, initial, system_fingerprint):
    return fingerprint(
        {
            "integrator": INTEGRATOR,
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
        }
    )


@dataclass(frozen=True)
class DynamicsOptions:
    timestep_fs: float
    steps: int
    max_evaluations: int
    max_frames: int
    recording_interval: int = 1
    max_energy_deviation: float = 1.0  # absolute kJ/mol, relative to initial total E

    def __post_init__(self):
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
                    "timestep_fs and requested duration must be finite and positive"
                )
            if not number(self.max_energy_deviation, nonnegative=True):
                raise DynamicsInputError(
                    "max_energy_deviation must be finite and nonnegative in kJ/mol"
                )
            object.__setattr__(self, "timestep_fs", float(self.timestep_fs))
            object.__setattr__(
                self, "max_energy_deviation", float(self.max_energy_deviation)
            )
        except (TypeError, ValueError, OverflowError) as error:
            raise DynamicsInputError(str(error)) from error


@dataclass(frozen=True)
class DynamicsFrame:
    step: int
    time_ps: float
    coordinates: Mapping[int, tuple[float, float, float]]
    velocities: Mapping[int, tuple[float, float, float]]
    evaluation: EvaluationResult | None
    kinetic_energy: float
    total_energy: float | None
    energy_deviation: float | None
    coordinate_fingerprint: str
    velocity_fingerprint: str
    coordinate_unit: str = field(default="angstrom", init=False)
    velocity_unit: str = field(default="angstrom/ps", init=False)
    energy_unit: str = field(default="kJ/mol", init=False)
    time_unit: str = field(default="ps", init=False)

    def __post_init__(self):
        object.__setattr__(self, "coordinates", owned_vectors(self.coordinates))
        object.__setattr__(self, "velocities", owned_vectors(self.velocities))

    def validate_integrity(self):
        from .integrity import validate_frame_content

        validate_frame_content(self)

    @property
    def potential_energy(self):
        self.validate_integrity()
        return None if self.evaluation is None else self.evaluation.potential_energy

    @property
    def fmax(self):
        self.validate_integrity()
        return force_metrics(self.evaluation)[0]

    def __deepcopy__(self, memo):
        # Validate standalone numerical/ownership content before sharing immutable data.
        from .integrity import validate_frame_content

        validate_frame_content(self)
        return self


@dataclass(frozen=True)
class DynamicsResult:
    options: DynamicsOptions
    masses: Mapping[int, float]
    frames: tuple[DynamicsFrame, ...]
    completed_steps: int
    evaluations: int
    failed_trial_evaluations: int
    final_evaluation_attempted: bool
    final_evaluation_verified: bool
    final_evaluation: EvaluationResult | None
    completed: bool
    termination_reason: str
    attempted_step: int | None
    failure_stage: str | None
    message: str
    failure_details: Mapping[str, str | float | int | bool]
    final_check_error: str | None
    max_abs_energy_deviation: float | None
    system_fingerprint: str
    dynamics_fingerprint: str
    stereochemistry: str
    integrator: str = field(default=INTEGRATOR, init=False)
    mass_unit: str = field(default="dalton", init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)

    def __post_init__(self):
        try:
            if not self.masses or any(
                type(site) is not int or not number(mass) or mass <= 0
                for site, mass in self.masses.items()
            ):
                raise InvalidDynamicsResultError(
                    "Masses require integer IDs and positive finite values"
                )
            object.__setattr__(
                self,
                "masses",
                MappingProxyType(
                    {site: float(mass) for site, mass in self.masses.items()}
                ),
            )
            object.__setattr__(self, "frames", tuple(self.frames))
            object.__setattr__(
                self, "failure_details", MappingProxyType(dict(self.failure_details))
            )
        except (AttributeError, TypeError, ValueError, OverflowError) as error:
            raise InvalidDynamicsResultError(str(error)) from error

    def validate_integrity(self):
        from .integrity import validate_result

        try:
            validate_result(self)
        except InvalidDynamicsResultError:
            raise
        except (
            IslandError,
            AttributeError,
            KeyError,
            TypeError,
            ValueError,
            OverflowError,
            IndexError,
            FloatingPointError,
        ) as error:
            raise InvalidDynamicsResultError(
                f"Malformed dynamics record: {error}"
            ) from error

    def __deepcopy__(self, memo):
        self.validate_integrity()
        return self

    @property
    def initial_state(self):
        self.validate_integrity()
        return self.frames[0]

    @property
    def final_state(self):
        self.validate_integrity()
        return self.frames[-1]

    @property
    def requested_steps(self):
        self.validate_integrity()
        return self.options.steps

    @property
    def model_fingerprint(self):
        self.validate_integrity()
        record = self.frames[0].evaluation
        return None if record is None else record.model_fingerprint

    @property
    def parameter_fingerprint(self):
        self.validate_integrity()
        record = self.frames[0].evaluation
        return None if record is None else record.parameter_fingerprint

    def to_system(self, system, *, allow_incomplete=False):
        """Apply coordinates only; velocities remain explicit state on this result."""
        from island.chemistry.coordinate_stereo import (
            assigned_cip_labels,
            validate_coordinate_stereochemistry,
        )

        self.validate_integrity()
        if type(allow_incomplete) is not bool:
            raise DynamicsInputError("allow_incomplete must be an explicit boolean")
        try:
            if not isinstance(system, MolecularSystem):
                raise InvalidDynamicsResultError("Application requires MolecularSystem")
            system.validate()
            if (
                set(system.topology.sites) != set(self.masses)
                or system_identity(system) != self.system_fingerprint
            ):
                raise InvalidDynamicsResultError(
                    "Target graph, masses, IDs or provenance differ"
                )
            if {site: system.topology.sites[site].mass for site in self.masses} != dict(
                self.masses
            ):
                raise InvalidDynamicsResultError("Target masses differ")
            if not self.completed and not allow_incomplete:
                raise DynamicsInputError(
                    "Incomplete dynamics requires allow_incomplete=True"
                )
            if self.frames[-1].evaluation is None:
                raise DynamicsInputError("No valid evaluated state can be applied")
            expected = assigned_cip_labels(system)
            for frame in self.frames:
                validate_coordinate_stereochemistry(
                    system, frame.coordinates, expected, stage="apply NVE coordinates"
                )
            if self.stereochemistry != ("passed" if expected else "not_assigned"):
                raise InvalidDynamicsResultError(
                    "Stereochemistry declaration contradicts target"
                )
            copied = deepcopy(system)
            frame = self.frames[-1]
            copied.coordinates = Coordinates(frame.coordinates)
            source = "nve_dynamics" if self.completed else "nve_dynamics_diagnostic"
            provenance = {
                "coordinate_source": source,
                "original_coordinate_source": previous_coordinate_source(
                    system.metadata
                ),
                "initial_coordinate_fingerprint": self.frames[0].coordinate_fingerprint,
                "coordinate_fingerprint": frame.coordinate_fingerprint,
                "velocity_fingerprint": frame.velocity_fingerprint,
                "dynamics_fingerprint": self.dynamics_fingerprint,
                "completed_steps": self.completed_steps,
                "time_ps": frame.time_ps,
                "completed": self.completed,
                "termination_reason": self.termination_reason,
                "molecular_system_is_complete_restart": False,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
            previous = coordinate_hash(
                {site: tuple(system.coordinates.get(site)) for site in self.masses}
            )
            copied.metadata = updated_coordinate_metadata(
                system, provenance, previous_fingerprint=previous, dynamics=provenance
            )
            copied.metadata.update(
                production_validated=False, simulation_readiness="not_established"
            )
            copied.validate()
            return copied
        except ImportError as error:
            raise DynamicsUnavailableError(
                "RDKit is required for assigned tetrahedral stereochemistry"
            ) from error
        except DynamicsInputError:
            raise
        except (IslandError, AttributeError, KeyError, TypeError, ValueError) as error:
            raise InvalidDynamicsResultError(str(error)) from error
