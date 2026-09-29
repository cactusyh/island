"""Owned local-minimization options and diagnostic results, in angstrom units."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType

import numpy as np

from island.core import Coordinates, MolecularSystem
from island.core.coordinate_provenance import coordinate_hash
from island.evaluation.models import EvaluationResult, fingerprint
from island.exceptions import (
    InvalidMinimizationResultError,
    IslandError,
    MinimizationInputError,
    MinimizationUnavailableError,
)


def system_identity(system):
    payload = system.to_dict()
    payload.pop("coordinates")
    payload["sites"].sort(key=lambda site: site["id"])
    for bond in payload["bonds"]:
        bond["site1"], bond["site2"] = sorted((bond["site1"], bond["site2"]))
    payload["bonds"].sort(key=lambda bond: (bond["site1"], bond["site2"]))
    payload["impropers"] = sorted(
        (i.site1, i.site2, i.site3, i.site4) for i in system.topology.impropers
    )
    return fingerprint(payload)


def force_metrics(result):
    if result is None:
        return None, None
    norms = np.hypot.reduce(np.array(list(result.forces.values()), dtype=float), axis=1)
    maximum = float(np.max(norms))
    rms = (
        0.0
        if maximum == 0
        else maximum * float(np.sqrt(np.mean((norms / maximum) ** 2)))
    )
    return maximum, rms


@dataclass(frozen=True)
class MinimizationOptions:
    """Force tolerance is max atomic vector norm in kJ/(mol*angstrom)."""

    force_tolerance: float = 0.1
    max_iterations: int = 500
    max_evaluations: int = 2000
    max_line_search_steps: int = 20
    energy_change_tolerance: float = 1e-12
    energy_increase_tolerance: float = 1e-8  # absolute kJ/mol

    def __post_init__(self):
        for key in ("max_iterations", "max_evaluations", "max_line_search_steps"):
            if type(getattr(self, key)) is not int or getattr(self, key) <= 0:
                raise MinimizationInputError(f"{key} must be a finite positive integer")
        for key in (
            "force_tolerance",
            "energy_change_tolerance",
            "energy_increase_tolerance",
        ):
            value = getattr(self, key)
            if (
                not isinstance(value, (int, float))
                or isinstance(value, bool)
                or not isfinite(value)
                or value < 0
                or (key == "force_tolerance" and value == 0)
            ):
                raise MinimizationInputError(
                    f"{key} must be finite and nonnegative (force tolerance positive)"
                )


@dataclass(frozen=True)
class MinimizationStep:
    iteration: int
    evaluations: int
    energy: float
    fmax: float
    rms_force: float
    coordinate_fingerprint: str


@dataclass(frozen=True)
class MinimizationResult:
    initial_coordinates: Mapping[int, tuple[float, float, float]]
    coordinates: Mapping[int, tuple[float, float, float]]
    initial_evaluation: EvaluationResult | None
    final_evaluation: EvaluationResult | None
    converged: bool
    termination_reason: str
    iterations: int
    evaluations: int
    optimizer_message: str
    optimizer_success: bool | None
    optimizer_version: str
    options: MinimizationOptions
    history: tuple[MinimizationStep, ...]
    system_fingerprint: str
    stereochemistry: str
    returned_state: str
    final_evaluation_verified: bool
    optimizer_name: str = field(default="SciPy L-BFGS-B", init=False)
    coordinate_unit: str = field(default="angstrom", init=False)
    energy_unit: str = field(default="kJ/mol", init=False)
    force_unit: str = field(default="kJ/(mol*angstrom)", init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)

    def __post_init__(self):
        from .integrity import number, require

        def owned_vector(xyz):
            vector = tuple(xyz)
            require(
                len(vector) == 3 and all(number(v) for v in vector),
                "Coordinates require finite three-component numeric vectors",
            )
            return tuple(float(v) for v in vector)

        try:
            for key in ("initial_coordinates", "coordinates"):
                object.__setattr__(
                    self,
                    key,
                    MappingProxyType(
                        {
                            site: owned_vector(xyz)
                            for site, xyz in getattr(self, key).items()
                        }
                    ),
                )
            object.__setattr__(self, "history", tuple(self.history))
        except (AttributeError, TypeError, ValueError, OverflowError) as error:
            raise InvalidMinimizationResultError(
                f"Malformed result contents: {error}"
            ) from error

    def validate_integrity(self):
        """Validate owned record content, not authenticity of caller-supplied energies."""
        from .integrity import validate

        try:
            validate(self)
        except InvalidMinimizationResultError:
            raise
        except (
            IslandError,
            AttributeError,
            KeyError,
            TypeError,
            ValueError,
            OverflowError,
            IndexError,
        ) as error:
            raise InvalidMinimizationResultError(
                f"Malformed minimization result: {error}"
            ) from error

    def __deepcopy__(self, memo):
        self.validate_integrity()
        return (
            self  # validation establishes that all supported nested data are immutable
        )

    @property
    def initial_energy(self):
        self.validate_integrity()
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.potential_energy
        )

    @property
    def final_energy(self):
        self.validate_integrity()
        return (
            None
            if self.final_evaluation is None
            else self.final_evaluation.potential_energy
        )

    @property
    def initial_fmax(self):
        self.validate_integrity()
        return force_metrics(self.initial_evaluation)[0]

    @property
    def final_fmax(self):
        self.validate_integrity()
        return force_metrics(self.final_evaluation)[0]

    @property
    def initial_rms_force(self):
        self.validate_integrity()
        return force_metrics(self.initial_evaluation)[1]

    @property
    def final_rms_force(self):
        self.validate_integrity()
        return force_metrics(self.final_evaluation)[1]

    @property
    def forces(self):
        self.validate_integrity()
        return (
            MappingProxyType({})
            if self.final_evaluation is None
            else self.final_evaluation.forces
        )

    @property
    def energy_components(self):
        self.validate_integrity()
        return (
            MappingProxyType({})
            if self.final_evaluation is None
            else self.final_evaluation.energy_components
        )

    @property
    def model_fingerprint(self):
        self.validate_integrity()
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.model_fingerprint
        )

    @property
    def parameter_fingerprint(self):
        self.validate_integrity()
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.parameter_fingerprint
        )

    @property
    def initial_coordinate_fingerprint(self):
        self.validate_integrity()
        return coordinate_hash(self.initial_coordinates)

    @property
    def coordinate_fingerprint(self):
        self.validate_integrity()
        return coordinate_hash(self.coordinates)

    def to_system(self, system: MolecularSystem, *, allow_unconverged=False):
        """Validate records and geometry, then copy; diagnostics require explicit opt-in."""
        from island.chemistry.coordinate_stereo import (
            assigned_cip_labels,
            validate_coordinate_stereochemistry,
        )
        from island.core.coordinate_provenance import (
            previous_coordinate_source,
            updated_coordinate_metadata,
        )

        self.validate_integrity()
        if type(allow_unconverged) is not bool:
            raise MinimizationInputError(
                "allow_unconverged must be an explicit boolean"
            )
        try:
            if not isinstance(system, MolecularSystem):
                raise InvalidMinimizationResultError(
                    "Application requires a MolecularSystem"
                )
            system.validate()
            if set(self.coordinates) != set(system.topology.sites):
                raise InvalidMinimizationResultError(
                    "Coordinate coverage differs from target topology"
                )
            if system_identity(system) != self.system_fingerprint:
                raise InvalidMinimizationResultError(
                    "System chemistry/provenance differs from the minimization input"
                )
        except (IslandError, AttributeError, KeyError, TypeError, ValueError) as error:
            raise InvalidMinimizationResultError(str(error)) from error
        if not self.converged and not allow_unconverged:
            raise MinimizationInputError(
                "Unconverged result requires allow_unconverged=True"
            )
        if self.final_evaluation is None:
            raise MinimizationInputError(
                "No valid compatible evaluated coordinates are available to apply"
            )
        try:
            expected = assigned_cip_labels(system)
            for label, coordinates in (
                ("initial", self.initial_coordinates),
                ("returned", self.coordinates),
            ):
                validate_coordinate_stereochemistry(
                    system, coordinates, expected, stage=f"apply minimization {label}"
                )
            if self.stereochemistry != ("passed" if expected else "not_assigned"):
                raise InvalidMinimizationResultError(
                    "Stereochemistry outcome contradicts target graph"
                )
        except ImportError as error:
            raise MinimizationUnavailableError(
                "RDKit is required to validate assigned tetrahedral stereochemistry"
            ) from error
        except IslandError as error:
            raise InvalidMinimizationResultError(str(error)) from error
        copied = deepcopy(system)
        copied.coordinates = Coordinates(self.coordinates)
        source = (
            "local_minimization" if self.converged else "local_minimization_diagnostic"
        )
        record = {
            "coordinate_source": source,
            "original_coordinate_source": previous_coordinate_source(system.metadata),
            "initial_coordinate_fingerprint": self.initial_coordinate_fingerprint,
            "coordinate_fingerprint": self.coordinate_fingerprint,
            "model_fingerprint": self.model_fingerprint,
            "parameter_fingerprint": self.parameter_fingerprint,
            "termination_reason": self.termination_reason,
            "converged": self.converged,
            "final_evaluation_verified": self.final_evaluation_verified,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        previous = coordinate_hash(
            {
                site: tuple(system.coordinates.get(site))
                for site in system.topology.sites
            }
        )
        copied.metadata = updated_coordinate_metadata(
            system, record, previous_fingerprint=previous, minimization=record
        )
        copied.metadata["production_validated"] = False
        copied.metadata["simulation_readiness"] = "not_established"
        copied.validate()
        return copied
