"""Owned local-minimization options and diagnostic results, in angstrom units."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType

import numpy as np

from island.core import Coordinates, MolecularSystem
from island.evaluation.models import EvaluationResult, fingerprint
from island.exceptions import MinimizationInputError


def coordinate_hash(coordinates):
    return fingerprint(
        {
            "unit": "angstrom",
            "sites": [(site, list(coordinates[site])) for site in sorted(coordinates)],
        }
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
        for key in ("initial_coordinates", "coordinates"):
            object.__setattr__(
                self,
                key,
                MappingProxyType(
                    {
                        site: tuple(float(v) for v in xyz)
                        for site, xyz in getattr(self, key).items()
                    }
                ),
            )
        object.__setattr__(self, "history", tuple(self.history))

    def __deepcopy__(self, memo):
        return self  # all owned contents are immutable

    @property
    def initial_energy(self):
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.potential_energy
        )

    @property
    def final_energy(self):
        return (
            None
            if self.final_evaluation is None
            else self.final_evaluation.potential_energy
        )

    @property
    def initial_fmax(self):
        return force_metrics(self.initial_evaluation)[0]

    @property
    def final_fmax(self):
        return force_metrics(self.final_evaluation)[0]

    @property
    def initial_rms_force(self):
        return force_metrics(self.initial_evaluation)[1]

    @property
    def final_rms_force(self):
        return force_metrics(self.final_evaluation)[1]

    @property
    def forces(self):
        return (
            MappingProxyType({})
            if self.final_evaluation is None
            else self.final_evaluation.forces
        )

    @property
    def energy_components(self):
        return (
            MappingProxyType({})
            if self.final_evaluation is None
            else self.final_evaluation.energy_components
        )

    @property
    def model_fingerprint(self):
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.model_fingerprint
        )

    @property
    def parameter_fingerprint(self):
        return (
            None
            if self.initial_evaluation is None
            else self.initial_evaluation.parameter_fingerprint
        )

    @property
    def initial_coordinate_fingerprint(self):
        return coordinate_hash(self.initial_coordinates)

    @property
    def coordinate_fingerprint(self):
        return coordinate_hash(self.coordinates)

    def to_system(self, system: MolecularSystem, *, allow_unconverged=False):
        """Copy compatible chemistry; diagnostic coordinate application is opt-in."""
        if (
            not isinstance(system, MolecularSystem)
            or system_identity(system) != self.system_fingerprint
        ):
            raise MinimizationInputError(
                "System chemistry/provenance differs from the minimization input"
            )
        if not self.converged and not allow_unconverged:
            raise MinimizationInputError(
                "Unconverged result requires allow_unconverged=True"
            )
        if (
            self.final_evaluation is None
            or self.stereochemistry not in ("passed", "not_assigned")
            or self.final_evaluation.coordinate_fingerprint
            != self.coordinate_fingerprint
            or self.final_evaluation.model_fingerprint != self.model_fingerprint
            or self.final_evaluation.parameter_fingerprint != self.parameter_fingerprint
        ):
            raise MinimizationInputError(
                "No valid compatible evaluated coordinates are available to apply"
            )
        if self.converged and (
            not self.final_evaluation_verified
            or self.final_fmax > self.options.force_tolerance
            or self.final_energy
            > self.initial_energy + self.options.energy_increase_tolerance
        ):
            raise MinimizationInputError(
                "Result does not satisfy the public convergence criterion"
            )
        copied = deepcopy(system)
        copied.coordinates = Coordinates(self.coordinates)
        copied.metadata["local_minimization"] = {
            "original_coordinate_source": deepcopy(
                system.metadata.get("coordinate_source")
            ),
            "initial_coordinate_fingerprint": coordinate_hash(self.initial_coordinates),
            "coordinate_fingerprint": self.coordinate_fingerprint,
            "model_fingerprint": self.model_fingerprint,
            "termination_reason": self.termination_reason,
            "converged": self.converged,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        copied.metadata["coordinate_source"] = (
            "local_minimization" if self.converged else "local_minimization_diagnostic"
        )
        copied.metadata["production_validated"] = False
        copied.metadata["simulation_readiness"] = "not_established"
        return copied
