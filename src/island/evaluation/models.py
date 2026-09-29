"""Backend-independent single-point contract; all public coordinates are angstroms."""

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType
from typing import Protocol, runtime_checkable

from island.exceptions import EvaluationError


def fingerprint(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode()
    ).hexdigest()


@dataclass(frozen=True)
class EvaluationResult:
    """Owned immutable single point; forces are F = -dE/dR in angstrom units."""

    potential_energy: float
    energy_components: Mapping[str, float]
    forces: Mapping[int, tuple[float, float, float]]
    coordinate_fingerprint: str
    parameter_fingerprint: str
    model_fingerprint: str
    evaluation_fingerprint: str
    backend_name: str
    backend_version: str
    platform: str
    settings: Mapping[str, str | bool | float]
    energy_unit: str = field(default="kJ/mol", init=False)
    force_unit: str = field(default="kJ/(mol*angstrom)", init=False)
    coordinate_unit: str = field(default="angstrom", init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)

    def __post_init__(self) -> None:
        if any(
            not isinstance(value, (str, bool, int, float))
            for value in self.settings.values()
        ):
            raise EvaluationError(
                "Calculation settings must contain immutable scalar values"
            )
        components = dict(self.energy_components)
        forces = {site: tuple(vector) for site, vector in self.forces.items()}
        if (
            not isfinite(self.potential_energy)
            or any(not isfinite(value) for value in components.values())
            or any(
                len(vector) != 3 or not all(isfinite(x) for x in vector)
                for vector in forces.values()
            )
        ):
            raise EvaluationError(
                "Backend returned nonfinite energy or malformed forces"
            )
        object.__setattr__(self, "energy_components", MappingProxyType(components))
        object.__setattr__(self, "forces", MappingProxyType(forces))
        object.__setattr__(self, "settings", MappingProxyType(dict(self.settings)))

    def __deepcopy__(self, memo: dict) -> "EvaluationResult":
        # Every value is immutable; callers cannot mutate an earlier frame.
        return self


@runtime_checkable
class PotentialEvaluator(Protocol):
    """A bound potential that evaluates explicitly labelled Cartesian frames."""

    def evaluate(
        self,
        coordinates: Mapping[int, object] | None = None,
        *,
        coordinate_unit: str = "angstrom",
    ) -> EvaluationResult: ...
