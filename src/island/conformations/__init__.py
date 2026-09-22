"""Coordinate-only molecular conformation generation."""

from typing import TYPE_CHECKING, Any

from island.conformations.base import ConformationGenerator
from island.conformations.random_walk import (
    RetryPolicy,
    SelfAvoidingRandomWalkGenerator,
)
from island.conformations.result import ConformationResult
from island.conformations.sterics import (
    FixedDistanceStericPolicy,
    StericPolicy,
    VanDerWaalsStericPolicy,
)
from island.conformations.torsions import TorsionSampler, UniformTorsionSampler
from island.core import MolecularSystem
from island.exceptions import UnsupportedConformationError

if TYPE_CHECKING:
    from island.conformations.local_templates import (
        ConnectionBondLengthPolicy,
        CovalentRadiiBondLengthPolicy,
        LocalTemplate,
        LocalTemplateConformationGenerator,
    )

_LOCAL_TEMPLATE_EXPORTS = {
    "ConnectionBondLengthPolicy",
    "CovalentRadiiBondLengthPolicy",
    "LocalTemplate",
    "LocalTemplateConformationGenerator",
}

__all__ = [
    "ConformationGenerator",
    "ConformationResult",
    "ConnectionBondLengthPolicy",
    "CovalentRadiiBondLengthPolicy",
    "FixedDistanceStericPolicy",
    "LocalTemplate",
    "LocalTemplateConformationGenerator",
    "RetryPolicy",
    "SelfAvoidingRandomWalkGenerator",
    "StericPolicy",
    "TorsionSampler",
    "UniformTorsionSampler",
    "VanDerWaalsStericPolicy",
    "generate_polymer_conformation",
]


def __getattr__(name: str) -> Any:
    if name in _LOCAL_TEMPLATE_EXPORTS:
        from island.conformations import local_templates

        return getattr(local_templates, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def generate_polymer_conformation(
    system: MolecularSystem,
    *,
    method: str = "random_walk",
    seed: int = 2026,
    **generator_settings: Any,
) -> ConformationResult:
    """Generate coordinates non-destructively using the requested method."""
    if method != "random_walk":
        raise UnsupportedConformationError(
            f"Unsupported conformation generation method: {method!r}"
        )
    generator = SelfAvoidingRandomWalkGenerator(seed=seed, **generator_settings)
    return generator.generate(system)
