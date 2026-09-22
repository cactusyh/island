"""Coordinate-only molecular conformation generation."""

from typing import Any

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

__all__ = [
    "ConformationGenerator",
    "ConformationResult",
    "FixedDistanceStericPolicy",
    "RetryPolicy",
    "SelfAvoidingRandomWalkGenerator",
    "StericPolicy",
    "TorsionSampler",
    "UniformTorsionSampler",
    "VanDerWaalsStericPolicy",
    "generate_polymer_conformation",
]


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
