"""Torsion sampling strategies independent of polymer chemistry."""

import random
from abc import ABC, abstractmethod


class TorsionSampler(ABC):
    """Sample a rotation increment about an inter-repeat bond, in radians."""

    @abstractmethod
    def sample(self, rng: random.Random, *, repeat_index: int, trial: int) -> float:
        """Return a torsion using the provided local RNG."""


class UniformTorsionSampler(TorsionSampler):
    """Uniformly sample rotation increments over ``[-pi, pi)``."""

    def sample(self, rng: random.Random, *, repeat_index: int, trial: int) -> float:
        del repeat_index, trial
        return rng.uniform(-3.141592653589793, 3.141592653589793)
