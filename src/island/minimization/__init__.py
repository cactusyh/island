"""Bounded local geometry minimization; SciPy is imported only when invoked."""

from .engine import minimize_geometry
from .models import MinimizationOptions, MinimizationResult, MinimizationStep

__all__ = [
    "MinimizationOptions",
    "MinimizationResult",
    "MinimizationStep",
    "minimize_geometry",
]
