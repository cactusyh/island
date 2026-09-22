"""Simulation-cell abstractions."""

from dataclasses import dataclass


@dataclass
class SimulationBox:
    """Orthorhombic simulation box, extensible to other cell geometries later."""

    lx: float
    ly: float
    lz: float
    periodic: tuple[bool, bool, bool] = (True, True, True)

    def __post_init__(self) -> None:
        if min(self.lx, self.ly, self.lz) <= 0:
            raise ValueError("Simulation box lengths must be positive")
        if len(self.periodic) != 3:
            raise ValueError("Periodicity must contain three flags")

    @property
    def lengths(self) -> tuple[float, float, float]:
        return (self.lx, self.ly, self.lz)
