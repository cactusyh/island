"""Coordinates keyed by stable site ID."""

from collections.abc import Iterable, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from island.exceptions import CoordinateError, ValidationError


class Coordinates:
    """A mapping from site IDs to three-dimensional Cartesian positions."""

    def __init__(self, positions: Mapping[int, ArrayLike] | None = None) -> None:
        self._positions: dict[int, NDArray[np.float64]] = {}
        for site_id, position in (positions or {}).items():
            self.set(site_id, position)

    def get(self, site_id: int) -> NDArray[np.float64]:
        try:
            return self._positions[site_id].copy()
        except KeyError as error:
            raise CoordinateError(f"No coordinates for site id: {site_id}") from error

    def set(self, site_id: int, position: ArrayLike) -> None:
        array = np.asarray(position, dtype=float)
        if array.shape != (3,):
            raise CoordinateError("A coordinate must have exactly three components")
        if not np.all(np.isfinite(array)):
            raise CoordinateError("Coordinate components must be finite")
        self._positions[site_id] = array.copy()

    def translate(
        self, displacement: ArrayLike, site_ids: Iterable[int] | None = None
    ) -> None:
        vector = np.asarray(displacement, dtype=float)
        if vector.shape != (3,):
            raise CoordinateError("A translation must have exactly three components")
        ids = list(self._positions) if site_ids is None else list(site_ids)
        for site_id in ids:
            if site_id not in self._positions:
                raise CoordinateError(f"No coordinates for site id: {site_id}")
            self._positions[site_id] += vector

    def copy(self) -> "Coordinates":
        return Coordinates(self._positions)

    def validate(self, topology: object) -> None:
        topology_ids = set(topology.sites)
        coordinate_ids = set(self._positions)
        missing = topology_ids - coordinate_ids
        extra = coordinate_ids - topology_ids
        if missing or extra:
            details = []
            if missing:
                details.append(f"missing coordinates for {sorted(missing)}")
            if extra:
                details.append(f"coordinates for unknown sites {sorted(extra)}")
            raise ValidationError("; ".join(details))

    def __len__(self) -> int:
        return len(self._positions)
