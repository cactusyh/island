"""NumPy-only geometry with explicit stable IDs, masses and length units."""

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from numbers import Real

import numpy as np

from island.exceptions import AnalysisError


@dataclass(frozen=True)
class GeometryOptions:
    selection: str = "all"
    site_ids: tuple[int, ...] | None = None
    weighting: str = "mass"
    endpoints: tuple[int, int] | None = None

    def __post_init__(self):
        if type(self.selection) is not str or self.selection not in (
            "all",
            "heavy",
            "explicit",
        ):
            raise AnalysisError("selection must be all, heavy, or explicit")
        if type(self.weighting) is not str or self.weighting not in ("mass", "uniform"):
            raise AnalysisError("weighting must be mass or uniform")
        if (self.selection == "explicit") != (self.site_ids is not None):
            raise AnalysisError("Only explicit selection requires site_ids")
        for name in ("site_ids", "endpoints"):
            value = getattr(self, name)
            if value is None:
                continue
            try:
                values = tuple(value)
            except TypeError as error:
                raise AnalysisError(f"Invalid {name}") from error
            if not values or any(type(s) is not int for s in values):
                raise AnalysisError(f"{name} must contain integer stable IDs")
            if name == "site_ids" and len(set(values)) != len(values):
                raise AnalysisError("Duplicate selected stable IDs")
            if name == "endpoints" and len(values) != 2:
                raise AnalysisError(
                    "endpoints must be a pair (equal IDs are permitted)"
                )
            object.__setattr__(
                self, name, tuple(sorted(values)) if name == "site_ids" else values
            )


@dataclass(frozen=True)
class GeometryResult:
    selected_ids: tuple[int, ...]
    selected_masses: tuple[float, ...]
    weighting: str
    normalized_weights: tuple[float, ...]
    center: tuple[float, float, float]
    gyration_tensor: tuple[tuple[float, ...], ...]
    eigenvalues: tuple[float, float, float]
    radius_of_gyration: float
    kappa_squared: float | None
    endpoint_ids: tuple[int, int] | None
    endpoint_convention: str
    end_to_end_distance: float | None
    re_over_rg: float | None
    diagnostics: tuple[str, ...]
    coordinate_unit: str = "angstrom"
    tensor_unit: str = "angstrom^2"
    mass_unit: str = "dalton"
    shape_unit: str = "dimensionless"

    def to_dict(self):
        return asdict(self)


def _inventory(mapping, label):
    if (
        not isinstance(mapping, Mapping)
        or not mapping
        or any(type(s) is not int for s in mapping)
    ):
        raise AnalysisError(f"{label} requires nonempty integer stable-ID mapping")
    return sorted(mapping)


def _number(x):
    return (
        isinstance(x, Real) and not isinstance(x, (bool, np.bool_)) and np.isfinite(x)
    )


def geometry_metrics(
    coordinates,
    masses,
    *,
    atomic_numbers=None,
    options=None,
    coordinate_unit="angstrom",
    mass_unit="dalton",
):
    """Compute selected-site geometry; coordinates/masses cover the full inventory.

    Tiny negative eigenvalues within 64*eps*max(abs(G)) are roundoff; anything
    more negative is an error. Eigenvalues are ascending, in angstrom squared.
    """
    options = GeometryOptions() if options is None else options
    if type(options) is not GeometryOptions:
        raise AnalysisError("Expected GeometryOptions")
    if coordinate_unit != "angstrom" or mass_unit != "dalton":
        raise AnalysisError("Explicit units must be angstrom and dalton")
    try:
        ids = _inventory(masses, "Masses")
        if _inventory(coordinates, "Coordinates") != ids:
            raise AnalysisError("Exact coordinate/mass site coverage required")
        if any(not _number(masses[s]) or masses[s] <= 0 for s in ids):
            raise AnalysisError(
                "Masses must be positive and finite, even for uniform weighting"
            )
        if any(isinstance(x, (bool, np.bool_)) for s in ids for x in coordinates[s]):
            raise AnalysisError("Boolean coordinate scalars are not supported")
        xyz = np.asarray([coordinates[s] for s in ids])
        if xyz.shape != (len(ids), 3) or xyz.dtype.kind not in "ifu":
            raise AnalysisError("Coordinates require finite numeric N x 3 vectors")
        xyz = xyz.astype(float)
        if not np.isfinite(xyz).all():
            raise AnalysisError("Coordinates must be finite")
        if atomic_numbers is not None and (
            _inventory(atomic_numbers, "Atomic numbers") != ids
            or any(
                type(z) is not int or not 1 <= z <= 118 for z in atomic_numbers.values()
            )
        ):
            raise AnalysisError("Atomic numbers must cover the authoritative inventory")
        if options.selection == "heavy":
            if atomic_numbers is None:
                raise AnalysisError("Heavy selection requires explicit atomic_numbers")
            selected = [s for s in ids if atomic_numbers[s] > 1]
        elif options.selection == "explicit":
            selected = list(options.site_ids)
            if not set(selected) <= set(ids):
                raise AnalysisError("Unknown selected stable IDs")
        else:
            selected = ids
        if not selected:
            raise AnalysisError("Empty atom selection")
        if options.endpoints is not None and not set(options.endpoints) <= set(ids):
            raise AnalysisError("Unknown endpoint stable ID")
        positions = xyz[[ids.index(s) for s in selected]]
        selected_masses = np.array([masses[s] for s in selected], dtype=float)
        diagnostics = []
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            weights = (
                selected_masses
                if options.weighting == "mass"
                else np.ones(len(selected))
            )
            # Scale before summation; avoid overflowing a finite mass inventory.
            weights = weights / weights.max()
            if np.any(weights == 0):
                raise AnalysisError("Mass dynamic range cannot be resolved safely")
            weights = weights / weights.sum()
            # Reference shift improves accuracy for translated configurations.
            relative = positions - positions[0]
            shift = np.sum(weights[:, None] * relative, axis=0)
            center = positions[0] + shift
            centered = relative - shift
            scale = float(np.max(np.abs(centered)))
            if scale == 0:
                tensor = np.zeros((3, 3))
                values = np.zeros(3)
                rg = 0.0
                kappa = None
                diagnostics.append("zero_spread: kappa_squared and Re/Rg undefined")
            else:
                reduced = centered / scale
                normalized = (reduced * weights[:, None]).T @ reduced
                normalized = 0.5 * (normalized + normalized.T)
                values = np.linalg.eigvalsh(normalized)
                tolerance = 64 * np.finfo(float).eps * float(np.max(np.abs(normalized)))
                if values[0] < -tolerance:
                    raise AnalysisError("Substantial negative gyration eigenvalue")
                if np.any(values < 0):
                    values = np.maximum(values, 0)
                    diagnostics.append(
                        "tiny_negative_eigenvalue_corrected_within_64eps_tensor_scale"
                    )
                fractions = values / values.sum()
                kappa = float(1.5 * np.dot(fractions, fractions) - 0.5)
                shape_tol = 64 * np.finfo(float).eps
                if not -shape_tol <= kappa <= 1 + shape_tol:
                    raise AnalysisError("Relative anisotropy outside physical range")
                if kappa < 0 or kappa > 1:
                    diagnostics.append("anisotropy_roundoff_corrected_within_64eps")
                    kappa = min(1.0, max(0.0, kappa))
                tensor = (normalized * scale) * scale
                values = (values * scale) * scale
                if not np.any(tensor) or np.any(
                    (np.diag(normalized) > 0) & (np.diag(tensor) == 0)
                ):
                    raise AnalysisError("Gyration tensor underflow")
                rg = float(np.sqrt(np.trace(tensor)))
            re = ratio = None
            if options.endpoints is None:
                diagnostics.append(
                    "endpoints_unavailable: no explicit endpoint pair supplied"
                )
            else:
                a, b = options.endpoints
                delta = xyz[ids.index(b)] - xyz[ids.index(a)]
                extent = float(np.max(np.abs(delta)))
                re = float(extent * np.linalg.norm(delta / extent)) if extent else 0.0
                if a == b:
                    diagnostics.append("same_site_endpoints: Re=0")
                elif re == 0:
                    diagnostics.append("coincident_endpoints: Re=0")
                if rg > 0:
                    ratio = float(np.float64(re) / rg)
            numeric = [*center, *tensor.flat, *values, rg]
            numeric += [x for x in (kappa, re, ratio) if x is not None]
            if not np.isfinite(numeric).all():
                raise AnalysisError("Nonfinite geometry result")
        return GeometryResult(
            tuple(selected),
            tuple(map(float, selected_masses)),
            options.weighting,
            tuple(map(float, weights)),
            tuple(map(float, center)),
            tuple(tuple(map(float, row)) for row in tensor),
            tuple(map(float, values)),
            rg,
            kappa,
            options.endpoints,
            "explicit_stable_ids" if options.endpoints else "unavailable",
            re,
            ratio,
            tuple(diagnostics),
        )
    except AnalysisError:
        raise
    except (
        TypeError,
        ValueError,
        OverflowError,
        FloatingPointError,
        np.linalg.LinAlgError,
    ) as error:
        raise AnalysisError(f"Geometry cannot be evaluated safely: {error}") from error
