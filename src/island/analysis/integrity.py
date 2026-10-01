"""Internal v1 report consistency, not computational authenticity.

Geometry identities use 256*float64 epsilon relative to tensor scale (no
absolute angstrom floor). Energy/temperature arithmetic retains dynamics'
1e-12 relative / 1e-10 absolute policy. No source files are opened.
"""

import json
import math
from dataclasses import fields
from pathlib import PurePosixPath

import numpy as np

from island.exceptions import AnalysisError

from .geometry import GeometryOptions, GeometryResult

EPS = 256 * np.finfo(float).eps


def require(condition, message):
    if not condition:
        raise AnalysisError(message)


def keys(value, expected):
    require(
        type(value) is dict and set(value) == set(expected),
        "Missing or unexpected report fields",
    )


def number(x, *, nonnegative=False, positive=False):
    require(
        type(x) in (int, float) and math.isfinite(x), "Expected finite numerical scalar"
    )
    require(not nonnegative or x >= 0, "Expected nonnegative value")
    require(not positive or x > 0, "Expected positive value")
    return x


def integer(x):
    require(type(x) is int and x >= 0, "Expected nonnegative integer")


def text(x):
    require(type(x) is str and bool(x), "Expected nonempty string")


def digest(x):
    require(
        type(x) is str and len(x) == 64 and all(c in "0123456789abcdef" for c in x),
        "Invalid content identity",
    )


def path(x):
    text(x)
    require(
        not PurePosixPath(x).is_absolute()
        and all(p not in ("", ".", "..") for p in x.split("/")),
        "Invalid relative artifact path",
    )


def arithmetic(a, b):
    require(
        math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-10),
        "Inconsistent energy/temperature arithmetic",
    )


def _pairs(rows):
    result = {}
    for key, value in rows:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _finite_tree(value):
    if type(value) is float:
        require(math.isfinite(value), "Nonfinite JSON number")
    elif type(value) is dict:
        for item in value.values():
            _finite_tree(item)
    elif type(value) is list:
        for item in value:
            _finite_tree(item)


def parse_report(raw):
    try:
        require(type(raw) is str, "AnalysisReport requires JSON text")
        p = json.loads(raw, object_pairs_hook=_pairs)
        _finite_tree(p)  # Also catches 1e999, not just NaN/Infinity constants.
        return p
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Invalid analysis JSON: {error}") from error


def _vector(value, n):
    require(type(value) is list and len(value) == n, "Invalid geometry dimensions")
    for x in value:
        number(x)
    return np.asarray(value, dtype=float)


def _geometry(g, options, inventory):
    keys(g, [f.name for f in fields(GeometryResult)])
    ids = g["selected_ids"]
    require(
        type(ids) is list and bool(ids) and all(type(s) is int for s in ids),
        "Invalid selected IDs",
    )
    require(
        ids == sorted(set(ids)) and set(ids) <= inventory,
        "Duplicate/unknown selected IDs",
    )
    if options.selection == "all":
        require(set(ids) == inventory, "All-site selection lacks full inventory")
    elif options.selection == "explicit":
        require(ids == list(options.site_ids), "Explicit selection mismatch")
    # v1 contains no element inventory: heavy membership cannot be re-derived.
    masses = _vector(g["selected_masses"], len(ids))
    weights = _vector(g["normalized_weights"], len(ids))
    require(
        np.all(masses > 0) and np.all(weights > 0), "Masses/weights must be positive"
    )
    require(g["weighting"] == options.weighting, "Weighting mismatch")
    expected = (
        masses / masses.max() if options.weighting == "mass" else np.ones(len(ids))
    )
    expected /= expected.sum()
    require(
        np.allclose(weights, expected, rtol=EPS, atol=0)
        and abs(weights.sum() - 1) <= EPS,
        "Incorrect normalized weights",
    )
    require(
        (g["coordinate_unit"], g["tensor_unit"], g["mass_unit"], g["shape_unit"])
        == ("angstrom", "angstrom^2", "dalton", "dimensionless"),
        "Unsupported geometry units",
    )
    _vector(g["center"], 3)
    require(
        type(g["gyration_tensor"]) is list and len(g["gyration_tensor"]) == 3,
        "Invalid tensor dimensions",
    )
    tensor = np.array([_vector(row, 3) for row in g["gyration_tensor"]])
    eig = _vector(g["eigenvalues"], 3)
    require(
        np.all(eig >= 0) and np.all(np.diff(eig) >= 0), "Invalid eigenvalue order/sign"
    )
    rg = number(g["radius_of_gyration"], nonnegative=True)
    diagnostics = g["diagnostics"]
    require(
        type(diagnostics) is list and all(type(s) is str and s for s in diagnostics),
        "Invalid diagnostics",
    )
    scale = float(np.max(np.abs(tensor)))
    require(len(ids) != 1 or scale == 0, "Single-site selection must have zero spread")
    if scale == 0:
        require(
            rg == 0
            and np.all(eig == 0)
            and g["kappa_squared"] is None
            and g["re_over_rg"] is None,
            "Invalid zero-spread metrics",
        )
        require(
            any(s.startswith("zero_spread:") for s in diagnostics),
            "Missing zero-spread diagnostic",
        )
    else:
        normalized = tensor / scale
        require(
            np.max(np.abs(normalized - normalized.T)) <= EPS,
            "Asymmetric gyration tensor",
        )
        computed = np.linalg.eigvalsh(normalized)
        require(computed[0] >= -64 * np.finfo(float).eps, "Indefinite gyration tensor")
        require(
            np.max(np.abs(eig / scale - np.maximum(computed, 0))) <= EPS,
            "Tensor/eigenvalue mismatch",
        )
        require(rg > 0, "Nonzero tensor requires positive Rg")
        require(
            math.isclose(
                rg / math.sqrt(scale),
                math.sqrt(float(np.trace(normalized))),
                rel_tol=EPS,
                abs_tol=0,
            ),
            "Tensor/Rg mismatch",
        )
        fractions = (eig / scale) / np.sum(eig / scale)
        kappa = number(g["kappa_squared"], nonnegative=True)
        require(
            kappa <= 1
            and abs(kappa - (1.5 * float(fractions @ fractions) - 0.5)) <= EPS,
            "Invalid shape anisotropy",
        )
    endpoints = g["endpoint_ids"]
    if endpoints is None:
        require(
            options.endpoints is None
            and g["endpoint_convention"] == "unavailable"
            and g["end_to_end_distance"] is None
            and g["re_over_rg"] is None,
            "Unavailable endpoint mismatch",
        )
        require(
            any(s.startswith("endpoints_unavailable:") for s in diagnostics),
            "Missing endpoint diagnostic",
        )
    else:
        require(
            type(endpoints) is list
            and len(endpoints) == 2
            and all(type(s) is int and s in inventory for s in endpoints),
            "Invalid endpoint IDs",
        )
        require(
            g["endpoint_convention"]
            == (
                "explicit_stable_ids"
                if options.endpoints is not None
                else "builder_head_tail"
            ),
            "Endpoint convention mismatch",
        )
        if options.endpoints is not None:
            require(endpoints == list(options.endpoints), "Explicit endpoints differ")
        re = number(g["end_to_end_distance"], nonnegative=True)
        if endpoints[0] == endpoints[1]:
            require(
                re == 0 and "same_site_endpoints: Re=0" in diagnostics,
                "Invalid equal endpoints",
            )
        elif re == 0:
            require(
                "coincident_endpoints: Re=0" in diagnostics,
                "Missing coincident endpoint diagnostic",
            )
        if rg > 0:
            ratio = number(g["re_over_rg"], nonnegative=True)
            require(
                math.isclose(ratio, re / rg, rel_tol=EPS, abs_tol=0), "Re/Rg mismatch"
            )
    return (
        tuple(ids),
        tuple(masses),
        tuple(weights),
        endpoints,
        g["endpoint_convention"],
    )


def validate_report(raw):
    """Return an owned, internally checked v1 payload, without source access."""
    from island.dynamics._checkpoint_data import frame_from
    from island.dynamics.integrity import validate_frame_content
    from island.dynamics.thermal import GAS_CONSTANT

    from .workflow import AnalysisOptions

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            p = parse_report(raw)
            keys(
                p,
                (
                    "schema",
                    "implementation_version",
                    "source",
                    "options",
                    "frames",
                    "sample_count",
                    "sampling",
                    "boundary_policy",
                    "conventions",
                    "production_validated",
                    "simulation_readiness",
                    "interpretation",
                ),
            )
            require(
                p["schema"] == "island.trajectory-analysis.v1"
                and p["implementation_version"] == "1",
                "Unsupported analysis schema/version",
            )
            require(
                p["production_validated"] is False
                and p["simulation_readiness"] == "not_established",
                "Analysis cannot promote scientific readiness",
            )
            require(
                p["sampling"]
                == "inclusive windows then stride over retained samples; no interpolation"
                and p["boundary_policy"]
                == "retain earlier accepted frame; original records unchanged"
                and p["interpretation"]
                == "retained samples only; no equilibrium or independence claim",
                "Unsupported analysis interpretation",
            )
            keys(
                p["options"],
                (
                    "geometry",
                    "start_step",
                    "end_step",
                    "start_time_ps",
                    "end_time_ps",
                    "stride",
                ),
            )
            keys(
                p["options"]["geometry"],
                ("selection", "site_ids", "weighting", "endpoints"),
            )
            options = AnalysisOptions(
                **{
                    **p["options"],
                    "geometry": GeometryOptions(**p["options"]["geometry"]),
                }
            )
            require(
                p["conventions"]
                == {
                    "weighting": options.geometry.weighting,
                    "center": "sum(w*r)/sum(w)",
                    "tensor": "sum(w*(r-center)*(r-center)^T)/sum(w)",
                    "Rg": "sqrt(trace(G))",
                    "eigenvalues": "ascending angstrom^2",
                    "kappa_squared": "1.5*sum(lambda^2)/sum(lambda)^2-0.5",
                    "temperature_dof": "3N; translation and rotation included",
                    "missing_values": "JSON null; CSV empty field; see diagnostics",
                },
                "Unsupported conventions",
            )
            source = p["source"]
            keys(
                source,
                (
                    "directory",
                    "manifest_checksum",
                    "manifest_file_sha256",
                    "status",
                    "accepted_step",
                    "run_id",
                    "evidence",
                    "production_validated",
                    "simulation_readiness",
                    "model_fingerprint",
                    "parameter_fingerprint",
                    "backend",
                    "artifacts",
                    "checkpoint",
                    "segments",
                    "origin",
                ),
            )
            text(source["directory"])
            text(source["run_id"])
            integer(source["accepted_step"])
            require(
                source["status"]
                in (
                    "paused",
                    "completed",
                    "budget_exhausted",
                    "stage_failed",
                    "starting",
                    "ready",
                ),
                "Invalid source status",
            )
            require(
                source["evidence"]
                in (
                    "live_parameterization",
                    "archived_parameters",
                    "synthetic_software_test",
                ),
                "Invalid source evidence",
            )
            require(
                source["production_validated"] is False
                and source["simulation_readiness"] == "not_established",
                "Invalid source readiness",
            )
            for key in (
                "manifest_checksum",
                "manifest_file_sha256",
                "model_fingerprint",
                "parameter_fingerprint",
            ):
                digest(source[key])
            require(
                type(source["artifacts"]) is dict and bool(source["artifacts"]),
                "Missing source artifacts",
            )
            for name, descriptor in source["artifacts"].items():
                path(name)
                keys(descriptor, ("bytes", "sha256"))
                integer(descriptor["bytes"])
                digest(descriptor["sha256"])
            path(source["checkpoint"])
            require(source["checkpoint"] in source["artifacts"], "Unlisted checkpoint")
            require(
                type(source["segments"]) is list and bool(source["segments"]),
                "Missing segments",
            )
            segments = []
            for row in source["segments"]:
                keys(row, ("file", "checksum"))
                path(row["file"])
                digest(row["checksum"])
                require(
                    row["file"] in source["artifacts"] and row["file"] not in segments,
                    "Unlisted/duplicate segment",
                )
                segments.append(row["file"])
            origin = source["origin"]
            keys(
                origin,
                ("initial_state", "initial_total_energy", "trajectory_fingerprint"),
            )
            digest(origin["trajectory_fingerprint"])
            initial = frame_from(origin["initial_state"], True)
            validate_frame_content(initial)
            require(
                initial.step == 0
                and initial.time_ps == 0
                and initial.evaluation is not None,
                "Invalid trajectory origin",
            )
            arithmetic(number(origin["initial_total_energy"]), initial.total_energy)
            ev = initial.evaluation
            require(
                source["model_fingerprint"] == ev.model_fingerprint
                and source["parameter_fingerprint"] == ev.parameter_fingerprint,
                "Source potential identity mismatch",
            )
            require(
                source["backend"]
                == {
                    "name": ev.backend_name,
                    "version": ev.backend_version,
                    "platform": ev.platform,
                    "settings": dict(ev.settings),
                },
                "Source backend mismatch",
            )
            frames = p["frames"]
            require(
                type(frames) is list and bool(frames), "Report requires nonempty frames"
            )
            integer(p["sample_count"])
            require(p["sample_count"] == len(frames), "Sample count mismatch")
            last_step, last_time, last_segment = -1, -1.0, -1
            selection = None
            for row in frames:
                keys(
                    row,
                    (
                        "step",
                        "time_ps",
                        "coordinate_fingerprint",
                        "source_segment",
                        "evaluation_fingerprint",
                        "potential_energy_kj_mol",
                        "kinetic_energy_kj_mol",
                        "total_energy_kj_mol",
                        "temperature_kelvin",
                        "geometry",
                    ),
                )
                integer(row["step"])
                number(row["time_ps"], nonnegative=True)
                require(
                    last_step < row["step"] <= source["accepted_step"]
                    and last_time < row["time_ps"],
                    "Unordered frames or invalid progress",
                )
                require(
                    (row["step"] == 0) == (row["time_ps"] == 0),
                    "Step/time zero mismatch",
                )
                if last_step >= 0:
                    require(
                        row["step"] - last_step >= options.stride,
                        "Impossible retained-sample stride",
                    )
                for low, high, x in (
                    (options.start_step, options.end_step, row["step"]),
                    (options.start_time_ps, options.end_time_ps, row["time_ps"]),
                ):
                    require(
                        (low is None or x >= low) and (high is None or x <= high),
                        "Frame outside selected window",
                    )
                require(row["source_segment"] in segments, "Unknown source segment")
                index = segments.index(row["source_segment"])
                require(index >= last_segment, "Reversed source segments")
                digest(row["coordinate_fingerprint"])
                digest(row["evaluation_fingerprint"])
                u = number(row["potential_energy_kj_mol"])
                k = number(row["kinetic_energy_kj_mol"], nonnegative=True)
                arithmetic(number(row["total_energy_kj_mol"]), u + k)
                arithmetic(
                    number(row["temperature_kelvin"], nonnegative=True),
                    2 * k / (3 * len(initial.coordinates) * GAS_CONSTANT),
                )
                current = _geometry(
                    row["geometry"], options.geometry, set(initial.coordinates)
                )
                require(
                    selection is None or current == selection,
                    "Selection/masses/endpoints changed across frames",
                )
                selection = current
                if row["step"] == 0:
                    require(
                        row["coordinate_fingerprint"] == initial.coordinate_fingerprint
                        and row["evaluation_fingerprint"] == ev.evaluation_fingerprint,
                        "Origin frame identity mismatch",
                    )
                    arithmetic(u, initial.potential_energy)
                    arithmetic(k, initial.kinetic_energy)
                last_step, last_time, last_segment = row["step"], row["time_ps"], index
            if set(selection[0]) == set(initial.coordinates):
                # All masses needed for this additional origin consistency check
                # are present; subset selections cannot supply missing masses.
                from island.dynamics.models import kinetic_energy

                arithmetic(
                    initial.kinetic_energy,
                    kinetic_energy(
                        dict(zip(selection[0], selection[1], strict=True)),
                        initial.velocities,
                    ),
                )
            return p
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Invalid analysis report: {error}") from error
