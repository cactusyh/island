"""Bounded, independently reconstructed completion evidence; never native rows.

No Class-II coupling formula is authorized here. A fit is an external candidate,
not an execution profile. Neither adapter edits any historical parameter source.
"""

import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from hashlib import sha256
from math import isclose, isfinite, sqrt
from pathlib import Path

from .inventory import canonical, digest
from .registry import inspect_profile

MODES = (
    "strict_source",
    "derived_from_declared_formula",
    "fit_from_external_reference",
)
VERSION = "island_completion_evidence_v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pinned(root, entry):
    root = Path(root).resolve()
    name = Path(entry["file"])
    require(not name.is_absolute() and ".." not in name.parts, "Unsafe evidence path")
    path = (root / name).resolve()
    require(path.is_relative_to(root), "Escaping evidence symlink")
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == entry["sha256"], "Evidence source hash changed")
    return raw


def cited_lj_row(root, row):
    raw = pinned(root, row).decode().splitlines()
    require(
        type(row["line"]) is int and 1 <= row["line"] <= len(raw), "Invalid source line"
    )
    require(raw[row["line"] - 1] == row["raw"], "Raw source citation changed")
    require(row["units"] == ["angstrom", "kcal/mol"], "Unresolved LJ units")
    require(bool(row["label"]) and "*" not in row["label"], "Ambiguous atom label")
    values = [Decimal(row["raw"].split()[i]) for i in row["columns"]]
    require(
        len(values) == 2 and all(v.is_finite() for v in values), "Invalid LJ values"
    )
    require(values[0] > 0 and values[1] >= 0, "Invalid sigma/epsilon")
    require(
        [str(v) for v in values] == row["values"],
        "Extracted values differ from raw row",
    )
    require(
        bool(row["citation"]) and bool(row["source_identity"]),
        "Missing citation/source",
    )
    return values


def _geometric(root, profile, inputs):
    declared = inspect_profile(profile)
    require(
        declared["mixing_rule"] == "geometric_sigma_and_epsilon",
        "No declared geometric formula",
    )
    require(len(inputs) == 2, "Two explicit self-pair rows required")
    require(
        all(r["source_identity"] == profile for r in inputs), "Cross-source borrowing"
    )
    a, b = [cited_lj_row(root, r) for r in inputs]
    with localcontext() as context:
        context.prec = 50
        output = [str((x * y).sqrt()) for x, y in zip(a, b, strict=True)]
    independent = [sqrt(float(x) * float(y)) for x, y in zip(a, b, strict=True)]
    require(
        all(
            isclose(float(x), y, rel_tol=1e-13, abs_tol=1e-15)
            for x, y in zip(output, independent, strict=True)
        ),
        "Independent reconstruction failed",
    )
    return {
        "schema": VERSION,
        "kind": "derived",
        "source_identity": profile,
        "adapter": "geometric_lj_v1",
        "equation": "sigma_ij=sqrt(sigma_i*sigma_j); epsilon_ij=sqrt(epsilon_i*epsilon_j)",
        "formula_source": {
            "profile_identity": profile,
            "declaration_sha256": digest(declared),
            "mixing_rule": declared["mixing_rule"],
        },
        "inputs": inputs,
        "units": ["angstrom", "kcal/mol"],
        "values": output,
        "conversion": "identity; no unit conversion",
        "code_version": VERSION,
        "code_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "independent_reconstruction": "binary64 sqrt of separately parsed inputs",
        "validation_status": "reconstructed_data_only",
        "execution_authority": False,
    }


def _fit(root, profile, source):
    import numpy as np

    require(
        inspect_profile(profile)["family"] == "opls",
        "Only the declared harmonic OPLS form is supported",
    )
    data = json.loads(pinned(root, source))
    require(data["schema"] == "island_isolated_bond_scan_v1", "Unsupported fit dataset")
    require(data["source_identity"] == profile, "Cross-family fit")
    require(
        data["reference_kind"] in ("qm", "experiment", "reference_engine"),
        "Missing external reference",
    )
    require(
        bool(data["reference_provenance"]) and bool(data["citation"]),
        "Missing reference provenance",
    )
    require(
        data["functional_form"] == "K*(r-r0)^2+offset",
        "Incompatible fit functional form",
    )
    require(
        data["units"] == {"r": "angstrom", "energy": "kcal/mol"}, "Unresolved fit units"
    )
    require(
        data["isolated_bond_contribution"] is True,
        "Total energies cannot be silently treated as isolated bond data",
    )
    require(
        len(data["labels"]) == 2
        and all(type(x) is str and x and "*" not in x for x in data["labels"]),
        "Ambiguous fit labels",
    )
    x, y, w = [np.array(data[k], dtype=float) for k in ("r", "energy", "weights")]
    require(
        x.ndim == y.ndim == w.ndim == 1 and len(x) == len(y) == len(w) >= 4,
        "Incomplete fit arrays",
    )
    require(
        all(np.isfinite(v).all() for v in (x, y, w))
        and (x > 0).all()
        and (w > 0).all(),
        "Invalid fit samples/weights",
    )
    design = np.column_stack((x * x, x, np.ones_like(x)))
    a = design * np.sqrt(w)[:, None]
    b = y * np.sqrt(w)
    coeff, _, rank, _ = np.linalg.lstsq(a, b, rcond=None)
    require(rank == 3, "Underdetermined fit")
    k, linear, constant = coeff
    r0 = -linear / (2 * k)
    offset = constant - k * r0 * r0
    require(
        k > 0 and r0 > 0 and all(isfinite(float(v)) for v in (k, r0, offset)),
        "Fit violates positive K/r0 constraints",
    )
    # Separate normal-equation reconstruction checks the QR/SVD result.
    second = np.linalg.solve(a.T @ a, a.T @ b)
    require(
        np.allclose(coeff, second, rtol=1e-9, atol=1e-9),
        "Independent fit reconstruction failed",
    )
    residual = design @ coeff - y
    objective = float(np.sum(w * residual**2))
    covariance = np.linalg.inv(a.T @ a) * objective / (len(x) - 3)
    return {
        "schema": VERSION,
        "kind": "fitted_external",
        "source_identity": profile,
        "adapter": "weighted_harmonic_bond_v1",
        "reference": source,
        "labels": data["labels"],
        "reference_kind": data["reference_kind"],
        "reference_provenance": data["reference_provenance"],
        "citation": data["citation"],
        "equation": data["functional_form"],
        "units": data["units"],
        "objective_function": "sum_i weights_i*(K*(r_i-r0)^2+offset-E_i)^2",
        "weights": data["weights"],
        "constraints": ["K>0", "r0>0", "isolated bond contribution"],
        "values": {
            "K": repr(float(k)),
            "r0": repr(float(r0)),
            "offset": repr(float(offset)),
        },
        "objective_value": objective,
        "uncertainty": {
            "kind": "residual-based quadratic-coefficient covariance; not physical-model uncertainty",
            "covariance": covariance.tolist(),
            "degrees_of_freedom": len(x) - 3,
        },
        "code_version": VERSION,
        "code_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "independent_reconstruction": "normal equations compared with weighted lstsq",
        "validation_status": "reconstructed_data_only",
        "execution_authority": False,
    }


@dataclass(frozen=True)
class CompletionRecord:
    json_text: str

    @property
    def payload(self):
        p = json.loads(self.json_text)
        require(
            p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
            "Completion checksum changed",
        )
        return p

    def validate_integrity(self, root):
        p = self.payload
        if p["adapter"] == "geometric_lj_v1":
            expected = _geometric(root, p["source_identity"], p["inputs"])
        elif p["adapter"] == "weighted_harmonic_bond_v1":
            expected = _fit(root, p["source_identity"], p["reference"])
        else:
            raise ValueError("No defensible analytical/fitting adapter for this term")
        require(
            {k: v for k, v in p.items() if k != "identity"} == expected,
            "Contradictory completion evidence",
        )
        return p


def derive_geometric_pair(root, profile, left, right):
    p = _geometric(root, profile, [left, right])
    return CompletionRecord(canonical({**p, "identity": digest(p)}))


def fit_harmonic_bond(root, profile, *, reference_file, expected_sha256):
    p = _fit(root, profile, {"file": reference_file, "sha256": expected_sha256})
    return CompletionRecord(canonical({**p, "identity": digest(p)}))
