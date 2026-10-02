"""Experimental PE templates; no automatic force-field or charge assignment hook."""

from dataclasses import dataclass
from math import fsum, sqrt
from pathlib import Path

import numpy as np

from island.exceptions import ChargeReferenceError
from island.workflows import storage
from island.workflows.bundle import system_data, system_from

from .correspondence import _correspondence, require
from .projection import POLICY, ChargeProjection
from .records import pack, unpack

TRAINING = ((3, 2026), (5, 2026), (7, 2026), (5, 80317))
HELD_OUT = ((5, 314159), (9, 314159), (9, 271828), (11, 314159))
CLASSES = (
    "terminal_C",
    "terminal_H",
    "near_end_C",
    "near_end_H",
    "interior_C",
    "interior_H",
)
CONSTRAINTS = [[1.0, 3.0, 1.0, 2.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 2.0, 4.0]]
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
MODEL_SCHEMA = "island_experimental_pe_template_v1"
PREDICTION_SCHEMA = "island_experimental_pe_prediction_v1"
VALIDATION_SCHEMA = "island_experimental_pe_validation_v1"


def _identity(p):
    return storage.checksum(storage.json_bytes(storage.encode(p)))


def pe_target_correspondence(system):
    """Strict graph/provenance mapping, even or odd PE DP>=3, <=1000 sites."""
    try:
        return _correspondence(system, reference=False)
    except Exception as error:
        raise ChargeReferenceError(
            f"Invalid experimental PE target: {error}"
        ) from error


def _classes(corr):
    assignments = {}
    for repeat in corr["repeats"]:
        for group in repeat["groups"]:
            j = group["source_repeat_atom_index"]
            role = repeat["role"]
            # Verified oriented CC path reversal: (r,j) -> (DP-1-r,3-j).
            base = (
                4
                if role == "interior"
                else 0
                if (role, j) in {("head", 1), ("tail", 2)}
                else 2
            )
            require(
                len(group["hydrogen_ids"]) == (3 if base == 0 else 2),
                "PE class multiplicity mismatch",
            )
            assignments[group["site_id"]] = base
            for h in group["hydrogen_ids"]:
                require(h not in assignments, "Duplicate class site")
                assignments[h] = base + 1
    return dict(sorted(assignments.items()))


def _chemical_identity(system, corr):
    # No coordinate or construction-seed content; IDs remain authoritative.
    return _identity(
        {
            "correspondence": corr,
            "masses": {
                s: system.topology.sites[s].mass for s in sorted(system.topology.sites)
            },
        }
    )


def _solve(design, target, constraints):
    """Full-rank homogeneous equality constrained LS via SVD null space."""
    a = np.asarray(design, dtype=float)
    b = np.asarray(target, dtype=float)
    c = np.asarray(constraints, dtype=float).reshape(-1, a.shape[1])
    require(
        a.ndim == 2
        and b.shape == (a.shape[0],)
        and a.shape[1] > 0
        and all(np.isfinite(x).all() for x in (a, b, c)),
        "Invalid least-squares arrays",
    )
    rank_c = 0
    z = np.eye(a.shape[1])
    if len(c):
        _, singular, vh = np.linalg.svd(c, full_matrices=True)
        rank_c = int(
            np.sum(singular > np.finfo(float).eps * max(c.shape) * singular[0])
        )
        require(rank_c == len(c), "Dependent constraint rows")
        z = vh[rank_c:].T
    reduced = a @ z
    solution, _, rank, singular = np.linalg.lstsq(reduced, b, rcond=None)
    require(rank == z.shape[1] and z.shape[1] > 0, "Underdetermined PE fit")
    coefficients = z @ solution
    return coefficients, {
        "algorithm": "numpy_svd_nullspace_lstsq_v1",
        "constraint_rank": rank_c,
        "reduced_rank": int(rank),
        "parameter_count": a.shape[1],
        "singular_values": singular.tolist(),
        "objective_e2": float(np.sum((a @ coefficients - b) ** 2)),
        "constraint_residuals_e": (c @ coefficients).tolist(),
    }


def _fit(training, mode):
    require(mode in {"baseline", "conserving"}, "Unknown PE model")
    require(len(training) == 4, "Require exactly four frozen training cases")
    rows, targets, refs, seen = [], [], [], []
    for payload in training:
        p = ChargeProjection(pack(payload)).payload
        d = p["observation"]["data"]
        source = system_from(d["system"])
        corr = pe_target_correspondence(source)
        seed = d["seeds"]["template_seed"]
        require(d["seeds"]["assembly_seed"] == seed, "Training seed mismatch")
        key = (corr["dp"], seed)
        require(
            key in TRAINING and key not in seen, "Training/held-out split violation"
        )
        require(p["projection"]["policy"] == POLICY, "Wrong training target policy")
        seen.append(key)
        classes = _classes(corr)
        q = p["projection"]["projected_charges"]
        w = sqrt(1.0 / (4 * len(classes)))
        for site, cls in classes.items():
            row = [0.0] * 6
            row[cls] = w
            rows.append(row)
            targets.append(w * q[site])
        refs.append(
            {
                "case": d["case"],
                "observation_identity": p["observation_identity"],
                "projection_identity": _identity(p),
                "system_sha256": d["system_sha256"],
                "historical_status": d["historical_status"],
                "historical_failure": d["historical_failure"],
            }
        )
    require(tuple(seen) == TRAINING, "Training cases must use declared order")
    constraints = CONSTRAINTS if mode == "conserving" else []
    coefficients, diagnostics = _solve(rows, targets, constraints)
    return {
        "mode": mode,
        "classes": list(CLASSES),
        "hydrogen_multiplicities": [3, 2, 2],
        "reversal": "(r,j)->(DP-1-r,3-j); head1=tail2, head2=tail1, interior1=interior2",
        "objective": "sum_cases(sum_atoms(error_e**2)/N_atoms)/4",
        "target": "post-BCC uniform_molecular_l2_v1 projected charges",
        "constraints": constraints,
        "coefficients_e": coefficients.tolist(),
        "solver": diagnostics,
        "training": refs,
        "held_out": [list(x) for x in HELD_OUT],
        "numerical_conservation_tolerance_e": 1e-12,
        **FLAGS,
    }


@dataclass(frozen=True)
class PETemplateModel:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p) == {"schema", "training_projections", "fit"}
                and p["schema"] == MODEL_SCHEMA,
                "Invalid PE model schema",
            )
            require(
                pack(p["fit"])
                == pack(_fit(p["training_projections"], p["fit"]["mode"])),
                "PE model differs from frozen training/definition",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid PE model: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return _identity(self.payload)


def fit_pe_template(projections, *, mode):
    try:
        training = [p.payload for p in projections]
        result = PETemplateModel(
            pack(
                {
                    "schema": MODEL_SCHEMA,
                    "training_projections": training,
                    "fit": _fit(training, mode),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot fit PE model: {error}") from error


def _prediction(model, source):
    system = system_from(source)
    corr = pe_target_correspondence(system)
    classes = _classes(corr)
    coefficients = model["fit"]["coefficients_e"]
    q = {s: coefficients[k] for s, k in classes.items()}
    total = fsum(q.values())
    if model["fit"]["mode"] == "conserving":
        require(abs(total) <= 1e-12, "Constrained model conservation failed")
    repeats = [
        {
            "repeat_index": r["repeat_index"],
            "role": r["role"],
            "total_e": fsum(
                q[s] for g in r["groups"] for s in [g["site_id"], *g["hydrogen_ids"]]
            ),
        }
        for r in corr["repeats"]
    ]
    return {
        "model_identity": _identity(model),
        "chemical_identity": _chemical_identity(system, corr),
        "assignments_e": q,
        "class_indices": classes,
        "correspondence": corr,
        "repeat_totals": repeats,
        "total_e": total,
        "unit": "elementary_charge",
        "interpretation": "experimental projected-post-BCC PE template; not AM1-BCC assignment",
        **FLAGS,
    }


@dataclass(frozen=True)
class PEChargePrediction:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p) == {"schema", "model", "source_system", "prediction"}
                and p["schema"] == PREDICTION_SCHEMA,
                "Invalid prediction schema",
            )
            model = PETemplateModel(pack(p["model"])).payload
            require(
                pack(p["prediction"]) == pack(_prediction(model, p["source_system"])),
                "Prediction contradicts model/source",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid PE prediction: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return _identity(self.payload)


def predict_pe_charges(model, system):
    try:
        m, source = model.payload, system_data(system)
        result = PEChargePrediction(
            pack(
                {
                    "schema": PREDICTION_SCHEMA,
                    "model": m,
                    "source_system": source,
                    "prediction": _prediction(m, source),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot predict PE charges: {error}") from error


def load_pe_record(path):
    try:
        raw = Path(path).read_text()
        kind = unpack(raw)["schema"]
        result = {
            MODEL_SCHEMA: PETemplateModel,
            PREDICTION_SCHEMA: PEChargePrediction,
            VALIDATION_SCHEMA: PETemplateValidation,
        }[kind](raw)
        result.validate_integrity()
        return result
    except Exception as error:
        raise ChargeReferenceError(f"Cannot load PE record: {error}") from error


def _errors(values):
    return {
        "maximum_absolute_e": max(abs(v) for v in values),
        "rmse_e": sqrt(fsum(v * v for v in values) / len(values)),
    }


def _comparison(prediction, q):
    assigned = prediction["assignments_e"]
    corr = prediction["correspondence"]
    groups, repeats = [], []
    for r in corr["repeats"]:
        sites = []
        for g in r["groups"]:
            site, hs = g["site_id"], g["hydrogen_ids"]
            sites.extend([site, *hs])
            herror = fsum(assigned[h] - q[h] for h in hs)
            groups.append(
                {
                    "repeat_index": r["repeat_index"],
                    "role": r["role"],
                    "distance_from_nearest_end_repeats": min(
                        r["repeat_index"], corr["dp"] - 1 - r["repeat_index"]
                    ),
                    "local_heavy_index": g["source_repeat_atom_index"],
                    "heavy_error_e": assigned[site] - q[site],
                    "hydrogen_count": len(hs),
                    "hydrogen_mean_error_e": herror / len(hs),
                    "hydrogen_sum_error_e": herror,
                }
            )
        repeats.append(
            {
                "repeat_index": r["repeat_index"],
                "role": r["role"],
                "target_total_e": fsum(q[s] for s in sites),
                "predicted_total_e": fsum(assigned[s] for s in sites),
            }
        )
    profiles = {}
    for category in ("role", "distance_from_nearest_end_repeats"):
        profiles[category] = {
            str(key): _errors(
                [g["heavy_error_e"] for g in groups if g[category] == key]
            )
            for key in sorted({g[category] for g in groups})
        }
    return {
        "heavy_atoms": _errors([g["heavy_error_e"] for g in groups]),
        "parent_H_means": _errors([g["hydrogen_mean_error_e"] for g in groups]),
        "parent_H_sums": _errors([g["hydrogen_sum_error_e"] for g in groups]),
        "groups": groups,
        "profiles": profiles,
        "repeat_totals": repeats,
        "target_total_e": fsum(q.values()),
        "predicted_total_e": prediction["total_e"],
    }


def _validation(models, targets, frozen):
    require(
        len(models) == 2
        and [m["fit"]["mode"] for m in models] == ["baseline", "conserving"],
        "Require both frozen models",
    )
    require([_identity(m) for m in models] == frozen, "Frozen model identities changed")
    require(
        models[0]["training_projections"] == models[1]["training_projections"],
        "Different training data",
    )
    require(0 < len(targets) <= 4, "Require declared held-out targets")
    seen, cases, observations = [], [], []
    for payload in targets:
        p = ChargeProjection(pack(payload)).payload
        d = p["observation"]["data"]
        corr = d["correspondence"]
        seed = d["seeds"]["template_seed"]
        key = (corr["dp"], seed)
        require(
            key in HELD_OUT and key not in seen and d["seeds"]["assembly_seed"] == seed,
            "Held-out split violation",
        )
        seen.append(key)
        observations.append(d)
        comparisons = {}
        for m in models:
            pred = _prediction(m, d["system"])
            comparisons[m["fit"]["mode"]] = {
                "raw_post_BCC": _comparison(pred, d["raw_charges"]),
                "projected_post_BCC": _comparison(
                    pred, p["projection"]["projected_charges"]
                ),
            }
        cases.append(
            {
                "case": d["case"],
                "projection_identity": _identity(p),
                "observation_identity": p["observation_identity"],
                "raw_import_status": d["historical_status"],
                "raw_failure": d["historical_failure"],
                "input_coordinate_fingerprint": d["input_coordinate_fingerprint"],
                "post_sqm_coordinate_fingerprint": d["post_sqm_coordinate_fingerprint"],
                "comparisons": comparisons,
            }
        )
    pair = [d for d in observations if d["correspondence"]["dp"] == 9]
    geometry = None
    if len(pair) == 2:

        def distances(d, field):
            ids = [
                g["site_id"]
                for r in d["correspondence"]["repeats"]
                for g in r["groups"]
            ]
            if field == "input":
                system = system_from(d["system"])
                x = np.array([system.coordinates.get(s) for s in ids])
            else:
                x = np.array([d["post_sqm_coordinates_angstrom"][s] for s in ids])
            return np.array(
                [np.linalg.norm(x[i] - x[j]) for i in range(len(x)) for j in range(i)]
            )

        geometry = {
            "cases": [d["case"] for d in pair],
            "convention": "verified oriented heavy atom pair-distance RMS; no individual-H match or basin claim",
            **{
                field + "_heavy_pair_distance_rms_angstrom": float(
                    np.sqrt(
                        np.mean(
                            (distances(pair[0], field) - distances(pair[1], field)) ** 2
                        )
                    )
                )
                for field in ("input", "post_sqm")
            },
        }
    return {
        "frozen_model_identities": frozen,
        "cases": cases,
        "dp9_geometry": geometry,
        "complete_held_out_targets": len(cases) == 4,
        "scientific_accuracy_threshold": None,
        "unit": "elementary_charge",
        **FLAGS,
    }


@dataclass(frozen=True)
class PETemplateValidation:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p)
                == {
                    "schema",
                    "models",
                    "targets",
                    "frozen_model_identities",
                    "validation",
                }
                and p["schema"] == VALIDATION_SCHEMA,
                "Invalid PE validation schema",
            )
            models = [PETemplateModel(pack(m)).payload for m in p["models"]]
            require(
                pack(p["validation"])
                == pack(
                    _validation(models, p["targets"], p["frozen_model_identities"])
                ),
                "Validation differs from frozen evidence",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid PE validation: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return _identity(self.payload)


def validate_pe_templates(models, targets, *, frozen_model_identities):
    try:
        ms, ts = [m.payload for m in models], [t.payload for t in targets]
        result = PETemplateValidation(
            pack(
                {
                    "schema": VALIDATION_SCHEMA,
                    "models": ms,
                    "targets": ts,
                    "frozen_model_identities": list(frozen_model_identities),
                    "validation": _validation(ms, ts, list(frozen_model_identities)),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot validate PE templates: {error}") from error
