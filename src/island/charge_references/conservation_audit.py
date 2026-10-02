"""Conditional raw/projected comparisons; no target atom charge assignment."""

from dataclasses import dataclass
from math import fsum, sqrt
from pathlib import Path

import numpy as np

from island.exceptions import ChargeReferenceError

from .audit import _selected, audit_summaries, summarize_charge_mapping
from .correspondence import require
from .projection import ChargeProjection
from .records import pack, unpack

SCHEMA = "island_conditional_conservation_audit_v1"


def _audit(projections):
    summaries = {"raw": [], "projected": []}
    cases = []
    data = []
    for projection in projections:
        payload = projection.payload
        d = payload["observation"]["data"]
        q = payload["projection"]
        data.append(d)
        cases.append(
            {
                "case": d["case"],
                "observation_identity": payload["observation_identity"],
                "projection_identity": projection.identity,
                "historical_status": d["historical_status"],
                "historical_failure": d["historical_failure"],
                "historical_reference_identity": d["historical_reference_identity"],
                "projection": q,
                "input_coordinate_fingerprint": d["input_coordinate_fingerprint"],
                "post_sqm_coordinate_fingerprint": d["post_sqm_coordinate_fingerprint"],
                "typed_coordinate_fingerprint": d["typed_coordinate_fingerprint"],
            }
        )
        for mode, values in summaries.items():
            values.append(
                summarize_charge_mapping(
                    d["correspondence"],
                    q[mode + "_charges"],
                    d["case"],
                    d["seeds"],
                    q[mode + "_residual_e"],
                )
            )
    require(len({d["case"] for d in data}) == len(data), "Duplicate audit cases")
    audits = {mode: audit_summaries(values) for mode, values in summaries.items()}
    for mode, a in audits.items():
        a["convention"] = (
            mode
            + " post-BCC observation; heavy atoms and parent H groups; diagnostic extrapolation only"
        )
    for mode, audit in audits.items():
        by_case = {s["reference"]: s for s in summaries[mode]}
        for comparison in audit["comparisons"]:
            left, right = (by_case[comparison[k]] for k in ("left", "right"))
            comparison["repeat_net_differences_e"] = {
                role: _selected(right, role)["net_charge"]
                - _selected(left, role)["net_charge"]
                for role in ("head", "central", "tail")
            }
    contributions = []
    for raw, projected in zip(
        audits["raw"]["comparisons"], audits["projected"]["comparisons"], strict=True
    ):
        require(
            (raw["left"], raw["right"], raw["kind"])
            == (projected["left"], projected["right"], projected["kind"]),
            "Comparison mismatch",
        )
        rows = []
        for r, p in zip(
            raw["matched_groups"], projected["matched_groups"], strict=True
        ):
            rows.append(
                {
                    "role": r["role"],
                    "source_repeat_atom_index": r["source_repeat_atom_index"],
                    **{
                        k: None if r[k] is None else p[k] - r[k]
                        for k in (
                            "heavy_difference",
                            "hydrogen_sum_difference",
                            "hydrogen_mean_difference",
                            "hydrogen_spread_difference",
                        )
                    },
                }
            )
        contributions.append(
            {
                "left": raw["left"],
                "right": raw["right"],
                "kind": raw["kind"],
                "convention": "projected difference minus raw difference",
                "matched_groups": rows,
                "repeat_net_differences_e": {
                    role: projected["repeat_net_differences_e"][role]
                    - raw["repeat_net_differences_e"][role]
                    for role in ("head", "central", "tail")
                },
            }
        )
    geometry = []
    for definition in ("[*:1]CC[*:2]", "[*:1]CCO[*:2]"):
        selected = [
            d
            for d in data
            if d["correspondence"]["definition"] == definition
            and d["correspondence"]["dp"] == 5
        ]
        if len(selected) != 2:
            continue
        left, right = selected

        def distances(d, key):
            ids = [
                g["site_id"]
                for r in d["correspondence"]["repeats"]
                for g in r["groups"]
            ]
            if key == "system":
                from island.workflows.bundle import system_from

                system = system_from(d["system"])
                xyz = np.array([system.coordinates.get(s) for s in ids])
            else:
                xyz = np.array([d[key][s] for s in ids])
            return np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)[
                np.triu_indices(len(ids), 1)
            ]

        metrics = {}
        for label, key in [
            ("input", "system"),
            ("post_sqm", "post_sqm_coordinates_angstrom"),
            ("typed", "typed_coordinates_angstrom"),
        ]:
            delta = distances(right, key) - distances(left, key)
            metrics[label + "_heavy_pair_distance_rms_difference_angstrom"] = sqrt(
                fsum(float(x * x) for x in delta) / len(delta)
            )
        geometry.append(
            {
                "left": left["case"],
                "right": right["case"],
                **metrics,
                "convention": "All heavy-atom pairs mapped by oriented repeat/source index; excludes ambiguous H permutations; finite printed geometry, not a basin certificate",
            }
        )
    return {
        "cases": cases,
        "raw": audits["raw"],
        "projected": audits["projected"],
        "projection_contributions": contributions,
        "geometry_comparisons": geometry,
        "historical_valid_reference_count": sum(
            d["historical_status"] == "passed" for d in data
        ),
        "raw_observation_count": len(data),
        "projected_case_count": len(data),
        "original_raw_conformation_acceptance": "unmet; original failures are not repaired",
        "interpretation": "Conditional projected comparison, not recovered unrounded QM or scientific validation",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


@dataclass(frozen=True)
class ConservationAudit:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p) == {"schema", "projections", "audit"} and p["schema"] == SCHEMA,
                "Invalid conservation audit schema",
            )
            require(bool(p["projections"]), "Empty audit")
            projections = [ChargeProjection(pack(v)) for v in p["projections"]]
            require(
                pack(p["audit"]) == pack(_audit(projections)),
                "Audit metrics or identities changed",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(
                f"Invalid conservation audit: {error}"
            ) from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)


def audit_charge_projections(projections):
    try:
        projections = list(projections)
        require(
            all(isinstance(p, ChargeProjection) for p in projections),
            "Expected charge projections",
        )
        result = ConservationAudit(
            pack(
                {
                    "schema": SCHEMA,
                    "projections": [p.payload for p in projections],
                    "audit": _audit(projections),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot audit projections: {error}") from error


def load_conservation_audit(path):
    try:
        result = ConservationAudit(Path(path).read_text())
        result.validate_integrity()
        return result
    except Exception as error:
        raise ChargeReferenceError(
            f"Cannot load conservation audit: {error}"
        ) from error
