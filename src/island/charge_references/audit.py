"""Continuous charge differences and explicitly hypothetical naive-transfer drift."""

from dataclasses import dataclass
from itertools import combinations
from math import fsum, sqrt

from island.exceptions import ChargeReferenceError

from .correspondence import require
from .records import ChargeReference, pack, unpack

AUDIT_SCHEMA = "island_oligomer_charge_audit_v1"


def summarize(reference):
    p = reference.payload
    charges = {s: a["charge"] for s, a in p["charge_result"]["assignments"].items()}
    return summarize_charge_mapping(
        p["correspondence"],
        charges,
        reference.identity,
        p["seeds"],
        p["charge_result"]["total_charge_residual"],
    )


def summarize_charge_mapping(correspondence, charges, identity, seeds, residual):
    """Shared comparison kernel; does not certify the input as a ChargeReference."""
    repeats = []
    for repeat in correspondence["repeats"]:
        groups = []
        for g in repeat["groups"]:
            values = [charges[s] for s in g["hydrogen_ids"]]
            count = len(values)
            total = fsum(values)
            groups.append(
                {
                    **g,
                    "heavy_charge": charges[g["site_id"]],
                    "hydrogen_raw": dict(zip(g["hydrogen_ids"], values)),
                    "hydrogen_count": count,
                    "hydrogen_sum": total,
                    "hydrogen_mean": total / count if count else None,
                    "hydrogen_spread_max_minus_min": max(values) - min(values)
                    if count
                    else None,
                }
            )
        repeats.append(
            {
                **repeat,
                "groups": groups,
                "net_charge": fsum(
                    g["heavy_charge"] + g["hydrogen_sum"] for g in groups
                ),
            }
        )
    return {
        "reference": identity,
        "definition": correspondence["definition"],
        "dp": correspondence["dp"],
        "seeds": seeds,
        "repeats": repeats,
        "molecular_charge_residual": residual,
    }


def _selected(summary, role):
    return next(
        r
        for r in summary["repeats"]
        if (r["central"] if role == "central" else r["role"] == role)
    )


def _difference(left, right, kind):
    rows = []
    for role in ("head", "central", "tail"):
        a, b = _selected(left, role), _selected(right, role)
        for ga, gb in zip(a["groups"], b["groups"], strict=True):
            for key in (
                "source_repeat_atom_index",
                "element",
                "heavy_environment",
                "hydrogen_count",
            ):
                require(
                    ga[key] == gb[key], f"Incompatible comparison environment: {key}"
                )
            rows.append(
                {
                    "role": role,
                    "source_repeat_atom_index": ga["source_repeat_atom_index"],
                    "heavy_difference": gb["heavy_charge"] - ga["heavy_charge"],
                    "hydrogen_count": ga["hydrogen_count"],
                    "hydrogen_sum_difference": gb["hydrogen_sum"] - ga["hydrogen_sum"],
                    "hydrogen_mean_difference": None
                    if ga["hydrogen_mean"] is None
                    else gb["hydrogen_mean"] - ga["hydrogen_mean"],
                    "hydrogen_spread_difference": None
                    if ga["hydrogen_mean"] is None
                    else gb["hydrogen_spread_max_minus_min"]
                    - ga["hydrogen_spread_max_minus_min"],
                }
            )
    metrics = {}
    for key in (
        "heavy_difference",
        "hydrogen_sum_difference",
        "hydrogen_mean_difference",
    ):
        values = [r[key] for r in rows if r[key] is not None]
        metrics[key] = {
            "maximum_absolute": max(map(abs, values)) if values else None,
            "rms": sqrt(fsum(v * v for v in values) / len(values)) if values else None,
            "group_count": len(values),
        }
    return {
        "kind": kind,
        "left": left["reference"],
        "right": right["reference"],
        "sign": "right_minus_left",
        "matched_groups": rows,
        "metrics": metrics,
    }


def audit_data(references):
    return audit_summaries([summarize(r) for r in references])


def audit_summaries(summaries):
    """Compare already validated, explicitly identified charge summaries."""
    require(
        len({r["reference"] for r in summaries}) == len(summaries),
        "Duplicate references",
    )
    comparisons = []
    diagnostics = []
    for definition in ("[*:1]CC[*:2]", "[*:1]CCO[*:2]"):
        matching = [s for s in summaries if s["definition"] == definition]
        for dp in (3, 5, 7):
            if not any(s["dp"] == dp for s in matching):
                diagnostics.append(
                    {
                        "definition": definition,
                        "dp": dp,
                        "reason": "No validated reference for this length",
                    }
                )
        seeds = {tuple(sorted(s["seeds"].items())) for s in matching if s["dp"] == 5}
        if len(seeds) < 2:
            diagnostics.append(
                {
                    "definition": definition,
                    "reason": "DP5 conformation comparison unavailable: fewer than two validated seeds",
                }
            )

    for a, b in combinations(summaries, 2):
        if a["definition"] != b["definition"]:
            continue
        if a["seeds"] == b["seeds"] and a["dp"] != b["dp"]:
            kind = "length_variation_at_declared_seeds"
        elif a["dp"] == b["dp"] == 5 and a["seeds"] != b["seeds"]:
            kind = "conformation_seed_variation"
        else:
            continue
        try:
            comparisons.append(_difference(a, b, kind))
        except ChargeReferenceError as error:
            diagnostics.append(
                {"left": a["reference"], "right": b["reference"], "reason": str(error)}
            )
    extrapolations = []
    for s in summaries:
        head, interior, tail = (
            _selected(s, r)["net_charge"] for r in ("head", "central", "tail")
        )
        extrapolations.append(
            {
                "reference": s["reference"],
                "hypothesis": "head + (N-2)*central + tail; NOT target atom charges",
                "head": head,
                "interior": interior,
                "tail": tail,
                "predicted_total_e": {
                    n: head + (n - 2) * interior + tail for n in (20, 50, 100)
                },
            }
        )
    return {
        "summaries": summaries,
        "comparisons": comparisons,
        "diagnostic_extrapolations": extrapolations,
        "diagnostics": diagnostics,
        "units": "elementary_charge",
        "convention": "Raw heavy atoms and parent hydrogen groups; no averaging assignment or neutrality correction",
        "limitation": "Length and conformation effects are not perfectly isolated; no scientific pass threshold",
    }


@dataclass(frozen=True)
class ChargeAudit:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p)
                == {
                    "schema",
                    "references",
                    "audit",
                    "production_validated",
                    "simulation_readiness",
                },
                "Unexpected audit fields",
            )
            require(
                p["schema"] == AUDIT_SCHEMA
                and p["production_validated"] is False
                and p["simulation_readiness"] == "not_established",
                "Invalid audit schema/readiness",
            )
            refs = [ChargeReference(pack(r)) for r in p["references"]]
            require(bool(refs), "Empty audit")
            require(
                audit_data(refs) == p["audit"], "Audit does not match its references"
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid charge audit: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)


def audit_charge_references(references):
    try:
        refs = list(references)
        result = ChargeAudit(
            pack(
                {
                    "schema": AUDIT_SCHEMA,
                    "references": [r.payload for r in refs],
                    "audit": audit_data(refs),
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot audit references: {error}") from error
