"""J6 independent, pinned-source authority controls; never a runtime resolver."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from pcff_j5_reference import read_source
from validate_pcff_fallbacks import compare, compiled, lammps
from validate_pcff_operational_support import write


def external_matcher(path, expected):
    """Explicit trusted validation tool path, not a persisted scientific loader."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError("External matcher checksum mismatch")
    spec = importlib.util.spec_from_file_location("pinned_lunar_matcher", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matcher_controls(raw, matcher, queries):
    if not queries or len({name for name, _ in queries}) != len(queries):
        raise ValueError("Expected nonempty distinct declared angle controls")
    rows = read_source(raw)
    table = {tuple(r["types"]): r for r in rows["quadratic_angle", "cff91_auto"]}
    auto = {}
    for r in rows["auto_equivalence", "cff91_auto"]:
        if r["type"] not in auto or r["version"] > auto[r["type"]]["version"]:
            auto[r["type"]] = r
    equiv = {k: SimpleNamespace(**r["map"]) for k, r in auto.items()}
    messages = []
    logger = SimpleNamespace(
        warn=messages.append, error=lambda m: (_ for _ in ()).throw(ValueError(m))
    )
    results = []
    for case, supplied in queries:
        checks = []
        for reverse_table in (False, True):
            mapping = dict(reversed(list(table.items()))) if reverse_table else table
            for types in (supplied, list(reversed(supplied))):
                found, key, order, transformed = matcher.match_3_body(
                    *types, logger, mapping, equiv, True, "auto-equiv"
                )
                row = mapping[key] if found else None
                checks.append(
                    {
                        "supplied_types": types,
                        "reverse_table": reverse_table,
                        "found": found,
                        "selected_pattern": list(key),
                        "ordered_types": list(order),
                        "transformed": list(transformed),
                        "source_row": row["id"] if row else None,
                        "raw_coefficients": [str(x) for x in row["values"]]
                        if row
                        else None,
                    }
                )
        results.append(
            {
                "case": case,
                "checks": checks,
                "reversal_invariant": checks[0]["raw_coefficients"]
                == checks[1]["raw_coefficients"],
                "row_iteration_invariant": checks[0]["raw_coefficients"]
                == checks[2]["raw_coefficients"],
            }
        )
    # Synthetic suffix-renaming control: same source coefficients, explicitly not a source file.
    renamed = {
        tuple("*99" if x.startswith("*") else x for x in key): value
        for key, value in table.items()
    }
    suffix = matcher.match_3_body(
        "s'_", "c=_", "h_", logger, renamed, False, True, "auto-equiv"
    )
    exact = matcher.match_3_body(
        "h_", "c'_", "h_", logger, table, False, True, "auto-equiv"
    )
    wrong_center = matcher.match_3_body(
        "h_", "undeclared_center", "h_", logger, table, False, True, "auto-equiv"
    )
    original = matcher.match_3_body(
        "s'_", "c=_", "h_", logger, table, False, True, "auto-equiv"
    )
    if (
        not exact[0]
        or wrong_center[0]
        or renamed[suffix[1]]["values"] != table[original[1]]["values"]
    ):
        raise AssertionError("Declared external-matcher controls failed")
    return {
        "negative_controls": {
            "exact_row": table[exact[1]]["id"],
            "wrong_center_rejected": not wrong_center[0],
            "suffix_renaming_preserves_coefficients": True,
        },
        "cases": results,
        "messages": messages,
        "synthetic_suffix_control": {
            "selected_pattern": list(suffix[1]),
            "description": "Only wildcard names replaced by *99; not a source modification or physical policy",
        },
        "authoritative_orientation_independent_precedence_established": all(
            c["reversal_invariant"] for c in results
        ),
    }


def run(args):
    declaration = json.loads(args.declaration.read_text())
    raw = args.source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != declaration["source"]["sha256"]:
        raise ValueError("FRC checksum mismatch")
    expected = declaration["references"]["external_files"][str(args.matcher)]
    matcher = external_matcher(args.matcher, expected)
    rows = read_source(raw)
    queries = []
    for case in declaration["cases"]:
        from island.workflows.storage import decode

        report = decode(
            json.loads(
                (args.reproduction / case["name"] / "inspection.json").read_text()
            )
        )
        seen = set()
        for q in report["interaction_queries"]:
            if (
                q["family"] == "quartic_angle"
                and q["resolution"]["status"] == "ambiguous"
            ):
                k = min(
                    tuple(q["supplied_types"]), tuple(reversed(q["supplied_types"]))
                )
                if k not in seen:
                    queries.append((case["name"], list(k)))
                    seen.add(k)
    if {name for name, _ in queries} != {
        "thioformaldehyde",
        "thioacetone",
        "alpha_lactam",
        "beta_lactam",
    }:
        raise ValueError("Missing or unexpected required ambiguity controls")
    controls = matcher_controls(raw, matcher, queries)
    # No unmatched/ambiguous angle is promoted to an operational source assignment.
    write(args.output / "matcher.json", controls)
    if args.matcher_only:
        return {
            "matcher": controls,
            "scope": "matching controls only",
            "new_model_gate": False,
            "full_source_gate": False,
        }
    term_rows = {r["line"]: r for r in rows["quadratic_angle", "cff91_auto"]}
    xyz = np.array(declaration["wildcard_controls"]["geometry_angstrom"])
    write(
        args.output / "numerical-declaration.json",
        {
            "rows": declaration["wildcard_controls"]["source_rows"],
            "coordinates_angstrom": xyz.tolist(),
            "tolerances": declaration["tolerances"],
            "scope": "Isolated alternative coefficients, not molecule models or source-policy acceptance",
        },
    )
    numerics = []
    for line in declaration["wildcard_controls"]["source_rows"]:
        row = term_rows[line]
        theta, k = map(float, row["values"])
        coeff = [theta * np.pi / 180, k * 4.184]
        energy, force = compiled("quadratic_angle", coeff, xyz)
        ref, rf, meta = lammps(
            args.lammps.resolve(),
            args.output / f"angle-{line}",
            xyz,
            "angle",
            "harmonic",
            f"{k:.17g} {theta:.17g}",
        )

        def energy_function(x, coeff=coeff):
            a, b = x[0] - x[1], x[2] - x[1]
            angle = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
            return coeff[1] * (angle - coeff[0]) ** 2

        fd = []
        for h in declaration["tolerances"]["finite_difference_angstrom"]:
            f = np.zeros_like(xyz)
            for i in range(3):
                for j in range(3):
                    plus, minus = xyz.copy(), xyz.copy()
                    plus[i, j] += h
                    minus[i, j] -= h
                    f[i, j] = -(energy_function(plus) - energy_function(minus)) / (
                        2 * h
                    )
            fd.append(
                {
                    "step_angstrom": h,
                    "maximum_force_error": float(np.max(np.abs(f - force))),
                }
            )
        if (
            fd[-1]["maximum_force_error"]
            > declaration["tolerances"]["force_atol_kj_mol_angstrom"]
        ):
            raise AssertionError("Declared final finite-difference criterion failed")
        numerics.append(
            {
                "source_row": row["id"],
                "pattern": row["types"],
                "source_values": [str(v) for v in row["values"]],
                "energy": energy,
                "force": force.tolist(),
                "comparison": compare(energy, force, ref, rf),
                "finite_difference": fd,
                **meta,
            }
        )
    return {
        "matcher": controls,
        "isolated_angles": numerics,
        "new_model_gate": False,
        "full_source_gate": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "source",
        "matcher",
        "declaration",
        "reproduction",
        "lammps",
        "output",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--matcher-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        result = run(args)
    except Exception as error:
        write(
            args.output / "failure.json",
            {"error": f"{type(error).__name__}: {error}", "new_model_gate": False},
        )
        raise
    write(args.output / "report.json", result)


if __name__ == "__main__":
    main()
