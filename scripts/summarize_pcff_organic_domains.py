"""Summarize J4 evidence without promoting partial source coverage to acceptance."""

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from island.workflows.storage import checksum, json_bytes, publish

MANDATORY = (
    "hydrogen_fluoride",
    "hydrogen_chloride",
    "hydrogen_bromide",
    "hydrogen_iodide",
)


def numerical_gate(report, names):
    if Counter(names) != Counter(MANDATORY):
        return False
    cases = report.get("cases", [])
    if (
        Counter(c.get("case") for c in cases) != Counter(names)
        or report.get("source_unchanged") is not True
    ):
        return False
    for case in cases:
        if (
            case.get("passed") is not True
            or "error" in case
            or not case.get("bundle_hashes")
            or not case.get("prepared_identity")
            or not case.get("model_identity")
            or not case.get("parameter_inventory_verified")
            or case.get("offline_loading") is not True
            or case.get("child_pid") == case.get("parent_pid")
            or not case.get("child_pid")
        ):
            return False
        delta = case.get("reconstruction_max_error")
        if (
            type(delta) not in (int, float)
            or not math.isfinite(delta)
            or not 0 <= delta <= 1e-10
        ):
            return False
        frames = case.get("frames", [])
        if len(frames) != 2:
            return False
        for f in frames:
            for key in ("energy_error", "force_max_error"):
                value = f.get(key)
                if (
                    type(value) not in (int, float)
                    or not math.isfinite(value)
                    or value < 0
                    or value > 1e-5
                ):
                    return False
            fd = f.get("finite_difference_errors", [])
            if (
                len(fd) != 3
                or not all(
                    type(x) in (int, float) and math.isfinite(x) and x >= 0 for x in fd
                )
                or fd[-1] > 1e-5
            ):
                return False
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    declaration = json.loads((a.root / "declaration.json").read_text())
    matrix = json.loads((a.root / "matrix/report.json").read_text())
    numeric = json.loads((a.root / "numerics-corrected/report.json").read_text())
    ledger = json.loads((a.root / "matrix/coverage.json").read_text())
    assert matrix["source"] == declaration["source"] == ledger["source"]
    assert Counter(c["name"] for c in matrix["cases"]) == Counter(
        c["name"] for c in declaration["cases"]
    )
    names = declaration["numerics"]["cases"]
    models = [c["name"] for c in matrix["cases"] if c["model_complete"]]
    gate = numerical_gate(numeric, names) and Counter(models) == Counter(names)
    for r in numeric["cases"]:
        m = next(c for c in matrix["cases"] if c["name"] == r["case"])
        # Fresh preparation and matrix assignment are deterministic in this path.
        gate = gate and r.get("model_identity") == m.get("model_identity")
        if r.get("passed"):
            child = json.loads(
                (a.root / "numerics-corrected" / r["case"] / "child.json").read_text()
            )
            gate = (
                gate
                and child["prepared"] == r["prepared_identity"]
                and child["native"] == r["model_identity"]
            )
    for row in ledger["type_rows"]:
        row["fixture_observations"] = [
            {
                "case": c["name"],
                "typed": c["typed"],
                "charges_complete": c["charges_complete"],
                "model_complete": c["model_complete"],
                "independently_verified": gate and c["name"] in names,
            }
            for c in matrix["cases"]
            if row["type"] in c.get("actual_counts", {})
        ]
    gaps = {
        "cases": [c for c in matrix["cases"] if not c["model_complete"]],
        "unresolved_labels": [
            r["type"] for r in ledger["type_rows"] if r["declared_predicate"] is None
        ],
        "full_source_complete": False,
    }
    summary = {
        "source": declaration["source"],
        "profile": declaration["profile"],
        "counts": matrix["counts"],
        "fixture_denominator": matrix["denominator"],
        "rule_labels": ledger["rule_labels"],
        "source_label_denominator": 133,
        "numerical_bundle_gate": bool(gate),
        "new_complete_models": len(models),
        "independently_verified_new_models": len(names) if gate else 0,
        "numerical": numeric,
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    for name, data in [
        ("summary.json", summary),
        ("coverage.json", ledger),
        ("gaps.json", gaps),
    ]:
        publish(a.output / name, json_bytes(data))
    publish(
        a.output / "input_hashes.json",
        json_bytes(
            {
                str(p.relative_to(a.root)): checksum(p.read_bytes())
                for p in a.root.rglob("*")
                if p.is_file()
                and a.output not in p.parents
                and p.relative_to(a.root).parts[0]
                in (
                    "matrix",
                    "numerics",
                    "numerics-corrected",
                    "controls",
                    "declaration.json",
                    "source-audit.json",
                    "historical-typing.json",
                    "historical-offline.json",
                )
            }
        ),
    )
    print("Numerical/bundle gate:", bool(gate), "; full source: False")
    return 1 if a.require_full_source or not gate else 0


if __name__ == "__main__":
    raise SystemExit(main())
