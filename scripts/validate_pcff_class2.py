"""Pinned-source Class II assignment audit. No energy or simulation acceptance.

Default gate: inventories and independent row/conversion checks pass.
--require-complete additionally requires every applicable family to be assigned.
Missing rows remain durable diagnostics under either gate.
"""

import argparse
import json
import math
import time
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    load_pcff_source,
    save_pcff_parameters,
    type_pcff_atoms,
)
from island.forcefields.pcff.source import PIN
from island.workflows.storage import json_bytes, publish

# Frozen before execution: (sites, bonds, angles, propers, angle-angle couplings).
CASES = [
    {"name": "butane", "smiles": "CCCC", "counts": [14, 13, 24, 27, 48]},
    {"name": "neopentane", "smiles": "CC(C)(C)C", "counts": [17, 16, 30, 36, 60]},
    {"name": "ethanol", "smiles": "CCO", "counts": [9, 8, 13, 12, 24]},
    {"name": "dimethyl_ether", "smiles": "COC", "counts": [9, 8, 13, 6, 24]},
    {
        "name": "PE_DP3",
        "psmiles": "[*:1]CC[*:2]",
        "dp": 3,
        "counts": [20, 19, 36, 45, 72],
    },
    {
        "name": "PEO_DP3",
        "psmiles": "[*:1]CCO[*:2]",
        "dp": 3,
        "counts": [23, 22, 39, 42, 72],
    },
    {
        "name": "PE_DP50",
        "psmiles": "[*:1]CC[*:2]",
        "dp": 50,
        "counts": [302, 301, 600, 891, 1200],
    },
]
# Independent oracle specification: source arity, ordinary-equivalence column,
# indices with degree units, and indices with unconverted angstrom units.
ORACLE = {
    "quartic_bond": (2, 2, [], [0]),
    "quartic_angle": (3, 3, [0], []),
    "torsion_3": (4, 4, [1, 3, 5], []),
    "wilson_out_of_plane": (4, 5, [1], []),
    "nonbond(9-6)": (1, 1, [], [0]),
    "bond-bond": (3, 3, [], []),
    "bond-angle": (3, 3, [], []),
    "bond-bond_1_3": (4, 4, [], []),
    "end_bond-torsion_3": (4, 4, [], []),
    "middle_bond-torsion_3": (4, 4, [], []),
    "angle-torsion_3": (4, 4, [], []),
    "angle-angle-torsion_1": (4, 4, [], []),
    "angle-angle": (4, 5, [], []),
}
TOLERANCES = {"absolute": 1e-10, "relative": 1e-13}


def independent_rows(raw):
    """Separate line scanner; no production parser or selection helpers."""
    table, equiv = {}, {}
    section = namespace = ""
    for lineno, line in enumerate(raw.decode().splitlines(), 1):
        words = line.split("!", 1)[0].split()
        if not words:
            continue
        if words[0].startswith("#"):
            section = words[0][1:]
            namespace = words[1] if len(words) > 1 else ""
            continue
        if namespace != "cff91" or not words[0][0].isdigit():
            continue
        if section == "equivalence":
            key = words[2]
            if key not in equiv or Decimal(words[0]) > Decimal(equiv[key][0]):
                equiv[key] = words
        if section in ORACLE:
            arity = ORACLE[section][0]
            table.setdefault(section, []).append(
                (lineno, Decimal(words[0]), words[2 : 2 + arity], words[2 + arity :])
            )
    return table, equiv


def oracle_check(payload, raw):
    table, equiv = independent_rows(raw)
    maximum = 0.0
    checked = 0
    for assignment in payload["assignments"]:
        family = assignment["family"]
        if assignment["status"] == "not_applicable":
            assert family == "wilson_out_of_plane"
            continue
        _, column, degrees, lengths = ORACLE[family]
        supplied = assignment["supplied_types"]
        queries = [supplied, [equiv[t][2 + column] for t in supplied]]
        candidates = []
        for labels in queries:
            for line, version, types, strings in table.get(family, []):
                mirrored = (
                    [labels[3], labels[1], labels[2], labels[0]]
                    if family == "angle-angle"
                    else labels[::-1]
                )
                for reversed_roles, match in [(False, labels), (True, mirrored)]:
                    if types != match:
                        continue
                    numbers = list(map(Decimal, strings))
                    if family == "bond-angle" and len(numbers) == 1:
                        numbers *= 2
                    if (
                        family in ("end_bond-torsion_3", "angle-torsion_3")
                        and len(numbers) == 3
                    ):
                        numbers *= 2
                    converted = [
                        float(v) * math.pi / 180
                        if i in degrees
                        else float(v)
                        if i in lengths
                        else float(v * Decimal("4.184"))
                        for i, v in enumerate(numbers)
                    ]
                    if reversed_roles and family == "bond-angle":
                        converted.reverse()
                    if reversed_roles and family in (
                        "end_bond-torsion_3",
                        "angle-torsion_3",
                    ):
                        converted = converted[3:] + converted[:3]
                    candidates.append((line, version, converted))
            if candidates:
                break
        if not candidates:
            assert assignment["status"] == "missing"
            continue
        newest = max(c[1] for c in candidates)
        chosen = [c for c in candidates if c[1] == newest]
        if any(c[2] != chosen[0][2] for c in chosen):
            assert assignment["status"] == "ambiguous"
            continue
        assert assignment["status"] == "assigned", assignment["id"]
        selected_lines = {
            int(c["record_id"].rsplit(":", 1)[1]) for c in assignment["selected"]
        }
        assert selected_lines == {c[0] for c in chosen}
        for actual, expected in zip(
            assignment["normalized_values"], chosen[0][2], strict=True
        ):
            maximum = max(maximum, abs(actual - expected))
            assert math.isclose(
                actual,
                expected,
                abs_tol=TOLERANCES["absolute"],
                rel_tol=TOLERANCES["relative"],
            )
        checked += 1
    return {
        "independently_checked_assignments": checked,
        "maximum_coefficient_error": maximum,
        "oracle": "separate source-line scanner, Decimal kcal conversion, independent orientation/selection",
    }


def build(case):
    if "smiles" in case:
        from island.chemistry import from_smiles

        return from_smiles(case["smiles"], add_hydrogens=True, random_seed=2026)
    from island.builders import build_linear_polymer

    return build_linear_polymer(case["psmiles"], dp=case["dp"], generate_3d=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args(argv)
    source = load_pcff_source(args.source)
    out = Path(args.output)
    out.mkdir(exist_ok=False, parents=True)
    declaration = {
        "schema": "island_pcff_class2_acceptance_v1",
        "source": source.identity,
        "cases": CASES,
        "coordinate_seed": 2026,
        "tolerances": TOLERANCES,
        "selection": "exact then ordinary family equivalence, highest version; conflicts diagnostic",
        "gates": [
            "inventory counts",
            "independent source rows/conversions",
            "integrity",
        ],
        "require_complete": args.require_complete,
        "scientific_validation": False,
        "reference_revision": PIN["commit"],
        "no_external_zero_fill": True,
    }
    publish(out / "declaration.json", json_bytes(declaration))
    outcomes = []
    for case in CASES:
        start = time.perf_counter()
        try:
            system = build(case)
            typing = type_pcff_atoms(system, source)
            charges = assign_automatic_pcff_charges(system, typing)
            assignment_start = time.perf_counter()
            result = assign_pcff_parameters(system, typing, charges)
            assignment_seconds = time.perf_counter() - assignment_start
            p = result.payload
            counts = [
                len(system.topology.sites),
                len(p["inventories"]["bonds"]),
                len(p["inventories"]["angles"]),
                len(p["inventories"]["proper_torsions"]),
                sum(p["coverage"]["angle-angle"].values()),
            ]
            assert counts == case["counts"], (counts, case["counts"])
            comparisons = oracle_check(p, source.raw)
            path = out / (case["name"] + ".json")
            save_pcff_parameters(result, path)
            outcome = {
                "name": case["name"],
                "software_gate": "passed",
                "counts": counts,
                "coverage": p["coverage"],
                "parameter_coverage_complete": p["parameter_coverage_complete"],
                "physical_model_complete": p["physical_model_complete"],
                "record_identity": result.identity,
                "file_sha256": sha256(path.read_bytes()).hexdigest(),
                "bytes": path.stat().st_size,
                "assignment_with_validation_seconds": assignment_seconds,
                **comparisons,
            }
        except Exception as error:  # noqa: BLE001 -- durable per-case failure reporting
            outcome = {
                "name": case["name"],
                "software_gate": "failed",
                "error": f"{type(error).__name__}: {error}",
            }
        outcome["total_seconds"] = time.perf_counter() - start
        outcomes.append(outcome)
        print(json.dumps(outcome), flush=True)
        publish(
            out / "report.json",
            json_bytes(
                {
                    "declaration_sha256": sha256(
                        (out / "declaration.json").read_bytes()
                    ).hexdigest(),
                    "outcomes": outcomes,
                }
            ),
            replace=True,
        )
    return int(
        any(
            o["software_gate"] != "passed"
            or (args.require_complete and not o.get("parameter_coverage_complete"))
            for o in outcomes
        )
    )


if __name__ == "__main__":
    raise SystemExit(main())
