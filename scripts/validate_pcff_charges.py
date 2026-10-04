"""Pinned real-source PCFF acceptance with illustrative manual typing.

No embedding, automatic typing, QM, parameterization, or energy evaluation.
The independent oracle reads fixed audited source rows and sums Decimal values.
"""

import argparse
import json
import sys
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_pcff_charges,
    assign_pcff_types,
    load_pcff_source,
    save_pcff_record,
)
from island.forcefields.pcff.source import PIN
from island.workflows.storage import json_bytes, publish

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from pcff_charges import molecule

CASES = ["ethane", "ethanol", "dimethyl_ether"]
# Audited fixed source lines, not production parser/selection output.
EXPECTED = {
    493: ["1.0", "1", "c", "c", "0.0000", "0.0000"],
    505: ["1.0", "1", "c", "h", "-0.0530", "0.0530"],
    519: ["2.1", "8", "c", "o", "0.1330", "-0.1330"],
    810: ["1.0", "1", "h*", "o", "0.4241", "-0.4241"],
}
EQ_LINES = {
    "c3": (214, "c"),
    "c2": (213, "c"),
    "hc": (244, "h"),
    "ho": (249, "h*"),
    "oh": (298, "o"),
    "oc": (296, "o"),
}
TOLERANCE = 1e-12


def reference(raw, system, types):
    lines = raw.decode().splitlines()
    table = {}
    for line, fields in EXPECTED.items():
        if lines[line - 1].split() != fields:
            raise ValueError(f"Audited source row changed: {line}")
        table[tuple(fields[2:4])] = (line, [Decimal(t) for t in fields[4:]])
    for label, (line, equivalent) in EQ_LINES.items():
        tokens = lines[line - 1].split()
        if tokens[2] != label or tokens[4] != equivalent:
            raise ValueError(f"Audited bond-equivalence row changed: {line}")
    q = {sid: Decimal(0) for sid in types}
    selections = {}
    for bond in system.topology.bonds.values():
        a, b = sorted((bond.site1, bond.site2))
        left, right = EQ_LINES[types[a]][1], EQ_LINES[types[b]][1]
        if (left, right) in table:
            line, inc = table[left, right]
        else:
            line, flipped = table[right, left]
            inc = flipped[::-1]
        q[a] += inc[0]
        q[b] += inc[1]
        selections[a, b] = line
    return {k: float(v) for k, v in q.items()}, selections


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    root = Path(args.output)
    # Exclusive experiment, declaration persisted before assignment.
    try:
        root.mkdir()
    except OSError as error:
        print(str(error), file=sys.stderr)
        return 1
    declaration = {
        "schema": "island_pcff_acceptance_v1",
        "source_pin": PIN,
        "cases": CASES,
        "charge_tolerance_e": TOLERANCE,
        "typing_provenance": "Illustrative manual labels from source comments; not upstream validated typing",
        "expected_increment_rows": EXPECTED,
        "expected_equivalence_rows": EQ_LINES,
        "gates": [
            "complete",
            "independent_row_selection",
            "independent_per_site_charges",
            "component_and_molecular_charge",
            "reversal_and_nonmutation",
            "record_roundtrip",
        ],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(root / "declaration.json", json_bytes(declaration))
    outcomes = []
    for case in CASES:
        try:
            raw = Path(args.source).read_bytes()
            if sha256(raw).hexdigest() != PIN["sha256"]:
                raise ValueError(
                    "Source checksum mismatch before reading scientific data"
                )
            source = load_pcff_source(args.source)
            system, types = molecule(case)
            before = system.copy()
            typing = assign_pcff_types(
                system, source, types, provenance=declaration["typing_provenance"]
            )
            result = assign_pcff_charges(system, typing)
            q, selected = reference(raw, system, types)
            error = max(abs(result.charges[sid] - value) for sid, value in q.items())
            if error > TOLERANCE or abs(result.payload["total_charge"]) > TOLERANCE:
                raise ValueError("Numerical acceptance failed")
            for c in result.payload["contributions"]:
                if (
                    c["source_matches"][0]["selection"]["record"]["line"]
                    != selected[tuple(c["sites"])]
                ):
                    raise ValueError("Selected record mismatch")
            reversed_system, reversed_types = molecule(case, reverse=True)
            reverse = assign_pcff_charges(
                reversed_system,
                assign_pcff_types(
                    reversed_system,
                    source,
                    reversed_types,
                    provenance=declaration["typing_provenance"],
                ),
            )
            if reverse.identity != result.identity:
                raise ValueError("Orientation/insertion dependence")
            result.validate_integrity(system)
            if (
                system.topology.sites != before.topology.sites
                or system.topology.bonds != before.topology.bonds
                or system.metadata != before.metadata
                or len(system.coordinates) != 0
            ):
                raise ValueError("Input mutated")
            save_pcff_record(result, root / f"{case}.json")
            from island.forcefields.pcff import load_pcff_record

            if (
                load_pcff_record(root / f"{case}.json", source, system=system).identity
                != result.identity
            ):
                raise ValueError("Roundtrip changed record")
            outcomes.append(
                {
                    "case": case,
                    "status": "passed",
                    "sites": len(q),
                    "bonds": len(selected),
                    "maximum_charge_error_e": error,
                    "total_charge_e": result.payload["total_charge"],
                    "charges_e": result.charges,
                    "result_identity": result.identity,
                    "typing_identity": typing.identity,
                    "source_sha256": sha256(raw).hexdigest(),
                }
            )
        except (PCFFError, OSError, ValueError) as error:
            outcomes.append({"case": case, "status": "failed", "error": str(error)})
    complete = len(outcomes) == 3 and all(o["status"] == "passed" for o in outcomes)
    try:
        unchanged = sha256(Path(args.source).read_bytes()).hexdigest() == PIN["sha256"]
    except OSError:
        unchanged = False
    report = {
        "schema": "island_pcff_acceptance_report_v1",
        "passed": complete,
        "outcomes": outcomes,
        "source_unchanged": unchanged,
        "upstream_typed_fixture": "unavailable; manual examples only",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(root / "report.json", json_bytes(report))
    print(json.dumps(report, indent=2))
    return 0 if complete and report["source_unchanged"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
