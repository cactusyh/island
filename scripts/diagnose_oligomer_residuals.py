"""Read-only, checksummed charge-stage diagnosis; never executes AmberTools.

AC has no independent element field: its generated name prefix is checked.
SQM indices are resolved through the named sqm.in inventory, never guessed.
Decimal retains printed values; prmtop uses the Amber 18.2223 charge scale.
"""

import argparse
import itertools
import re
from decimal import Decimal
from pathlib import Path

from island.workflows import storage
from island.workflows.bundle import system_from


def charge(token):
    value = Decimal(token)
    if not value.is_finite():
        raise ValueError("Nonfinite charge")
    return value


def inventory(rows, expected):
    result = {}
    for name, token in rows:
        if name in result or name not in expected:
            raise ValueError("Duplicate or unknown atom name")
        result[name] = charge(token)
    if set(result) != set(expected):
        raise ValueError("Incomplete atom inventory")
    return result


def ac(text, expected):
    rows = [r.split() for r in text.splitlines() if r.startswith("ATOM ")]
    for r in rows:
        if r[2] != r[2][:4] + "MOL":
            raise ValueError("Unexpected AC generated name/residue convention")
        name, atom_type = r[2][:4], r[-1]
        # BCC AC atom types may be numeric; these are not atomic numbers.
        if name not in expected or (
            not atom_type.isdigit() and atom_type[0].upper() != expected[name]
        ):
            raise ValueError("AC atom type/element mismatch")
    return inventory(((r[2][:4], r[-2]) for r in rows), expected)


def mol2(text, expected):
    atoms = text.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
    rows = [r.split() for r in atoms.splitlines() if r.strip()]
    for r in rows:
        if r[1] not in expected or r[5][0].upper() != expected[r[1]]:
            raise ValueError("MOL2 element/type mismatch")
    return inventory(((r[1], r[8]) for r in rows), expected)


def sqm(text, inputs, expected):
    named = [r.split() for r in inputs.splitlines() if re.match(r"^\s*\d+\s+\w+\s", r)]
    elements = {1: "H", 6: "C", 8: "O"}
    names = [r[1] for r in named]
    if len(set(names)) != len(names) or set(names) != set(expected):
        raise ValueError("SQM input lineage mismatch")
    for r in named:
        if elements[int(r[0])] != expected[r[1]]:
            raise ValueError("SQM input element mismatch")
    table = text.rsplit("Element       Mulliken Charge", 1)[1].split(
        "Total Mulliken Charge"
    )[0]
    rows = []
    seen = set()
    for line in table.splitlines():
        r = line.split()
        if not r:
            continue
        i = int(r[0])
        if i in seen or not 1 <= i <= len(names) or r[1] != expected[names[i - 1]]:
            raise ValueError("SQM output index/element mismatch")
        seen.add(i)
        rows.append((names[i - 1], r[2]))
    total = text.rsplit("Total Mulliken Charge =", 1)[1].split()[0]
    return inventory(rows, expected), total


def prmtop(text, expected):
    def field(name, width):
        body = (
            re.split(r"%FLAG " + name + r"\s*\n", text)[1]
            .split("\n", 1)[1]
            .split("%FLAG")[0]
        )
        return [
            line[i : i + width].strip()
            for line in body.splitlines()
            for i in range(0, len(line), width)
            if line[i : i + width].strip()
        ]

    names = field("ATOM_NAME", 4)
    numbers = field("ATOMIC_NUMBER", 8)
    values = field("CHARGE", 16)
    if len(names) != len(numbers) or len(names) != len(values):
        raise ValueError("prmtop array mismatch")
    for n, z in zip(names, numbers, strict=True):
        if {"H": 1, "C": 6, "O": 8}[expected[n]] != int(z):
            raise ValueError("prmtop element mismatch")
    return {
        n: q / Decimal("18.2223")
        for n, q in inventory(zip(names, values, strict=True), expected).items()
    }


def diagnose(root, manifest):
    # Validate every original checksum before parsing any scientific artifact.
    files = [
        (root / r["case"] / p, h)
        for r in manifest["cases"]
        for p, h in r["source_checksums"].items()
    ]
    for path, h in files:
        if storage.checksum(path.read_bytes()) != h:
            raise ValueError(f"Original checksum mismatch: {path}")
    results = []
    for row in manifest["cases"]:
        case = root / row["case"]
        a = next(case.glob("island-ambertools-*"))
        system = system_from(
            storage.decode(storage.read_json(case / "input-system.json"))
        )
        lineage = storage.read_json(a / "lineage.json")["input_name_to_site_id"]
        expected = {n: system.topology.sites[s].element for n, s in lineage.items()}
        if (
            set(lineage.values()) != set(system.topology.sites)
            or len(lineage) != system.number_of_sites
        ):
            raise ValueError("Source lineage coverage mismatch")
        if any(not n.startswith(e) for n, e in expected.items()):
            raise ValueError("Generated element prefix mismatch")
        stages = {}
        stages["input.mol2"] = mol2((a / "input.mol2").read_text(), expected)
        for name in [
            "ANTECHAMBER_BOND_TYPE.AC0",
            "ANTECHAMBER_BOND_TYPE.AC",
            "ANTECHAMBER_AC.AC0",
            "ANTECHAMBER_AC.AC",
        ]:
            stages[name] = ac((a / name).read_text(), expected)
        stages["sqm.out"], printed_total = sqm(
            (a / "sqm.out").read_text(), (a / "sqm.in").read_text(), expected
        )
        for name in ["ANTECHAMBER_AM1BCC_PRE.AC", "ANTECHAMBER_AM1BCC.AC"]:
            stages[name] = ac((a / name).read_text(), expected)
        stages["typed.mol2"] = mol2((a / "typed.mol2").read_text(), expected)
        stages["result.prmtop"] = prmtop((a / "result.prmtop").read_text(), expected)
        path = [
            "sqm.out",
            "ANTECHAMBER_AM1BCC_PRE.AC",
            "ANTECHAMBER_AM1BCC.AC",
            "typed.mol2",
            "result.prmtop",
        ]
        changes = []
        for before, after in itertools.pairwise(path):
            delta = {
                str(lineage[n]): stages[after][n] - stages[before][n] for n in expected
            }
            changes.append(
                {
                    "before": before,
                    "after": after,
                    "total_change": str(sum(delta.values())),
                    "max_absolute_change": str(max(map(abs, delta.values()))),
                    "per_site_change": {s: str(q) for s, q in delta.items()},
                }
            )
        results.append(
            {
                "case": row["case"],
                "historical_status": row["status"],
                "sites": len(expected),
                "source_checksums": row["source_checksums"],
                "sqm_printed_total": printed_total,
                "stages": {
                    k: {
                        "sum_e": str(sum(v.values())),
                        "printed_decimal_places": sorted(
                            {-q.as_tuple().exponent for q in v.values()}
                        )
                        if k != "result.prmtop"
                        else None,
                        "per_site_e": {str(lineage[n]): str(q) for n, q in v.items()},
                    }
                    for k, v in stages.items()
                },
                "adjacent_changes": changes,
            }
        )
    if any(storage.checksum(p.read_bytes()) != h for p, h in files):
        raise ValueError("Source changed during inspection")
    return {
        "schema": "island.charge-stage-diagnosis.v1",
        "cases": results,
        "verified_files": len(files),
        "production_validated": False,
        "simulation_readiness": "not_established",
        "precision": {
            "sqm": "3 fractional decimals",
            "AC/mol2": "6 fractional decimals",
            "prmtop": "E16.8 / 18.2223",
        },
        "limitation": "No unrounded per-atom SQM charge record is retained; AC elements rely on verified generated names.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument(
        "--manifest", type=Path, default=Path("docs/phase_4f3_acceptance.json")
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = diagnose(args.source, storage.read_json(args.manifest))
    storage.publish(args.output, storage.json_bytes(report))


if __name__ == "__main__":
    main()
