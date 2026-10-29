"""Independent ASCII, charge and graph audit; imports no production resolver."""

import argparse
import json
from decimal import Decimal
from hashlib import sha256
from itertools import permutations
from math import pi
from pathlib import Path

from island.charge_references.records import unpack
from island.workflows.bundle import system_from
from scripts.pcff_j5_reference import (
    ARITY,
    ROLE,
    charge_raw,
    equivalents,
    oriented_values,
    raw_inventory,
    read_source,
)


def strict_resolution(rows, family, types):
    """Independent literal/ordinary lookup, keeping center/shared-angle roles."""
    eq, evidence = equivalents(rows, types, "equivalence", [ROLE[family]] * len(types))
    for query in [types] + ([eq] if eq else []):
        if family == "wilson_out_of_plane":
            perms = [(p[0], 1, p[1], p[2]) for p in permutations((0, 2, 3))]
        elif family == "angle-angle":
            perms = [(0, 1, 2, 3), (3, 1, 2, 0)]
        else:
            perms = [tuple(range(len(query))), tuple(reversed(range(len(query))))]
        matches = [
            (r, p)
            for r in rows[family, "cff91"]
            for p in perms
            if r["types"] == [query[i] for i in p]
        ]
        if not matches:
            continue
        latest = max(r["version"] for r, _ in matches)
        matches = [(r, p) for r, p in matches if r["version"] == latest]
        values = []
        for r, perm in matches:
            v = oriented_values(r, family, perm)
            # Exact IEEE comparison requires the declared multiplication order:
            # degree * (pi/180), rather than (degree*pi)/180. Raw Decimal rows
            # above remain the authoritative evidence; no tolerance is changed.
            angular = (
                [0]
                if family == "quartic_angle"
                else [1, 3, 5]
                if family == "torsion_3"
                else [1]
                if family == "wilson_out_of_plane"
                else []
            )
            for i in angular:
                v[i] = float(r["values"][i]) * (pi / 180)
            values.append(v)
        return {
            "status": "assigned"
            if all(v == values[0] for v in values)
            else "ambiguous",
            "row_ids": sorted({r["id"] for r, _ in matches}),
            "values": values[0],
            "ordinary_evidence": evidence,
        }
    return {"status": "missing", "row_ids": []}


def audit(study, distribution):
    study, distribution = Path(study), Path(distribution)
    report = json.loads((study / "report.json").read_text())
    for filename, pin in report["source_hashes"].items():
        assert sha256((distribution / filename).read_bytes()).hexdigest() == pin
    raw = (distribution / "pcff.frc").read_bytes()
    lines = raw.decode().splitlines()
    rows = read_source(raw)
    indexed = {r["id"]: r for group in rows.values() for r in group if "types" in r}
    numerical = json.loads((study / "source-rows.json").read_text())
    explicit_zero_columns = 0
    for rid, r in numerical.items():
        assert r["raw"] == lines[r["line"] - 1]
        assert r["source_sha256"] == report["source_hashes"]["pcff.frc"]
        assert r["source_file"] == "pcff.frc"
        word = r["raw"].split("!", 1)[0].split()
        k = ARITY[r["family"]]
        assert word[2 : 2 + k] == r["types"]
        assert word[2 + k :] == r["original_values"] == r["converted_values"]
        assert list(map(Decimal, r["original_values"])) == indexed[rid]["values"]
        zeros = [i for i, v in enumerate(r["original_values"]) if Decimal(v) == 0]
        assert zeros == r["explicit_zero_columns"]
        explicit_zero_columns += len(zeros)
    fixtures = []
    for fixture in report["fixtures"]:
        name = fixture["fixture"]
        system = system_from(unpack((study / f"{name}-graph.json").read_text()))
        p = unpack((study / f"{name}-parameters.json").read_text())
        labels = p["typing"]["assignments"]
        adj = {i: set() for i in system.topology.sites}
        for b in system.topology.bonds.values():
            adj[b.site1].add(b.site2)
            adj[b.site2].add(b.site1)
        for i, a in system.topology.sites.items():
            neighbors = [system.topology.sites[j].element for j in adj[i]]
            expected = {"Si": "sio", "O": "osi", "C": "c3"}.get(a.element)
            if a.element == "H":
                assert len(neighbors) == 1
                expected = {"Si": "hsi", "C": "hc"}[neighbors[0]]
            if a.element == "O":
                assert sorted(neighbors) == ["Si", "Si"]
            if a.element == "C":
                assert sorted(neighbors) == ["H", "H", "H", "Si"]
            assert labels[i] == expected
            atom = p["typing"]["entries"][i]["source_atom_record"]["record"]
            assert atom["raw"] == lines[atom["line"] - 1]
            assert atom["raw"].split()[2] == expected
            assert Decimal(atom["raw"].split()[3]) == Decimal(str(atom["data"]["mass"]))
        expected_requests = {
            f + ":" + ",".join(map(str, sites))
            for f, sites in raw_inventory(system)
            if f != "wilson_out_of_plane"
        }
        assert expected_requests == set(p["interactions"])
        for key, q in p["queries"].items():
            family, csv = key.split(":")
            independent = strict_resolution(rows, family, csv.split(","))
            assert independent["status"] == q["status"], key
            if q["status"] == "assigned":
                assert independent["row_ids"] == sorted(
                    {r["record_id"] for r in q["selected"]}
                )
                assert independent["values"] == q["normalized_values"]
        native_charge = charge_raw(rows, system, labels)
        assert native_charge["complete"] == p["native_charge_record"]["complete"]
        assert {i: float(q) for i, q in native_charge["partial_charges"].items()} == p[
            "exact_native_charge_vector"
        ]
        assert all(Decimal(c["residual"]) == 0 for c in native_charge["components"])
        for r in p["source_rows"].values():
            assert numerical[r["id"]] == r
        fixtures.append(
            {
                "fixture": name,
                "raw_requests": len(expected_requests),
                "unique_queries": len(p["queries"]),
                "charge_vector_decimal_e": native_charge["partial_charges"],
                "charges": native_charge,
                "missing_requests": sum(
                    r["status"] != "assigned" for r in p["interactions"].values()
                ),
            }
        )
    return {
        "status": "independent_source_reconstruction_passed",
        "numerical_rows": len(numerical),
        "explicit_zero_columns_source_evidenced": explicit_zero_columns,
        "fixtures": fixtures,
        "execution_profiles": [],
        "acceptance_gate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.study, args.distribution)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "fixtures"}))


if __name__ == "__main__":
    main()
