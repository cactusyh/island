"""Independent pinned-source obstructions to full-source acceptance (not a model).

Uses the separately authored Decimal reader, never native selected coefficients.
Missing source information cannot be cured by marking a coverage ledger complete.
"""

import argparse
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from pcff_j5_reference import charge_raw, equivalents, read_source, resolve_raw

from island.forcefields.pcff import (
    load_pcff_automatic_record,
    load_pcff_parameters,
    load_pcff_source,
)
from island.workflows import storage
from island.workflows.bundle import system_from

PIN = "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"


def increment_searches(rows, types):
    """All declared charge paths, preserving endpoint roles and candidates."""
    paths = [("direct", list(types), [])]
    for family, column in (
        ("equivalence", "bond"),
        ("auto_equivalence", "bond_increment"),
    ):
        mapped, ids = equivalents(rows, types, family, [column, column])
        if mapped is not None:
            paths.append((f"{family}.{column}", mapped, ids))
    result = []
    for path, query, evidence in paths:
        matches = []
        for row in rows["bond_increments", "cff91_auto"]:
            for order in ((0, 1), (1, 0)):
                if row["types"] == [query[i] for i in order]:
                    matches.append(
                        {
                            "row": row["id"],
                            "version": str(row["version"]),
                            "orientation": list(order),
                            "endpoint_values": [
                                str(row["values"][order.index(i)]) for i in range(2)
                            ],
                        }
                    )
        result.append(
            {
                "path": path,
                "types": query,
                "equivalence_rows": evidence,
                "candidates": matches,
            }
        )
    return result


def source_obligations(raw, unresolved_labels):
    if sha256(raw).hexdigest() != PIN:
        raise ValueError("Completion audit requires the exact declared PCFF source")
    rows = read_source(raw)
    atoms = {}
    section = None
    for line, text in enumerate(raw.decode().splitlines(), 1):
        fields = text.split("!", 1)[0].split()
        if not fields:
            continue
        if fields[0].startswith("#"):
            section = fields[0][1:]
        elif section == "atom_types" and fields[0][0].isdigit():
            atoms[fields[2]] = {
                "row": f"atom_types:cff91:{line}",
                "line": line,
                "element": fields[4],
                "connections": int(fields[5]),
                "description": " ".join(fields[6:]),
            }
    if len(atoms) != 133 or len(rows["bond_increments", "cff91_auto"]) != 564:
        raise ValueError("Independent source inventory disagrees with declared pin")
    if (
        len(set(unresolved_labels)) != len(unresolved_labels)
        or not set(unresolved_labels) <= atoms.keys()
    ):
        raise ValueError("Invalid unresolved label inventory")
    types = []
    for label in unresolved_labels:
        atom = atoms[label]
        domain = (
            "zeolite_surface"
            if label
            in {"az", "sz", "oah", "oas", "ob", "osh", "oss", "hb", "hoa", "hos"}
            else "metal"
            if "metal" in atom["description"]
            else "ionic"
            if label in {"Br", "Cl", "ca+"}
            else "organic_alias_or_specialized_environment"
        )
        evidence = [
            {"row": r["id"], "columns": r["map"], "version": str(r["version"])}
            for family in ("equivalence", "auto_equivalence")
            for (f, _), values in rows.items()
            if f == family
            for r in values
            if r["type"] == label
        ]
        types.append(
            {
                "label": label,
                **atom,
                "domain": domain,
                "equivalences": evidence,
                "direct_increment_rows": [
                    r["id"]
                    for r in rows["bond_increments", "cff91_auto"]
                    if label in r["types"]
                ],
                "nonbond_search": resolve_raw(rows, "nonbond(9-6)", [label]),
                "implemented_and_independently_verified": False,
                "required_evidence": (
                    "Metal/surface phase and coordination model; element and zero graph degree alone do not establish its domain"
                    if domain == "metal"
                    else "Explicit oxidation/base-charge and coordination convention; graph bonds cannot imply isotope or surface phase"
                    if domain in ("ionic", "zeolite_surface")
                    else "Disambiguate source alias from existing specific types, then verify graph predicate, charge state, all active terms and forces"
                ),
            }
        )
    guan_rows = [
        r
        for r in rows["bond_increments", "cff91_auto"]
        if r["types"] in (["c+", "nr"], ["h*", "nr"])
    ]
    if len(guan_rows) != 2:
        raise ValueError("Unexpected guanidinium source candidates")
    total = sum(
        (
            sum(r["values"]) * (3 if r["types"] == ["c+", "nr"] else 6)
            for r in guan_rows
        ),
        Decimal(0),
    )
    return {
        "source_sha256": PIN,
        "atom_type_count": len(atoms),
        "increment_count": 564,
        "unresolved_labels": types,
        "required_missing_increment_searches": {
            " / ".join(pair): increment_searches(rows, pair)
            for pair in (("nh+", "cp"), ("nh+", "hn"), ("na", "hn2"))
        },
        "guanidinium": {
            "rows": [
                {
                    "id": r["id"],
                    "types": r["types"],
                    "values": list(map(str, r["values"])),
                }
                for r in guan_rows
            ],
            "printed_total_e": str(total),
            "formal_total_e": 1,
            "residual_e": str(total - 1),
            "tolerance_e": "1e-12",
            "source_native_charge_gate": abs(total - 1) <= Decimal("1e-12"),
            "correction_authorized": False,
        },
        "nonzero_wilson_rows": [
            {"id": r["id"], "types": r["types"], "values": list(map(str, r["values"]))}
            for (f, _), values in rows.items()
            if f == "wilson_out_of_plane"
            for r in values
            if r["values"][1] != 0
        ],
        "full_source_acceptance": False,
        "interpretation": "Necessary source obligations, not global certification or a new charge/model policy",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--declaration", required=True, type=Path)
    p.add_argument("--retained-matrix", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    declaration = storage.read_json(a.declaration)
    if declaration["source_sha256"] != PIN:
        raise ValueError("Declaration source mismatch")
    result = source_obligations(a.source.read_bytes(), declaration["unresolved_labels"])
    rows = read_source(a.source.read_bytes())
    source = load_pcff_source(a.source)
    cases = []
    for name in ("pyridinium", "guanidinium", "nitrogen_amine"):
        directory = a.retained_matrix / name
        paths = [
            directory / filename
            for filename in ("system.json", "typing.json", "charge.json")
        ]
        hashes = {p.name: sha256(p.read_bytes()).hexdigest() for p in paths}
        system = system_from(storage.decode(storage.read_json(paths[0])))
        typing = load_pcff_automatic_record(paths[1], source, system=system)
        charge = load_pcff_automatic_record(paths[2], source, system=system)
        independent = charge_raw(rows, system, typing.assignments)
        if independent["complete"] != charge.complete:
            raise ValueError(f"Independent/native charge status differs: {name}")
        native_partial = charge.payload["native_charge_record"]["partial_charges"]
        if set(native_partial) != set(independent["partial_charges"]):
            raise ValueError(f"Independent/native site inventory differs: {name}")
        error = max(
            abs(Decimal(str(native_partial[i])) - Decimal(value))
            for i, value in independent["partial_charges"].items()
        )
        if error > Decimal("1e-12"):
            raise ValueError(f"Independent/native charge values differ: {name}")
        cases.append(
            {
                "case": name,
                "input_hashes": hashes,
                "typing_identity": typing.identity,
                "charge_identity": charge.identity,
                "native_complete": charge.complete,
                "independent": independent,
                "max_charge_difference_e": str(error),
            }
        )
        if hashes != {p.name: sha256(p.read_bytes()).hexdigest() for p in paths}:
            raise ValueError("Historical input changed during audit")
    result["retained_charge_cases"] = cases
    gaps = []
    for name in ("hetero_peo", "sulfur_thioether", "halogenated"):
        path = a.retained_matrix / name / "assignment.json"
        before = sha256(path.read_bytes()).hexdigest()
        assignment = load_pcff_parameters(path, source)
        for row in assignment.payload["assignments"]:
            if row["status"] in ("assigned", "not_applicable"):
                continue
            checked = resolve_raw(rows, row["requested_family"], row["supplied_types"])
            if row["status"] == "missing" and checked["status"] != "missing":
                raise ValueError(
                    f"Native missing row disagrees with raw-source audit: {name}/{row['id']}"
                )
            gaps.append(
                {
                    "case": name,
                    "assignment_identity": assignment.identity,
                    "file_sha256": before,
                    "native_interaction": row,
                    "independent_source_search": checked,
                }
            )
        if sha256(path.read_bytes()).hexdigest() != before:
            raise ValueError("Historical assignment changed during inspection")
    result["interaction_obligations"] = gaps
    result["unresolved_domain_counts"] = dict(
        Counter(r["domain"] for r in result["unresolved_labels"])
    )
    storage.publish(a.output, storage.json_bytes(result))
    print(
        {
            "full_source_acceptance": False,
            "unresolved_labels": len(result["unresolved_labels"]),
            "charge_cases": {c["case"]: c["native_complete"] for c in cases},
        }
    )
    return 1 if a.require_full_source else 0


if __name__ == "__main__":
    raise SystemExit(main())
