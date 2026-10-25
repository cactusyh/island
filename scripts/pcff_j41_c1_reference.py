"""Independent source inventory for mixed silicon; never a chemical typer."""

import ctypes
import re
from collections import Counter, defaultdict
from hashlib import sha256

from scripts.pcff_j5_reference import ARITY, match, read_source
from scripts.pcff_j41_reference import require

NEIGHBORS = ("C", "H", "O", "N", "S", "F", "Cl", "Br", "I", "Si", "P")


def mixed_inventory(raw):
    rows, section, namespace = [], None, ""
    for n, line in enumerate(raw.decode("ascii").splitlines(), 1):
        words = line.split("!", 1)[0].split()
        if not words:
            continue
        if words[0].startswith("#"):
            section, namespace = words[0][1:], " ".join(words[1:])
        elif section and re.fullmatch(r"\d+(?:\.\d+)?", words[0]) and len(words) >= 3:
            rows.append(
                {
                    "section": section,
                    "namespace": namespace,
                    "line": n,
                    "raw": line,
                    "fields": words,
                }
            )
    elements = defaultdict(set)
    for r in rows:
        if r["section"] == "atom_types":
            elements[r["fields"][2]].add(r["fields"][4])
    # One-step aliases are family-dependent candidate inventory only, not rules.
    for r in rows:
        if r["section"] in ("equivalence", "auto_equivalence"):
            for token in r["fields"][3:]:
                elements[token].update(elements[r["fields"][2]])
    silicon = {t for t, es in elements.items() if "Si" in es}
    out = {}
    families = sorted(
        set(ARITY)
        | {"atom_types", "equivalence", "auto_equivalence", "torsion-torsion"}
    )
    for partner in NEIGHBORS:
        buckets = {f: [] for f in families}
        for r in rows:
            family = r["section"]
            n = ARITY.get(family, 1)
            ts = r["fields"][2 : 2 + n]
            common = family in (
                "atom_types",
                "equivalence",
                "auto_equivalence",
                "nonbond(9-6)",
            )
            has_si = bool(set(ts) & silicon)
            has_partner = any(partner in elements[t] for t in ts)
            if (has_si and has_partner) or (common and (has_si or has_partner)):
                buckets.setdefault(family, []).append(
                    f"{family}:{r['namespace']}:{r['line']}"
                )
        out["Si-" + partner] = {
            "families": buckets,
            "counts": {k: len(v) for k, v in buckets.items()},
            "interpretation": "candidate source inventory, not chemical or completeness authority; empty families remain empty",
        }
    return {
        "schema": "island.j41-c1-independent-inventory.v1",
        "source_sha256": sha256(raw).hexdigest(),
        "environments": out,
        "rows": rows,
        "wildcard_candidates": [
            r
            for r in rows
            if any(
                t.startswith("*")
                for t in r["fields"][2 : 2 + ARITY.get(r["section"], 1)]
            )
        ],
        "source_section_counts": dict(
            sorted(Counter(r["section"] for r in rows).items())
        ),
        "interpretation": "one-step ordinary/automatic aliases only; no inferred source row, wildcard authorization, missing zero or P type",
    }


def compiled_source_audit(source, library):
    """Compare raw Si/P rows and compiled C leaf selection with native catalog."""
    from island.forcefields.pcff import inspect_pcff_full_source
    from island.forcefields.pcff.fallbacks import MSI_POLICY, lookup

    raw = read_source(source.raw)
    inv = mixed_inventory(source.raw)
    indexed = {r["id"]: r for r in inspect_pcff_full_source(source)["records"]}
    checked = []
    queries = set()
    for row in inv["rows"]:
        key = f"{row['section']}:{row['namespace']}:{row['line']}"
        if key in indexed:
            native = indexed[key]
            require(native["raw"] == row["raw"], "Raw source mismatch")
            require(
                list(map(float, row["fields"][2 + len(native["types"]) :]))
                == native["original_values"],
                "Coefficient mismatch",
            )
            checked.append(key)
            if set(native["types"]) & {"si", "sio", "sz", "p", "p="} and not any(
                t.startswith("*") for t in native["types"]
            ):
                queries.add((row["section"], row["namespace"], tuple(native["types"])))
    for f in (
        "bond-bond_1_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
    ):
        for ts in (("hc", "c3", "si", "hsi"), ("hc", "c3", "sio", "c3")):
            queries.add((f, "cff91", ts))
    atom = ctypes.c_char * 5

    class Row(ctypes.Structure):
        _fields_ = [
            ("version", ctypes.c_float),
            ("ref", ctypes.c_int),
            ("types", atom * 6),
            ("parameters", ctypes.c_double * 8),
        ]

    lib = ctypes.CDLL(str(library))
    lib.select_rows.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(atom),
        ctypes.POINTER(Row),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.select_rows.restype = ctypes.c_int
    results = []
    for family, ns, query in sorted(queries):
        slots = {}
        for r in raw[family, ns]:
            key = tuple(r["types"])
            if key not in slots or r["version"] > slots[key]["version"]:
                slots[key] = r
        table = list(slots.values())
        data = (Row * len(table))()
        for dst, r in zip(data, table):
            dst.version = float(r["version"])
            for i, t in enumerate(r["types"]):
                require(len(t) <= 4, "C matcher label overflow")
                dst.types[i].value = t.encode()
        q = (atom * len(query))(*(atom(*t.encode()) for t in query))
        back = ctypes.c_int(-1 if family == "angle-angle" else 0)
        ix = lib.select_rows(len(query), len(table), q, data, ctypes.byref(back))
        selected = table[ix] if ix >= 0 else None
        native = lookup(family, list(query), ns, indexed, [], policy=MSI_POLICY)
        agrees = (selected is None and native["status"] == "missing") or (
            selected is not None
            and selected["id"] in {r["record_id"] for r in native.get("candidates", [])}
        )
        require(agrees, (query, selected, native))
        independent = match(raw, family, ns, list(query), guarded_msi=True)
        require(independent["status"] == native["status"], (query, independent, native))
        if native["status"] == "assigned":
            import numpy as np

            require(
                selected["id"] in {r["record_id"] for r in native["selected"]},
                "Converter selected a superseded row",
            )
            np.testing.assert_allclose(
                independent["values"],
                native["normalized_values"],
                atol=1e-10,
                rtol=1e-13,
            )
        results.append(
            {
                "family": family,
                "namespace": ns,
                "query": query,
                "compiled_row": selected["id"] if selected else None,
                "compiled_source_values": list(map(str, selected["values"]))
                if selected
                else None,
                "native": native,
                "independent_raw_match": independent,
                "agrees": agrees,
            }
        )
    return {
        "raw_rows_checked": len(checked),
        "raw_record_ids": checked,
        "compiled_queries": len(results),
        "queries": results,
        "native_status_counts": dict(Counter(r["native"]["status"] for r in results)),
        "scope": "leaf selection only; no chemical typing, charge assignment or complete assembly implied",
    }
