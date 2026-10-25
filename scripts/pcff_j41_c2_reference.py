"""C2 independent full-table searches and phosphorus connectivity controls.

The ASCII parser/oracle is separate from the production source parser/resolver.
Compiled leaf selection is a reference observation, never typing authority.
"""

import ctypes
from hashlib import sha256
from types import SimpleNamespace

from scripts.pcff_j5_reference import match, read_source
from scripts.pcff_j41_c1_reference import mixed_inventory
from scripts.pcff_j41_reference import audit_queries, require


def p_audit(system, decision):
    p = decision.payload
    inv = mixed_inventory(decision.distribution.source.raw)
    declarations = [
        r for r in inv["rows"] if r["section"] == "atom_types" and r["fields"][4] == "P"
    ]
    require(
        [r["raw"] for r in declarations]
        == [r["raw"] for r in p["source_evidence"]["atom_rows"]],
        "Changed P declarations",
    )
    # Independent template block scanner, with no production pattern matcher.
    lines = decision.distribution.templates.decode("ascii").splitlines()
    blocks = []
    for i, line in enumerate(lines):
        if line.split()[:1] == ["type:"] and line.split()[1] in ("p", "p="):
            end = next(
                j for j in range(i + 1, len(lines)) if lines[j].strip() == "end_type"
            )
            blocks.append(
                {
                    "type": line.split()[1],
                    "line": i + 1,
                    "raw": "\n".join(lines[i : end + 1]),
                }
            )
    require(
        blocks == p["source_evidence"]["template_blocks"], "Changed P template evidence"
    )
    result = []
    atoms = system.topology.sites
    for row in p["phosphorus_decisions"]:
        site = row["chemical_state"]["site"]
        bonds = [
            b for b in system.topology.bonds.values() if site in (b.site1, b.site2)
        ]
        neighbors = [b.site2 if b.site1 == site else b.site1 for b in bonds]
        require(row["degree"] == len(neighbors), "Changed P degree")
        require(
            row["hydrogen_neighbors"]
            == sum(atoms[i].element == "H" for i in neighbors),
            "Changed P hydrogen substitution",
        )
        require(
            row["selected_type"] is None and row["status"] == "typing_blocked",
            "Unjustified P selection",
        )
        candidates = []
        for decl, candidate in zip(declarations, row["candidates"], strict=True):
            w = decl["fields"]
            require(
                candidate["type"] == w[2]
                and candidate["source_line"] == decl["line"]
                and candidate["source_raw"] == decl["raw"],
                "Wrong P source candidate",
            )
            degree_ok = len(neighbors) == int(w[5])
            require(
                candidate["degree_matches"] == degree_ok
                and not candidate["authorized"],
                "Unjustified candidate authorization",
            )
            candidates.append(
                {
                    "type": w[2],
                    "source_line": decl["line"],
                    "declared_connections": int(w[5]),
                    "degree_matches": degree_ok,
                }
            )
        orders = sorted(b.order for b in bonds)
        oxo_match = orders == [1, 1, 1, 2] and any(
            b.order == 2
            and atoms[b.site2 if b.site1 == site else b.site1].element
            in ("O", "S", "N")
            for b in bonds
        )
        result.append(
            {
                "site": site,
                "degree": len(neighbors),
                "bond_order_sum": sum(orders),
                "formal_charge": atoms[site].formal_charge,
                "candidates": candidates,
                "oxo_template_matches": oxo_match,
                "oxo_template_connections_conflict": oxo_match
                and all(
                    int(r["fields"][5]) != len(neighbors)
                    for r in declarations
                    if r["fields"][2] == "p="
                ),
                "selected_type": None,
            }
        )
    return {
        "declarations": declarations,
        "templates": blocks,
        "sites": result,
        "converter_does_not_perceive_types": True,
    }


def source_audit(system, decision):
    p = decision.payload
    phosphorus = p_audit(system, decision)
    legacy = p["c1_assessment"]
    if legacy["inspection"] is None:
        return {
            "phosphorus": phosphorus,
            "status": "blocked_at_typing",
            "queries": [],
            "search_tables": {},
        }
    graph_audit = audit_queries(
        system, SimpleNamespace(payload=legacy), decision.distribution.source.raw
    )
    raw = read_source(decision.distribution.source.raw)
    tables, searches = {}, {}
    for q in graph_audit["queries"]:
        for s in q["independent_resolution"].get("searches", []):
            f, ns = s["family"], s["namespace"]
            key = f + ":" + ns
            table = raw[f, ns]
            tables[key] = {
                "row_ids": [r["id"] for r in table],
                "raw_table_sha256": sha256(
                    "\n".join(r["raw"] for r in table).encode("ascii")
                ).hexdigest(),
                "row_count": len(table),
            }
            query_key = (f, ns, tuple(s["types"]))
            searches[query_key] = {
                "family": f,
                "namespace": ns,
                "types": s["types"],
                "result": match(raw, f, ns, s["types"], guarded_msi=True),
            }
    return {
        "phosphorus": phosphorus,
        "status": "complete" if p["source_complete"] else "source_incomplete",
        "queries": graph_audit["queries"],
        "charges": graph_audit["charges"],
        "search_tables": tables,
        "unique_searches": [searches[k] for k in sorted(searches)],
        "interpretation": "all rows of each declared family/namespace searched; no lower-order row supplies an absent coupling",
    }


def compiled_searches(raw_bytes, library, queries):
    rows = read_source(raw_bytes)
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
    for q in queries:
        f, ns, query = q["family"], q["namespace"], q["types"]
        slots = {}
        for r in rows[f, ns]:
            key = tuple(r["types"])
            if key not in slots or r["version"] > slots[key]["version"]:
                slots[key] = r
        table = list(slots.values())
        data = (Row * len(table))()
        for dst, r in zip(data, table, strict=True):
            dst.version = float(r["version"])
            for i, t in enumerate(r["types"]):
                require(len(t) <= 4, "Unsafe C label")
                dst.types[i].value = t.encode("ascii")
        require(all(len(t) <= 4 for t in query), "Unsafe query label")
        query_data = (atom * len(query))(*(atom(*t.encode("ascii")) for t in query))
        direction = ctypes.c_int(-1 if f == "angle-angle" else 0)
        index = lib.select_rows(
            len(query), len(table), query_data, data, ctypes.byref(direction)
        )
        selected = table[index] if index >= 0 else None
        oracle = match(rows, f, ns, query, guarded_msi=True)
        require(oracle == q["result"], "Changed raw query")
        require(
            (selected is None and oracle["status"] == "missing")
            or (
                selected is not None
                and selected["id"]
                in oracle[
                    "selected_ids"
                    if oracle["status"] == "assigned"
                    else "candidate_ids"
                ]
            ),
            "Compiled full-table search disagrees",
        )
        results.append(
            {
                "family": f,
                "namespace": ns,
                "types": query,
                "compiled_row": selected["id"] if selected else None,
                "compiled_values": list(map(str, selected["values"]))
                if selected
                else None,
                "independent": oracle,
                "reverse": direction.value,
                "scope": "selection observation; ambiguous candidates remain blocked",
            }
        )
    return results
