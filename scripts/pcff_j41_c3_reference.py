"""Independent ASCII oracle for authorized single-source Si/P comparisons.

No alternate-source chemical typing or model publication is authorized here.
Trial labels test explicit source environments, not a transferable type system.
"""

from decimal import Decimal
from hashlib import sha256

from island.workflows import storage
from island.workflows.bundle import system_data
from scripts.pcff_j5_reference import ORDINARY, ROLE, match, raw_inventory, read_source
from scripts.pcff_j41_c1_reference import mixed_inventory
from scripts.pcff_j41_reference import require


def audit_inventory(variant):
    p = variant.payload
    independently = mixed_inventory(variant.raw)
    require(
        [(r["line"], r["raw"], r["fields"]) for r in p["rows"]]
        == [(r["line"], r["raw"], r["fields"]) for r in independently["rows"]],
        "Independent raw inventory differs",
    )
    numeric = read_source(variant.raw)
    indexed = {r["id"]: r for r in p["rows"]}
    count = 0
    for rows in numeric.values():
        for row in rows:
            if "values" not in row:
                continue
            require(
                row["types"] == indexed[row["id"]]["types"]
                and row["values"] == list(map(Decimal, indexed[row["id"]]["values"])),
                "Independent coefficient differs",
            )
            count += 1
    return {
        "variant_identity": variant.identity,
        "raw_rows_checked": len(p["rows"]),
        "numeric_rows_checked": count,
        "mixed_silicon_inventory": independently["environments"],
        "family_counts": independently["source_section_counts"],
        "source_sha256": p["sha256"],
    }


def compass_labels(system, variant):
    """Narrow diagnostic mappings supported by the source's printed templates.

    The P probe is neutral cyclic H-substituted phosphazene, not phosphate.
    Oxo P is deliberately not mapped to the broad four-neighbor P template.
    """
    atoms = system.topology.sites
    bonds = list(system.topology.bonds.values())
    adj = {i: [] for i in atoms}
    for b in bonds:
        adj[b.site1].append((b.site2, b.order))
        adj[b.site2].append((b.site1, b.order))
    labels, evidence, blocked = {}, [], []
    templates = variant.payload["templates"][variant.filename]
    declarations = {r["types"][0]: r for r in variant.payload["atom_types"]}
    for i, atom in atoms.items():
        ns = adj[i]
        elements = [atoms[j].element for j, _ in ns]
        orders = sorted(o for _, o in ns)
        label = None
        if atom.formal_charge:
            pass
        elif atom.element == "Si" and orders == [1, 1, 1, 1]:
            if set(elements) <= {"H", "C", "O"}:
                label = "si4c" if "O" in elements and "H" not in elements else "si4"
        elif atom.element == "C" and orders == [1, 1, 1, 1]:
            if set(elements) <= {"H", "C", "Si"}:
                label = "c4"
        elif atom.element == "O" and orders == [1, 1] and "Si" in elements:
            if set(elements) <= {"Si", "H"}:
                label = "o2z"
        elif atom.element == "H" and len(ns) == 1:
            if set(elements) <= {"C", "Si", "P"}:
                label = "h1"
        elif atom.element == "P" and orders == [1, 1, 1, 2]:
            if sorted(elements) == ["H", "H", "N", "N"]:
                label = "p4="
        elif atom.element == "N" and orders == [1, 2] and elements == ["P", "P"]:
            label = "n2="
        if label is None or label not in declarations:
            blocked.append(
                {
                    "site": i,
                    "element": atom.element,
                    "formal_charge": atom.formal_charge,
                    "neighbor_elements": elements,
                    "bond_orders": orders,
                    "reason": "no_unambiguous_state_specific_source_mapping",
                }
            )
        else:
            labels[i] = label
            evidence.append(
                {
                    "site": i,
                    "type": label,
                    "atom_row": declarations[label]["id"],
                    "template_lines": [
                        t["line"] for t in templates if t["type"] == label
                    ],
                }
            )
    return {
        "labels": labels,
        "evidence": evidence,
        "blockers": blocked,
        "status": "diagnostic_candidate_only" if not blocked else "typing_blocked",
        "chemical_typing_qualification": False,
        "phosphorus_scope": "phosphazene source reference 4 only; no phosphate/phosphonate oxidation-state transfer",
    }


def resolve_compass(rows, family, types):
    stages = [(list(types), [])]
    eq = rows["equivalence", "compass"]
    mapped, evidence = [], []
    for label in types:
        found = [r for r in eq if r["type"] == label]
        if not found:
            break
        top = max(r["version"] for r in found)
        winners = [r for r in found if r["version"] == top]
        require(
            len({tuple(r["map"].items()) for r in winners}) == 1,
            "Ambiguous COMPASS equivalence",
        )
        role = "bond" if family == "bond_increments" else ROLE[family]
        require(role in ORDINARY, "Unknown equivalence role")
        mapped.append(winners[0]["map"][role])
        evidence.extend(r["id"] for r in winners)
    if len(mapped) == len(types) and mapped != list(types):
        stages.append((mapped, evidence))
    searches = []
    for labels, ids in stages:
        result = match(rows, family, "compass", labels, guarded_msi=True)
        table = rows[family, "compass"]
        searches.append(
            {
                "family": family,
                "namespace": "compass",
                "types": labels,
                "equivalence_ids": ids,
                "result": result,
                "row_ids": [r["id"] for r in table],
                "raw_table_sha256": sha256(
                    "\n".join(r["raw"] for r in table).encode()
                ).hexdigest(),
            }
        )
        if result["status"] != "missing":
            return {**result, "searches": searches}
    return {"status": "missing", "candidate_ids": [], "searches": searches}


def compass_trial(system, variant):
    mapping = compass_labels(system, variant)
    out = {
        "source_sha256": variant.expected_sha256,
        "system_sha256": sha256(
            storage.json_bytes(storage.encode(system_data(system)))
        ).hexdigest(),
        "mapping": mapping,
        "interactions": [],
        "charges": None,
        "source_complete": False,
        "execution_authorized": False,
        "execution_profile": None,
    }
    if mapping["blockers"]:
        return out
    labels = mapping["labels"]
    rows = read_source(variant.raw)
    degrees = {i: 0 for i in labels}
    for b in system.topology.bonds.values():
        degrees[b.site1] += 1
        degrees[b.site2] += 1
    for family, sites in raw_inventory(system):
        ts = [labels[i] for i in sites]
        result = (
            {
                "status": "not_applicable",
                "reason": "Wilson term requires degree-three center",
                "searches": [],
            }
            if family == "wilson_out_of_plane" and degrees[sites[1]] != 3
            else resolve_compass(rows, family, ts)
        )
        out["interactions"].append(
            {"family": family, "sites": list(sites), "types": ts, "resolution": result}
        )
    charge = {i: Decimal(0) for i in labels}
    increments, missing = [], []
    indexed = {r["id"]: r for rs in rows.values() for r in rs}
    for b in system.topology.bonds.values():
        sites = [b.site1, b.site2]
        ts = [labels[i] for i in sites]
        q = resolve_compass(rows, "bond_increments", ts)
        if q["status"] != "assigned":
            missing.append({"sites": sites, "types": ts, "resolution": q})
            continue
        # Match's oriented values use energy conversions; increments are charges.
        # Read original decimal endpoints and orient explicitly instead.
        choices = []
        for rid in q["selected_ids"]:
            row = indexed[rid]
            query = q["searches"][-1]["types"]
            choices.append(
                row["values"] if row["types"] == query else row["values"][::-1]
            )
        require(all(v == choices[0] for v in choices), "Conflicting charge increments")
        for i, value in zip(sites, choices[0], strict=True):
            charge[i] += value
        increments.append(
            {
                "sites": sites,
                "types": ts,
                "resolution": q,
                "increments_e": list(map(str, choices[0])),
            }
        )
    formal = sum(a.formal_charge for a in system.topology.sites.values())
    total = sum(charge.values())
    adj = {i: set() for i in labels}
    for b in system.topology.bonds.values():
        adj[b.site1].add(b.site2)
        adj[b.site2].add(b.site1)
    remaining, components = set(labels), []
    while remaining:
        todo, group = [min(remaining)], set()
        while todo:
            i = todo.pop()
            if i not in group:
                group.add(i)
                todo.extend(adj[i] - group)
        remaining -= group
        component_formal = sum(system.topology.sites[i].formal_charge for i in group)
        component_total = sum(charge[i] for i in group)
        components.append(
            {
                "sites": sorted(group),
                "formal_charge_e": component_formal,
                "total_e": str(component_total),
                "complete": not any(set(q["sites"]) & group for q in missing)
                and abs(component_total - component_formal) <= Decimal("1e-12"),
            }
        )
    out["charges"] = {
        "vector_e": {i: str(v) for i, v in charge.items()},
        "increments": increments,
        "missing": missing,
        "total_e": str(total),
        "formal_charge_e": formal,
        "components": components,
        "unit": "elementary_charge",
        "tolerance_e": "1e-12",
        "complete": not missing
        and all(c["complete"] for c in components)
        and abs(total - formal) <= Decimal("1e-12"),
        "origin": "COMPASS/source_native_bond_increments",
        "normalized": False,
    }
    out["source_complete"] = out["charges"]["complete"] and all(
        q["resolution"]["status"] in ("assigned", "not_applicable")
        for q in out["interactions"]
    )
    return out


def validate_compass_reconstruction(system, variant, value):
    require(
        value == compass_trial(system, variant),
        "Contradictory source, graph, ordering, charge, cross-term or qualification evidence",
    )
    return True


def trial_searches(trial):
    queries = {}
    resolutions = [q["resolution"] for q in trial["interactions"]]
    if trial["charges"]:
        resolutions += [q["resolution"] for q in trial["charges"]["increments"]]
        resolutions += [q["resolution"] for q in trial["charges"]["missing"]]
    for q in resolutions:
        for s in q["searches"]:
            queries[s["family"], s["namespace"], tuple(s["types"])] = {
                k: s[k] for k in ("family", "namespace", "types", "result")
            }
    return [queries[k] for k in sorted(queries)]
