"""Independent, bounded C/H/O/S reference labels and converter BB13 checks.

These declared acceptance rules do not import the production PCFF typer/resolver.
They are not a new runtime typing profile.
"""

from collections import Counter

import numpy as np
from pcff_j5_reference import charge_raw, raw_inventory, read_source, resolve_raw

from island.charge_references.records import unpack


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def labels_for(system):
    atoms = system.topology.sites
    adj = {i: set() for i in atoms}
    for b in system.topology.bonds.values():
        require(b.order == 1, "Reference domain requires saturated single bonds")
        i, j = b.site1, b.site2
        adj[i].add(j)
        adj[j].add(i)
    labels = {}
    for i, a in atoms.items():
        e = Counter(atoms[j].element for j in adj[i])
        require(not a.formal_charge, "Reference domain is neutral")
        if a.element == "C" and len(adj[i]) == 4 and set(e) <= {"C", "H", "O", "S"}:
            labels[i] = "c" + str(e["H"]) if e["H"] else "c"
        elif a.element in {"O", "S"} and len(adj[i]) == 2:
            require(
                e in ({"C": 2}, {"C": 1, "H": 1}), "Unexpected heteroatom environment"
            )
            labels[i] = ("o" if a.element == "O" else "s") + ("h" if e["H"] else "c")
        elif a.element == "H" and len(adj[i]) == 1:
            p = atoms[next(iter(adj[i]))].element
            require(p in {"C", "O", "S"}, "Unexpected H parent")
            labels[i] = {"C": "hc", "O": "ho", "S": "hs"}[p]
        else:
            raise ValueError("Outside declared independent reference domain")
    return labels, adj


def selection_check(system, model, source, *, reference_labels=None):
    if reference_labels is None:
        labels, adj = labels_for(system)
    else:
        labels = dict(reference_labels)
        adj = {i: set() for i in system.topology.sites}
        for b in system.topology.bonds.values():
            adj[b.site1].add(b.site2)
            adj[b.site2].add(b.site1)
        require(set(labels) == set(adj), "Independent label coverage mismatch")
    raw = read_source(source.read_bytes())
    assignment = unpack(model.assignment.json_text)

    def key(f, s):
        s = tuple(s)
        if f == "angle-angle":
            s = min(s, (s[3], s[1], s[2], s[0]))
        elif f == "wilson_out_of_plane":
            p = sorted((s[0], s[2], s[3]))
            s = (p[0], s[1], *p[1:])
        else:
            s = min(s, s[::-1])
        return f, s

    expected = Counter(key(f, s) for f, s in raw_inventory(system))
    actual = Counter(key(a["family"], a["sites"]) for a in assignment["assignments"])
    require(expected == actual, "Independent interaction inventory mismatch")
    rows, zeros = [], []
    terms = {t["assignment_id"]: t for t in model.payload["terms"]}
    for a in assignment["assignments"]:
        family, sites = a["family"], a["sites"]
        require(
            a["supplied_types"] == [labels[i] for i in sites],
            "Independent labels differ",
        )
        if family == "wilson_out_of_plane" and len(adj[sites[1]]) != 3:
            require(a["status"] == "not_applicable", "Invalid Wilson applicability")
            continue
        q = resolve_raw(raw, family, [labels[i] for i in sites], guarded_msi=True)
        if (
            family == "bond-bond_1_3"
            and q["status"] == "missing"
            and "cp" not in labels.values()
        ):
            require(a["status"] == "missing", "Missing row rewritten as assigned")
            term = terms[a["id"]]
            require(
                term["origin"] == "converter_derived_zero"
                and term["source_rows"] == []
                and term["coefficients"] == [0.0],
                "False BB13 provenance",
            )
            eq = [
                resolve_raw(
                    raw, "quartic_bond", [labels[i] for i in pair], guarded_msi=True
                )
                for pair in (sites[:2], sites[2:])
            ]
            require(
                all(e["status"] == "assigned" for e in eq),
                "Missing reference BB13 equilibrium",
            )
            np.testing.assert_allclose(
                term["equilibria"], [e["values"][0] for e in eq], atol=1e-12, rtol=0
            )
            zeros.append(
                {
                    "sites": sites,
                    "origin": "converter_initialization",
                    "raw_search": q,
                    "equilibrium_rows": [e["selected_ids"] for e in eq],
                }
            )
            continue
        require(
            q["status"] == a["status"] == "assigned",
            "Independent source coverage failure",
        )
        require(
            q["selected_ids"] == sorted({r["record_id"] for r in a["selected"]}),
            "Selected-row mismatch",
        )
        np.testing.assert_allclose(
            q["values"], a["normalized_values"], atol=1e-12, rtol=0
        )
        rows.append(
            {
                "family": family,
                "sites": sites,
                "source_rows": q["selected_ids"],
                "independent_values": q["values"],
            }
        )
    q = charge_raw(raw, system, labels)
    require(q["complete"], "Independent native charges incomplete")
    charges = {i: float(v) for i, v in q["partial_charges"].items()}
    np.testing.assert_allclose(
        [charges[n["site"]] for n in model.payload["nonbonded"]],
        [n["charge"] for n in model.payload["nonbonded"]],
        atol=1e-12,
        rtol=0,
    )
    return (
        labels,
        charges,
        {
            "inventory": dict(Counter(f for f, s in expected.elements())),
            "rows": rows,
            "converter_zeros": zeros,
            "independent_charge": q,
        },
    )


def converter_zero_check(sections):
    """Check actual compiled converter's coefficients and its own bond inventory."""
    bonds = {int(r[0]): list(map(float, r[1:])) for r in sections["Bond Coeffs"]}
    by_pair = {
        tuple(sorted((int(r[2]), int(r[3])))): bonds[int(r[1])][0]
        for r in sections["Bonds"]
    }
    bb13 = {int(r[0]): list(map(float, r[1:])) for r in sections["BondBond13 Coeffs"]}
    for row in sections["Dihedrals"]:
        a, b, c, d = map(int, row[2:])
        expected = [0.0, by_pair[tuple(sorted((a, b)))], by_pair[tuple(sorted((c, d)))]]
        np.testing.assert_allclose(bb13[int(row[1])], expected, atol=1e-12, rtol=0)
    return {
        "dihedrals": len(sections["Dihedrals"]),
        "bb13_types": len(bb13),
        "all_zero_with_reference_equilibria": True,
    }
