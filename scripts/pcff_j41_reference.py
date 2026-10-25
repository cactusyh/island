"""Independent raw-line Si/P inventory and graph oracle; no runtime selection imports."""

import re
from collections import Counter
from hashlib import sha256

from scripts.pcff_j5_reference import (
    ARITY,
    charge_raw,
    raw_inventory,
    read_source,
    resolve_raw,
)


def inventory(raw):
    rows, section, namespace = ([], None, "")
    for n, line in enumerate(raw.decode("ascii").splitlines(), 1):
        words = line.split("!", 1)[0].split()
        if not words:
            continue
        if words[0].startswith("#"):
            section, namespace = (words[0][1:], " ".join(words[1:]))
        elif (
            section and re.fullmatch("\\d+(?:\\.\\d+)?", words[0]) and (len(words) >= 3)
        ):
            rows.append(
                {
                    "section": section,
                    "namespace": namespace,
                    "line": n,
                    "raw": line,
                    "fields": words,
                }
            )
    primary = {
        r["fields"][2]
        for r in rows
        if r["section"] == "atom_types" and r["fields"][4] in ("Si", "P")
    }
    descriptive = {
        r["fields"][2]
        for r in rows
        if r["section"] == "atom_types"
        and re.search("silicon|silane|siloxane|phosph", r["raw"], re.IGNORECASE)
    }
    labels = primary | descriptive
    aliases = set()
    for r in rows:
        if (
            r["section"] in ("equivalence", "auto_equivalence")
            and r["fields"][2] in labels
        ):
            aliases.update(r["fields"][3:])
    selected = []
    for r in rows:
        family = r["section"]
        arity = ARITY.get(family)
        tokens = r["fields"][2 : 2 + arity] if arity else r["fields"][2:3]
        if family in ("equivalence", "auto_equivalence"):
            tokens = r["fields"][2:]
        reasons = []
        if set(tokens) & labels:
            reasons.append("explicit_Si_P_or_described_companion_label")
        if set(tokens) & aliases:
            reasons.append("ordinary_or_automatic_equivalence_reachable")
        if any(t.startswith("*") for t in tokens):
            reasons.append("wildcard_candidate_not_coverage_authority")
        if reasons:
            selected.append({**r, "selection_reasons": reasons})
    return {
        "schema": "island.j41-independent-source-inventory.v1",
        "source_sha256": sha256(raw).hexdigest(),
        "primary_labels": sorted(primary),
        "descriptive_labels": sorted(descriptive),
        "one_step_equivalence_tokens": sorted(aliases),
        "closure_policy": "one-step family/position mappings; generic h/o aliases retain unrelated candidates without authorizing them; all wildcard rows retained",
        "section_counts": dict(sorted(Counter(r["section"] for r in rows).items())),
        "selected_counts": dict(
            sorted(Counter(r["section"] for r in selected).items())
        ),
        "rows": selected,
    }


def labels_for(system):
    """Independent neutral finite saturated Si fixture oracle, never a public typer."""
    atoms = system.topology.sites
    adj = {i: [] for i in atoms}
    for b in system.topology.bonds.values():
        require(b.order == 1 and (not b.aromatic), "Independent J41 control failed")
        adj[b.site1].append(b.site2)
        adj[b.site2].append(b.site1)
    out = {}
    for i, a in atoms.items():
        require(a.formal_charge == 0, "Independent J41 control failed")
        ns = Counter(atoms[j].element for j in adj[i])
        if a.element == "Si":
            require(
                len(adj[i]) == 4 and set(ns) <= {"C", "H", "O"},
                "Independent J41 control failed",
            )
            out[i] = "sio" if ns["O"] else "si"
        elif a.element == "H":
            require(len(adj[i]) == 1, "Independent J41 control failed")
            out[i] = {"Si": "hsi", "C": "hc"}[atoms[adj[i][0]].element]
        elif a.element == "C":
            require(len(adj[i]) == 4, "Independent J41 control failed")
            out[i] = "c" + str(ns["H"]) if ns["H"] else "c"
        elif a.element == "O":
            require(ns == {"Si": 2}, "Independent J41 control failed")
            out[i] = "osi"
        else:
            raise ValueError("Outside independent neutral saturated Si fixture")
    return out


def audit_queries(system, assessment, raw):
    """Independently enumerate the graph and resolve every applicable source row."""
    p = assessment.payload
    labels = labels_for(system)
    require(labels == p["typing"]["assignments"], "Independent J41 control failed")
    queries = p["inspection"]["interaction_queries"]
    expected = set(raw_inventory(system))
    actual = {(q["family"], tuple(q["sites"])) for q in queries}
    require(
        actual == expected and len(actual) == len(queries),
        "Independent J41 control failed",
    )
    rows, checked = (read_source(raw), [])
    for q in queries:
        if q["resolution"]["status"] == "not_applicable":
            require(
                q["family"] == "wilson_out_of_plane", "Independent J41 control failed"
            )
            require(
                sum(
                    b.site1 == q["sites"][1] or b.site2 == q["sites"][1]
                    for b in system.topology.bonds.values()
                )
                == 4,
                "Independent J41 control failed",
            )
            continue
        found = resolve_raw(
            rows, q["family"], [labels[i] for i in q["sites"]], guarded_msi=True
        )
        require(found["status"] == q["resolution"]["status"], (q, found))
        if found["status"] == "assigned":
            import numpy as np

            np.testing.assert_allclose(
                found["values"],
                q["resolution"]["normalized_values"],
                atol=1e-10,
                rtol=1e-13,
            )
        checked.append(
            {
                "family": q["family"],
                "sites": q["sites"],
                "independent_resolution": found,
            }
        )
    charges = charge_raw(rows, system, labels)
    require(
        charges["complete"] == p["inspection"]["charges_complete"],
        "Independent J41 control failed",
    )
    native = p["inspection"]["native_charge_record"]["native_charge_record"]
    for i, value in charges["partial_charges"].items():
        require(
            abs(float(value) - native["partial_charges"][i]) <= 1e-12,
            "Independent J41 control failed",
        )
    return {
        "queries": checked,
        "charges": charges,
        "inventory_equal": True,
        "typing_equal": True,
    }


def require(ok, message):
    if not ok:
        raise ValueError(message)
