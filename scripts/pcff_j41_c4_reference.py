"""C4 independent full-converter lookup audit, never a runtime profile.

This oracle intentionally differs from the historical ISLAND automatic fallback
resolver: the pinned msi2lmp does not read auto-equivalences or bond increments.
No coefficient from an unsuccessful conversion gains execution authority.
"""

import struct
from decimal import Decimal
from functools import cache
from hashlib import sha256

from island.workflows import storage
from island.workflows.bundle import system_data
from scripts.pcff_j5_reference import ROLE, charge_raw, raw_inventory, read_source
from scripts.pcff_j41_c3_reference import compass_trial
from scripts.pcff_j41_reference import require


def table(rows):
    """SearchAndFill: identical ordered key, greatest C float version, first tie.

    The table slot stays at the first occurrence. Ordinary equivalence keys
    include all six labels, not just the input label.
    """
    slots = {}
    for row in rows:
        key = tuple(row.get("types", [row.get("type"), *row.get("map", {}).values()]))
        if key not in slots or f32(row["version"]) > f32(slots[key]["version"]):
            slots[key] = row
    return list(slots.values())


def f32(value):
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def select(rows, query, *, reverse=True):
    """Exact pass, then wildcard pass, file order; forward before reverse.

    File-order selection is observable converter behavior, not authority to
    resolve conflicting physical wildcard coefficients.
    """
    orientations = [list(query)] + ([list(query)[::-1]] if reverse else [])
    for wildcard in (False, True):
        matches = []
        for index, row in enumerate(table(rows)):
            for direction, labels in enumerate(orientations):
                if all(
                    a == b or (wildcard and a.startswith("*"))
                    for a, b in zip(row["types"], labels, strict=True)
                ):
                    matches.append((index, row, direction))
                    break
        if matches:
            index, row, direction = matches[0]
            return {
                "status": "source_row_selected",
                "row_id": row["id"],
                "table_index": index,
                "reverse": direction,
                "wildcard": wildcard,
                "candidate_ids": [r["id"] for _, r, _ in matches],
                "values": list(map(str, row["values"])),
            }
    return {
        "status": "source_row_absent",
        "row_id": None,
        "table_index": -1,
        "reverse": 0 if reverse else -1,
        "candidate_ids": [],
    }


def equivalent(rows, ns, family, labels):
    mapped, evidence = [], []
    for label in labels:
        # find_equiv_type returns the first label; no recursion and no auto map.
        found = next(
            (r for r in table(rows["equivalence", ns]) if r["type"] == label), None
        )
        if found is None:
            return None, evidence
        mapped.append(found["map"][ROLE[family]])
        evidence.append(found["id"])
    return mapped, evidence


def permutations_for(family, labels):
    if family == "wilson_out_of_plane":
        return [
            [labels[i] for i in p]
            for p in (
                (0, 1, 2, 3),
                (0, 1, 3, 2),
                (3, 1, 0, 2),
                (3, 1, 2, 0),
                (2, 1, 0, 3),
                (2, 1, 3, 0),
            )
        ]
    if family == "angle-angle":
        return [list(labels), [labels[i] for i in (3, 1, 2, 0)]]
    return [list(labels)]


def resolve(rows, ns, family, labels):
    """Caller permutations surround find_match; ordinary fallback comes last."""
    if family == "angle-angle":
        return angle_angle(rows, ns, labels)
    mapped, ids = equivalent(rows, ns, family, labels)
    paths = [("direct", list(labels), [])]
    if mapped is not None:
        paths.append(("ordinary_equivalence", mapped, ids))
    searches = []
    for path, query, evidence in paths:
        for oriented in permutations_for(family, query):
            result = select(
                rows[family, ns],
                oriented,
                reverse=family not in ("angle-angle", "wilson_out_of_plane"),
            )
            searches.append(
                {
                    "path": path,
                    "family": family,
                    "namespace": ns,
                    "types": oriented,
                    "equivalence_ids": evidence,
                    "result": result,
                }
            )
            if result["row_id"]:
                # Strict source guard: conflicting wildcard values do not
                # become an authorized physical selection by file order.
                candidates = {r["id"]: r for r in rows[family, ns]}
                ambiguous = (
                    result.get("wildcard", False)
                    and len(
                        {
                            tuple(candidates[i]["values"])
                            for i in result["candidate_ids"]
                        }
                    )
                    > 1
                )
                return {
                    **result,
                    "status": "source_row_ambiguous" if ambiguous else result["status"],
                    "searches": searches,
                }
    return {**searches[-1]["result"], "searches": searches}


def angle_angle(rows, ns, labels):
    """Reproduce the caller's any-of-three success, without promoting its zeros."""
    mapped, ids = equivalent(rows, ns, "angle-angle", labels)
    paths = [("direct", list(labels), [])]
    if mapped is not None:
        paths.append(("ordinary_equivalence", mapped, ids))
    searches = []
    for path, query, evidence in paths:
        results = []
        direction = -1
        for p in ((0, 1, 2, 3), (3, 1, 0, 2), (0, 1, 3, 2)):
            oriented = [query[i] for i in p]
            for trial in permutations_for("angle-angle", oriented):
                # The C caller shares backwards across all three slots. A
                # successful first lookup changes -1 to 0 and enables reversal
                # on subsequent slots. Record it, but never grant authority to
                # a reversed key that moves the central/shared-angle roles.
                result = select(rows["angle-angle", ns], trial, reverse=direction != -1)
                if result["row_id"]:
                    direction = result["reverse"]
                    if direction:
                        result = {
                            **result,
                            "status": "source_row_ambiguous",
                            "reason": "converter_stateful_reverse_changes_angle_roles",
                        }
                else:
                    result["reverse"] = direction
                searches.append(
                    {
                        "path": path,
                        "family": "angle-angle",
                        "namespace": ns,
                        "types": trial,
                        "equivalence_ids": evidence,
                        "result": result,
                    }
                )
                if result["row_id"]:
                    break
            results.append(result)
        if any(r["row_id"] for r in results):
            return {
                **results[0],
                "status": "source_row_ambiguous"
                if any(r["status"] == "source_row_ambiguous" for r in results)
                else results[0]["status"],
                "searches": searches,
                "caller_slots": results,
                "caller_success_is_not_completeness": True,
            }
    return {
        **results[0],
        "searches": searches,
        "caller_slots": results,
        "caller_success_is_not_completeness": True,
    }


def dependencies(family, sites):
    a = list(sites)
    if family in ("bond-bond", "bond-angle"):
        return [("quartic_bond", a[:2]), ("quartic_bond", a[1:])]
    if family in ("end_bond-torsion_3", "bond-bond_1_3"):
        return [("quartic_bond", a[:2]), ("quartic_bond", a[2:])]
    if family == "middle_bond-torsion_3":
        return [("quartic_bond", a[1:3])]
    if family in ("angle-torsion_3", "angle-angle-torsion_1"):
        return [("quartic_angle", a[:3]), ("quartic_angle", a[1:])]
    if family == "angle-angle":
        return [
            ("quartic_angle", [a[i] for i in p])
            for p in ((0, 1, 2), (2, 1, 3), (0, 1, 3))
        ]
    return []


def reconstruct(system, variant, labels, typing):
    """One source and exact final graph; blocked P never receives trial charges."""
    rows = read_source(variant.raw)
    ns = "cff91" if variant.filename == "pcff.frc" else "compass"
    payload = {
        "schema": "island.j41-c4-converter-audit.v1",
        "source_sha256": variant.expected_sha256,
        "companion_hashes": variant.payload["companion_hashes"],
        "variant_identity": variant.identity,
        "system_sha256": sha256(
            storage.json_bytes(storage.encode(system_data(system)))
        ).hexdigest(),
        "atom_order": sorted(system.topology.sites),
        "labels": labels,
        "typing": typing,
        "status": "typing_blocked" if labels is None else "source_incomplete",
        "interactions": [],
        "dependencies": [],
        "charges": None,
        "execution_profile": None,
        "execution_authorized": False,
        "automatic_equivalence_used_by_converter": False,
        "converter_computes_charge_increments": False,
    }
    if labels is None:
        return payload

    @cache
    def lookup(family, query):
        return resolve(rows, ns, family, query)

    require(set(labels) == set(system.topology.sites), "Incomplete supplied labels")
    degrees = {i: 0 for i in labels}
    for b in system.topology.bonds.values():
        degrees[b.site1] += 1
        degrees[b.site2] += 1
    for family, sites in raw_inventory(system):
        types = [labels[i] for i in sites]
        if family == "wilson_out_of_plane" and degrees[sites[1]] != 3:
            result = {
                "status": "not_applicable",
                "reason": "MakeLists degree-three Wilson center only",
            }
        elif family == "bond-bond_1_3" and not any(t.startswith("cp") for t in types):
            # GetParameters.c initializes this slot to zero and does not look it
            # up outside cp. This is a converter suppression, NOT a coefficient.
            result = {
                "status": "converter_zero_initialization",
                "reason": "quo_cp guard skips lookup; no physical coefficient authorized",
                "converter_key": types,
                "independent_physical_lookup": lookup(family, tuple(types)),
            }
        else:
            result = lookup(family, tuple(types))
        payload["interactions"].append(
            {
                "family": family,
                "sites": list(sites),
                "types": types,
                "resolution": result,
            }
        )
        for lower, ids in dependencies(family, sites):
            query = [labels[i] for i in ids]
            resolution = lookup(lower, tuple(query))
            payload["dependencies"].append(
                {
                    "consumer_family": family,
                    "consumer_sites": list(sites),
                    "family": lower,
                    "sites": ids,
                    "types": query,
                    "status": "converter_derived_term"
                    if resolution["status"] == "source_row_selected"
                    else "missing_coupling_dependency",
                    "algorithm": "get_r0/get_t0: first matching graph bond/angle type, equilibrium coefficient slot 0; source units Angstrom/degrees",
                    "resolution": resolution,
                    "value": resolution.get("values", [None])[0],
                }
            )
    if ns == "cff91":
        charge = charge_raw(rows, system, labels)
        payload["charges"] = {
            **charge,
            "origin": "PCFF/source_native_bond_increments",
            "unit": "elementary_charge",
            "normalized": False,
        }
    else:
        trial = compass_trial(system, variant)
        require(trial["mapping"]["labels"] == labels, "Changed COMPASS trial labels")
        payload["charges"] = trial["charges"]
    payload["missing_charge_increment"] = not payload["charges"]["complete"]
    payload["missing_couplings"] = [
        q
        for q in payload["interactions"]
        if q["resolution"]["status"] not in ("source_row_selected", "not_applicable")
    ]
    complete = (
        not payload["missing_couplings"] and not payload["missing_charge_increment"]
    )
    # Even source-complete controls require numerical/workflow qualification.
    payload["status"] = (
        "scientifically_unqualified" if complete else "source_incomplete"
    )
    return payload


def trace_audit(text, variant):
    """Reconstruct real full-converter lookup calls from an instrumented build.

    The original converter is run separately and must agree on exit/output.
    Tracing changes only diagnostic writes, not lookup or initialization.
    """
    rows = read_source(variant.raw)
    ns = "cff91" if variant.filename == "pcff.frc" else "compass"
    checked = []
    for line in text.splitlines():
        if not line.startswith("C4_MATCH "):
            continue
        w = line.split()
        family, n, index, direction, initial = w[1].lstrip("#"), *map(int, w[2:6])
        query = w[6 : 6 + n]
        require(len(w) == 6 + n, "Corrupt converter trace")
        if family == "atom_types":
            # Independent source declarations retain mass as the first value.
            candidates = [
                {
                    "id": r["id"],
                    "types": r["types"],
                    "version": Decimal(r["fields"][0]),
                    "values": [Decimal(r["mass"])],
                }
                for r in variant.payload["atom_types"]
            ]
        else:
            if family == "nonbond":
                family = "nonbond(9-6)"
            candidates = rows[family, ns]
        selected = select(candidates, query, reverse=initial != -1)
        require(
            selected["table_index"] == index, "Raw row disagrees with full converter"
        )
        if index >= 0:
            require(selected["reverse"] == direction, "Converter orientation mismatch")
        checked.append(
            {
                "family": family,
                "types": query,
                "initial_reverse": initial,
                "converter_table_index": index,
                "converter_reverse": direction,
                "independent": selected,
            }
        )
    require(bool(checked), "No executed full-converter trace")
    return checked
