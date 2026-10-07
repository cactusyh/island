"""J16 reference boundary: independent graph labels and raw source requests."""

from pcff_j5_reference import resolve_raw
from pcff_j13_reference import labels_for as cho_labels
from pcff_j15_reference import labels_for as urethane_labels

# Exact pinned GetParameters.c missing-record messages and exit codes. Converter
# matching uses ordinary equivalences; it never supplies automatic heuristics.
MISSING = {
    "bond": (12, "quartic_bond"),
    "angle": (13, "quartic_angle"),
    "bondbond": (14, "bond-bond"),
    "bondangle": (15, "bond-angle"),
    "torsion": (16, "torsion_3"),
    "endbonddihedral": (17, "end_bond-torsion_3"),
    "midbonddihedral": (18, "middle_bond-torsion_3"),
    "angledihedral": (19, "angle-torsion_3"),
    "angleangledihedral": (20, "angle-angle-torsion_1"),
    "bond13": (21, "bond-bond_1_3"),
    "oop": ((22, 23), "wilson_out_of_plane"),
    "angleangle": ((24, 25), "angle-angle"),
}


def labels_for(system, case):
    if case["role"] == "historical_positive":
        return cho_labels(system)[0]
    if case["name"] != "amide":
        return urethane_labels(system)
    # Independently specified neutral acetamide control, never a runtime typer.
    atoms = system.topology.sites
    adj = {i: set() for i in atoms}
    orders = {i: [] for i in atoms}
    for b in system.topology.bonds.values():
        for i, j in ((b.site1, b.site2), (b.site2, b.site1)):
            adj[i].add(j)
            orders[i].append(b.order)
    labels = {}
    for i, a in atoms.items():
        if (
            a.formal_charge
            or a.metadata.get("isotope")
            or a.metadata.get("radical_electrons")
        ):
            raise ValueError("Outside independent amide control")
        pattern = (a.element, tuple(sorted(orders[i])))
        if pattern == ("C", (1, 1, 1, 1)):
            labels[i] = "c3"
        elif pattern == ("C", (1, 1, 2)):
            labels[i] = "c_1"
        elif pattern == ("O", (2,)):
            labels[i] = "o_1"
        elif pattern == ("N", (1, 1, 1)):
            labels[i] = "n"
        elif pattern == ("H", (1,)):
            labels[i] = "hn" if atoms[next(iter(adj[i]))].element == "N" else "hc"
        else:
            raise ValueError("Wrong independent amide environment")
    return labels


def missing_request(raw, code, text, *, system=None, reference_labels=None):
    import re

    matches = re.findall(r"Unable to find (\w+) data for ([^\n]+)", text)
    if len(matches) != 1 or matches[0][0] not in MISSING:
        raise ValueError("Unclassified converter failure; not a chemical rejection")
    kind, pattern = matches[0]
    exits, family = MISSING[kind]
    if code not in (exits if isinstance(exits, tuple) else (exits,)):
        raise ValueError("Converter failure exit/stage mismatch")
    labels = pattern.split()
    if system is not None:
        from itertools import permutations

        from pcff_j5_reference import raw_inventory

        possible = set()
        for f, sites in raw_inventory(system):
            if f != family:
                continue
            values = tuple(reference_labels[i] for i in sites)
            if f in {"angle-angle", "wilson_out_of_plane"}:
                possible.update(
                    (p[0], values[1], p[1], p[2])
                    for p in permutations((values[0], values[2], values[3]))
                )
            else:
                possible.update((values, values[::-1]))
        if tuple(labels) not in possible:
            raise ValueError(
                "Converter missing request does not belong to final graph roles"
            )
    independent = resolve_raw(raw, family, labels, guarded_msi=True)
    return {
        "family": family,
        "types": labels,
        "independent_search": independent,
        "classification": "source_parameter_missing"
        if independent["status"] == "missing"
        else "wildcard_conflict"
        if independent["status"] == "ambiguous"
        else "automatic_fallback_not_implemented_by_converter"
        if independent.get("family") != family
        else "converter_lookup_or_permutation_difference",
    }
