"""Independent declared aliphatic C/H/N reference; never a runtime typer.

FRC ref1 descriptors at 99/100/105/118/123/130 and equivalences 247/265.
The reference does not import ISLAND's chemical predicates or source resolver.
It only accepts saturated C/H/N graphs declared in the J14 experiment.
"""

from collections import Counter

from pcff_j13_reference import require
from pcff_j13_reference import selection_check as check_rows


def labels_for(system):
    atoms = system.topology.sites
    neighbors = {i: [] for i in atoms}
    for b in system.topology.bonds.values():
        require(b.order == 1, "Reference requires saturated graph")
        neighbors[b.site1].append(b.site2)
        neighbors[b.site2].append(b.site1)
    labels = {}
    for i, atom in atoms.items():
        require(
            not atom.metadata.get("isotope")
            and not atom.metadata.get("radical_electrons")
            and not atom.metadata.get("aromatic"),
            "Reference chemical state mismatch",
        )
        env = Counter(atoms[j].element for j in neighbors[i])
        if atom.element == "C":
            require(
                atom.formal_charge == 0
                and sum(env.values()) == 4
                and set(env) <= {"C", "H", "N"},
                "Reference carbon mismatch",
            )
            labels[i] = "c" + str(env["H"]) if env["H"] else "c"
        elif atom.element == "N":
            require(
                atom.formal_charge in (0, 1)
                and sum(env.values()) == 3 + atom.formal_charge
                and set(env) <= {"C", "H"}
                and env["C"] >= 1,
                "Reference nitrogen mismatch",
            )
            labels[i] = "n4" if atom.formal_charge else "na"
        elif atom.element == "H":
            require(
                atom.formal_charge == 0 and len(neighbors[i]) == 1,
                "Reference hydrogen mismatch",
            )
            parent = atoms[neighbors[i][0]]
            require(parent.element in ("C", "N"), "Reference hydrogen parent mismatch")
            labels[i] = (
                "hc"
                if parent.element == "C"
                else "h+"
                if parent.formal_charge
                else "hn"
            )
        else:
            raise ValueError("Outside independently declared C/H/N reference")
    return labels


def selection_check(system, model, source):
    return check_rows(system, model, source, reference_labels=labels_for(system))


def aa_raw_reference_coefficients(sections, labels, source):
    """Separate corrected oracle: preserve B and the shared arm in all K lookups.

    GetParameters.c find_match reverses ABCD to DCBA, moving the central atom.
    Its find_angleangle_data calls can therefore select a different-center row.
    Raw converter output is retained. Rebuild K and theta from independent raw
    rows and converter inventories; no native model/assignment is accepted here.
    LAMMPS improper_class2.cpp uses K1 ABC*CBD, K2 ABC*ABD, K3 ABD*CBD.
    """
    from pcff_j5_reference import read_source, resolve_raw

    raw = read_source(source.read_bytes())
    angles = {
        min(tuple(r[2:]), tuple(reversed(r[2:]))): r[1] for r in sections["Angles"]
    }
    angle_coeff = {r[0]: float(r[1]) for r in sections["Angle Coeffs"]}
    old = {r[0]: list(map(float, r[1:])) for r in sections["AngleAngle Coeffs"]}
    type_labels = {
        str(index): labels[sid] for index, sid in enumerate(sorted(labels), 1)
    }
    values, evidence = {}, []
    for row in sections["Impropers"]:
        a, b, c, d = row[2:]
        permutations = ((a, b, c, d), (d, b, a, c), (a, b, d, c))
        matches = [
            resolve_raw(
                raw, "angle-angle", [type_labels[i] for i in p], guarded_msi=True
            )
            for p in permutations
        ]
        require(
            all(m["status"] == "assigned" for m in matches),
            "Raw AA source coverage missing/ambiguous",
        )
        theta = [
            angle_coeff[angles[min(p, p[::-1])]]
            for p in ((a, b, c), (a, b, d), (c, b, d))
        ]
        coefficients = [m["values"][0] / 4.184 for m in matches] + theta
        require(
            row[1] not in values or values[row[1]] == coefficients,
            "Conflicting reference AA type roles",
        )
        values[row[1]] = coefficients
        evidence.append(
            {
                "improper": row,
                "types": [type_labels[i] for i in row[2:]],
                "raw_converter": old[row[1]],
                "corrected": coefficients,
                "role_patterns": permutations,
                "source_matches": matches,
            }
        )
    commands = [
        "improper_coeff " + k + " aa " + " ".join(format(v, ".17g") for v in v)
        for k, v in values.items()
    ]
    return commands, evidence
