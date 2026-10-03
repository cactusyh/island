"""Owned data-only typing/native-charge records; not parameterized systems."""

from dataclasses import dataclass
from math import fsum, isfinite

from island import AtomSite
from island.charge_references.records import pack, unpack
from island.forcefields.oplsaa.source import (
    CAPABILITIES,
    FLAGS,
    PIN,
    boundary,
    digest,
    require,
)


def identity(payload):
    return digest(pack(payload).encode())


@boundary
def chemical_graph(system):
    require(
        system.representation == "atomistic" and system.box is None,
        "Only nonperiodic atomistic systems supported",
    )
    topology = system.topology
    topology.validate_bond_graph()
    require(bool(topology.sites), "Empty system")
    sites = []
    for key in sorted(topology.sites):
        atom = topology.sites[key]
        require(
            type(key) is int and isinstance(atom, AtomSite) and atom.id == key,
            "Require integer stable atom IDs",
        )
        require(
            type(atom.formal_charge) is int
            and atom.formal_charge == 0
            and not atom.metadata.get("isotope")
            and not atom.metadata.get("radical_electrons", 0),
            "Charged, isotopic or radical sites unsupported in this adapter",
        )
        require(
            type(atom.mass) in (float, int) and isfinite(atom.mass) and atom.mass > 0,
            "Invalid atomic mass",
        )
        numbers = {
            "H": 1,
            "C": 6,
            "N": 7,
            "O": 8,
            "F": 9,
            "S": 16,
            "Cl": 17,
            "Br": 35,
            "I": 53,
        }
        require(
            type(atom.atomic_number) is int
            and atom.atomic_number == numbers.get(atom.element),
            "Unsupported or inconsistent element/atomic number",
        )
        bonds = [b for b in topology.bonds.values() if key in (b.site1, b.site2)]
        valence = fsum(b.order for b in bonds)
        allowed = {
            "H": (1,),
            "C": (4,),
            "N": (3,),
            "O": (2,),
            "F": (1,),
            "Cl": (1,),
            "Br": (1,),
            "I": (1,),
            "S": (2, 4, 6),
        }
        if atom.metadata.get("aromatic"):
            valid = (
                atom.element == "C"
                and valence in (4, 4.5)
                or atom.element == "N"
                and valence in (3, 4, 4.5)
                or atom.element in ("O", "S")
                and len(bonds) == 2
                and all(b.aromatic and b.order == 1.5 for b in bonds)
            )
        else:
            valid = valence in allowed.get(atom.element, ())
        require(valid, f"Invalid valence or missing explicit H at site {key}")
        sites.append(
            {
                "id": key,
                "element": atom.element,
                "atomic_number": atom.atomic_number,
                "mass": atom.mass,
                "formal_charge": atom.formal_charge,
                "chemical_metadata": {
                    k: atom.metadata.get(k)
                    for k in (
                        "aromatic",
                        "chiral_tag",
                        "cip_label",
                        "isotope",
                        "radical_electrons",
                        "explicit_hydrogen_count",
                        "no_implicit_hydrogens",
                    )
                },
            }
        )
    bonds = []
    for pair, b in sorted(topology.bonds.items()):
        require(
            b.order in (1, 1.5, 2, 3) and (not b.aromatic or b.order == 1.5),
            "Unsupported bond semantics",
        )
        bonds.append(
            {"sites": tuple(sorted(pair)), "order": b.order, "aromatic": b.aromatic}
        )
    return {
        "sites": sites,
        "bonds": bonds,
        "components": [
            sorted(c) for c in sorted(topology.connected_components(), key=min)
        ],
    }


def typing_data(p, system, source):
    require(
        set(p)
        == {
            "schema",
            "source",
            "graph",
            "index_to_site_id",
            "matches",
            "environment",
            "data",
        },
        "Invalid typing fields",
    )
    require(p["schema"] == "island_foyer_typing_v1", "Unsupported typing schema")
    require(p["source"] == source.identity, "Source identity mismatch")
    graph = chemical_graph(system)
    require(pack(p["graph"]) == pack(graph), "Chemical graph/stereo identity mismatch")
    ids = [a["id"] for a in graph["sites"]]
    require(
        pack(p["index_to_site_id"]) == pack(dict(enumerate(ids))),
        "Index lineage mismatch",
    )
    require(
        set(p["matches"]) == set(ids) and all(type(k) is int for k in p["matches"]),
        "Typing coverage mismatch",
    )
    types, _ = source.entries
    assignments = {}
    for sid, row in p["matches"].items():
        require(
            set(row) == {"whitelist", "blacklist", "atomtype"}, "Invalid Foyer match"
        )
        for key in ("whitelist", "blacklist"):
            require(
                type(row[key]) is list
                and row[key] == sorted(set(row[key]))
                and all(t in types for t in row[key]),
                "Invalid matched types",
            )
        selected = set(row["whitelist"]) - set(row["blacklist"])
        require(
            len(selected) == 1 and row["atomtype"] in selected,
            f"Untyped/ambiguous site {sid}: {sorted(selected)}",
        )
        entry = types[row["atomtype"]]
        require(
            entry.get("def") and entry["element"] == system.topology.sites[sid].element,
            "Selected type lacks definition or has wrong element",
        )
        expected_black = set()
        for t in row["whitelist"]:
            expected_black.update(
                x.strip() for x in types[t].get("overrides", "").split(",") if x.strip()
            )
        require(
            set(row["blacklist"]) == expected_black, "Override evidence inconsistent"
        )
        assignments[sid] = row["atomtype"]
    env = p["environment"]
    require(
        type(env) is dict
        and env.get("foyer") == PIN["foyer_version"]
        and all(type(k) is str and type(v) is str and v for k, v in env.items()),
        "Invalid dependency provenance",
    )
    return {
        "assignments": assignments,
        "graph_identity": identity(graph),
        "capabilities": CAPABILITIES,
        **FLAGS,
    }


@dataclass(frozen=True)
class OPLSTypingResult:
    json_text: str

    @boundary
    def validate_integrity(self, system, source):
        p = unpack(self.json_text)
        require(
            pack(p["data"]) == pack(typing_data(p, system, source)),
            "Typing record inconsistent",
        )

    @property
    @boundary
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


def charge_data(p, system, source):
    require(
        set(p) == {"schema", "typing", "typing_identity", "charges", "data"},
        "Invalid charge fields",
    )
    require(p["schema"] == "island_opls_native_charges_v1", "Unsupported charge schema")
    typing = OPLSTypingResult(pack(p["typing"]))
    typing.validate_integrity(system, source)
    require(p["typing_identity"] == typing.identity, "Typing identity mismatch")
    assignments = typing.payload["data"]["assignments"]
    require(
        set(p["charges"]) == set(assignments)
        and all(type(k) is int for k in p["charges"]),
        "Charge coverage mismatch",
    )
    _, native = source.entries
    for sid, row in p["charges"].items():
        t = assignments[sid]
        require(t in native, f"Missing native charge for {t}")
        require(
            set(row) == {"type", "charge", "location", "unit"}
            and row["type"] == t
            and row["unit"] == "elementary_charge"
            and row["location"] == f'NonbondedForce/Atom[@type="{t}"]'
            and type(row["charge"]) in (float, int)
            and isfinite(row["charge"])
            and row["charge"] == native[t],
            "Native charge/source mismatch",
        )
    totals = []
    for component in p["typing"]["graph"]["components"]:
        total = fsum(p["charges"][s]["charge"] for s in component)
        # These are exact decimal library entries, not noisy QM output. A
        # floating summation allowance scales with inventory, not rounding bins.
        tolerance = len(component) * 1e-12
        require(
            abs(total) <= tolerance, f"Native component neutrality failed: {total} e"
        )
        totals.append({"sites": component, "total_e": total, "tolerance_e": tolerance})
    total = fsum(r["charge"] for r in p["charges"].values())
    require(
        abs(total) <= len(assignments) * 1e-12, "Native molecular neutrality failed"
    )
    return {
        "components": totals,
        "total_e": total,
        "unit": "elementary_charge",
        "molecular_tolerance_e": len(assignments) * 1e-12,
        "source": source.identity,
        "capabilities": CAPABILITIES,
        **FLAGS,
    }


@dataclass(frozen=True)
class OPLSChargeResult:
    json_text: str

    @boundary
    def validate_integrity(self, system, source):
        p = unpack(self.json_text)
        require(
            pack(p["data"]) == pack(charge_data(p, system, source)),
            "Native charge record inconsistent",
        )

    @property
    @boundary
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def assign_native_charges(system, typing, source):
    require(isinstance(typing, OPLSTypingResult), "Expected Foyer typing result")
    typing.validate_integrity(system, source)
    _, native = source.entries
    rows = {}
    for sid, t in typing.payload["data"]["assignments"].items():
        require(t in native, f"Missing native charge for {t}")
        rows[sid] = {
            "type": t,
            "charge": native[t],
            "unit": "elementary_charge",
            "location": f'NonbondedForce/Atom[@type="{t}"]',
        }
    p = {
        "schema": "island_opls_native_charges_v1",
        "typing": typing.payload,
        "typing_identity": typing.identity,
        "charges": rows,
        "data": {},
    }
    p["data"] = charge_data(p, system, source)
    return OPLSChargeResult(pack(p))


@boundary
def save_opls_result(result, system, source, path):
    """Validate before exclusive data-only publication; not an MD restart."""
    from pathlib import Path

    from island.workflows.storage import publish

    require(
        isinstance(result, (OPLSTypingResult, OPLSChargeResult)),
        "Invalid OPLS result type",
    )
    result.validate_integrity(system, source)
    publish(Path(path), result.json_text.encode())


@boundary
def load_opls_result(path, system, source):
    """Offline validation requires the authoritative graph and pinned XML."""
    from pathlib import Path

    text = Path(path).read_text()
    schema = unpack(text)["schema"]
    require(
        schema in ("island_foyer_typing_v1", "island_opls_native_charges_v1"),
        "Unsupported OPLS result schema",
    )
    cls = OPLSTypingResult if schema == "island_foyer_typing_v1" else OPLSChargeResult
    result = cls(text)
    result.validate_integrity(system, source)
    return result
