"""Versioned source-resolved OPLS parameters; historical 4G1 records unchanged."""

import json
from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations
from math import degrees, isfinite
from pathlib import Path
from xml.etree import ElementTree as ET

from island.charge_references.records import pack, unpack
from island.forcefields.parameterized import ParameterizedSystem

from .adapter import type_atoms
from .models import OPLSChargeResult, assign_native_charges, chemical_graph, identity
from .source import FLAGS, boundary, require

VERIFICATION = {
    "policy": "force_inventory_and_numerical_content_v1",
    **json.loads(Path(__file__).with_name("resolution_pin.json").read_text()),
}

SECTIONS = {
    "bonds": ("HarmonicBondForce", "Bond", ("length", "k")),
    "angles": ("HarmonicAngleForce", "Angle", ("angle", "k")),
    "rb": ("RBTorsionForce", "Proper", tuple(f"c{i}" for i in range(6))),
}
UNITS = {
    "bonds": {"length": "nm", "k": "kJ/(mol*nm^2)"},
    "angles": {"angle": "radian", "k": "kJ/(mol*radian^2)"},
    "rb": {f"c{i}": "kJ/mol" for i in range(6)},
}
POLICY = {
    "mixing": "geometric_sigma_and_epsilon",
    "lj": "12-6",
    "coulomb_14": 0.5,
    "lj_14": 0.5,
    "excluded_graph_distances": [1, 2],
    "scaled_graph_distance": 3,
    "method": "NoCutoff",
    "impropers": "no_explicit_source_section",
}


def inventory(system):
    """Regenerate unique simple paths and shortest-distance pairs from bonds."""
    ids = sorted(system.topology.sites)
    neighbours = {i: set() for i in ids}
    for a, b in system.topology.bonds:
        neighbours[a].add(b)
        neighbours[b].add(a)
    angles = sorted(
        (a, b, c) for b in ids for a, c in combinations(sorted(neighbours[b]), 2)
    )
    propers = set()
    for b, c in system.topology.bonds:
        for a in neighbours[b] - {c}:
            for d in neighbours[c] - {a, b}:
                v = (a, b, c, d)
                propers.add(min(v, v[::-1]))
    pairs = []
    for a in ids:
        seen = {a}
        frontier = {a}
        for distance in (1, 2, 3):
            frontier = (
                set().union(*(neighbours[i] for i in frontier)) - seen
                if frontier
                else set()
            )
            seen |= frontier
            for b in sorted(frontier):
                if a < b:
                    pairs.append(
                        {
                            "sites": (a, b),
                            "distance": distance,
                            "lj_scale": 0.5 if distance == 3 else 0.0,
                            "coulomb_scale": 0.5 if distance == 3 else 0.0,
                        }
                    )
    return {
        "bonds": sorted(system.topology.bonds),
        "angles": angles,
        "rb": sorted(propers),
    }, sorted(pairs, key=lambda r: r["sites"])


def selected_record(root, types, names, kind):
    """Pinned OpenMM matching: reversible, first angle/specific RB, first wildcard RB.

    Bond candidate indices are a set upstream. Reject differing simultaneous
    bond matches rather than relying on Python set order. Numerically identical
    duplicate rows retain all candidate IDs and a deterministic primary ID.
    """
    section, tag, values = SECTIONS[kind]
    matches = []
    classes = {}
    for t, a in types.items():
        classes.setdefault(a["class"], set()).add(t)
    for index, row in enumerate(root.find(section)):
        require(row.tag == tag, f"Unsupported {section} entry {row.tag}")
        attrs = row.attrib
        groups = []
        wildcard = False
        for n in range(1, len(names) + 1):
            key = f"type{n}" if f"type{n}" in attrs else f"class{n}"
            require(key in attrs, "Malformed source matching fields")
            value = attrs[key]
            wildcard |= value == ""
            groups.append(
                set(types)
                if value == ""
                else {value}
                if key.startswith("type")
                else classes.get(value, set())
            )
        if any(
            all(t in g for t, g in zip(order, groups, strict=True))
            for order in (names, names[::-1])
        ):
            matches.append((index, attrs, wildcard))
    require(bool(matches), f"Missing {section} parameters for {names}")
    if kind == "bonds":
        require(
            len({tuple(float(r[k]) for k in values) for _, r, _ in matches}) == 1,
            f"Ambiguous differing bond selections for {names}",
        )
        match = matches[0]
    elif kind == "angles":
        match = matches[0]
    else:
        match = next((r for r in matches if not r[2]), matches[0])
    index, attrs, _ = match
    return index, dict(attrs), [r[0] for r in matches]


def resolved_data(system, source, charges):
    source.validate_integrity()
    charges.validate_integrity(system, source)
    root = ET.fromstring(source.xml)
    require(
        [e.tag for e in root]
        == [
            "AtomTypes",
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "RBTorsionForce",
            "NonbondedForce",
        ],
        "Unsupported source forms",
    )
    require(
        root.attrib["combining_rule"] == "geometric"
        and root.find("NonbondedForce").attrib
        == {"coulomb14scale": "0.5", "lj14scale": "0.5"},
        "Unsupported nonbonded policy",
    )
    types, _ = source.entries
    assigned = charges.payload["typing"]["data"]["assignments"]
    paths, pairs = inventory(system)
    result = {}
    selections = {}
    native_rows = charges.payload["charges"]
    for kind, paths_for_kind in paths.items():
        section, tag, values = SECTIONS[kind]
        rows = []
        for sites in paths_for_kind:
            names = tuple(assigned[s] for s in sites)
            key = (kind, names)
            if key not in selections:
                selections[key] = selected_record(root, types, names, kind)
            index, attrs, candidates = selections[key]
            raw = {k: float(attrs[k]) for k in values}
            require(
                all(isfinite(v) for v in raw.values()), "Nonfinite source parameter"
            )
            if kind == "bonds":
                require(raw["k"] >= 0 and raw["length"] > 0, "Invalid bond parameter")
                converted = {
                    "length_angstrom": raw["length"] * 10,
                    "k_kj_mol_angstrom2": raw["k"] / 100,
                }
            elif kind == "angles":
                require(raw["k"] >= 0, "Invalid angle parameter")
                converted = {
                    "angle_degrees": degrees(raw["angle"]),
                    "k_kj_mol_radian2": raw["k"],
                }
            else:
                converted = {
                    "coefficients_kj_mol": [raw[f"c{i}"] for i in range(6)],
                    "convention": "sum_c_n_cos(theta_minus_pi)^n",
                }
            rows.append(
                {
                    "sites": sites,
                    "types": names,
                    "classes": tuple(types[t]["class"] for t in names),
                    "source_record": f"{section}/{tag}[{index}]",
                    "matching_record_indices": candidates,
                    "source_attributes": attrs,
                    "original": raw,
                    "original_units": UNITS[kind],
                    "converted": converted,
                }
            )
        result[kind] = rows
    nb = {
        r.attrib["type"]: (i, r.attrib)
        for i, r in enumerate(root.find("NonbondedForce"))
    }
    result["sites"] = {}
    for sid, t in assigned.items():
        require(t in nb, f"Missing LJ record {t}")
        index, attrs = nb[t]
        q, sigma, epsilon = (float(attrs[k]) for k in ("charge", "sigma", "epsilon"))
        require(
            all(isfinite(v) for v in (q, sigma, epsilon))
            and sigma >= 0
            and epsilon >= 0,
            "Invalid nonbonded values",
        )
        require(q == native_rows[sid]["charge"], "Native charge mismatch")
        result["sites"][sid] = {
            "type": t,
            "class": types[t]["class"],
            "source_record": f"NonbondedForce/Atom[{index}]",
            "source_attributes": dict(attrs),
            "original_units": {
                "charge": "elementary_charge",
                "sigma": "nm",
                "epsilon": "kJ/mol",
            },
            "mass_dalton": system.topology.sites[sid].mass,
            "charge_e": q,
            "sigma_angstrom": sigma * 10,
            "epsilon_kj_mol": epsilon,
        }
    result.update(
        pairs=pairs,
        policy=POLICY,
        impropers=[],
        coverage={k: len(v) for k, v in paths.items()},
        source=source.identity,
        graph_identity=identity(chemical_graph(system)),
        charge_identity=charges.identity,
        typing_identity=charges.payload["typing_identity"],
        conversions={
            "bond": "r_A=10*r_nm; k_A=k_nm/100; harmonic_one_half",
            "angle": "degrees=180*rad/pi; k unchanged; harmonic_one_half",
            "rb": "six coefficients unchanged, including c0",
        },
        **FLAGS,
    )
    return result


@dataclass(frozen=True)
class OPLSParameterizationResult:
    json_text: str

    @property
    @boundary
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @boundary
    def validate_integrity(self, system, source):
        p = self.payload
        require(
            set(p) == {"schema", "charges", "resolved", "upstream_verification"},
            "Malformed parameter result",
        )
        require(
            p["schema"] == "island_opls_parameters_v1", "Unsupported parameter schema"
        )
        charges = OPLSChargeResult(pack(p["charges"]))
        require(
            pack(p["resolved"]) == pack(resolved_data(system, source, charges)),
            "Resolved parameters contradict source or graph",
        )
        require(
            p["upstream_verification"] == VERIFICATION,
            "Missing upstream verification declaration",
        )

    @boundary
    def to_parameterized_system(self, system, source):
        self.validate_integrity(system, source)
        p = self.payload
        from island.forcefields.charges.models import ChargeAssignment
        from island.forcefields.nonbonded.models import NonbondedPolicy

        return ParameterizedSystem(
            deepcopy(system),
            backend_name="Foyer OPLS-AA",
            site_assignments=deepcopy(p["resolved"]["sites"]),
            interaction_assignments={
                k: {r["sites"]: deepcopy(r) for r in p["resolved"][k]} for k in SECTIONS
            },
            metadata={"opls_parameters_v1": p, **FLAGS},
            charge_assignments={
                sid: ChargeAssignment(
                    sid,
                    r["charge_e"],
                    "source_native",
                    identity(source.identity),
                    entry_id=r["source_record"],
                    atom_type=r["type"],
                )
                for sid, r in p["resolved"]["sites"].items()
            },
            nonbonded_policy=NonbondedPolicy(
                "Foyer OPLS-AA",
                "1",
                "geometric",
                0.0,
                0.0,
                0.0,
                0.0,
                0.5,
                0.5,
                identity(source.identity),
            ),
            aggregate_signature=self.identity,
        )


@boundary
def parameterize_oplsaa(system, source):
    charges = assign_native_charges(system, type_atoms(system, source), source)
    resolved = resolved_data(system, source, charges)
    from .upstream import verify_resolved

    verify_resolved(
        system, source, charges.payload["typing"]["data"]["assignments"], resolved
    )
    result = OPLSParameterizationResult(
        pack(
            {
                "schema": "island_opls_parameters_v1",
                "charges": charges.payload,
                "resolved": resolved,
                "upstream_verification": VERIFICATION,
            }
        )
    )
    result.validate_integrity(system, source)
    return result
