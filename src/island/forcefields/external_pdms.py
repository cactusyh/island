"""External PDMS Class-II data contract, separate from native PCFF authority.

No parameter values or execution profiles are shipped here. A checksum proves
consistency, not scientific authority. Even a complete import is a data record;
it cannot grant itself a backend qualification or bypass periodic gates.
"""

from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.exceptions import IslandError
from island.workflows import storage

from .pcff.automatic import chemical_graph
from .pcff.class2 import FAMILIES, interaction_requests

SCHEMA = "island_external_pdms_class2_bundle_v1"
PROFILE = "island_pcff_external_pdms_class2_v1"
ROWS_SCHEMA = "island_external_class2_source_rows_v1"
UNITS = {
    "energy": "kcal/mol",
    "length": "angstrom",
    "angle_equilibrium": "degree",
    "angle_displacement": "radian",
    "mass": "dalton",
    "charge": "elementary_charge",
    "conversion": "identity_source_units_v1; kcal_to_kJ=4.184; degree_to_rad=pi/180",
}
SPECIAL = {"lj": [0, 0, 1], "coulomb": [0, 0, 1]}
MIXING = "class2_sixth_power_rmin_epsilon_v1"
LIMITS = {
    "production_validated": False,
    "simulation_readiness": "not_established",
    "ewald": "ewald_unvalidated",
    "pressure": "pressure_unavailable",
    "cutoff_recovery": "blocked",
    "long_npt": "blocked",
}
AUTHORITY = {
    "classification": "external_parameter_bundle_candidate",
    "native_pcff_coverage": False,
    "msi2lmp_parity": False,
    "execution_authorized": False,
}


def digest(value):
    return sha256(storage.json_bytes(value)).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(value):
    """Preserve source decimals; floats, NaN and implicit charge edits are refused."""
    require(type(value) is str and value.strip() == value, "decimal_string_required")
    try:
        result = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("invalid_decimal") from exc
    require(result.is_finite(), "nonfinite_decimal")
    return result


def pdms_final_graph(system):
    """Explicit-H neutral CHO/Si graph, including molecule and component identity.

    The historical generic final-graph validator does not support Si. This
    separately versioned contract leaves that validator and its hashes intact.
    Coordinates are not chemical identity and remain in the system record.
    """
    system.validate()
    graph = chemical_graph(system)
    require(graph["representation"] == "atomistic", "atomistic_graph_required")
    adjacency = {a["id"]: set() for a in graph["sites"]}
    for b in graph["bonds"]:
        require(b["order"] == 1 and not b["aromatic"], "single_bonds_required")
        a, z = b["sites"]
        adjacency[a].add(z)
        adjacency[z].add(a)
    numbers = {"H": 1, "C": 6, "O": 8, "Si": 14}
    valences = {"H": 1, "C": 4, "O": 2, "Si": 4}
    elements = {a["element"] for a in graph["sites"]}
    require({"Si", "O"} <= elements <= set(numbers), "outside_pdms_domain")
    element_of = {a["id"]: a["element"] for a in graph["sites"]}
    require(
        any({element_of[i] for i in b["sites"]} == {"Si", "O"} for b in graph["bonds"]),
        "silicon_oxygen_bond_required",
    )
    for a in graph["sites"]:
        e, i = a["element"], a["id"]
        require(a["atomic_number"] == numbers[e], f"wrong_element:{i}")
        require(len(adjacency[i]) == valences[e], f"wrong_valence_or_missing_H:{i}")
        require(
            a["formal_charge"] == 0
            and not any(
                a["metadata"][k] for k in ("isotope", "radical_electrons", "aromatic")
            ),
            f"unsupported_chemical_state:{i}",
        )
    components = sorted(
        [sorted(c) for c in system.topology.connected_components()], key=lambda c: c[0]
    )
    graph.update(
        schema="island_external_pdms_final_graph_v1",
        components=components,
        component_membership={str(i): n for n, c in enumerate(components) for i in c},
        molecule_membership={
            str(i): str(system.topology.sites[i].metadata.get("chain_id", n))
            for n, c in enumerate(components)
            for i in c
        },
        polymer_provenance=storage.encode(system.metadata.get("polymer", {})),
    )
    return graph


def interaction_key(family, sites):
    return family + ":" + ",".join(map(str, sites))


def requirements(system):
    """Enumerate every applicable coefficient and equilibrium dependency."""
    graph = pdms_final_graph(system)
    requests, _ = interaction_requests(graph, expanded=True)
    result, inapplicable = {}, []
    for family, sites, dependencies, reason in requests:
        key = interaction_key(family, sites)
        if reason:
            inapplicable.append({"key": key, "reason": reason})
            continue
        result[key] = {
            "family": family,
            "sites": list(sites),
            "dependencies": [
                # Bond/angle reverse orientations share the same equilibrium.
                interaction_key(f, min(s, s[::-1]))
                for f, s in dependencies
            ],
        }
    return {
        "graph": graph,
        "graph_identity": digest(graph),
        "atom_order": [a["id"] for a in graph["sites"]],
        "interactions": result,
        "inapplicable": inapplicable,
    }


def bundle_template(system):
    """An explicitly incomplete input template: no charges or coefficients."""
    r = requirements(system)
    return {
        "schema": SCHEMA,
        "profile": PROFILE,
        "profile_version": 1,
        "graph": r["graph"],
        "graph_identity": r["graph_identity"],
        "atom_order": r["atom_order"],
        "atoms": [
            {
                "id": a["id"],
                "element": a["element"],
                "mass": str(a["mass"]),
                "type": None,
                "evidence_row": None,
            }
            for a in r["graph"]["sites"]
        ],
        "source_family": None,
        "primary_source": None,
        "files": {},
        "source_rows": {},
        "provenance": None,
        "evidence_references": [],
        "executable_provenance": [],
        "charge_origin": None,
        "charges": [],
        "charge_increments": {},
        "charge_transformations": [],
        "interactions": {},
        "coverage": r["interactions"],
        "inapplicable": r["inapplicable"],
        "nonbonded_policy": {
            "mixing": MIXING,
            "special_bonds": deepcopy(SPECIAL),
            "evidence_row": None,
        },
        "units": deepcopy(UNITS),
        "authority": deepcopy(AUTHORITY),
        "limitations": deepcopy(LIMITS),
    }


def source_row(row, document, assets):
    """Reconstruct one row from its independent, hash-checked source document.

    External suppliers provide a versioned JSON row table plus raw evidence.
    No FRC fallback, wildcard, implicit zero, or unit inference is performed.
    Rows must include raw-line citations whose numeric columns reproduce every
    supplied coefficient, including explicitly evidenced zeros.
    """
    require(
        type(row) is dict
        and set(row)
        == {"family", "types", "values", "units", "provenance", "citations"},
        "malformed_source_row",
    )
    require(
        row["provenance"] and type(row["provenance"]) is str, "row_provenance_missing"
    )
    require(
        type(row["types"]) is list and all(type(t) is str and t for t in row["types"]),
        "row_types_invalid",
    )
    require(not any("*" in t for t in row["types"]), "wildcard_not_authorized")
    require(type(row["values"]) is list and row["values"], "row_values_missing")
    values = [number(v) for v in row["values"]]
    require(
        type(row["units"]) is list and len(row["units"]) == len(values),
        "row_units_missing",
    )
    if row["family"] in FAMILIES:
        require(
            row["units"] == FAMILIES[row["family"]][2]
            and len(values) == len(row["units"]),
            "source_row_units_or_shape_invalid",
        )
    require(row in document["rows"].values(), "row_not_in_primary_source")
    require(type(row["citations"]) is list and row["citations"], "raw_evidence_missing")
    reconstructed = []
    for citation in row["citations"]:
        require(
            set(citation) == {"file", "line", "raw", "columns"}, "malformed_citation"
        )
        raw = assets[citation["file"]].decode("utf-8").splitlines()
        n = citation["line"]
        require(type(n) is int and 1 <= n <= len(raw), "source_line_invalid")
        require(raw[n - 1] == citation["raw"], "source_line_changed")
        columns = citation["columns"]
        require(
            type(columns) is list and all(type(c) is int and c >= 0 for c in columns),
            "source_columns_invalid",
        )
        words = raw[n - 1].split("!", 1)[0].split()
        reconstructed.extend(number(words[c]) for c in columns)
    require(reconstructed == values, "source_coefficient_changed_or_inserted_zero")
    return values


@dataclass(frozen=True)
class ExternalBundleValidation:
    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def complete(self):
        return self.payload["status"] == "parameter_bundle_validated"


@dataclass(frozen=True)
class ExternalParameterBundle:
    """Owned immutable candidate, never a PCFF model or execution authority."""

    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return sha256(self.json_text.encode()).hexdigest()

    def validate(self, system, *, assets, expected_sha256):
        return inspect_external_bundle(
            self, system, assets=assets, expected_sha256=expected_sha256
        )


def inspect_external_bundle(bundle, system, *, assets, expected_sha256):
    """Collect precise blocked stages, including every missing interaction.

    expected_sha256 is an out-of-band user/reviewer pin, never read from the
    candidate itself. Files are relocated by relative name; their bytes are
    checked on every validation. A complete record establishes data integrity
    only. C5 has no qualified external execution adapter.
    """
    errors = []

    def check(ok, stage, context=None):
        if not ok:
            errors.append({"stage": stage, "context": context})
        return bool(ok)

    try:
        p = bundle.payload
        check(bundle.identity == expected_sha256, "bundle_pin_mismatch")
        template = bundle_template(system)
        check(set(p) == set(template), "bundle_fields_invalid")
        for key in (
            "schema",
            "profile",
            "profile_version",
            "graph",
            "graph_identity",
            "atom_order",
            "coverage",
            "inapplicable",
            "units",
            "authority",
            "limitations",
        ):
            check(p.get(key) == template[key], key + "_mismatch")
        family = p.get("source_family")
        check(
            type(family) is str and family.startswith("external:") and len(family) > 9,
            "external_source_family_required",
        )
        check(
            type(p.get("provenance")) is str and bool(p["provenance"]),
            "provenance_missing",
        )
        check(p.get("charge_transformations") == [], "charge_redistribution_forbidden")
        files, raw_assets = p.get("files", {}), {}
        for name, record in files.items():
            check(
                set(record) == {"sha256", "family", "role", "provenance"},
                "file_fields_invalid",
                name,
            )
            check(record.get("family") == family, "mixed_source_family", name)
            check(bool(record.get("provenance")), "file_provenance_missing", name)
            path = storage.child(assets, name)
            if not check(path.is_file(), "asset_unavailable", name):
                continue
            raw = path.read_bytes()
            if check(
                sha256(raw).hexdigest() == record.get("sha256"),
                "source_hash_mismatch",
                name,
            ):
                raw_assets[name] = raw
        check(bool(files), "source_dependencies_missing")
        primary = p.get("primary_source")
        document = None
        if check(primary in raw_assets, "primary_source_unavailable", primary):
            from island.dynamics._checkpoint_data import strict_load

            document = strict_load(raw_assets[primary].decode())
            check(
                set(document) == {"schema", "source_family", "rows", "provenance"},
                "primary_source_fields_invalid",
            )
            check(
                document.get("schema") == ROWS_SCHEMA
                and document.get("source_family") == family,
                "primary_source_identity_mismatch",
            )
            check(document.get("rows") == p.get("source_rows"), "source_rows_mismatch")
        check(
            bool(p.get("evidence_references"))
            and all(n in raw_assets for n in p.get("evidence_references", [])),
            "evidence_unavailable",
        )
        executables = p.get("executable_provenance", [])
        check(bool(executables), "executable_provenance_missing")
        for e in executables:
            check(
                set(e) == {"file", "version", "build", "revision", "role"}
                and e.get("file") in raw_assets
                and all(type(v) is str and v for v in e.values()),
                "executable_provenance_invalid",
            )
        rows = p.get("source_rows", {})
        valid_rows = {}
        if document:
            for rid, row in rows.items():
                try:
                    valid_rows[rid] = source_row(row, document, raw_assets)
                except (
                    ValueError,
                    KeyError,
                    TypeError,
                    IndexError,
                    UnicodeError,
                ) as exc:
                    check(False, "source_row_invalid", {"row": rid, "reason": str(exc)})
        atoms = p.get("atoms", [])
        check(len(atoms) == len(template["atoms"]), "atom_coverage_mismatch")
        labels = {}
        for actual, expected in zip(atoms, template["atoms"]):
            check(set(actual) == set(expected), "atom_fields_invalid", expected["id"])
            for key in ("id", "element", "mass"):
                check(
                    actual.get(key) == expected[key],
                    "atom_" + key + "_mismatch",
                    expected["id"],
                )
            t = actual.get("type")
            check(
                type(t) is str and bool(t) and "*" not in t,
                "explicit_type_missing",
                expected["id"],
            )
            labels[expected["id"]] = t
            row = rows.get(actual.get("evidence_row"), {})
            check(
                row.get("family") == "atom_type"
                and row.get("types") == [t]
                and row.get("values") == [expected["mass"]]
                and row.get("units") == ["dalton"]
                and actual.get("evidence_row") in valid_rows,
                "typing_mass_evidence_missing",
                expected["id"],
            )
        expected_interactions = template["coverage"]
        interactions = p.get("interactions", {})
        check(set(interactions) <= set(expected_interactions), "unexpected_interaction")
        for key, required in expected_interactions.items():
            item = interactions.get(key)
            if not check(
                item is not None, "missing_interaction", {"key": key, **required}
            ):
                continue
            check(
                set(item) == {"row", "types", "dependencies"},
                "interaction_fields_invalid",
                key,
            )
            check(
                item.get("dependencies") == required["dependencies"],
                "equilibrium_dependencies_mismatch",
                key,
            )
            types = [labels.get(i) for i in required["sites"]]
            check(item.get("types") == types, "interaction_type_order_mismatch", key)
            row = rows.get(item.get("row"), {})
            f = required["family"]
            units = FAMILIES[f][2]
            check(
                row.get("family") == f
                and row.get("types") == types
                and row.get("units") == units
                and len(row.get("values", [])) == len(units)
                and item.get("row") in valid_rows,
                "coefficient_evidence_missing",
                key,
            )
            for dep in required["dependencies"]:
                check(
                    dep in interactions,
                    "equilibrium_dependency_missing",
                    {"consumer": key, "dependency": dep},
                )
        policy = p.get("nonbonded_policy", {})
        check(policy.get("mixing") == MIXING, "lj_mixing_mismatch")
        check(policy.get("special_bonds") == SPECIAL, "special_bonds_mismatch")
        policy_row = rows.get(policy.get("evidence_row"), {})
        check(
            policy_row.get("family") == "nonbonded_policy"
            and policy_row.get("types") == [MIXING]
            and policy_row.get("values") == ["0", "0", "1", "0", "0", "1"]
            and policy_row.get("units") == ["dimensionless"] * 6
            and policy.get("evidence_row") in valid_rows,
            "nonbonded_policy_evidence_missing",
        )
        origin = p.get("charge_origin")
        check(
            origin in ("external_explicit_vector", "external_source_bond_increments"),
            "charge_origin_unsupported",
        )
        charges = p.get("charges", [])
        check(
            [q.get("id") for q in charges] == template["atom_order"],
            "charge_order_coverage_mismatch",
        )
        qvalues = {}
        for q in charges:
            check(set(q) == {"id", "value", "evidence_row"}, "charge_fields_invalid")
            try:
                qvalues[q["id"]] = number(q["value"])
            except ValueError as exc:
                check(False, "charge_value_invalid", str(exc))
            row = rows.get(q.get("evidence_row"), {})
            check(
                row.get("family") == "charge"
                and row.get("types") == [str(q.get("id"))]
                and row.get("values") == [q.get("value")]
                and row.get("units") == ["elementary_charge"]
                and q.get("evidence_row") in valid_rows,
                "charge_vector_evidence_mismatch",
                q.get("id"),
            )
        if len(qvalues) == len(template["atoms"]):
            check(sum(qvalues.values(), Decimal(0)) == 0, "source_charge_total_invalid")
            for component in template["graph"]["components"]:
                check(
                    all(i in qvalues for i in component)
                    and sum((qvalues[i] for i in component), Decimal(0)) == 0,
                    "source_component_charge_total_invalid",
                    component,
                )
        increments = p.get("charge_increments", {})
        if origin == "external_source_bond_increments":
            reconstructed = dict.fromkeys(template["atom_order"], Decimal(0))
            expected_bonds = {
                interaction_key("increment", b["sites"]): b["sites"]
                for b in template["graph"]["bonds"]
            }
            check(set(increments) == set(expected_bonds), "charge_increments_missing")
            for key, (i, j) in expected_bonds.items():
                rid = increments.get(key)
                row = rows.get(rid, {})
                good = check(
                    row.get("family") == "bond_increment"
                    and row.get("types") == [labels[i], labels[j]]
                    and row.get("units") == ["elementary_charge"]
                    and len(row.get("values", [])) == 1
                    and rid in valid_rows,
                    "charge_increment_evidence_missing",
                    key,
                )
                if good:
                    v = valid_rows[rid][0]
                    reconstructed[i] += v
                    reconstructed[j] -= v
            check(reconstructed == qvalues, "increment_charge_vector_mismatch")
        else:
            check(increments == {}, "unexpected_charge_increments")
    except (
        ValueError,
        TypeError,
        KeyError,
        AttributeError,
        IndexError,
        OSError,
        UnicodeError,
        IslandError,
        ArithmeticError,
    ) as exc:
        # Invalid external input must produce a classified result, not escape a CLI.
        check(False, "malformed_bundle", str(exc))
    status = "parameter_incomplete" if errors else "parameter_bundle_validated"
    return ExternalBundleValidation(
        pack(
            {
                "schema": "island_external_bundle_validation_v1",
                "status": status,
                "bundle_identity": bundle.identity,
                "diagnostics": errors,
                "execution_authorized": False,
                "execution_profile": None,
                "backend_status": "scientifically_unqualified",
                "limitations": LIMITS,
            }
        )
    )


def load_external_pdms_bundle(path, system, *, expected_sha256, assets=None):
    """Explicit path + out-of-band pin. No CAR/MDF or implicit source discovery."""
    path = Path(path)
    if not path.is_file():
        return None, ExternalBundleValidation(
            pack(
                {
                    "schema": "island_external_bundle_validation_v1",
                    "status": "asset_unavailable",
                    "bundle_identity": None,
                    "diagnostics": [
                        {"stage": "asset_unavailable", "context": path.name}
                    ],
                    "execution_authorized": False,
                    "execution_profile": None,
                    "backend_status": "scientifically_unqualified",
                    "limitations": LIMITS,
                }
            )
        )
    bundle = ExternalParameterBundle(path.read_text())
    return bundle, bundle.validate(
        system,
        assets=path.parent if assets is None else assets,
        expected_sha256=expected_sha256,
    )
