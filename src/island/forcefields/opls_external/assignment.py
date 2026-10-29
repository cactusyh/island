"""Replayable graph/type/charge/row assignment, isolated from native records."""

from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import degrees

from island.charge_references.records import pack, unpack
from island.dynamics._checkpoint_data import strict_load
from island.exceptions import OPLSDependencyError
from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.oplsaa.parameters import inventory

from .graph import extension_graph
from .source import FORMS, ExternalOPLSError, decimal, inspect_xml, load_source

FAMILY = "external_opls"


@dataclass(frozen=True)
class ExternalOPLSOptions:
    source_path: object = None
    expected_sha256: str | None = None


@lru_cache(maxsize=64)
def _typing(graph_text, xml):
    """Pure cache; its key includes the entire source and graph, never coordinates."""
    from xml.etree import ElementTree as ET

    from island.forcefields.oplsaa.adapter import check_foyer_installation

    Provider, find, Graph, _ = check_foyer_installation()
    graph, root = unpack(graph_text), ET.fromstring(xml)
    types = {a.attrib["name"]: dict(a.attrib) for a in root.find("AtomTypes")}
    provider = Provider(
        {n: a["def"] for n, a in types.items() if a.get("def")},
        {
            n: {s.strip() for s in a.get("overrides", "").split(",") if s.strip()}
            for n, a in types.items()
        },
        set(),
    )
    external = Graph()
    ids = [a["id"] for a in graph["sites"]]
    index = {i: n for n, i in enumerate(ids)}
    for n, atom in enumerate(graph["sites"]):
        external.add_atom(
            n,
            atom["element"],
            atomic_number=atom["atomic_number"],
            symbol=atom["element"],
        )
    for row in graph["bonds"]:
        a, b = row["sites"]
        external.add_bond(index[a], index[b], bond_order=row["order"])
    try:
        matches = find(external, provider, max_iter=len(types) * len(ids) + 1)
        require(set(matches) == set(range(len(ids))), "typing_coverage")
        require(
            all(row["atomtype"] in types for row in matches.values()), "unknown_type"
        )
    except Exception as error:
        raise ExternalOPLSError(
            "typing_blocked",
            [
                {
                    "stage": "untyped_or_ambiguous_environment",
                    "reason": str(error),
                    "index_to_site_id": ids,
                }
            ],
        ) from error
    return pack(
        {
            i: {
                "atomtype": matches[n]["atomtype"],
                "whitelist": sorted(matches[n]["whitelist"]),
                "blacklist": sorted(matches[n]["blacklist"]),
            }
            for n, i in enumerate(ids)
        }
    )


def _select(root, types, names, section):
    tag, numeric = FORMS[section]
    matches = []
    for index, element in enumerate(root.findall(section + "/" + tag)):
        a = element.attrib
        arity = len(names)
        roles, wildcard = [], False
        for j in range(1, arity + 1):
            require(
                (f"type{j}" in a) != (f"class{j}" in a), "source_role_label_required"
            )
            key = f"type{j}" if f"type{j}" in a else f"class{j}"
            value = a[key]
            wildcard |= value == ""
            roles.append((key.startswith("type"), value))
        require(
            set(a)
            == numeric
            | {
                f"type{j}" if f"type{j}" in a else f"class{j}"
                for j in range(1, arity + 1)
            },
            "unknown_parameter_columns",
        )
        if any(
            all(
                not value or (t if is_type else types[t]["class"]) == value
                for t, (is_type, value) in zip(order, roles, strict=True)
            )
            for order in (names, names[::-1])
        ):
            matches.append((index, dict(a), wildcard))
    require(bool(matches), "missing_source_row")
    specific = [r for r in matches if not r[2]]
    chosen = specific or matches
    require(
        len({tuple(decimal(a[k]) for k in sorted(numeric)) for _, a, _ in chosen}) == 1,
        "ambiguous_conflicting_source_rows",
    )
    index, a, wildcard = chosen[0]
    raw = {k: a[k] for k in sorted(numeric)}
    if section == "HarmonicBondForce":
        require(
            decimal(raw["length"]) > 0 and decimal(raw["k"]) >= 0, "invalid_bond_values"
        )
        converted = {
            "length_angstrom": float(Fraction(decimal(raw["length"])) * 10),
            "k_kj_mol_angstrom2": float(Fraction(decimal(raw["k"])) / 100),
        }
    elif section == "HarmonicAngleForce":
        require(
            0 < decimal(raw["angle"]) < Decimal("3.141592653589794")
            and decimal(raw["k"]) >= 0,
            "invalid_angle_values",
        )
        converted = {
            "angle_degrees": degrees(float(raw["angle"])),
            "k_kj_mol_radian2": float(raw["k"]),
        }
    else:
        prefix, start, count = (
            ("c", 0, 6) if section == "RBTorsionForce" else ("k", 1, 4)
        )
        converted = {
            "coefficients_kj_mol": [
                float(raw[f"{prefix}{n}"]) for n in range(start, start + count)
            ]
        }
    row = {
        "source_record": f"{section}/{tag}[{index}]",
        "types": list(names),
        "source_attributes": a,
        "original_values": raw,
        "original_units": "openmm_XML",
        "converted": converted,
        "conversion_code": "external_opls_openmm_to_A_v1",
        "candidate_indices": [r[0] for r in matches],
        "selected_wildcard": wildcard,
    }
    return {**row, "evidence_checksum": digest(row)}


def resolve(system, graph, source):
    from island.forcefields.oplsaa.adapter import check_foyer_installation

    source.validate_integrity()
    try:
        environment = check_foyer_installation()[
            3
        ]  # Recheck implementation pins on cache hits.
    except OPLSDependencyError as error:
        raise ExternalOPLSError(
            "asset_unavailable",
            [{"stage": "pinned_typing_engine_unavailable", "reason": str(error)}],
        ) from error
    return unpack(
        _resolved(
            pack(graph),
            pack(inventory(system)),
            source.xml,
            source.manifest_text,
            pack(environment),
        )
    )


@lru_cache(maxsize=32)
def _resolved(graph_text, inventory_text, xml, manifest_text, environment_text):
    graph, document = unpack(graph_text), strict_load(manifest_text)
    root, types, nb, evidence = inspect_xml(xml, document)
    missing_elements = {a["element"] for a in graph["sites"]} - set(
        document["model"]["elements"]
    )
    if missing_elements:
        raise ExternalOPLSError(
            "typing_blocked",
            [
                {
                    "stage": "source_has_no_declared_types_for_elements",
                    "missing_elements": sorted(missing_elements),
                }
            ],
        )
    if graph["impropers"]:
        raise ExternalOPLSError(
            "parameter_incomplete",
            [
                {
                    "stage": "improper_form_not_supported_v1",
                    "required": graph["impropers"],
                }
            ],
        )
    matches = unpack(_typing(graph_text, xml))
    assigned = {i: row["atomtype"] for i, row in matches.items()}
    # Fraction preserves every supplied decimal digit even beyond Decimal's
    # ambient precision. A small residual must never round itself to neutrality.
    sites, component_charges = {}, defaultdict(Fraction)
    component_of = {i: n for n, c in enumerate(graph["components"]) for i in c}
    for atom in graph["sites"]:
        i, t = atom["id"], assigned[atom["id"]]
        a, n = types[t], nb[t]
        require(a["element"] == atom["element"], "assigned_type_element_mismatch")
        require(
            decimal(a["mass"]) == Decimal(str(atom["mass"])),
            "source_graph_mass_mismatch",
        )
        component_charges[component_of[i]] += Fraction(decimal(n["charge"]))
        sites[i] = {
            "type": t,
            "class": a["class"],
            "mass_dalton": float(a["mass"]),
            "charge_e": float(n["charge"]),
            "exact_charge_e": n["charge"],
            "sigma_angstrom": float(Fraction(decimal(n["sigma"])) * 10),
            "epsilon_kj_mol": float(n["epsilon"]),
            "source_attributes": dict(n),
            "mass_source_attributes": dict(a),
        }
    require(
        all(q == 0 for q in component_charges.values()),
        "source_component_charge_not_exactly_neutral_no_normalization_permitted",
    )
    paths, pairs = unpack(inventory_text)
    special = root.find("NonbondedForce").attrib
    for pair in pairs:
        pair["lj_scale"] = float(special["lj14scale"]) if pair["distance"] == 3 else 0.0
        pair["coulomb_scale"] = (
            float(special["coulomb14scale"]) if pair["distance"] == 3 else 0.0
        )
    result, missing = {}, []
    torsion = (
        "RBTorsionForce"
        if root.find("RBTorsionForce") is not None
        else "OPLSTorsionForce"
    )
    for kind, section in (
        ("bonds", "HarmonicBondForce"),
        ("angles", "HarmonicAngleForce"),
        ("rb", torsion),
    ):
        rows = []
        for sites_for_row in paths[kind]:
            names = tuple(assigned[i] for i in sites_for_row)
            try:
                row = _select(root, types, names, section)
                rows.append({**row, "sites": list(sites_for_row)})
            except (KeyError, ValueError) as error:
                missing.append(
                    {
                        "stage": "parameter_row_missing_or_ambiguous",
                        "family": kind,
                        "sites": list(sites_for_row),
                        "types": list(names),
                        "reason": str(error),
                    }
                )
        result[kind] = rows
    if missing:
        raise ExternalOPLSError("parameter_incomplete", missing)
    used = {
        r["source_record"] for kind in ("bonds", "angles", "rb") for r in result[kind]
    }
    used_types = set(assigned.values())
    for section, tag, label in (
        ("AtomTypes", "Type", "name"),
        ("NonbondedForce", "Atom", "type"),
    ):
        used.update(
            f"{section}/{tag}[{j}]"
            for j, r in enumerate(root.find(section))
            if r.attrib[label] in used_types
        )
    evidence = [r for r in evidence if r["source_record"] in used]
    require(
        used == {r["source_record"] for r in evidence},
        "selected_source_evidence_incomplete",
    )
    result.update(
        sites=sites,
        pairs=pairs,
        torsion_form=torsion,
        typing_matches=matches,
        exact_component_charges={str(i): str(q) for i, q in component_charges.items()},
        source_identity=digest(
            {
                "manifest_sha256": sha256(manifest_text.encode()).hexdigest(),
                "xml_sha256": sha256(xml).hexdigest(),
            }
        ),
        graph_identity=graph["identity"],
        environment=unpack(environment_text),
        source_row_evidence=evidence,
        model=document["model"],
        special_pairs={k: special[k] for k in sorted(special)},
    )
    return pack(result)


def _make(system, source):
    graph, data, diagnostics, status = None, None, [], "asset_unavailable"
    source_error = None
    if type(source) is ExternalOPLSOptions:
        try:
            source = load_source(source.source_path, source.expected_sha256)
        except ExternalOPLSError as error:
            source, source_error = None, error
    try:
        graph = extension_graph(system)
    except Exception as error:  # noqa: BLE001 -- classified final-graph boundary
        status, diagnostics = (
            "typing_blocked",
            [{"stage": "final_graph_invalid", "reason": str(error)}],
        )
    if graph is not None:
        if source is None:
            if source_error is not None:
                status, diagnostics = source_error.status, source_error.diagnostics
            else:
                diagnostics = [
                    {"stage": "explicit_external_source_and_independent_pin_required"}
                ]
        else:
            try:
                source.validate_integrity()
                diagnostics = source.execution_diagnostics()
                if diagnostics:
                    status = "scientifically_unqualified"
                else:
                    data = resolve(system, graph, source)
                    status = "prepared"
            except ExternalOPLSError as error:
                status, diagnostics = error.status, error.diagnostics
            except Exception as error:  # noqa: BLE001 -- classified assignment boundary
                status, diagnostics = (
                    "parameter_incomplete",
                    [{"stage": "assignment_validation_failed", "reason": str(error)}],
                )
    p = {
        "schema": "island_external_opls_preparation_v1",
        "resolver": "island_external_opls_resolver_v1",
        "family": FAMILY,
        "input_binding": graph["identity"]
        if graph is not None
        else _raw_binding(system),
        "status": status,
        "graph": graph,
        "assignment": data,
        "diagnostics": diagnostics,
        "profile_identity": source.document["profile_id"]
        if status == "prepared"
        else None,
        "requested_profile_identity": None
        if source is None
        else source.document["profile_id"],
        "source_identity": None if source is None else source.identity,
        "native_opls_coverage": False,
        "scientifically_validated": False,
        "backend_scope": "finite_nonperiodic_NoCutoff_OpenMM",
        "execution_permitted": status == "prepared",
        "limitations": LIMITS,
    }
    return {**p, "identity": digest(p)}


def _raw_binding(system):
    try:
        payload = deepcopy(system.to_dict())
        payload.pop("coordinates")
        return digest(payload)
    except (AttributeError, TypeError, ValueError):
        return None


@dataclass(frozen=True)
class ExternalOPLSPreparation:
    json_text: str
    _system: object
    _source: object

    def __post_init__(self):
        object.__setattr__(self, "_system", deepcopy(self._system))
        self.validate_integrity()

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]

    def validate_integrity(self, system=None):
        try:
            p = self.payload
            require(
                pack(p) == pack(_make(self._system, self._source)),
                "preparation_record_or_source_changed",
            )
            if system is not None:
                if p["graph"] is None:
                    require(
                        _raw_binding(system) == p["input_binding"],
                        "stale_invalid_input",
                    )
                else:
                    require(
                        pack(extension_graph(system)) == pack(p["graph"]),
                        "stale_final_graph_or_atom_order",
                    )
        except Exception as error:
            raise ExternalOPLSError(
                "parameter_incomplete",
                [{"stage": "preparation_integrity_failed", "reason": str(error)}],
                result=self,
            ) from error

    def require_executable(self, system=None):
        self.validate_integrity(system)
        p = self.payload
        if p["status"] != "prepared":
            raise ExternalOPLSError(p["status"], p["diagnostics"], result=self)


def prepare_external_opls(system, options):
    require(type(options) is ExternalOPLSOptions, "ExternalOPLSOptions_required")
    if options.source_path is None:
        source = None
    else:
        try:
            source = load_source(options.source_path, options.expected_sha256)
        except ExternalOPLSError:
            source = options  # Replay failed input too; it confers no authority.
    return ExternalOPLSPreparation(pack(_make(system, source)), system, source)
