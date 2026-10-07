"""Unified, force-field-neutral typed graph and assignment contracts.

These records are adapters around family-owned records. They never infer types,
borrow parameters, normalize charges, or alter historical serialized payloads.
"""

from copy import deepcopy
from dataclasses import dataclass
from math import fsum, isfinite

from island.exceptions import ValidationError

from .final import _digest, pack, unpack

TYPED_SCHEMA = "island_unified_typed_graph_v1"
CHARGES_SCHEMA = "island_unified_graph_charges_v1"
SOURCE_SCHEMA = "island_unified_force_field_source_v1"
ASSIGNMENT_SCHEMA = "island_unified_parameter_assignment_v1"
DIAGNOSTICS_SCHEMA = "island_unified_assignment_diagnostics_v1"
FAMILIES = {"PCFF", "OPLS-AA", "GAFF", "GAFF2"}
UNITS = {"elementary_charge"}


def _require(ok, message):
    if not ok:
        raise ValidationError(message)


def _finite(value, label):
    _require(
        type(value) in (int, float) and not isinstance(value, bool) and isfinite(value),
        f"{label} must be finite",
    )


def _source(source, family=None):
    _require(
        type(source) is dict
        and type(source.get("sha256")) is str
        and len(source["sha256"]) == 64,
        "Source identity must contain SHA-256",
    )
    if family is not None:
        _require(
            source.get("family", family)
            in (family, "PCFF", "OPLS-AA", "GAFF", "GAFF2"),
            "Source family mismatch",
        )


def _ids(graph):
    return {s["id"] for s in graph.payload["sites"]}


def _record_identity(payload):
    return _digest({k: v for k, v in payload.items() if k != "identity"})


def _json_safe(value):
    if isinstance(value, dict):
        return {
            str(k) if isinstance(k, int) else k: _json_safe(v) for k, v in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


@dataclass(frozen=True)
class UnifiedForceFieldSource:
    json_text: str

    def validate_integrity(self):
        p = unpack(self.json_text)
        _require(p.get("schema") == SOURCE_SCHEMA, "Unsupported unified source schema")
        _require(
            p.get("identity") == _record_identity(p), "Unified source identity mismatch"
        )
        _require(
            p.get("force_field") in FAMILIES, "Unsupported unified force-field family"
        )
        _source(p["source"], p["force_field"])
        _require(
            type(p["provenance"]) is str and p["provenance"].strip(),
            "Source provenance required",
        )
        _require(
            type(p["evidence"]) is list
            and all(type(x) is str and x.strip() for x in p["evidence"]),
            "Source evidence references required",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]


def unified_force_field_source(
    *,
    force_field,
    source,
    provenance,
    evidence=None,
    executable=None,
    source_files=None,
):
    _require(force_field in FAMILIES, "Unsupported unified force-field family")
    _source(source, force_field)
    _require(
        type(provenance) is str and provenance.strip(), "Source provenance required"
    )
    p = {
        "schema": SOURCE_SCHEMA,
        "force_field": force_field,
        "source": deepcopy(source),
        "source_files": deepcopy(source_files or {}),
        "executable": deepcopy(executable),
        "provenance": provenance,
        "evidence": list(evidence or []),
    }
    p["identity"] = _record_identity(p)
    result = UnifiedForceFieldSource(pack(p))
    result.validate_integrity()
    return result


@dataclass(frozen=True)
class UnifiedTypedGraph:
    json_text: str

    def validate_integrity(self, graph=None):
        p = unpack(self.json_text)
        _require(
            p.get("schema") == TYPED_SCHEMA, "Unsupported unified typed-graph schema"
        )
        _require(
            p.get("identity") == _record_identity(p),
            "Unified typed-graph identity mismatch",
        )
        _require(
            p.get("force_field") in FAMILIES, "Unsupported unified force-field family"
        )
        _require(
            type(p.get("graph_identity")) is str
            and p["graph_identity"] == p.get("final_graph_identity"),
            "Unified graph identity mismatch",
        )
        _source(p["source"], p["force_field"])
        _require(type(p["atom_types"]) is dict, "Unified atom types required")
        if graph is not None:
            graph.validate_integrity()
            _require(
                graph.identity == p["graph_identity"],
                "Unified typed graph/final graph mismatch",
            )
            ids = _ids(graph)
        else:
            ids = {int(i) for i in p["atom_types"]}
        _require(
            set(p["atom_types"]) == {str(i) for i in ids},
            "Unified atom type coverage mismatch",
        )
        _require(
            all(type(v) is str and v for v in p["atom_types"].values()),
            "Unified atom labels required",
        )
        for key in ("typing_profile", "typing_method", "provenance"):
            _require(type(p[key]) is str and p[key].strip(), f"Unified {key} required")
        _require(
            type(p["evidence"]) is list
            and all(type(x) is str and x.strip() for x in p["evidence"]),
            "Unified typing evidence required",
        )
        _require(
            type(p["automatic_perception"]) is str
            and p["automatic_perception"]
            in {"not_performed", "performed", "unavailable"},
            "Invalid automatic perception status",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]

    @property
    def atom_types(self):
        return {int(i): v for i, v in self.payload["atom_types"].items()}


def unified_typed_graph(
    graph,
    *,
    force_field,
    atom_types,
    source,
    typing_profile,
    typing_method,
    provenance,
    evidence=None,
    automatic_perception="not_performed",
):
    graph.validate_integrity()
    _require(force_field in FAMILIES, "Unsupported unified force-field family")
    _source(source, force_field)
    ids = _ids(graph)
    _require(
        type(atom_types) is dict and set(atom_types) == ids,
        "Exact unified atom type coverage required",
    )
    p = {
        "schema": TYPED_SCHEMA,
        "final_graph_identity": graph.identity,
        "graph_identity": graph.identity,
        "component_identity": _digest(graph.payload["components"]),
        "molecule_identity": _digest(graph.payload["molecule_membership"]),
        "force_field": force_field,
        "source": deepcopy(source),
        "atom_types": {str(i): atom_types[i] for i in sorted(ids)},
        "typing_profile": typing_profile,
        "typing_method": typing_method,
        "provenance": provenance,
        "evidence": list(evidence or []),
        "automatic_perception": automatic_perception,
    }
    p["identity"] = _record_identity(p)
    result = UnifiedTypedGraph(pack(p))
    result.validate_integrity(graph)
    return result


@dataclass(frozen=True)
class UnifiedGraphCharges:
    json_text: str

    def validate_integrity(self, graph=None, typed_graph=None):
        p = unpack(self.json_text)
        _require(p.get("schema") == CHARGES_SCHEMA, "Unsupported unified charge schema")
        _require(
            p.get("identity") == _record_identity(p), "Unified charge identity mismatch"
        )
        _require(
            p.get("force_field") in FAMILIES, "Unsupported charge force-field family"
        )
        _source(p["source"], p["force_field"])
        if graph is not None:
            graph.validate_integrity()
            _require(
                graph.identity == p["graph_identity"],
                "Unified charge/final graph mismatch",
            )
        if typed_graph is not None:
            typed_graph.validate_integrity(graph)
            _require(
                typed_graph.identity == p["typed_graph_identity"],
                "Unified charge/typed graph mismatch",
            )
        ids = _ids(graph) if graph is not None else {int(i) for i in p["charges"]}
        _require(
            set(p["charges"]) == {str(i) for i in ids},
            "Unified charge coverage mismatch",
        )
        _require(p["unit"] in UNITS, "Unsupported unified charge units")
        for value in p["charges"].values():
            _finite(value, "Unified charge")
        _require(
            type(p["component_totals"]) is dict, "Unified component totals required"
        )
        _finite(p["total_charge"], "Unified total charge")
        if graph is not None:
            components = graph.payload["components"]
            _require(
                set(p["component_totals"]) == {str(c[0]) for c in components},
                "Unified component total coverage mismatch",
            )
            charge_values = {int(i): value for i, value in p["charges"].items()}
            for component in components:
                declared = p["component_totals"][str(component[0])]
                _finite(declared, "Unified component total")
                _require(
                    abs(fsum(charge_values[i] for i in component) - declared) <= 1e-12,
                    f"Unified component total mismatch: {component}",
                )
        _require(
            abs(fsum(p["charges"].values()) - p["total_charge"]) <= 1e-12,
            "Unified total charge mismatch",
        )
        _require(
            type(p["charge_method"]) is str and p["charge_method"],
            "Charge method required",
        )
        _require(
            type(p["charge_origin"]) is str and p["charge_origin"],
            "Charge origin required",
        )
        _require(
            type(p["provenance"]) is str and p["provenance"].strip(),
            "Charge provenance required",
        )
        _require(
            type(p["evidence"]) is list
            and all(type(x) is str and x.strip() for x in p["evidence"]),
            "Charge evidence required",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]

    @property
    def charges(self):
        return {int(i): v for i, v in self.payload["charges"].items()}


def unified_graph_charges(
    graph,
    typed_graph,
    *,
    force_field,
    charges,
    charge_method,
    charge_origin,
    source,
    component_totals,
    total_charge,
    provenance,
    evidence=None,
    unit="elementary_charge",
):
    graph.validate_integrity()
    typed_graph.validate_integrity(graph)
    _require(
        force_field in FAMILIES and typed_graph.payload["force_field"] == force_field,
        "Charge family mismatch",
    )
    _source(source, force_field)
    ids = _ids(graph)
    _require(
        type(charges) is dict and set(charges) == ids,
        "Exact unified charge coverage required",
    )
    for value in charges.values():
        _finite(value, "Unified charge")
    p = {
        "schema": CHARGES_SCHEMA,
        "graph_identity": graph.identity,
        "typed_graph_identity": typed_graph.identity,
        "force_field": force_field,
        "source": deepcopy(source),
        "charges": {str(i): charges[i] for i in sorted(ids)},
        "charge_method": charge_method,
        "charge_origin": charge_origin,
        "unit": unit,
        "component_totals": deepcopy(component_totals),
        "total_charge": total_charge,
        "provenance": provenance,
        "evidence": list(evidence or []),
    }
    p["identity"] = _record_identity(p)
    result = UnifiedGraphCharges(pack(p))
    result.validate_integrity(graph, typed_graph)
    return result


@dataclass(frozen=True)
class UnifiedAssignmentDiagnostics:
    json_text: str

    def validate_integrity(self):
        p = unpack(self.json_text)
        _require(
            p.get("schema") == DIAGNOSTICS_SCHEMA
            and p.get("identity") == _record_identity(p),
            "Unified diagnostics identity mismatch",
        )
        for key in (
            "typing_complete",
            "charges_complete",
            "native_charge_available",
            "parameter_rows_complete",
            "executable_potential_complete",
            "periodic_backend_support",
        ):
            _require(
                type(p.get(key)) is bool, f"Unified diagnostic {key} must be boolean"
            )
        _require(
            p["scientific_verification_status"]
            in {"not_attested", "software_test", "independent_numeric"},
            "Invalid scientific verification status",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]


def unified_assignment_diagnostics(
    *,
    typing_complete,
    charges_complete,
    native_charge_available,
    parameter_rows_complete,
    executable_potential_complete,
    periodic_backend_support,
    scientific_verification_status="not_attested",
    reasons=None,
):
    p = {
        "schema": DIAGNOSTICS_SCHEMA,
        "typing_complete": typing_complete,
        "charges_complete": charges_complete,
        "native_charge_available": native_charge_available,
        "parameter_rows_complete": parameter_rows_complete,
        "executable_potential_complete": executable_potential_complete,
        "periodic_backend_support": periodic_backend_support,
        "scientific_verification_status": scientific_verification_status,
        "reasons": list(reasons or []),
    }
    p["identity"] = _record_identity(p)
    result = UnifiedAssignmentDiagnostics(pack(p))
    result.validate_integrity()
    return result


@dataclass(frozen=True)
class UnifiedParameterAssignment:
    json_text: str

    def validate_integrity(self, graph=None, typed_graph=None, charges=None):
        p = unpack(self.json_text)
        _require(
            p.get("schema") == ASSIGNMENT_SCHEMA
            and p.get("identity") == _record_identity(p),
            "Unified assignment identity mismatch",
        )
        _require(p.get("force_field") in FAMILIES, "Unsupported assignment family")
        if graph is not None:
            graph.validate_integrity()
            _require(
                graph.identity == p["graph_identity"],
                "Unified assignment graph mismatch",
            )
        if typed_graph is not None:
            typed_graph.validate_integrity(graph)
            _require(
                typed_graph.identity == p["typed_graph_identity"],
                "Unified assignment typed graph mismatch",
            )
        if charges is not None:
            charges.validate_integrity(graph, typed_graph)
            _require(
                charges.identity == p["charges_identity"],
                "Unified assignment charge mismatch",
            )
        UnifiedAssignmentDiagnostics(pack(p["diagnostics"])).validate_integrity()
        _source(p["source"], p["force_field"])

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]


def unified_parameter_assignment(
    graph, typed_graph, charges, *, force_field, source, native_payload, diagnostics
):
    graph.validate_integrity()
    typed_graph.validate_integrity(graph)
    charges.validate_integrity(graph, typed_graph)
    _require(
        force_field
        == typed_graph.payload["force_field"]
        == charges.payload["force_field"],
        "Assignment family mismatch",
    )
    _source(source, force_field)
    diagnostics.validate_integrity()
    p = {
        "schema": ASSIGNMENT_SCHEMA,
        "graph_identity": graph.identity,
        "typed_graph_identity": typed_graph.identity,
        "charges_identity": charges.identity,
        "force_field": force_field,
        "source": deepcopy(source),
        "native_payload": _json_safe(native_payload),
        "diagnostics": diagnostics.payload,
    }
    p["identity"] = _record_identity(p)
    result = UnifiedParameterAssignment(pack(p))
    result.validate_integrity(graph, typed_graph, charges)
    return result


def unify_pcff_typed_graph(graph, record, *, source=None):
    from island.forcefields.pcff.typed_graph import PCFFTypedGraph

    if type(record) is not PCFFTypedGraph:
        raise TypeError("PCFFTypedGraph required")
    p = record.payload
    src = source or p["source"]
    return unified_typed_graph(
        graph,
        force_field="PCFF",
        atom_types=p["assignments"],
        source=src,
        typing_profile=p["profile"]["name"],
        typing_method="external_source_bound",
        provenance=p["provenance"],
        evidence=p["evidence_references"],
        automatic_perception=p["validation"]["automatic_perception"],
    )


def unify_pcff_charges(graph, typed, record):
    from island.graph.charges import charge_record_from_pcff

    neutral = charge_record_from_pcff(graph, record)
    p = neutral.payload
    source = p["source"]
    return unified_graph_charges(
        graph,
        typed,
        force_field="PCFF",
        charges=neutral.charges,
        charge_method=p["method"],
        charge_origin=p["origin"],
        source=source,
        component_totals={
            str(component[0]): sum(neutral.charges[i] for i in component)
            for component in graph.payload["components"]
        },
        total_charge=p["total_charge"],
        provenance=p["provenance"],
        evidence=p["evidence"],
    )


def unify_pcff_assignment(graph, typed, charges, record):
    p = record.payload
    complete = bool(p.get("parameter_coverage_complete"))
    diagnostics = unified_assignment_diagnostics(
        typing_complete=True,
        charges_complete=True,
        native_charge_available=charges.payload["charge_origin"] == "native_increments",
        parameter_rows_complete=complete,
        executable_potential_complete=False,
        periodic_backend_support=not graph.periodic,
        reasons=p.get("assignments", []),
    )
    return unified_parameter_assignment(
        graph,
        typed,
        charges,
        force_field="PCFF",
        source=p["source"],
        native_payload=p,
        diagnostics=diagnostics,
    )


def unify_opls_typed_graph(graph, typing, source):
    p = typing.payload
    _require_opls_graph_binding(graph, p["graph"])
    assignments = p["data"]["assignments"]
    return unified_typed_graph(
        graph,
        force_field="OPLS-AA",
        atom_types=assignments,
        source=source.identity,
        typing_profile="foyer_" + p["environment"]["foyer"],
        typing_method="foyer_source_typing",
        provenance="Validated OPLS-AA typing record",
        evidence=[typing.identity],
        automatic_perception="performed",
    )


def _require_opls_graph_binding(graph, native_graph):
    neutral = graph.payload
    neutral_sites = [
        (
            s["id"],
            s["element"],
            s["formal_charge"],
            s["metadata"].get("aromatic", False),
        )
        for s in neutral["sites"]
    ]
    native_sites = [
        (
            s["id"],
            s["element"],
            s["formal_charge"],
            s["chemical_metadata"].get("aromatic", False),
        )
        for s in native_graph["sites"]
    ]
    _require(
        neutral_sites == native_sites and neutral["bonds"] == native_graph["bonds"],
        "OPLS record chemical graph differs from final graph",
    )


def unify_opls_charges(graph, typed, charges, source):
    p = charges.payload
    return unified_graph_charges(
        graph,
        typed,
        force_field="OPLS-AA",
        charges={i: row["charge"] for i, row in p["charges"].items()},
        charge_method="native_source",
        charge_origin="native_source",
        source=source.identity,
        component_totals={
            str(c["sites"][0]): c["total_e"] for c in p["data"]["components"]
        },
        total_charge=p["data"]["total_e"],
        provenance="Validated OPLS-AA charge record",
        evidence=[charges.identity],
    )


def unify_opls_assignment(graph, typed, charges, result, source):
    p = result.payload
    diagnostics = unified_assignment_diagnostics(
        typing_complete=True,
        charges_complete=True,
        native_charge_available=True,
        parameter_rows_complete=True,
        executable_potential_complete=True,
        periodic_backend_support=not graph.periodic,
        reasons=[],
    )
    return unified_parameter_assignment(
        graph,
        typed,
        charges,
        force_field="OPLS-AA",
        source=source.identity,
        native_payload=p,
        diagnostics=diagnostics,
    )


def unify_amber_preparation(graph, result, *, force_field):
    _require(force_field in {"GAFF", "GAFF2"}, "GAFF family required")
    imported = result.imported_result
    charge = imported.charge_result
    typed = unified_typed_graph(
        graph,
        force_field=force_field,
        atom_types=dict(imported.atom_types),
        source={
            "family": force_field,
            "sha256": imported.source_sha256,
            "source": imported.source,
        },
        typing_profile=imported.parser_version,
        typing_method="ambertools_imported_prmtop",
        provenance=str(imported.provenance),
        evidence=[imported.result_signature],
        automatic_perception="performed",
    )
    charges = unified_graph_charges(
        graph,
        typed,
        force_field=force_field,
        charges={i: a.charge for i, a in charge.assignments.items()},
        charge_method=result.record["charge_method"],
        charge_origin="provided"
        if result.record["charge_method"] == "provided"
        else "am1bcc",
        source={
            "family": force_field,
            "sha256": imported.source_sha256,
            "source": imported.source,
        },
        component_totals={
            str(c.site_ids[0]): c.assigned_charge for c in charge.component_diagnostics
        },
        total_charge=sum(a.charge for a in charge.assignments.values()),
        provenance=str(imported.provenance),
        evidence=[result.record_signature],
    )
    diagnostics = unified_assignment_diagnostics(
        typing_complete=True,
        charges_complete=charge.complete,
        native_charge_available=False,
        parameter_rows_complete=True,
        executable_potential_complete=True,
        periodic_backend_support=not graph.periodic,
        reasons=[],
    )
    assignment = unified_parameter_assignment(
        graph,
        typed,
        charges,
        force_field=force_field,
        source={
            "family": force_field,
            "sha256": imported.source_sha256,
            "source": imported.source,
        },
        native_payload={
            "record": dict(result.record),
            "imported_result_signature": imported.result_signature,
        },
        diagnostics=diagnostics,
    )
    return typed, charges, assignment
