"""Lossless chemical inputs and historically signed Amber reconstruction."""

from collections.abc import Mapping
from dataclasses import fields, is_dataclass, replace

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics.thermal import VelocityInitialization
from island.evaluation import EvaluationResult
from island.exceptions import WorkflowError
from island.forcefields import import_amber_prmtop
from island.forcefields.ambertools.models import AmberToolsPreparationResult
from island.minimization import MinimizationOptions, MinimizationResult
from island.minimization.models import MinimizationStep


def record(value):
    """Data-only dataclass fields; no type names or arbitrary-object loading."""
    if is_dataclass(value):
        return {f.name: record(getattr(value, f.name)) for f in fields(value) if f.init}
    if isinstance(value, Mapping):
        return {k: record(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return (
            tuple(record(v) for v in value)
            if isinstance(value, tuple)
            else [record(v) for v in value]
        )
    return value


def system_data(system):
    payload = system.to_dict()
    if system.topology.impropers:
        payload["impropers"] = [record(item) for item in system.topology.impropers]
    return payload


def system_from(payload):
    try:
        if set(payload) - {"impropers"} != {
            "representation",
            "box",
            "sites",
            "bonds",
            "coordinates",
            "metadata",
        }:
            raise ValueError("Unexpected system fields")
        if payload["representation"] != "atomistic" or payload["box"] is not None:
            raise ValueError("Only nonperiodic atomistic systems supported")
        graph = Topology()
        for row in payload["sites"]:
            if row["kind"] != "AtomSite" or type(row["id"]) is not int:
                raise ValueError("Expected integer atom site IDs")
            graph.add_site(AtomSite(**{k: v for k, v in row.items() if k != "kind"}))
        for row in payload["bonds"]:
            graph.add_bond(**row)
        if payload.get("impropers"):
            from island.core.topology import Improper

            graph.impropers = [Improper(**row) for row in payload["impropers"]]
        system = MolecularSystem(
            graph, Coordinates(payload["coordinates"]), metadata=payload["metadata"]
        )
        system.validate()
        return system
    except Exception as error:
        raise WorkflowError(f"Invalid saved chemical system: {error}") from error


def preparation_from(system, payload, prmtop):
    """Reparse source, then verify ORIGINAL signatures; never re-sign changes."""
    r = payload["record"]
    imported = import_amber_prmtop(
        system,
        prmtop,
        r["lineage"]["prmtop_index_to_site_id"],
        source=payload["import_source"],
        force_field=r["requested_force_field"],
        charge_method="provided" if r["charge_method"] == "provided" else "AM1-BCC",
        charge_tolerance=r["charge_validation_tolerance_e"],
    )
    imported = replace(
        imported,
        provenance=payload["import_provenance"],
        result_signature=r["imported_result_signature"],
    )
    result = AmberToolsPreparationResult(imported, r, payload["record_signature"])
    result.validate_integrity(system)
    return result


def minimum_from(payload):
    p = dict(payload)
    for name in ("initial_evaluation", "final_evaluation"):
        if p[name] is not None:
            p[name] = EvaluationResult(**p[name])
    p["options"] = MinimizationOptions(**p["options"])
    p["history"] = tuple(MinimizationStep(**r) for r in p["history"])
    result = MinimizationResult(**p)
    result.validate_integrity()
    if not result.converged or not result.final_evaluation_verified:
        raise WorkflowError("Saved minimum is not independently verified and converged")
    return result


def initialization_from(payload):
    result = VelocityInitialization(**payload)
    result.validate_integrity()
    return result
