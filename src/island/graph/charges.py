"""Explicit adapters into the neutral charge contract.

Adapters are named and source-bound; they never mutate or rewrite legacy records.
"""

from .final import force_neutral_charges


def charge_record_from_pcff(graph, result):
    from island.forcefields.pcff.automatic import PCFFAutomaticChargeResult
    from island.forcefields.pcff.typed_graph import PCFFGraphCharges

    if type(result) not in (PCFFGraphCharges, PCFFAutomaticChargeResult):
        raise TypeError("Validated PCFF charge result required")
    result.validate_integrity()
    p = result.payload
    typed = p.get("typed_graph", p.get("automatic_typing"))
    _require_pcff_graph_binding(graph, typed)
    if "charge_data" in p:
        data = p["charge_data"]
        origin = p["origin"]
        source = p["source"]
    else:
        data = p["native_charge_record"]
        origin = "native_increments"
        source = p["source"]
    return force_neutral_charges(
        graph,
        data["partial_charges"],
        force_field="PCFF",
        source=source,
        method="native_increment_resolver"
        if origin == "native_increments"
        else "provided",
        origin=origin,
        provenance=data.get("provenance", "PCFF validated charge result"),
        evidence=data.get("evidence_references", []),
    )


def _require_pcff_graph_binding(graph, typed):
    if typed is None:
        raise ValueError("PCFF record has no bound chemical graph")
    neutral = graph.payload
    if [
        (
            s["id"],
            s["element"],
            s["formal_charge"],
            s["metadata"].get("aromatic", False),
        )
        for s in neutral["sites"]
    ] != [
        (
            s["id"],
            s["element"],
            s["formal_charge"],
            s["metadata"].get("aromatic", False),
        )
        for s in typed["graph"]["sites"]
    ] or neutral["bonds"] != typed["graph"]["bonds"]:
        raise ValueError("PCFF record chemical graph differs from final graph")


def charge_record_from_opls(graph, result, *, source, method="native_source_charge"):
    result.validate_integrity()
    p = result.payload
    charges = {i: row["charge"] for i, row in p["charges"].items()}
    return force_neutral_charges(
        graph,
        charges,
        force_field="OPLS-AA",
        source=source.identity,
        method=method,
        origin="native_source",
        provenance="Validated OPLS native charge record",
        evidence=[result.identity],
    )


def charge_record_from_amber(graph, result, *, force_field, source, method):
    result.validate_integrity(graph)
    charges = {int(i): row["charge"] for i, row in result.record["charges"].items()}
    return force_neutral_charges(
        graph,
        charges,
        force_field=force_field,
        source=source,
        method=method,
        origin="provided" if method == "provided" else "am1bcc",
        provenance="Validated AmberTools preparation record",
        evidence=[result.record_signature],
    )
