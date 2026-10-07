"""Unified record contracts; family-native scientific tests remain separate."""

from math import inf, nan

import pytest

from island.exceptions import ValidationError
from island.forcefields import (
    ForceFieldRequest,
    OPLSOptions,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    assign_pcff_parameters,
    assign_typed_pcff_charges,
    bind_pcff_types,
    load_pcff_source,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
from island.graph import (
    UnifiedGraphCharges,
    UnifiedTypedGraph,
    build_psmiles_graph,
    unified_force_field_source,
    unified_graph_charges,
    unified_typed_graph,
    unify_pcff_assignment,
    unify_pcff_charges,
    unify_pcff_typed_graph,
)
from island.graph.final import pack

SOURCE = "../lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"


def records():
    system, graph, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)
    source = load_pcff_source(SOURCE)
    labels = {
        i: {"C": "c", "H": "h"}[a.element] for i, a in system.topology.sites.items()
    }
    old_typed = bind_pcff_types(
        system,
        source,
        labels,
        provenance="J19 fixture",
        evidence_references=["fixture"],
    )
    old_charges = assign_typed_pcff_charges(
        system, old_typed, resolution_policy=COMPATIBILITY_POLICY
    )
    typed = unify_pcff_typed_graph(graph, old_typed)
    charges = unify_pcff_charges(graph, typed, old_charges)
    assignment = assign_pcff_parameters(
        system, old_typed, old_charges, resolution_policy=COMPATIBILITY_POLICY
    )
    unified = unify_pcff_assignment(graph, typed, charges, assignment)
    return system, graph, source, old_typed, old_charges, typed, charges, unified


def test_pcff_unified_native_and_parameter_diagnostics():
    system, graph, _source, _old_typed, _old_charges, typed, charges, assignment = (
        records()
    )
    assert typed.payload["force_field"] == "PCFF"
    assert charges.payload["charge_origin"] == "native_increments"
    assert charges.payload.get("native_charge_available", True)
    assert assignment.payload["diagnostics"]["typing_complete"]
    assert assignment.payload["diagnostics"]["charges_complete"]
    assert assignment.payload["diagnostics"]["executable_potential_complete"] is False
    # The common request path consumes the unified records through named PCFF adapters.
    options = PCFFOptions(
        SOURCE,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        resolution_policy=COMPATIBILITY_POLICY,
    )
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff", options, final_graph=graph, typed_graph=typed, graph_charges=charges
        ),
    )
    assert prepared.metadata["production_validated"] is False


def test_provided_charge_origin_and_exact_vector():
    _system, graph, source, _ot, _oq, typed, _charges, _assignment = records()
    values = {i: 0.0 for i in typed.atom_types}
    provided = unified_graph_charges(
        graph,
        typed,
        force_field="PCFF",
        charges=values,
        charge_method="independent_vector",
        charge_origin="provided",
        source=source.identity,
        component_totals={str(c[0]): 0.0 for c in graph.payload["components"]},
        total_charge=0.0,
        provenance="independent software fixture",
        evidence=["fixture"],
    )
    assert provided.charges == values
    assert provided.payload["charge_origin"] == "provided"
    assert provided.payload["charge_method"] == "independent_vector"


@pytest.mark.parametrize(
    "mutation", ["id", "label", "charge", "unit", "source", "family", "graph"]
)
def test_unified_mutations_reject(mutation):
    _system, graph, _source, _ot, _oq, typed, charges, _assignment = records()
    if mutation == "id":
        p = typed.payload
        p["atom_types"]["999"] = p["atom_types"].pop("1")
        with pytest.raises(ValidationError):
            UnifiedTypedGraph(pack(p)).validate_integrity(graph)
    elif mutation == "label":
        p = typed.payload
        p["atom_types"]["1"] = "changed"
        with pytest.raises(ValidationError):
            UnifiedTypedGraph(pack(p)).validate_integrity(graph)
    elif mutation == "charge":
        p = charges.payload
        p["charges"]["1"] = 0.1
        with pytest.raises(ValidationError):
            UnifiedGraphCharges(pack(p)).validate_integrity(graph, typed)
    elif mutation == "unit":
        p = charges.payload
        p["unit"] = "coulomb"
        with pytest.raises(ValidationError):
            UnifiedGraphCharges(pack(p)).validate_integrity(graph, typed)
    elif mutation == "source":
        p = typed.payload
        p["source"]["sha256"] = "0" * 64
        with pytest.raises(ValidationError):
            UnifiedTypedGraph(pack(p)).validate_integrity(graph)
    elif mutation == "family":
        p = typed.payload
        p["force_field"] = "GAFF"
        with pytest.raises(ValidationError):
            UnifiedTypedGraph(pack(p)).validate_integrity(graph)
    else:
        p = typed.payload
        p["graph_identity"] = "0" * 64
        with pytest.raises(ValidationError):
            UnifiedTypedGraph(pack(p)).validate_integrity(graph)


@pytest.mark.parametrize("value", [nan, inf, -inf])
def test_unified_nonfinite_charges_reject(value):
    _system, graph, source, _ot, _oq, typed, _charges, _assignment = records()
    values = {i: 0.0 for i in typed.atom_types}
    values[1] = value
    with pytest.raises(ValidationError, match="finite"):
        unified_graph_charges(
            graph,
            typed,
            force_field="PCFF",
            charges=values,
            charge_method="provided",
            charge_origin="provided",
            source=source.identity,
            component_totals={str(c[0]): 0.0 for c in graph.payload["components"]},
            total_charge=0.0,
            provenance="fixture",
            evidence=["fixture"],
        )


def test_family_source_and_graph_identities_are_separate():
    _system, graph, source, _ot, _oq, typed, charges, _assignment = records()
    ff_source = unified_force_field_source(
        force_field="PCFF",
        source=source.identity,
        provenance="Pinned FRC",
        evidence=["pcff.pin.json"],
    )
    assert ff_source.payload["identity"] != typed.identity
    assert typed.payload["source"]["sha256"] == ff_source.payload["source"]["sha256"]
    assert typed.payload["graph_identity"] == graph.identity
    assert charges.payload["typed_graph_identity"] == typed.identity


def test_no_cross_family_fallback_from_unified_records():
    _system, graph, _source, _ot, _oq, typed, _charges, _assignment = records()
    wrong = unified_typed_graph(
        graph,
        force_field="OPLS-AA",
        atom_types=typed.atom_types,
        source={"family": "OPLS-AA", "sha256": "b" * 64},
        typing_profile="explicit_fixture",
        typing_method="external",
        provenance="fixture",
        evidence=["fixture"],
    )
    wrong_charges = unified_graph_charges(
        graph,
        wrong,
        force_field="OPLS-AA",
        charges={i: 0.0 for i in wrong.atom_types},
        charge_method="provided",
        charge_origin="provided",
        source={"family": "OPLS-AA", "sha256": "b" * 64},
        component_totals={str(c[0]): 0.0 for c in graph.payload["components"]},
        total_charge=0.0,
        provenance="fixture",
        evidence=["fixture"],
    )
    with pytest.raises(
        PreparedForceFieldError, match="oplsaa unified assignment adapter"
    ):
        prepare_forcefield(
            _system,
            ForceFieldRequest(
                "oplsaa",
                OPLSOptions("missing.xml"),
                final_graph=graph,
                typed_graph=wrong,
                graph_charges=wrong_charges,
            ),
        )


@pytest.mark.parametrize("family", ["PCFF", "OPLS-AA", "GAFF", "GAFF2"])
def test_source_family_binding_rejects_mismatch(family):
    wrong = "GAFF" if family == "PCFF" else "PCFF"
    with pytest.raises(ValidationError, match="family"):
        unified_force_field_source(
            force_field=family,
            source={"family": wrong, "sha256": "a" * 64},
            provenance="fixture",
            evidence=["fixture"],
        )


@pytest.mark.parametrize("value", ["a" * 63, "a" * 65, "g" * 64, "A" * 64, 7, None])
def test_source_hash_format_rejects(value):
    with pytest.raises(ValidationError, match="SHA"):
        unified_force_field_source(
            force_field="PCFF",
            source={"family": "PCFF", "sha256": value},
            provenance="fixture",
            evidence=["fixture"],
        )


def test_recomputed_outer_identity_cannot_hide_nested_graph_or_source_mutations():
    _system, graph, source, _ot, _oq, typed, charges, assignment = records()
    typed_payload = typed.payload
    typed_payload["component_identity"] = "0" * 64
    typed_payload["identity"] = __import__(
        "island.graph.unified", fromlist=["_record_identity"]
    )._record_identity(typed_payload)
    with pytest.raises(ValidationError, match="component"):
        UnifiedTypedGraph(pack(typed_payload)).validate_integrity(graph)
    typed_payload = typed.payload
    typed_payload["molecule_identity"] = "0" * 64
    typed_payload["identity"] = __import__(
        "island.graph.unified", fromlist=["_record_identity"]
    )._record_identity(typed_payload)
    with pytest.raises(ValidationError, match="molecule"):
        UnifiedTypedGraph(pack(typed_payload)).validate_integrity(graph)
    assignment_payload = assignment.payload
    assignment_payload["source"] = {"family": "GAFF", "sha256": "b" * 64}
    assignment_payload["identity"] = __import__(
        "island.graph.unified", fromlist=["_record_identity"]
    )._record_identity(assignment_payload)
    from island.graph import UnifiedParameterAssignment

    with pytest.raises(ValidationError, match="family/source"):
        UnifiedParameterAssignment(pack(assignment_payload)).validate_integrity(
            graph, typed, charges
        )
    assert source.identity["sha256"] == typed.payload["source"]["sha256"]


def test_pcff_common_path_rejects_unified_source_hash_not_loaded_frc():
    system, graph, _source, _ot, _oq, typed, charges, _assignment = records()
    from island.graph import UnifiedGraphCharges, UnifiedTypedGraph
    from island.graph.unified import _record_identity

    tp = typed.payload
    cp = charges.payload
    tp["source"]["sha256"] = "b" * 64
    cp["source"]["sha256"] = "b" * 64
    tp["identity"] = _record_identity(tp)
    cp["typed_graph_identity"] = tp["identity"]
    cp["identity"] = _record_identity(cp)
    bad_typed = UnifiedTypedGraph(pack(tp))
    bad_charges = UnifiedGraphCharges(pack(cp))
    options = PCFFOptions(
        SOURCE,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        resolution_policy=COMPATIBILITY_POLICY,
    )
    with pytest.raises(Exception, match="source hash"):
        prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                options,
                final_graph=graph,
                typed_graph=bad_typed,
                graph_charges=bad_charges,
            ),
        )


@pytest.fixture(scope="module")
def public_request_case():
    system, graph, _source, _ot, _oq, typed, charges, _assignment = records()
    options = PCFFOptions(
        SOURCE,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        resolution_policy=COMPATIBILITY_POLICY,
    )
    return system, graph, options, typed, charges


def test_unified_public_request_requires_final_graph(public_request_case):
    from island.forcefields import ForceFieldRequestError

    system, graph, options, typed, charges = public_request_case
    with pytest.raises(ForceFieldRequestError, match="require final_graph"):
        ForceFieldRequest("pcff", options, typed_graph=typed, graph_charges=charges)
    request = ForceFieldRequest("pcff", options, graph, typed, charges)
    # A deserialized or tampered frozen request must be checked again.
    object.__setattr__(request, "final_graph", None)
    with pytest.raises(ForceFieldRequestError, match="require final_graph"):
        prepare_forcefield(system, request)


def tampered_public_records(typed, charges, defect):
    from island.graph.unified import _record_identity

    tp, cp = typed.payload, charges.payload
    if defect in ("component_identity", "molecule_identity"):
        tp[defect] = "0" * 64
    elif defect == "component_totals":
        key = next(iter(cp["component_totals"]))
        cp["component_totals"][key] += 1.0
    elif defect == "component_keys":
        cp["component_totals"]["999"] = cp["component_totals"].pop(
            next(iter(cp["component_totals"]))
        )
    else:
        cp["source"]["sha256"] = "b" * 64
    tp["identity"] = _record_identity(tp)
    cp["typed_graph_identity"] = tp["identity"]
    cp["identity"] = _record_identity(cp)
    return UnifiedTypedGraph(pack(tp)), UnifiedGraphCharges(pack(cp))


@pytest.mark.parametrize(
    "defect",
    [
        "component_identity",
        "molecule_identity",
        "component_totals",
        "component_keys",
        "source",
    ],
)
def test_public_request_and_preparation_reject_rechecksummed_records(
    public_request_case,
    monkeypatch,
    defect,
):
    from island.forcefields import pcff

    system, graph, options, typed, charges = public_request_case
    bad_typed, bad_charges = tampered_public_records(typed, charges, defect)
    # Both envelopes and standalone identities remain valid.
    bad_typed.validate_integrity()
    bad_charges.validate_integrity()
    with pytest.raises(ValidationError) as direct:
        bad_typed.validate_integrity(graph)
        bad_charges.validate_integrity(graph, bad_typed)
    with pytest.raises(ValidationError) as public:
        ForceFieldRequest("pcff", options, graph, bad_typed, bad_charges)
    assert str(public.value) == str(direct.value)

    def forbidden(*args, **kwargs):
        pytest.fail("Tampered unified records reached legacy PCFF adaptation")

    monkeypatch.setattr(pcff, "bind_pcff_types", forbidden)
    request = ForceFieldRequest("pcff", options, graph, typed, charges)
    object.__setattr__(request, "typed_graph", bad_typed)
    object.__setattr__(request, "graph_charges", bad_charges)
    with pytest.raises(ValidationError) as preparation:
        prepare_forcefield(system, request)
    assert str(preparation.value) == str(direct.value)


@pytest.mark.parametrize("defect", ["molecule_identity", "component_totals", "source"])
def test_pcff_revalidates_after_request_before_legacy_adaptation(
    public_request_case,
    monkeypatch,
    defect,
):
    from island.forcefields import pcff

    system, graph, options, typed, charges = public_request_case
    # Own records so an injected late mutation cannot modify the shared fixture.
    typed = UnifiedTypedGraph(typed.json_text)
    charges = UnifiedGraphCharges(charges.json_text)
    request = ForceFieldRequest("pcff", options, graph, typed, charges)
    bad_typed, bad_charges = tampered_public_records(typed, charges, defect)
    load = pcff.load_pcff_source

    def load_and_mutate(*args, **kwargs):
        source = load(*args, **kwargs)
        object.__setattr__(typed, "json_text", bad_typed.json_text)
        object.__setattr__(charges, "json_text", bad_charges.json_text)
        return source

    def forbidden(*args, **kwargs):
        pytest.fail("Adapter-boundary validation was skipped")

    monkeypatch.setattr(pcff, "load_pcff_source", load_and_mutate)
    monkeypatch.setattr(pcff, "bind_pcff_types", forbidden)
    with pytest.raises(ValidationError):
        prepare_forcefield(system, request)
