"""Facade contracts using explicitly synthetic native records, not scientific data."""

import json
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace

import pytest
from test_oplsaa_parameters import synthetic as opls_fixture  # noqa: F401

# ruff: noqa: F811 -- pytest fixture injection
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.forcefields import (
    AmberToolsOptions,
    ForceFieldRequest,
    ForceFieldRequestError,
    OPLSOptions,
    PCFFOptions,
    PreparedForceFieldError,
    adopt_forcefield,
    create_evaluator,
    prepare_forcefield,
)


@pytest.fixture
def opls(opls_fixture):
    return opls_fixture


def test_lazy_import():
    subprocess.run(
        [
            sys.executable,
            "-c",
            """
import sys
from island.forcefields import ForceFieldRequest, OPLSOptions, PCFFOptions
ForceFieldRequest('oplsaa', OPLSOptions('/not/read/during/inspection'))
ForceFieldRequest('pcff', PCFFOptions('/not/read', (0,0,1), (0,0,1)))
assert not set(('rdkit','foyer','parmed','openmm','scipy')) & sys.modules.keys()
""",
        ],
        check=True,
    )


@pytest.mark.parametrize(
    "family,options",
    [
        ("unknown", OPLSOptions("missing")),
        ("gaff", OPLSOptions("missing")),
        ("gaff", AmberToolsOptions("gaff2", "provided", {})),
        ("pcff", AmberToolsOptions("gaff", "am1bcc")),
        ("oplsaa", PCFFOptions("missing", (0, 0, 1), (0, 0, 1))),
    ],
)
def test_wrong_options_before_execution(family, options):
    with pytest.raises(ForceFieldRequestError):
        ForceFieldRequest(family, options)


def test_explicit_policies_and_owned_request():
    with pytest.raises(TypeError):
        PCFFOptions("missing")
    with pytest.raises(ForceFieldRequestError):
        PCFFOptions("missing", (True, 0, 1), (0, 0, 1))
    charges = {3: 0.0}
    options = AmberToolsOptions("gaff", "provided", charges)
    request = ForceFieldRequest("gaff", options)
    charges[3] = 5
    assert request.options.provided_charges[3] == 0
    with pytest.raises(TypeError):
        AmberToolsOptions("gaff")  # no implicit QM


def agree(system, prepared, direct):
    facade = create_evaluator(system, prepared)
    a, b = facade.evaluate_fresh(), direct.evaluate_fresh()
    assert a == b
    with facade.open_session() as session:
        assert session.evaluate() == a
        assert session.evaluate_fresh() == a
    return facade


def test_opls_adoption_ownership_and_tampering(opls):
    from island.evaluation.oplsaa import OPLSSinglePointEvaluator

    system, source, native = opls
    before = deepcopy(system.to_dict())
    prepared = adopt_forcefield(system, "oplsaa", native, source=source)
    assert prepared.native_result.identity == native.identity
    evaluator = agree(
        system, prepared, OPLSSinglePointEvaluator(system, native, source)
    )
    metadata = prepared.metadata
    metadata["production_validated"] = True
    with pytest.raises(PreparedForceFieldError, match="Contradictory"):
        replace(prepared, _json_text=json.dumps(metadata))
    with pytest.raises(PreparedForceFieldError):
        replace(prepared, _native=replace(native, json_text="{}"))
    assert prepared.metadata["production_validated"] is False
    moved = deepcopy(system)
    moved.coordinates.set(39, (1.3, 0, 0))
    create_evaluator(moved, prepared).evaluate()
    assert system.to_dict() == before
    system.topology.sites[13].mass = 19
    with pytest.raises(PreparedForceFieldError, match="binding"):
        create_evaluator(system, prepared)
    assert evaluator.evaluate().potential_energy == pytest.approx(2)
    prepared.validate_integrity()  # caller mutation does not change accepted data


@pytest.mark.parametrize(
    "field,value", [("element", "F"), ("atomic_number", 9), ("mass", 19)]
)
def test_wrong_binding(opls, field, value):
    system, source, native = opls
    prepared = adopt_forcefield(system, "oplsaa", native, source=source)
    changed = deepcopy(system)
    setattr(changed.topology.sites[13], field, value)
    with pytest.raises(PreparedForceFieldError):
        create_evaluator(changed, prepared)


def test_opls_dispatch_and_native_failure(opls, monkeypatch):
    import island.forcefields.oplsaa as backend
    from island.exceptions import OPLSAssignmentError

    system, source, native = opls
    calls = []
    monkeypatch.setattr(backend, "load_oplsaa_source", lambda path: source)

    def parameterize(s, src):
        calls.append((s, src))
        return native

    monkeypatch.setattr(backend, "parameterize_oplsaa", parameterize)
    result = prepare_forcefield(
        system, ForceFieldRequest("oplsaa", OPLSOptions("explicit.xml"))
    )
    assert result.native_result.identity == native.identity and len(calls) == 1

    def fail(*args):
        raise OPLSAssignmentError("native-charge residual; sites 13,39")

    monkeypatch.setattr(backend, "parameterize_oplsaa", fail)
    with pytest.raises(OPLSAssignmentError, match="13,39"):
        prepare_forcefield(
            system, ForceFieldRequest("oplsaa", OPLSOptions("explicit.xml"))
        )


def test_missing_source(opls, tmp_path):
    from island.exceptions import OPLSAssignmentError

    with pytest.raises(OPLSAssignmentError) as error:
        prepare_forcefield(
            opls[0], ForceFieldRequest("oplsaa", OPLSOptions(tmp_path / "absent"))
        )
    assert error.value.__cause__ is not None
    path = tmp_path / "wrong.xml"
    path.write_text("<ForceField/>")
    with pytest.raises(OPLSAssignmentError, match="hash"):
        prepare_forcefield(opls[0], ForceFieldRequest("oplsaa", OPLSOptions(path)))


def test_pcff_adoption_and_fresh_dispatch(parameters, monkeypatch):
    import island.forcefields.pcff as backend
    from island.evaluation import PCFFSinglePointEvaluator
    from island.forcefields.pcff.automatic import graph_system

    # Coordinates are supplied for the synthetic authoritative graph.
    system = graph_system(
        parameters.payload["charge_record"]["automatic_typing"]["graph"]
    )
    from island import Coordinates

    system.coordinates = Coordinates(
        {
            i: (float(n % 3), float(n // 3), 0.2 * (n % 2))
            for n, i in enumerate(sorted(system.topology.sites))
        }
    )
    model = backend.define_pcff_model(
        parameters,
        special_pairs=backend.special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    prepared = adopt_forcefield(system, "pcff", model)
    assert prepared.native_result.identity == model.identity
    # Compare bound identities here; general Class II numeric geometry tests are retained.
    assert (
        create_evaluator(system, prepared).model_fingerprint
        == PCFFSinglePointEvaluator(system, model).model_fingerprint
    )
    monkeypatch.setattr(backend, "load_pcff_source", lambda _: parameters.source)
    fresh = prepare_forcefield(
        system,
        ForceFieldRequest("pcff", PCFFOptions("explicit.frc", (0, 0, 1), (0, 0, 1))),
    )
    assert fresh.identity == prepared.identity
    with pytest.raises(PreparedForceFieldError):
        adopt_forcefield(system, "oplsaa", model)
    with pytest.raises(PreparedForceFieldError, match="Contradictory"):
        replace(
            prepared,
            _native=backend.define_pcff_model(
                parameters,
                special_pairs=backend.special_pair_policy(
                    lj=(0, 0, 0.5), coulomb=(0, 0, 1)
                ),
            ),
        )


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
@pytest.mark.parametrize("method", ["provided", "am1bcc"])
def test_amber_routing_and_adoption(tmp_path, monkeypatch, family, method):
    from test_ambertools_integrity import _valid_result

    from island.evaluation import OpenMMSinglePointEvaluator
    from island.forcefields.ambertools import AmberToolsParameterizationEngine

    system, native = _valid_result(tmp_path, method)
    if family == "gaff":
        from island.forcefields.ambertools.models import digest

        inner = deepcopy(
            dict(native.imported_result.provenance["ambertools_preparation"])
        )
        inner["requested_force_field"] = "gaff"
        inner["force_field_data"]["path"] = inner["force_field_data"]["path"].replace(
            "gaff2", "gaff"
        )
        inner["leaprc"]["path"] = inner["leaprc"]["path"].replace("gaff2", "gaff")
        cmd = inner["stages"][0]["command"]
        cmd[cmd.index("-at") + 1] = "gaff"
        cmd = inner["stages"][1]["command"]
        cmd[cmd.index("-s") + 1] = "1"
        imported = replace(
            native.imported_result,
            provenance={
                **dict(native.imported_result.provenance),
                "force_field": "gaff",
                "ambertools_preparation": inner,
            },
            result_signature="",
        )
        imported = replace(imported, result_signature=imported.content_signature())
        outer = {**inner, "imported_result_signature": imported.result_signature}
        native = replace(
            native,
            imported_result=imported,
            record=outer,
            record_signature=digest(outer),
        )
    calls = []

    def run(self, s, options):
        calls.append(options)
        return native

    monkeypatch.setattr(AmberToolsParameterizationEngine, "parameterize", run)
    options = AmberToolsOptions(
        family,
        method,
        {i: 0.0 for i in system.topology.sites} if method == "provided" else None,
    )
    prepared = prepare_forcefield(system, ForceFieldRequest(family, options))
    assert (
        calls == [options]
        and prepared.native_result.record_signature == native.record_signature
    )
    agree(system, prepared, OpenMMSinglePointEvaluator(system, native.imported_result))
    moved = deepcopy(system)
    moved.coordinates.set(min(system.topology.sites), (0.01, 0, 0))
    create_evaluator(moved, prepared)
    with pytest.raises(PreparedForceFieldError, match="family"):
        adopt_forcefield(system, "gaff" if family == "gaff2" else "gaff2", native)


@pytest.mark.parametrize("change", ["stereo", "bond", "site"])
def test_graph_stereo_and_inventory_rejected(opls, change):
    system, source, native = opls
    prepared = adopt_forcefield(system, "oplsaa", native, source=source)
    changed = deepcopy(system)
    if change == "stereo":
        changed.topology.sites[13].metadata["cip_label"] = "R"
    elif change == "bond":
        changed.topology.bonds.clear()
    else:
        changed.topology.sites.pop(39)
    with pytest.raises(PreparedForceFieldError):
        create_evaluator(changed, prepared)


def test_changed_source_and_native_family(opls):
    system, source, native = opls
    with pytest.raises(PreparedForceFieldError):
        adopt_forcefield(
            system, "oplsaa", native, source=replace(source, xml=b"<bad/>")
        )
    for family in ("pcff", "gaff", "gaff2"):
        with pytest.raises(PreparedForceFieldError):
            adopt_forcefield(system, family, native)


def test_amber_unexpected_method_rejected(tmp_path, monkeypatch):
    from test_ambertools_integrity import _valid_result

    from island.forcefields.ambertools import AmberToolsParameterizationEngine

    system, native = _valid_result(tmp_path, "am1bcc")
    monkeypatch.setattr(
        AmberToolsParameterizationEngine, "parameterize", lambda *a: native
    )
    with pytest.raises(PreparedForceFieldError, match="method"):
        prepare_forcefield(
            system,
            ForceFieldRequest(
                "gaff2",
                AmberToolsOptions(
                    "gaff2", "provided", {i: 0 for i in system.topology.sites}
                ),
            ),
        )


def test_unsupported_pcff_preserves_native_diagnostic(synthetic, monkeypatch):
    from test_pcff_automatic import explicit

    import island.forcefields.pcff as backend
    from island.exceptions import PCFFError

    system = explicit([(0, 1)], [3, 3])
    system.topology.sites[min(system.topology.sites)].formal_charge = 1
    monkeypatch.setattr(backend, "load_pcff_source", lambda _: synthetic)
    with pytest.raises(PCFFError, match="Incomplete|unsupported|Unsupported"):
        prepare_forcefield(
            system,
            ForceFieldRequest("pcff", PCFFOptions("source", (0, 0, 1), (0, 0, 1))),
        )


def test_reconstructed_request_revalidates_before_source_access(monkeypatch):
    import island.forcefields.pcff as backend

    def forbidden(*args):
        pytest.fail("invalid request reached source access")

    monkeypatch.setattr(backend, "load_pcff_source", forbidden)
    options = PCFFOptions("not-read", (0, 0, 1), (0, 0, 1))
    object.__setattr__(options, "lj", (0, 0, True))
    with pytest.raises(ForceFieldRequestError):
        ForceFieldRequest("pcff", options)
