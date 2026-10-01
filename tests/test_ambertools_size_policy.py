"""Default preservation, explicit bounds and signed v2/v3 policy semantics."""

import copy
import hashlib
from dataclasses import replace

import pytest

from island.exceptions import (
    AmberToolsInputError,
    InvalidAmberImportResultError,
    WorkflowError,
)
from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine
from island.forcefields.ambertools.lineage import (
    generated_atom_name,
    parse_and_validate_mol2,
    prepare_input,
)
from island.forcefields.ambertools.models import (
    PREPARATION_SCHEMA,
    AmberToolsPreparationResult,
    digest,
)
from island.workflows import WorkflowConfig, inspect_chain


@pytest.mark.parametrize("value", [0, -1, True, 1.0, None, 1001, 10**400])
def test_invalid_size_budgets(value):
    with pytest.raises(AmberToolsInputError):
        AmberToolsOptions("gaff", "provided", {}, max_atoms=value)
    with pytest.raises(WorkflowError):
        WorkflowConfig(
            "[*]CC[*]",
            3,
            "run",
            "gaff",
            "provided",
            provided_charges={},
            max_atoms=value,
        )


def test_default_am1bcc_and_historical_config():
    assert AmberToolsOptions("gaff2", "am1bcc").max_atoms == 100
    assert AmberToolsOptions("gaff2", "provided", {}).max_atoms == 100
    with pytest.raises(AmberToolsInputError):
        AmberToolsOptions("gaff2", "am1bcc", max_atoms=101)
    cfg = WorkflowConfig(
        "[*]CC[*]",
        20,
        "run",
        "gaff2",
        "provided",
        provided_charges={},
        max_atoms=1000,
        charge_source="known file, unverified method",
    )
    assert WorkflowConfig.from_dict(cfg.to_dict()) == cfg
    assert cfg.amber_options().max_atoms == 1000
    assert copy.deepcopy(cfg.amber_options()).charge_source == cfg.charge_source
    legacy = replace(
        cfg, schema="island_single_chain_config_v1", max_atoms=100, charge_source=None
    )
    assert "max_atoms" not in legacy.to_dict()
    assert WorkflowConfig.from_dict(legacy.to_dict()).to_dict() == legacy.to_dict()
    with pytest.raises(WorkflowError):
        replace(legacy, max_atoms=1000)
    with pytest.raises(WorkflowError):
        WorkflowConfig.from_dict({**legacy.to_dict(), "max_atoms": 100})


@pytest.fixture
def large():
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer
    from island.workflows.bundle import system_data, system_from

    built = build_linear_polymer("[*]CC[*]", dp=20, coordinate_method="local_templates")
    data = system_data(built)
    mapping = {s["id"]: 11 * s["id"] + 7 for s in data["sites"]}
    for site in data["sites"]:
        site["id"] = mapping[site["id"]]
    for bond in data["bonds"]:
        bond["site1"], bond["site2"] = mapping[bond["site1"]], mapping[bond["site2"]]
    # Coordinates.to_dict uses stable string keys.
    data["coordinates"] = {
        mapping[int(s)]: xyz for s, xyz in data["coordinates"].items()
    }
    return system_from(data)


def test_above_100_roundtrip_and_preexecution_rejection(large, tmp_path, monkeypatch):
    from island.forcefields.ambertools import engine

    before = copy.deepcopy(large.to_dict())
    charges = {
        s: 0.01 if i % 2 == 0 else -0.01
        for i, s in enumerate(sorted(large.topology.sites))
    }

    def forbidden(*a):
        raise AssertionError("external discovery must not be reached")

    monkeypatch.setattr(engine, "_discover", forbidden)
    for opts in (
        AmberToolsOptions("gaff", "provided", charges),
        AmberToolsOptions("gaff2", "am1bcc"),
    ):
        with pytest.raises(AmberToolsInputError, match="actual count=122.*100 sites"):
            AmberToolsParameterizationEngine().parameterize(large, opts)
    with pytest.raises(AmberToolsInputError):
        prepare_input(large, None, max_atoms=1000)
    prepared = prepare_input(large, charges, max_atoms=1000)
    path = tmp_path / "roundtrip.mol2"
    path.write_text(prepared.mol2_text)
    mapping, values, _ = parse_and_validate_mol2(
        path, large, prepared.names, prepared.expected_cip, stage="test"
    )
    assert set(mapping.values()) == set(charges)
    assert values == charges
    assert large.to_dict() == before
    # Workflow inspection propagates the same size ceiling; no preparation occurs.
    cfg = WorkflowConfig(
        "[*]CC[*]",
        20,
        str(tmp_path),
        "gaff",
        "provided",
        provided_charges={},
        max_atoms=1000,
    )
    assert inspect_chain(cfg).number_of_sites == 122


@pytest.mark.parametrize("element", ["C", "O", "H", "Cl", "Br"])
def test_names_at_entire_1000_site_budget(element):
    names = [generated_atom_name(element, i) for i in range(1000)]
    assert len(set(names)) == 1000 and all(
        len(n) == 4 and n.startswith(element) for n in names
    )
    assert [int(n[len(element) :], 36) for n in names] == list(range(1000))


def upgrade(result, system, max_atoms=1000):
    from island.forcefields.ambertools.policy import policy_record

    r = copy.deepcopy(dict(result.record))
    r.update(
        schema=PREPARATION_SCHEMA,
        engine_version="3",
        size_policy=policy_record(max_atoms, "provided", len(system.topology.sites)),
    )
    charges = {
        str(s): a.charge
        for s, a in result.imported_result.charge_result.assignments.items()
    }
    sha = hashlib.sha256(
        (
            "\n".join(f"{charges[str(s)]:.10f}" for s in sorted(system.topology.sites))
            + "\n"
        ).encode()
    ).hexdigest()
    r["artifact_sha256"]["charges.txt"] = sha
    r["provided_charge_input"] = {
        "charges": charges,
        "sha256": digest(charges),
        "source": "synthetic test",
        "method_provenance": "unverified_user_supplied",
        "charges_file_sha256": sha,
    }
    return resign(result, r)


def resign(result, record):
    shared = {k: v for k, v in record.items() if k != "imported_result_signature"}
    imported = replace(
        result.imported_result,
        provenance={
            **dict(result.imported_result.provenance),
            "ambertools_preparation": shared,
        },
        result_signature="",
    )
    imported = replace(imported, result_signature=imported.content_signature())
    record = {**shared, "imported_result_signature": imported.result_signature}
    return AmberToolsPreparationResult(imported, record, digest(record))


def test_v3_policy_and_charge_integrity(tmp_path):
    from test_ambertools_integrity import _valid_result

    system, historical = _valid_result(tmp_path)
    historical.validate_integrity(system)
    assert historical.schema == "island_ambertools_preparation_v2"
    new = upgrade(historical, system)
    new.validate_integrity(system)
    copy.deepcopy(new).validate_integrity(system)
    for key, value in (
        ("actual_atoms", 999),
        ("charge_method", "am1bcc"),
        ("max_atoms", 3),
        ("max_atoms", True),
    ):
        r = copy.deepcopy(dict(new.record))
        r["size_policy"][key] = value
        with pytest.raises(InvalidAmberImportResultError):
            resign(new, r).validate_integrity(system)
    for change in ("checksum", "charges", "mode"):
        r = copy.deepcopy(dict(new.record))
        if change == "checksum":
            r["provided_charge_input"]["sha256"] = "0" * 64
        elif change == "charges":
            r["provided_charge_input"]["charges"]["101"] = 0.1
            r["provided_charge_input"]["sha256"] = digest(
                r["provided_charge_input"]["charges"]
            )
        else:
            r["charge_method"] = "am1bcc"
        with pytest.raises(InvalidAmberImportResultError):
            resign(new, r).validate_integrity(system)
    # A mere outer rehash cannot contradict the signed inner budget.
    r = copy.deepcopy(dict(new.record))
    r["size_policy"]["max_atoms"] = 999
    with pytest.raises(
        InvalidAmberImportResultError, match="signed imported provenance"
    ):
        replace(new, record=r, record_signature=digest(r)).validate_integrity(system)
    other = upgrade(historical, system, max_atoms=100)
    assert (
        other.imported_result.model_content_signature()
        == new.imported_result.model_content_signature()
    )
    assert (
        other.imported_result.content_signature()
        != new.imported_result.content_signature()
    )


def test_halogen_chain_actual_mol2_roundtrip(tmp_path):
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer

    system = build_linear_polymer(
        "[*]C(Cl)(Br)C[*]", dp=30, coordinate_method="local_templates"
    )
    charges = dict.fromkeys(system.topology.sites, 0.0)
    prepared = prepare_input(system, charges, max_atoms=1000)
    assert system.number_of_sites > 100
    assert any(n.startswith("Cl") for n in prepared.names)
    assert any(n.startswith("Br") and int(n[2:], 36) > 100 for n in prepared.names)
    path = tmp_path / "halogens.mol2"
    path.write_text(prepared.mol2_text)
    mapping, q, _ = parse_and_validate_mol2(
        path, system, prepared.names, prepared.expected_cip, stage="roundtrip"
    )
    assert set(mapping.values()) == set(system.topology.sites) and q == charges


def test_large_input_failure_preserves_caller(large, monkeypatch, tmp_path):
    from types import SimpleNamespace

    from island.exceptions import AmberToolsStageError
    from island.forcefields.ambertools import engine

    before = copy.deepcopy(large.to_dict())
    monkeypatch.setattr(
        engine,
        "_discover",
        lambda options: SimpleNamespace(
            amberhome=tmp_path, executables={"antechamber": "unexecuted-test-command"}
        ),
    )

    def fail(stage, *a, **kw):
        raise AmberToolsStageError(
            stage, "injected timeout", artifact_dir=str(tmp_path)
        )

    monkeypatch.setattr(engine, "_run_stage", fail)
    with pytest.raises(AmberToolsStageError, match="injected timeout"):
        AmberToolsParameterizationEngine().parameterize(
            large,
            AmberToolsOptions(
                "gaff2",
                "provided",
                dict.fromkeys(large.topology.sites, 0.0),
                max_atoms=1000,
                work_root=tmp_path,
            ),
        )
    assert large.to_dict() == before
    assert list(tmp_path.glob("island-ambertools-*/input.mol2"))


def test_v3_prepared_workflow_and_physical_identity(tmp_path):
    pytest.importorskip("openmm")
    from test_ambertools_integrity import _valid_result

    from island.evaluation import OpenMMSinglePointEvaluator
    from island.workflows.chain import _preparation_budget

    system, old = _valid_result(tmp_path)
    a, b = upgrade(old, system, 100), upgrade(old, system, 1000)
    ea, eb = (
        OpenMMSinglePointEvaluator(system, a.imported_result),
        OpenMMSinglePointEvaluator(system, b.imported_result),
    )
    assert ea.model_fingerprint == eb.model_fingerprint
    assert (
        ea.parameter_fingerprint != eb.parameter_fingerprint
    )  # exhaustive audit identity
    assert ea.evaluate().potential_energy == eb.evaluate().potential_energy
    cfg = WorkflowConfig(
        "[*]CC[*]",
        3,
        "unused",
        "gaff2",
        "provided",
        provided_charges={},
        max_atoms=1000,
    )
    _preparation_budget(cfg, system, b.record)
    with pytest.raises(WorkflowError, match="size budget"):
        _preparation_budget(replace(cfg, max_atoms=100), system, b.record)


def test_workflow_oversize_before_any_tool_probe(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    from island.forcefields.ambertools import engine
    from island.workflows import start_workflow

    probes = []

    def discovery(options, *, probe_versions=True):
        probes.append(probe_versions)
        assert probe_versions is False

    monkeypatch.setattr(engine, "_discover", discovery)
    cfg = WorkflowConfig(
        "[*]CC[*]",
        20,
        str(tmp_path / "workflow"),
        "gaff2",
        "provided",
        provided_charges={},
    )
    result = start_workflow(cfg)
    assert probes == [False]
    assert result["status"] == "stage_failed"
    assert result["stages"][-1]["stage"] == "construction"
