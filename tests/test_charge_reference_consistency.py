"""Rechecksummed contradictory evidence and injected matrix orchestration."""

from copy import deepcopy
from types import SimpleNamespace

import pytest
from test_oligomer_charge_references import synthetic_reference

from island.charge_references import ChargeReference, load_charge_reference, save_record
from island.charge_references.records import pack
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.models import digest
from island.workflows import storage
from island.workflows.bundle import system_from
from scripts import validate_oligomer_charges as cli


def resign(p):
    r = p["preparation"]["record"]
    p["import_content"]["provenance"]["ambertools_preparation"] = deepcopy(
        {k: v for k, v in r.items() if k != "imported_result_signature"}
    )
    r["imported_result_signature"] = digest(p["import_content"])
    p["preparation"]["record_signature"] = digest(r)
    return ChargeReference(pack(p))


@pytest.mark.parametrize(
    "mutation",
    [
        lambda c: c.__setitem__(c.index("-c") + 1, "rc"),
        lambda c: c.extend(["-cf", "charges.txt"]),
        lambda c: c.extend(["-c", "bcc"]),
        lambda c: c.extend(["-at", "gaff"]),
        lambda c: c.__delitem__(slice(c.index("-fi"), c.index("-fi") + 2)),
    ],
)
def test_contradictory_command_even_with_all_checksums_recomputed(tmp_path, mutation):
    p = synthetic_reference().payload
    mutation(p["preparation"]["record"]["stages"][0]["command"])
    with pytest.raises(ChargeReferenceError):
        save_record(resign(p), tmp_path / "bad.json")
    assert not (tmp_path / "bad.json").exists()
    storage.publish(tmp_path / "bad.json", resign(p).json_text.encode())
    with pytest.raises(ChargeReferenceError):
        load_charge_reference(tmp_path / "bad.json")


@pytest.mark.parametrize("change", ["tolerance", "sites", "residual"])
def test_matrix_declared_settings_and_metrics(tmp_path, monkeypatch, change):
    ref = (
        synthetic_reference(tolerance=0.01, residual=0.004)
        if change == "tolerance"
        else synthetic_reference()
    )
    p = ref.payload
    system = system_from(p["system"])
    row = dict(
        cli.MATRIX[0],
        case="pe-dp3-seed2026",
        status="passed",
        sites=system.number_of_sites,
        residual=p["charge_result"]["total_charge_residual"],
        record_signature=p["preparation"]["record_signature"],
        import_signature=p["preparation"]["record"]["imported_result_signature"],
        elapsed_seconds=1.0,
    )
    if change == "sites":
        row["sites"] += 1
    if change == "residual":
        row["residual"] += 0.001
    source = tmp_path / "source"
    source.mkdir()
    output = tmp_path / "out"
    output.mkdir()
    storage.publish(source / "declared-settings.json", storage.json_bytes(cli.SETTINGS))
    storage.publish(
        source / "matrix.json", storage.json_bytes({"cases": [row], "complete": False})
    )
    prep = SimpleNamespace(
        record=p["preparation"]["record"],
        record_signature=row["record_signature"],
        imported_result=SimpleNamespace(result_signature=row["import_signature"]),
    )
    monkeypatch.setattr(cli, "reconstruct", lambda *_: (system, prep))
    monkeypatch.setattr(cli, "create_charge_reference", lambda *_: ref)
    assert cli.audit_existing(source, output) is False
    outcomes = storage.read_json(output / "outcomes.json")
    assert outcomes["validated_reference_count"] == 0
    assert outcomes["cases"][0]["status"] == "reference_validation_failed"
    assert not (output / (row["case"] + ".json")).exists()


@pytest.mark.parametrize("complete", [False, True])
def test_valid_partial_and_complete_matrix_gate(tmp_path, monkeypatch, complete):
    """Synthetic signed evidence and injected reconstruction, not live QM."""
    source = tmp_path / "source"
    source.mkdir()
    refs = {}
    rows = []
    for c in cli.MATRIX[: 8 if complete else 3]:
        name = f"{c['chemistry']}-dp{c['dp']}-seed{c['seed']}"
        ref = synthetic_reference(c["dp"], c["seed"], oxygen=c["chemistry"] == "peo")
        p = ref.payload
        refs[name] = ref
        rows.append(
            dict(
                c,
                case=name,
                status="passed",
                sites=len(p["charge_result"]["assignments"]),
                residual=p["charge_result"]["total_charge_residual"],
                record_signature=p["preparation"]["record_signature"],
                import_signature=p["preparation"]["record"][
                    "imported_result_signature"
                ],
                elapsed_seconds=1.0,
            )
        )

    def reconstruct(path):
        p = refs[path.name].payload
        return system_from(p["system"]), SimpleNamespace(
            reference=refs[path.name],
            record=p["preparation"]["record"],
            record_signature=p["preparation"]["record_signature"],
            imported_result=SimpleNamespace(
                result_signature=p["preparation"]["record"]["imported_result_signature"]
            ),
        )

    monkeypatch.setattr(cli, "reconstruct", reconstruct)
    monkeypatch.setattr(
        cli,
        "create_charge_reference",
        lambda s, p: p.reference,
    )
    storage.publish(source / "declared-settings.json", storage.json_bytes(cli.SETTINGS))
    storage.publish(
        source / "matrix.json",
        storage.json_bytes({"cases": rows, "complete": complete}),
    )
    assert cli.main(
        ["--audit-existing", str(source), "--output", str(tmp_path / "out")]
    ) == (0 if complete else 1)
    result = storage.read_json(tmp_path / "out/outcomes.json")
    assert result["validated_reference_count"] == len(rows)
    assert result["complete"] == complete


@pytest.mark.parametrize(
    "field,value",
    [
        ("status", "unknown"),
        ("dp", True),
        ("sites", False),
        ("residual", "0.0"),
        ("elapsed_seconds", -1),
    ],
)
def test_matrix_rejects_malformed_row(tmp_path, field, value):
    source = tmp_path / "source"
    source.mkdir()
    row = dict(
        cli.MATRIX[0],
        case="pe-dp3-seed2026",
        status="failed",
        failure="preserved failure",
    )
    row[field] = value
    storage.publish(source / "declared-settings.json", storage.json_bytes(cli.SETTINGS))
    storage.publish(
        source / "matrix.json", storage.json_bytes({"cases": [row], "complete": False})
    )
    assert (
        cli.main(["--audit-existing", str(source), "--output", str(tmp_path / "out")])
        == 1
    )
    assert (tmp_path / "out/failure.json").exists()
