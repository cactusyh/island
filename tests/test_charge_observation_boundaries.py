"""Rechecksummed synthetic evidence: semantic boundaries, not QM acceptance."""

import json

import pytest
from test_charge_conservation import synthetic_observation

from island.charge_references import (
    RawChargeObservation,
    load_raw_observation,
    project_charge_observation,
    save_record,
)
from island.charge_references.observations import _observation_data
from island.charge_references.projection import POLICY
from island.charge_references.records import pack
from island.exceptions import ChargeReferenceError
from island.workflows import storage


def changed(output=None, identities=()):
    p = synthetic_observation().payload
    e = p["evidence"]
    m = json.loads(e["manifest_json"])
    row = m["cases"][0]
    if output is not None:
        path = "island-ambertools-test/sqm.out"
        e["files"][path] = output(e["files"][path])
        row["source_checksums"][path] = storage.checksum(e["files"][path].encode())
    for field in identities:
        row[field] = "a" * 64
    e["manifest_json"] = storage.json_bytes(m).decode()
    # Update all affected derived fields exactly as the reviewed implementation
    # did. The rejection must concern evidence semantics, not stale digests.
    p["data"]["artifact_sha256"] = row["source_checksums"]
    p["data"]["source_manifest_sha256"] = storage.checksum(e["manifest_json"].encode())
    for field, target in [
        ("reference_identity", "historical_reference_identity"),
        ("record_signature", "preparation_signature"),
        ("import_signature", "import_signature"),
    ]:
        p["data"][target] = row.get(field)
    return p


@pytest.mark.parametrize(
    "transform",
    [
        lambda s: s.replace("Final Structure", "Initial Structure"),
        lambda s: s.replace("Final Structure\n", ""),
        lambda s: s + s,
    ],
)
def test_final_evidence_required(transform):
    p = changed(output=transform)
    with pytest.raises(ChargeReferenceError):
        _observation_data(p["evidence"])


@pytest.mark.parametrize(
    "fields",
    [
        ("reference_identity",),
        ("record_signature",),
        ("import_signature",),
        ("reference_identity", "record_signature", "import_signature"),
    ],
)
def test_failed_identities_rejected(fields, tmp_path):
    p = changed(identities=fields)
    record = RawChargeObservation(pack(p))
    with pytest.raises(ChargeReferenceError, match="successful"):
        record.validate_integrity()
    with pytest.raises(ChargeReferenceError):
        project_charge_observation(record, policy=POLICY)
    with pytest.raises(ChargeReferenceError):
        save_record(record, tmp_path / "unpublished.json")
    assert not (tmp_path / "unpublished.json").exists()
    storage.publish(tmp_path / "forged.json", pack(p).encode())
    with pytest.raises(ChargeReferenceError):
        load_raw_observation(tmp_path / "forged.json")


@pytest.mark.parametrize(
    "kind", ["truncated", "duplicate", "index", "element", "nonfinite"]
)
def test_invalid_final_table_after_complete_initial(kind):
    def transform(s):
        initial = s.replace("Final Structure", "Initial Structure").replace(
            "Calculation Completed", ""
        )
        lines = s.splitlines()
        if kind == "truncated":
            del lines[-2]
        elif kind == "duplicate":
            lines[2] = lines[1]
        elif kind == "index":
            lines[1] = lines[1].replace("1 1", "1 2")
        elif kind == "element":
            lines[1] = lines[1].replace(" C ", " O ")
        else:
            lines[1] = lines[1].replace("0 0 0", "nan 0 0")
        return initial + "\n" + "\n".join(lines)

    with pytest.raises(ChargeReferenceError):
        _observation_data(changed(output=transform)["evidence"])


def test_initial_then_final_and_null_failed_identities():
    def transform(s):
        return s.replace("Final Structure", "Initial Structure").replace(
            "Calculation Completed", ""
        ) + s.replace("0 0 0", "1 2 3")

    e = changed(output=transform)["evidence"]
    m = json.loads(e["manifest_json"])
    for key in ("reference_identity", "record_signature", "import_signature"):
        m["cases"][0][key] = None
    e["manifest_json"] = storage.json_bytes(m).decode()
    data = _observation_data(e)
    assert set(data["post_sqm_coordinates_angstrom"].values()) == {(1.0, 2.0, 3.0)}
    assert data["historical_status"] == "failed"
    assert data["historical_reference_identity"] is None
    record = RawChargeObservation(
        pack(
            {
                "schema": synthetic_observation().payload["schema"],
                "evidence": e,
                "data": data,
            }
        )
    )
    record.validate_integrity()


@pytest.mark.parametrize("defect", ["geometry", "identity"])
def test_rechecksummed_nested_records_reject_before_publication(defect, tmp_path):
    from island.charge_references import ChargeProjection, ConservationAudit
    from island.charge_references.conservation_audit import SCHEMA

    original = synthetic_observation()
    projection = project_charge_observation(original, policy=POLICY).payload
    p = (
        changed(output=lambda s: s.replace("Final Structure", "Initial Structure"))
        if defect == "geometry"
        else changed(identities=("import_signature",))
    )
    projection["observation"] = p
    projection["observation_identity"] = storage.checksum(
        storage.json_bytes(storage.encode(p))
    )
    records = [
        RawChargeObservation(pack(p)),
        ChargeProjection(pack(projection)),
        ConservationAudit(
            pack({"schema": SCHEMA, "projections": [projection], "audit": {}})
        ),
    ]
    for i, record in enumerate(records):
        with pytest.raises(
            ChargeReferenceError, match="final SQM|successful historical"
        ):
            record.validate_integrity()
        with pytest.raises(ChargeReferenceError):
            save_record(record, tmp_path / f"{i}.json")
        assert not (tmp_path / f"{i}.json").exists()
    storage.publish(tmp_path / "bad.json", pack(p).encode())
    with pytest.raises(ChargeReferenceError, match="final SQM|successful historical"):
        load_raw_observation(tmp_path / "bad.json")
