"""J89 software controls; the real-source integration case requires explicit pins."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.forcefields.bundle_adapters.intake import _read
from island.forcefields.bundle_adapters.opls_authorized import ALLOWED_FIXTURES
from island.forcefields.bundle_adapters.schema import loads
from scripts import phase_4j89_opls as j89
from tests.test_phase_4j84 import put
from tests.test_phase_4j88 import SPECS, control

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j89.scope()


def software_inputs(path, **kwargs):
    return control(path, specs=SPECS, **kwargs)


def test_complete_authorized_input_contract_all_eight_graphs_software_only(tmp_path):
    external = software_inputs(tmp_path)
    _, reports = j89.reconstruct(REPO, **external)
    r = reports["result.json"]
    assert r["decision"] == j89.READY
    assert (
        r["exact_graph_fixtures"]
        == r["qualified_graph_reports"]
        == r["required_graph_count"]
        == 8
    )
    assert r["source_files"] > 0 and r["source_records_inspected"] > 0
    assert r["parameter_rows_inspected"] > 0
    assert r["audit_sha256"] == external["audit_sha256"]
    assert r["manifest_sha256"] == external["manifest_sha256"]
    assert (
        not r["scope"]["runtime_registration"]
        and r["scope"]["accepted_coefficients"] == 0
    )
    assert (
        r["chemical_family_transfer"]
        == r["scope"]["simulation_readiness"]
        == "not_established"
    )
    audit = reports["audit-references.json"]
    assert audit["reviewer"] and audit["authority_basis"] and audit["records"]
    for graph in reports["graph-coverage-report.json"]["reports"]:
        assert graph["decision"] == j89.READY and not graph["missing_terms"]
        assert set(graph["coverage"]) == {
            "atom_typing",
            "bond",
            "angle",
            "torsion",
            "improper",
            "charge",
            "nonbonded",
            "interaction_14",
            "topology",
            "class_ii_cross_terms",
        }
        assert all(v["covered"] == v["required"] for v in graph["coverage"].values())
        assert all(
            row["audit_reference"]["audit_sha256"] == external["audit_sha256"]
            and row["audit_reference"]["claim_sha256"] == row["claim_sha256"]
            for row in graph["matched_source_records"]
        )


def test_real_user_supplied_authorized_bundle(tmp_path):
    # No discovery, implicit fixture or manufactured replacement for this test.
    path = os.environ.get("ISLAND_J89_REAL_INPUTS")
    pin = os.environ.get("ISLAND_J89_REAL_INPUTS_SHA256")
    if path is None or pin is None:
        pytest.skip(
            "Real authorized user inputs and descriptor SHA-256 not supplied; no scientific coverage tested"
        )
    descriptor = loads(_read(path, pin))
    external = {k: descriptor[k] for n in j89.INPUTS for k in (n, n + "_sha256")}
    output = j89.emit(REPO, tmp_path / "real-source-metadata", **external)
    assert output["result"]["decision"] == j89.READY
    assert output["result"]["qualified_graph_reports"] == 8
    assert (
        j89.replay(
            tmp_path / "real-source-metadata",
            output["snapshot_sha256"],
            repo=REPO,
            **external,
        )
        == output["result"]
    )


@pytest.mark.parametrize(
    "missing", [n + suffix for n in j89.INPUTS for suffix in ("", "_sha256")]
)
def test_missing_required_input_blocks_before_any_external_open(
    tmp_path, monkeypatch, missing
):
    external = {
        k: (tmp_path / n if k == n else "0" * 64)
        for n in j89.INPUTS
        for k in (n, n + "_sha256")
    }
    del external[missing]

    def forbidden(*_a, **_kw):
        pytest.fail("Incomplete input set must not open supplied files")

    monkeypatch.setattr(j89, "exact_eight_request", forbidden)
    monkeypatch.setattr(j89, "audit_identity", forbidden)
    _, reports = j89.reconstruct(REPO, **external)
    assert reports["result.json"]["decision"] == "source_bundle_required"
    assert (
        reports["result.json"]["source_files"]
        == reports["result.json"]["source_records_inspected"]
        == 0
    )
    assert not reports["input-gate.json"]["external_source_inspection_performed"]


@pytest.mark.parametrize("missing", j89.INPUTS)
def test_missing_json_file_is_a_source_required_receipt(tmp_path, missing):
    external = software_inputs(tmp_path)
    external[missing].unlink()
    out = j89.emit(REPO, tmp_path / "blocked", **external)
    assert out["result"]["decision"] == "source_bundle_required"
    assert (
        "required_external_file_unavailable:" + missing
        in out["result"]["blocking_reasons"]
    )
    assert (tmp_path / "blocked/qualification-receipt.json").is_file()


def test_unauthorized_source_stays_source_required(tmp_path):
    external = software_inputs(tmp_path, authorized=False)
    _, reports = j89.reconstruct(REPO, **external)
    assert reports["result.json"]["decision"] == "source_bundle_required"
    assert reports["result.json"]["qualified_graph_reports"] == 0
    assert reports["audit-references.json"]["status"] == "source_authorization_blocked"


@pytest.mark.parametrize("field", ["reviewer", "authority_basis", "policy"])
def test_missing_trusted_audit_identity_or_provenance_does_not_open_parameters(
    tmp_path, monkeypatch, field
):
    external = software_inputs(tmp_path)
    audit = j89.read(external["audit"])
    audit[field] = "unknown"
    external["audit_sha256"] = put(external["audit"], audit)
    original = j89.inherited.reconstruct

    def empty_only(repo, **kwargs):
        assert not kwargs, "Untrusted audit must not trigger external source capture"
        return original(repo)

    monkeypatch.setattr(j89.inherited, "reconstruct", empty_only)
    _, out = j89.reconstruct(REPO, **external)
    assert out["result.json"]["decision"] == "source_bundle_required"
    assert not out["result.json"]["source_bytes_rehashed"]


@pytest.mark.parametrize(
    "which", ["request", "manifest", "audit", "record_source", "primary_artifact"]
)
def test_all_pinned_bytes_are_reopened_on_replay(tmp_path, which):
    external = software_inputs(tmp_path)
    root = tmp_path / "output"
    emitted = j89.emit(REPO, root, **external)
    if which in external:
        path = external[which]
    else:
        manifest = j89.read(external["manifest"])
        path = external["manifest"].parent / next(
            f["path"] for f in manifest["files"] if f["role"] == which
        )
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA256_mutation"):
        j89.replay(root, emitted["snapshot_sha256"], repo=REPO, **external)


@pytest.mark.parametrize(
    "mutation", ["fixture_bytes", "graph_id", "subset", "fixture_hash"]
)
def test_exact_all_eight_j87_scope_required(tmp_path, mutation):
    external = software_inputs(tmp_path)
    req = j89.read(external["request"])
    if mutation == "fixture_bytes":
        path = external["request"].parent / req["fixtures"][0]["path"]
        path.write_bytes(path.read_bytes() + b"\n")
    elif mutation == "subset":
        removed = req["fixtures"].pop()
        req["graph_identities"].remove(removed["graph_identity"])
    elif mutation == "graph_id":
        req["fixtures"][0]["graph_identity"] = "0" * 64
    else:
        req["fixtures"][0]["sha256"] = "0" * 64
    external["request_sha256"] = put(external["request"], req)
    with pytest.raises(ValueError, match="J87|SHA256_mutation"):
        j89.reconstruct(REPO, **external)


@pytest.mark.parametrize(
    "mutation",
    [
        "version",
        "energy_unit",
        "length_unit",
        "torsion_phase",
        "torsion_normalization",
        "scaling_14_lj",
        "scaling_14_coulomb",
        "mixing_rule",
        "long_range",
    ],
)
def test_exact_version_units_and_conventions_are_not_reinterpreted(tmp_path, mutation):
    external = software_inputs(tmp_path)
    request = j89.read(external["request"])
    if mutation == "version":
        request["expected_version"] = "different-release"
    else:
        request["expected_conventions"][mutation] = "incompatible-policy"
    external["request_sha256"] = put(external["request"], request)
    _, reports = j89.reconstruct(REPO, **external)
    assert reports["result.json"]["decision"] == "parameter_bundle_incomplete"
    assert reports["result.json"]["qualified_graph_reports"] == 0
    assert reports["convention-comparison.json"]["blocking_reasons"]


@pytest.mark.parametrize(
    "family", ["bond", "angle", "torsion", "improper", "charge", "nonbonded"]
)
def test_each_missing_class_i_term_family_remains_explicit(tmp_path, family):
    def change(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == family)
        ]

    external = software_inputs(tmp_path, change=change)
    _, reports = j89.reconstruct(REPO, **external)
    assert reports["result.json"]["decision"] == "parameter_bundle_incomplete"
    graphs = reports["graph-coverage-report.json"]["reports"]
    assert all(g["coverage"][family]["covered"] == 0 for g in graphs)
    assert any(
        g["coverage"][family]["required"] > 0 and g["missing_terms"] for g in graphs
    )
    assert reports["result.json"]["accepted_coefficients"] == 0


@pytest.mark.parametrize(
    "forms", [None, {}, "bad", [None], ["bond"], [{}], [{"family": "bond"}]]
)
def test_malformed_forms_produce_structured_blockers_before_source_inspection(
    tmp_path, monkeypatch, capsys, forms
):
    external = software_inputs(tmp_path)
    request = j89.read(external["request"])
    request["expected_forms"] = forms
    external["request_sha256"] = put(external["request"], request)
    flags = [
        flag
        for n, v in external.items()
        for flag in ("--" + n.replace("_", "-"), str(v))
    ]
    monkeypatch.setattr(
        sys, "argv", ["j89", "ingest", "--directory", str(tmp_path / "out"), *flags]
    )
    assert j89.main() == 2
    result = json.loads(capsys.readouterr().out)
    assert (
        result["status"] == "blocked"
        and "expected_form" in result["blocking_reasons"][0]
    )
    assert not (tmp_path / "out").exists()


def test_relocated_external_replay_requires_originals_not_metadata(tmp_path):
    source = tmp_path / "original"
    external = software_inputs(source)
    output = j89.emit(REPO, source / "evidence", **external)
    moved = tmp_path / "moved"
    shutil.copytree(source, moved)
    for name in j89.INPUTS:
        external[name] = moved / external[name].relative_to(source)
    shutil.rmtree(source)
    assert (
        j89.replay(moved / "evidence", output["snapshot_sha256"], repo=REPO, **external)
        == output["result"]
    )
    with pytest.raises(
        ValueError, match="explicit_external_replay_paths_and_pins_required"
    ):
        j89.replay(moved / "evidence", output["snapshot_sha256"], repo=REPO)
    external["manifest"].unlink()
    with pytest.raises(ValueError, match="external_originals_unavailable"):
        j89.replay(moved / "evidence", output["snapshot_sha256"], repo=REPO, **external)


def test_no_source_is_deterministic_and_creates_no_scientific_coverage(tmp_path):
    one = j89.emit(REPO, tmp_path / "one")
    two = j89.emit(REPO, tmp_path / "two")
    assert one == two
    assert j89.replay(tmp_path / "one", repo=REPO) == one["result"]
    r = one["result"]
    assert r["decision"] == "source_bundle_required"
    assert r["required_graph_count"] == r["exact_graph_fixtures"] == 8
    assert (
        r["source_records_inspected"]
        == r["parameter_rows_inspected"]
        == r["source_files"]
        == r["qualified_graph_reports"]
        == 0
    )
    assert r["audit_sha256"] is r["manifest_sha256"] is r["request_sha256"] is None


def test_persistence_contains_no_raw_tables_or_coefficient_arrays(tmp_path):
    external = software_inputs(tmp_path)
    root = tmp_path / "metadata"
    j89.emit(REPO, root, **external)

    def walk(v):
        if isinstance(v, dict):
            assert (
                not {
                    "raw_text",
                    "values",
                    "record",
                    "request_raw",
                    "manifest_raw",
                    "audit_raw",
                    "evidence_graph_json",
                    "graph_reports_json",
                }
                & v.keys()
            )
            for value in v.values():
                walk(value)
        elif isinstance(v, list):
            for value in v:
                walk(value)

    for path in root.iterdir():
        walk(j89.read(path))


@pytest.mark.parametrize(
    "mutation", ["decision", "scope", "coverage", "audit", "hashes", "inputs"]
)
def test_strict_receipt_rejects_resigned_claims(tmp_path, monkeypatch, mutation):
    from scripts import validate_phase_4j89_receipt as validator

    root = tmp_path / "sealed"
    j89.emit(REPO, root)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, root)
    receipt = j89.read(root / "receipt.json")
    if mutation == "decision":
        receipt["decision"] = j89.READY
    elif mutation == "scope":
        receipt["scope"]["automatic_assignment_enabled"] = True
    else:
        name = {
            "coverage": "graph-coverage-report.json",
            "audit": "audit-references.json",
            "hashes": "source-locations.json",
            "inputs": "input-gate.json",
        }[mutation]
        obj = j89.read(root / name)
        if mutation == "coverage":
            obj["reports"][0]["missing_terms"] = []
        elif mutation == "audit":
            obj["records"] = [{"audit_sha256": "0" * 64}]
        elif mutation == "hashes":
            obj["source_hashes"] = ["0" * 64]
        else:
            obj["blocking_reasons"] = []
        j89.write(
            root / name, j89.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        manifest = j89.read(root / "manifest.json")
        manifest["artifacts"][name] = j89.pin(root / name)
        j89.write(
            root / "manifest.json",
            j89.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j89.pin(root / "manifest.json")
    j89.write(
        root / "receipt.json",
        j89.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, root)


@pytest.mark.parametrize(
    "name", ["sqm", "lammps", "msi2lmp", "fit_parameters", "evaluate_forces"]
)
def test_prohibited_processes_never_execute(tmp_path, name):
    marker = tmp_path / "ran"
    fake = tmp_path / name
    fake.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    fake.chmod(0o700)
    run = subprocess.run(
        [
            sys.executable,
            "-c",
            "from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["
            + repr(str(fake))
            + "])",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode and not marker.exists()


def test_allowlist_has_exactly_the_historical_eight_graphs():
    assert set(ALLOWED_FIXTURES) == {
        "hydrocarbon",
        "aromatic",
        "ether",
        "ester",
        "carbonyl",
        "amide",
        "nitrile",
        "halogenated",
    }


def test_invalid_pin_does_not_open_any_external_input(tmp_path, monkeypatch):
    external = {
        k: (tmp_path / n if k == n else "0" * 64)
        for n in j89.INPUTS
        for k in (n, n + "_sha256")
    }
    external["manifest_sha256"] = "not-a-checksum"
    monkeypatch.setattr(
        j89, "exact_eight_request", lambda *_: pytest.fail("Must not open request")
    )
    _, reports = j89.reconstruct(REPO, **external)
    assert reports["result.json"]["decision"] == "source_bundle_required"
    assert (
        "required_external_SHA256_pin_missing_or_invalid:manifest"
        in reports["input-gate.json"]["blocking_reasons"]
    )


def test_missing_opaque_attachment_returns_source_required(
    tmp_path, monkeypatch, capsys
):
    external = software_inputs(tmp_path)
    manifest = j89.read(external["manifest"])
    file = next(f for f in manifest["files"] if f["role"] == "primary_artifact")
    (external["manifest"].parent / file["path"]).unlink()
    flags = [
        flag
        for n, v in external.items()
        for flag in ("--" + n.replace("_", "-"), str(v))
    ]
    monkeypatch.setattr(
        sys, "argv", ["j89", "ingest", "--directory", str(tmp_path / "out"), *flags]
    )
    assert j89.main() == 2
    report = json.loads(capsys.readouterr().out)
    assert report["decision"] == "source_bundle_required"
    assert "source_file_unavailable" in report["blocking_reasons"][0]
    assert not (tmp_path / "out").exists()


def test_external_strict_receipt_requires_sources_and_pins(tmp_path, monkeypatch):
    from scripts import validate_phase_4j89_receipt as validator

    external = software_inputs(tmp_path)
    root = tmp_path / "sealed"
    emitted = j89.emit(REPO, root, **external)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, root, snapshot_sha=emitted["snapshot_sha256"], **external)
    assert (
        validator.validate(
            REPO, root, snapshot_sha=emitted["snapshot_sha256"], **external
        )["status"]
        == "passed"
    )
    with pytest.raises(ValueError, match="external_snapshot_pin_required"):
        validator.validate(REPO, root)
