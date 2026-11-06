"""Manufactured software controls only; no scientific source supplied or qualified."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.forcefields.bundle_adapters.opls_authorized import (
    ALLOWED_FIXTURES,
    capture_request,
    request_fields,
)
from scripts import phase_4j88_opls as j88
from tests.test_phase_4j84 import put
from tests.test_phase_4j87 import SPECS, request_for, source_set

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j88.scope()


def control(tmp_path, *, change=None, authorized=True, specs=None):
    specs = specs or SPECS[:1]
    _, cap, args = source_set(tmp_path / "source", specs, change, authorized)
    path, _ = request_for(tmp_path / "request", cap, specs)
    request = j88.read(path)
    for f in request["fixtures"]:
        original = REPO / "tests/fixtures/phase_4j87/graphs" / (f["id"] + ".json")
        shutil.copyfile(original, path.parent / f["path"])
        f.update(ALLOWED_FIXTURES[f["id"]])
    digest = put(path, request)
    return {
        "request": path,
        "request_sha256": digest,
        "manifest": args[0],
        "manifest_sha256": args[1],
        "audit": args[2],
        "audit_sha256": args[3],
    }


def result(tmp_path, **kwargs):
    external = control(tmp_path, **kwargs)
    _, outputs = j88.reconstruct(REPO, **external)
    return outputs, external


def test_valid_authorized_opls_control_all_eight_exact_graphs(tmp_path):
    outputs, _ = result(tmp_path, specs=SPECS)
    r = outputs["result.json"]
    assert r["decision"] == j88.READY
    assert r["qualified_graph_reports"] == r["exact_graph_fixtures"] == 8
    assert r["scope"]["accepted_coefficients"] == 0
    assert r["chemical_family_transfer"] == "not_established"
    assert not r["scope"]["runtime_registration"]
    for graph in outputs["graph-coverage-report.json"]["reports"]:
        assert graph["decision"] == j88.READY
        assert not graph["missing_terms"] and not graph["blocking_reasons"]
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
        assert all(c["required"] == c["covered"] for c in graph["coverage"].values())
        assert all(
            row["location"]
            and row["pointer"]
            and row["claim_sha256"]
            and row["source_sha256"] in r["source_hashes"]
            for row in graph["matched_source_records"]
        )


@pytest.mark.parametrize("missing", ["authorization", "audit", "bundle"])
def test_missing_source_authorization_or_independent_audit(tmp_path, missing):
    external = control(tmp_path, authorized=missing != "authorization")
    if missing == "audit":
        external.update(audit=None, audit_sha256=None)
    if missing == "bundle":
        external = {}
    _, out = j88.reconstruct(REPO, **external)
    assert out["result.json"]["decision"] == "source_bundle_required"
    assert not out["result.json"]["source_authorized"]
    assert not out["result.json"]["qualified_graph_reports"]


@pytest.mark.parametrize(
    "mutation", ["version", "unit", "torsion_phase", "forms", "policy14"]
)
def test_request_source_convention_and_version_mismatch(tmp_path, mutation):
    external = control(tmp_path)
    req = j88.read(external["request"])
    if mutation == "version":
        req["expected_version"] = "other-release"
    elif mutation == "forms":
        req["expected_forms"][0]["equation"] += " altered"
    else:
        key = {"unit": "energy_unit", "policy14": "scaling_14_lj"}.get(
            mutation, mutation
        )
        req["expected_conventions"][key] = "incompatible"
    external["request_sha256"] = put(external["request"], req)
    _, outputs = j88.reconstruct(REPO, **external)
    assert outputs["result.json"]["decision"] == "parameter_bundle_incomplete"
    assert outputs["convention-comparison.json"]["blocking_reasons"]


@pytest.mark.parametrize(
    "family", ["bond", "angle", "torsion", "improper", "charge", "nonbonded"]
)
def test_missing_term_coverage_is_explicit_and_never_zero_filled(tmp_path, family):
    def change(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == family)
        ]

    outputs, _ = result(tmp_path, change=change)
    graph = outputs["graph-coverage-report.json"]["reports"][0]
    assert outputs["result.json"]["decision"] == "parameter_bundle_incomplete"
    assert graph["coverage"][family]["covered"] == 0
    assert graph["coverage"][family]["required"] > 0
    assert graph["missing_terms"]
    assert graph["coverage"]["class_ii_cross_terms"]["required"] == 0


@pytest.mark.parametrize(
    "mutation",
    ["graph_identity", "fixture_sha", "fixture_size", "fixture_id", "fixture_bytes"],
)
def test_immutable_j87_fixture_scope_rejects_mutations(tmp_path, mutation):
    external = control(tmp_path)
    req = j88.read(external["request"])
    f = req["fixtures"][0]
    if mutation == "fixture_bytes":
        path = external["request"].parent / f["path"]
        path.write_bytes(path.read_bytes() + b"\n")
    else:
        key, value = {
            "graph_identity": ("graph_identity", "0" * 64),
            "fixture_sha": ("sha256", "0" * 64),
            "fixture_size": ("size", 1),
            "fixture_id": ("id", "invented"),
        }[mutation]
        f[key] = value
        external["request_sha256"] = put(external["request"], req)
    with pytest.raises(ValueError, match="J87|immutable|SHA256_mutation"):
        j88.reconstruct(REPO, **external)


@pytest.mark.parametrize(
    "field", ["request", "manifest", "audit", "record_source", "primary_artifact"]
)
def test_source_bytes_must_be_rehashed_during_replay(tmp_path, field):
    external = control(tmp_path)
    target = tmp_path / "output"
    emitted = j88.emit(REPO, target, **external)
    if field in external:
        path = external[field]
    else:
        m = j88.read(external["manifest"])
        role = "record_source" if field == "record_source" else "primary_artifact"
        path = external["manifest"].parent / next(
            f["path"] for f in m["files"] if f["role"] == role
        )
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA256_mutation"):
        j88.replay(target, emitted["snapshot_sha256"], repo=REPO, **external)


@pytest.mark.parametrize(
    "forms", [None, {}, "bad", [None], [1], ["bad"], [[]], [{}], [{"family": "bond"}]]
)
def test_malformed_expected_forms_have_structured_cli_blockers(
    tmp_path, capsys, monkeypatch, forms
):
    # A pinned request with malformed forms must stop before opening any source.
    req = j88.read(REPO / j88.REQUEST_PATH)
    req["expected_forms"] = forms
    path = tmp_path / "request.json"
    digest = put(path, req)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "j88",
            "ingest",
            "--directory",
            str(tmp_path / "out"),
            "--request",
            str(path),
            "--request-sha256",
            digest,
        ],
    )
    assert j88.main() == 2
    blocked = json.loads(capsys.readouterr().out)
    assert blocked["status"] == "blocked" and blocked["blocking_reasons"]
    assert "expected_form" in blocked["blocking_reasons"][0]
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "field,value",
    [
        ("family", []),
        ("family", "cross_term"),
        ("equation", None),
        ("equation", ""),
        ("ordering", {}),
        ("fields", "k"),
        ("fields", [[]]),
        ("fields", []),
        ("units", None),
        ("units", []),
        ("units", [False]),
        ("equilibrium_fields", "r0"),
        ("equilibrium_fields", [{}]),
        ("equilibrium_fields", ["unknown"]),
    ],
)
def test_complete_form_field_schema_precedes_any_access(tmp_path, field, value):
    external = control(tmp_path)
    req = j88.read(external["request"])
    req["expected_forms"][0][field] = value
    with pytest.raises(ValueError):
        request_fields(req)


def test_duplicate_and_incomplete_form_collections_are_blocked(tmp_path):
    external = control(tmp_path)
    req = j88.read(external["request"])
    req["expected_forms"] *= 2
    with pytest.raises(ValueError, match="duplicate_expected_form"):
        request_fields(req)
    req["expected_forms"] = []
    external["request_sha256"] = put(external["request"], req)
    _, out = j88.reconstruct(REPO, **external)
    assert out["result.json"]["decision"] == "parameter_bundle_incomplete"
    assert (
        "complete_expected_class_i_forms_required"
        in out["result.json"]["blocking_reasons"]
    )


def test_persisted_reports_never_embed_raw_source_or_parameter_values(tmp_path):
    external = control(tmp_path)
    root = tmp_path / "result"
    j88.emit(REPO, root, **external)

    def visit(value):
        if isinstance(value, dict):
            assert (
                not {
                    "raw_text",
                    "manifest_raw",
                    "audit_raw",
                    "request_raw",
                    "values",
                    "evidence_graph_json",
                    "graph_reports_json",
                }
                & value.keys()
            )
            for v in value.values():
                visit(v)
        elif isinstance(value, list):
            for v in value:
                visit(v)

    for path in root.iterdir():
        visit(j88.read(path))
    assert len(list(root.iterdir())) == len(j88.REPORTS) + 2


def test_external_replay_requires_explicit_paths_pins_and_snapshot_pin(tmp_path):
    external = control(tmp_path)
    root = tmp_path / "output"
    emitted = j88.emit(REPO, root, **external)
    with pytest.raises(ValueError, match="external_snapshot_pin_required"):
        j88.replay(root, repo=REPO, **external)
    with pytest.raises(ValueError, match="external_replay_pin_mutation"):
        j88.replay(root, emitted["snapshot_sha256"], repo=REPO)
    assert (
        j88.replay(root, emitted["snapshot_sha256"], repo=REPO, **external)
        == emitted["result"]
    )


def test_relocated_source_request_and_report_replay(tmp_path):
    original = tmp_path / "original"
    external = control(original)
    emitted = j88.emit(REPO, original / "output", **external)
    moved = tmp_path / "moved"
    shutil.copytree(original, moved)
    for name in ("request", "manifest", "audit"):
        external[name] = moved / external[name].relative_to(original)
    shutil.rmtree(original)
    assert (
        j88.replay(moved / "output", emitted["snapshot_sha256"], repo=REPO, **external)
        == emitted["result"]
    )


def test_deterministic_no_source_replay_is_pinned_to_original_j87_graphs(tmp_path):
    one = j88.emit(REPO, tmp_path / "one")
    two = j88.emit(REPO, tmp_path / "two")
    assert one == two
    r = j88.replay(tmp_path / "one", repo=REPO)
    assert r == one["result"]
    assert r["decision"] == "source_bundle_required"
    assert r["exact_graph_fixtures"] == 8
    assert (
        r["source_records_inspected"]
        == r["source_files"]
        == r["qualified_graph_reports"]
        == 0
    )
    assert r["missing_term_instances"] > 0


@pytest.mark.parametrize(
    "mutation", ["decision", "scope", "coverage", "location", "forms", "snapshot"]
)
def test_strict_receipt_rejects_resigned_semantic_mutations(
    tmp_path, monkeypatch, mutation
):
    from scripts import validate_phase_4j88_receipt as validator

    root = tmp_path / "evidence"
    j88.emit(REPO, root)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, root)
    receipt = j88.read(root / "receipt.json")
    if mutation in ("decision", "scope"):
        if mutation == "decision":
            receipt["decision"] = j88.READY
        else:
            receipt["scope"]["runtime_registration"] = True
    else:
        name = {
            "coverage": "graph-coverage-report.json",
            "location": "source-locations.json",
            "forms": "convention-comparison.json",
            "snapshot": "captured-input.json",
        }[mutation]
        obj = j88.read(root / name)
        if mutation == "coverage":
            obj["reports"][0]["missing_terms"] = []
        elif mutation == "location":
            obj["records"] = [{"location": "fabricated"}]
        elif mutation == "forms":
            obj["performed"] = True
        else:
            obj["manifest_sha256"] = "0" * 64
        j88.write(
            root / name,
            j88.signed({k: v for k, v in obj.items() if k != "identity"})
            if mutation != "snapshot"
            else obj,
        )
        m = j88.read(root / "manifest.json")
        m["artifacts"][name] = j88.pin(root / name)
        j88.write(
            root / "manifest.json",
            j88.signed({k: v for k, v in m.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j88.pin(root / "manifest.json")
    j88.write(
        root / "receipt.json",
        j88.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation|source_claims"):
        validator.validate(REPO, root)


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "lammps",
        "msi2lmp",
        "fit_parameters",
        "generate_parameters",
        "reserved_validation",
        "evaluate_forces",
    ],
)
def test_prohibited_child_process_cannot_run(tmp_path, name):
    marker = tmp_path / "ran"
    fake = tmp_path / name
    fake.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    fake.chmod(0o700)
    child = subprocess.run(
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
    assert child.returncode and not marker.exists()


def test_no_implicit_family_switch_or_graph_scope_transfer(tmp_path):
    external = control(tmp_path)
    req = j88.read(external["request"])
    req["family"] = "GAFF2"
    with pytest.raises(ValueError, match="OPLS_AA_only"):
        request_fields(req)
    req["family"] = "OPLS-AA"
    req["fixtures"][0]["graph_identity"] = "0" * 64
    req["graph_identities"] = ["0" * 64]
    with pytest.raises(ValueError, match="J87_fixture"):
        request_fields(req)


def test_symlink_request_ancestor_and_missing_pin_rejected(tmp_path):
    external = control(tmp_path / "original")
    with pytest.raises(ValueError, match="external_source_pin_required"):
        capture_request(external["request"], None)
    link = tmp_path / "link"
    link.symlink_to(external["request"].parent, target_is_directory=True)
    external["request"] = link / external["request"].name
    with pytest.raises(ValueError, match="external_path_symlink_forbidden"):
        j88.reconstruct(REPO, **external)


def test_claim_or_source_location_mutation_rejected_by_independent_audit(tmp_path):
    external = control(tmp_path)
    m = j88.read(external["manifest"])
    f = next(f for f in m["files"] if f["role"] == "record_source")
    path = external["manifest"].parent / f["path"]
    src = j88.read(path)
    src["records"][0]["location"] = "unaudited location"
    put(path, src)
    with pytest.raises(ValueError, match="SHA256_mutation"):
        j88.reconstruct(REPO, **external)


def test_request_mapping_mutation_is_structured(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "j88",
            "ingest",
            "--request",
            str(tmp_path / "bad"),
            "--request-sha256",
            "0" * 64,
        ],
    )
    assert j88.main() == 2
    obj = json.loads(capsys.readouterr().out)
    assert (
        obj["status"] == "blocked" and not obj["scope"]["automatic_assignment_enabled"]
    )


def test_ambiguous_typing_is_explicit_in_term_and_blocker_reports(tmp_path):
    def change(_, rows):
        old = next(
            r
            for r in rows
            if r["kind"] == "atom_type" and r["value"]["atom_type"] == "control_C"
        )
        new = json.loads(json.dumps(old))
        new["record_id"] = "atom_type:competing"
        new["value"].update(id="competing", atom_type="competing_C")
        rows.append(new)

    outputs, _ = result(tmp_path, change=change)
    assert outputs["result.json"]["decision"] == "parameter_bundle_incomplete"
    diagnostic = outputs["missing-terms-and-blockers.json"]["graphs"][0]
    assert diagnostic["missing_terms"] and diagnostic["blocker_categories"]["ambiguous"]


def test_resigned_independent_audit_cannot_fabricate_claim_hash(tmp_path):
    external = control(tmp_path)
    audit = j88.read(external["audit"])
    review = next(r for r in audit["source_reviews"] if r["claims"])
    review["claims"][0]["claim_sha256"] = "0" * 64
    external["audit_sha256"] = put(external["audit"], audit)
    with pytest.raises(ValueError, match="audit|claim"):
        j88.reconstruct(REPO, **external)


def test_source_location_blockers_are_collected_without_raw_rows(tmp_path):
    from island.forcefields.bundle_adapters.opls_authorized import (
        InputBlockers,
        capture_opls,
    )
    from island.forcefields.bundle_adapters.schema import (
        manifest_identity,
        source_digest,
    )

    external = control(tmp_path)
    manifest = j88.read(external["manifest"])
    f = next(f for f in manifest["files"] if f["role"] == "record_source")
    path = external["manifest"].parent / f["path"]
    source = j88.read(path)
    for row in source["records"][:2]:
        row["location"] = ""
    f["sha256"] = put(path, source)
    f["size"] = path.stat().st_size
    manifest["source_sha256"] = source_digest(manifest["files"])
    manifest["identity"] = manifest_identity(manifest)
    external["manifest_sha256"] = put(external["manifest"], manifest)
    with pytest.raises(InputBlockers) as error:
        capture_opls(
            external["manifest"],
            external["manifest_sha256"],
            external["audit"],
            external["audit_sha256"],
            capture_timestamp=j88.WHEN,
        )
    assert len(error.value.blockers) == 2
    assert all("source_location_blocked:" in b for b in error.value.blockers)


def test_missing_graph_audit_cannot_be_replaced_by_matching_chemistry(tmp_path):
    def change(_, rows):
        rows[:] = [r for r in rows if r["kind"] != "graph_scope"]

    outputs, _ = result(tmp_path, change=change)
    assert outputs["result.json"]["decision"] == "parameter_bundle_incomplete"
    assert (
        "exact_graph_scope_not_evidenced" in outputs["result.json"]["blocking_reasons"]
    )


def test_no_source_ingest_timestamp_stays_replayable(tmp_path):
    emitted = j88.emit(REPO, tmp_path / "out", capture_timestamp="2026-11-07T03:04:05Z")
    assert j88.replay(tmp_path / "out", repo=REPO) == emitted["result"]


def test_cli_external_receipt_and_replay_reopen_explicit_inputs(tmp_path):
    external = control(tmp_path)
    root = tmp_path / "output"
    flags = []
    for name, value in external.items():
        flags.extend(["--" + name.replace("_", "-"), str(value)])
    command = [sys.executable, "-m", "scripts.phase_4j88_opls"]
    run = subprocess.run(
        [*command, "ingest", "--directory", str(root), *flags],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr + run.stdout
    emitted = json.loads(run.stdout)
    replay = subprocess.run(
        [
            *command,
            "replay",
            "--directory",
            str(root),
            "--snapshot-sha256",
            emitted["snapshot_sha256"],
            *flags,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert replay.returncode == 0, replay.stderr + replay.stdout
    assert json.loads(replay.stdout) == emitted["result"]
    assert emitted["result"]["decision"] == j88.READY
