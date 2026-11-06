"""Manufactured source/audit controls only; zero scientific OPLS qualification."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.schema import data
from island.forcefields.bundle_adapters.opls_qualification import (
    capture_opls,
    capture_request,
    qualify_opls,
)
from scripts import phase_4j87_opls as j87
from tests.test_phase_4j82 import control
from tests.test_phase_4j84 import WHEN, put
from tests.test_phase_4j85 import supplied

REPO = Path(__file__).resolve().parents[1]
SPECS = j87.read(REPO / "tests/fixtures/phase_4j85/organic-motifs.json")["motifs"]
READY = "family_scoped_qualification_ready_runtime_blocked"


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j87.scope()


def source_set(tmp_path, specs=None, change=None, authorized=True):
    specs = specs or SPECS[:1]
    graphs = [normalize(j87.spec_from_data(spec)) for spec in specs]

    def additional(m, rows):
        for index, g in enumerate(graphs[1:], 1):
            _, p = control(g)
            for r in p.rules:
                rid = "atom_type:" + r.id
                if not any(row["record_id"] == rid for row in rows):
                    rows.append(
                        {
                            "record_id": rid,
                            "kind": "atom_type",
                            "location": "Manufactured software control; typing row "
                            + r.id,
                            "origin": "directly_reported",
                            "value": {
                                k: v for k, v in data(r).items() if k != "evidence_id"
                            },
                        }
                    )
            for r in p.parameters:
                rid = "parameter:" + r.id
                if not any(row["record_id"] == rid for row in rows):
                    rows.append(
                        {
                            "record_id": rid,
                            "kind": "parameter",
                            "location": "Manufactured software control; parameter "
                            + r.id,
                            "origin": "directly_reported",
                            "value": {
                                k: v
                                for k, v in data(r).items()
                                if k
                                not in (
                                    "profile_family",
                                    "evidence_id",
                                    "convention_identity",
                                    "source_location",
                                )
                            },
                        }
                    )
            rows.append(
                {
                    "record_id": "graph_scope:" + str(index),
                    "kind": "graph_scope",
                    "location": "Manufactured software control; graph " + str(index),
                    "origin": "directly_reported",
                    "value": {"graph_identity": g.identity},
                }
            )
        if change:
            change(m, rows)

    _, cap, args = supplied(
        tmp_path, g=graphs[0], authorized=authorized, change=additional
    )
    return graphs, cap, args


def request_for(tmp_path, cap, specs=None):
    specs = specs or SPECS[:1]
    files = []
    for spec in specs:
        path = tmp_path / (spec["id"] + ".json")
        digest = put(
            path,
            {
                "schema": "island.j87.exact-graph-fixture.v1",
                "id": spec["id"],
                "spec": spec,
            },
        )
        files.append(
            {
                "id": spec["id"],
                "path": path.name,
                "sha256": digest,
                "size": path.stat().st_size,
                "graph_identity": normalize(j87.spec_from_data(spec)).identity,
            }
        )
    conventions = {}
    for f in cap["files"]:
        if f["raw_text"]:
            for row in json.loads(f["raw_text"])["records"]:
                if row["kind"] == "convention":
                    conventions[row["value"]["name"]] = row["value"]["value"]
    # Fake expected declarations belong to this software control, not physical OPLS.
    m = cap["manifest"]
    path = tmp_path / "request.json"
    digest = put(
        path,
        {
            "schema": "island.j87.opls-qualification-request.v1",
            "family": "OPLS-AA",
            "profile_id": m["profile_id"],
            "graph_identities": sorted(f["graph_identity"] for f in files),
            "fixtures": sorted(files, key=lambda f: f["id"]),
            "expected_version": m["version"],
            "expected_conventions": conventions,
            "expected_forms": m["functional_form_convention"],
        },
    )
    return path, digest


def run_request(tmp_path, *, specs=None, change=None, authorized=True):
    _graphs, cap, args = source_set(tmp_path / "source", specs, change, authorized)
    path, digest = request_for(tmp_path / "graphs", cap, specs)
    rc = capture_request(path, digest)
    return qualify_opls(rc, cap, j87.spec_from_data), rc, cap, args


def test_valid_audited_input_branch_all_eight_exact_graphs_software_only(tmp_path):
    result, _, cap, _ = run_request(tmp_path, specs=SPECS)
    assert result["decision"] == READY
    assert len(result["graph_reports"]) == 8
    assert result["chemical_family_transfer"] == "not_established"
    assert result["scope_of_qualification"] == "declared_exact_graphs_only"
    assert result["source_hashes"] and result["source_locations"]
    assert set(result["source_hashes"]) == {
        f["sha256"] for f in cap["manifest"]["files"]
    }
    for report in result["graph_reports"]:
        assert not report["missing_terms"] and not report["blocking_reasons"]
        assert report["decision"] == READY
        assert report["matched_source_records"]
        assert all(c["required"] == c["covered"] for c in report["coverage"].values())
        assert all(
            r["source_sha256"] in result["source_hashes"]
            and r["location"]
            and r["pointer"]
            for r in report["matched_source_records"]
        )
    assert (
        not result["automatic_assignment_enabled"]
        and not result["runtime_registration"]
    )
    assert result["accepted_coefficients"] == 0 and not result["production_validated"]
    assert result["simulation_readiness"] == "not_established"


def test_missing_bundle_has_no_scientific_coverage():
    rc = capture_request(REPO / j87.REQUEST_PATH, j87.REQUEST_SHA)
    result = qualify_opls(rc, capture_opls(capture_timestamp=WHEN), j87.spec_from_data)
    assert result["decision"] == "source_bundle_required"
    assert not result["convention_comparison"]["performed"]
    assert result["source_hashes"] == () and result["source_locations"] == []
    assert all(
        r["missing_terms"] and not r["matched_source_records"]
        for r in result["graph_reports"]
    )
    assert all(
        c["covered"] == 0
        for r in result["graph_reports"]
        for c in r["coverage"].values()
        if c["required"]
    )


@pytest.mark.parametrize(
    "why", ("authorization", "no_audit", "synthetic", "polyjarvis")
)
def test_unauthorized_sources_never_qualify(tmp_path, why):
    def change(m, _):
        if why == "authorization":
            m["provenance"]["authorization"] = "not_authorized"
        elif why == "polyjarvis":
            m["provenance"]["upstream_identity"] = "PolyJarvis numeric patches"

    _, cap, args = source_set(tmp_path, change=change, authorized=why != "synthetic")
    if why == "no_audit":
        cap = capture_opls(args[0], args[1], capture_timestamp=WHEN)
    path, digest = request_for(tmp_path / "graphs", cap)
    result = qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)
    assert result["decision"] == "source_bundle_required"
    assert result["accepted_coefficients"] == 0


@pytest.mark.parametrize("family", ("GAFF", "GAFF2", "PCFF", "COMPASS"))
def test_other_families_rejected_before_source_files_opened(
    tmp_path, monkeypatch, family
):
    import island.forcefields.bundle_adapters.opls_qualification as module

    _, _, args = supplied(tmp_path, family, authorized=True)
    monkeypatch.setattr(
        module, "capture", lambda *a, **k: pytest.fail("other-family source opened")
    )
    with pytest.raises(ValueError, match="J87_OPLS_AA_only"):
        capture_opls(*args, capture_timestamp=WHEN)


@pytest.mark.parametrize("target", ("source", "manifest", "audit"))
def test_source_hash_mutation_rejected(tmp_path, target):
    _, _, args = source_set(tmp_path)
    path = (
        args[0]
        if target == "manifest"
        else args[2]
        if target == "audit"
        else args[0].parent / "sources/records_0.json"
    )
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA256_mutation"):
        capture_opls(*args, capture_timestamp=WHEN)


@pytest.mark.parametrize(
    "field",
    (
        "energy_unit",
        "length_unit",
        "torsion_phase",
        "torsion_normalization",
        "scaling_14_lj",
        "scaling_14_coulomb",
        "exclusions_12",
        "exclusions_13",
        "charge_origin",
        "mixing_rule",
    ),
)
def test_pinned_expected_convention_mutation_cannot_qualify(tmp_path, field):
    _, cap, _ = source_set(tmp_path)
    path, _ = request_for(tmp_path / "graphs", cap)
    request = j87.read(path)
    request["expected_conventions"][field] += "_different"
    digest = put(path, request)
    result = qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert "request_source_convention_mismatch:" + field in result["blocking_reasons"]
    assert all(r["decision"] != READY for r in result["graph_reports"])


@pytest.mark.parametrize(
    "family", ("bond", "angle", "torsion", "improper", "charge", "nonbonded")
)
def test_missing_required_record_no_zero_or_analog_fallback(tmp_path, family):
    def remove(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == family)
        ]

    specs = (
        [next(s for s in SPECS if s["id"] == "aromatic")]
        if family == "improper"
        else None
    )
    result, _, _, _ = run_request(tmp_path, specs=specs, change=remove)
    assert result["decision"] == "parameter_bundle_incomplete"
    missing = [m for r in result["graph_reports"] for m in r["missing_terms"]]
    assert any(m["family"] == family for m in missing)
    assert result["accepted_coefficients"] == 0 and not result["runtime_registration"]


@pytest.mark.parametrize(
    "name", ("scaling_14_lj", "scaling_14_coulomb", "charge_origin")
)
def test_missing_policy_or_charge_provenance_blocked(tmp_path, name):
    def remove(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "convention" and r["value"]["name"] == name)
        ]

    result, _, _, _ = run_request(tmp_path, change=remove)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert any(name in b for b in result["blocking_reasons"])


def test_missing_terminal_term_preserved_explicitly(tmp_path):
    def remove(_, rows):
        terminal = next(
            r
            for r in rows
            if r["kind"] == "parameter"
            and r["value"]["family"] == "bond"
            and "control_H" in r["value"]["atom_types"]
        )
        rows.remove(terminal)

    result, _, _, _ = run_request(tmp_path, change=remove)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert result["graph_reports"][0]["missing_terms"]


def test_parameter_unit_mutation_rejected(tmp_path):
    def mutate(_, rows):
        r = next(
            r
            for r in rows
            if r["kind"] == "parameter" and r["value"]["family"] == "bond"
        )
        r["value"]["form"]["units"][0] = "wrong_unit"

    result, _, _, _ = run_request(tmp_path, change=mutate)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert any(
        "equation_units_or_ordering_mismatch" in b for b in result["blocking_reasons"]
    )


def test_unsupported_motif_stays_blocked_with_software_complete_records(tmp_path):
    specs = [
        {
            "id": "unsupported",
            "repeat_units": [["A", "*CB(O)C*"]],
            "degree_of_polymerization": 1,
        }
    ]
    result, _, _, _ = run_request(tmp_path, specs=specs)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert "unsupported_motif_not_in_catalog" in result["blocking_reasons"]


@pytest.mark.parametrize(
    "mutation", ("bytes", "graph", "profile", "family", "version", "forms", "contract")
)
def test_request_graph_and_profile_binding(tmp_path, mutation):
    _, cap, _ = source_set(tmp_path)
    path, digest = request_for(tmp_path / "graphs", cap)
    request = j87.read(path)
    if mutation == "bytes":
        f = path.parent / request["fixtures"][0]["path"]
        f.write_bytes(f.read_bytes() + b"\n")
        with pytest.raises(ValueError, match="SHA256_mutation"):
            capture_request(path, digest)
        return
    if mutation == "graph":
        request["fixtures"][0]["graph_identity"] = "0" * 64
        request["graph_identities"] = ["0" * 64]
    elif mutation == "profile":
        request["profile_id"] = "not_this_source_profile"
    elif mutation == "family":
        request["family"] = "GAFF2"
    elif mutation == "version":
        request["expected_version"] = "different-release"
    elif mutation == "forms":
        request["expected_forms"] = []
    else:
        request["expected_conventions"] = {}
    digest = put(path, request)
    if mutation in ("graph", "profile", "family"):
        with pytest.raises(
            ValueError,
            match="exact_graph_scope_mismatch|profile_identity_mismatch|J87_OPLS_AA_only",
        ):
            qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)
    else:
        assert (
            qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)[
                "decision"
            ]
            == "parameter_bundle_incomplete"
        )


def test_relocated_source_and_request_reingestion(tmp_path):
    _, cap, args = source_set(tmp_path / "original")
    path, digest = request_for(tmp_path / "request", cap)
    source = tmp_path / "moved-source"
    fixture = tmp_path / "moved-request"
    shutil.copytree(args[0].parent, source)
    shutil.copytree(path.parent, fixture)
    review = tmp_path / "moved-audit.json"
    shutil.copyfile(args[2], review)
    fresh = capture_opls(
        source / "manifest.json", args[1], review, args[3], capture_timestamp=WHEN
    )
    first = qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)
    assert first == qualify_opls(
        capture_request(fixture / path.name, digest), fresh, j87.spec_from_data
    )


def test_cli_authorized_control_and_relocated_deterministic_replay(tmp_path):
    _, cap, args = source_set(tmp_path / "source")
    path, digest = request_for(tmp_path / "graphs", cap)
    target = tmp_path / "intake"
    command = [
        sys.executable,
        "-m",
        "scripts.phase_4j87_opls",
        "ingest",
        "--directory",
        str(target),
        "--request",
        str(path),
        "--request-sha256",
        digest,
        "--manifest",
        str(args[0]),
        "--manifest-sha256",
        args[1],
        "--audit",
        str(args[2]),
        "--audit-sha256",
        args[3],
        "--capture-timestamp",
        WHEN,
    ]
    run = subprocess.run(command, capture_output=True, text=True, check=False)
    assert run.returncode == 0, run.stderr + run.stdout
    output = json.loads(run.stdout)
    assert output["result"]["decision"] == READY
    other = tmp_path / "relocated"
    shutil.copytree(target, other)
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j87_opls",
            "replay",
            "--directory",
            str(other),
            "--snapshot-sha256",
            output["snapshot_sha256"],
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr + child.stdout
    assert (
        json.loads(child.stdout)
        == output["result"]
        == j87.replay(target, output["snapshot_sha256"])
    )


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("J87") / "evidence"
    j87.emit(REPO, root)
    return root


def test_default_receipt_is_source_blocked_for_eight_graphs(artifacts):
    result = j87.replay(artifacts)
    assert result["decision"] == "source_bundle_required"
    assert (
        result["exact_graph_fixtures"] == 8 and result["qualified_graph_reports"] == 0
    )
    assert (
        result["source_hashes"] == [] and result["scope"]["accepted_coefficients"] == 0
    )
    assert not result["source_authorized"] and result["family"] == "OPLS-AA"


@pytest.mark.parametrize(
    "mutation",
    ("scope", "decision", "missing", "locations", "conventions", "fixture", "receipt"),
)
def test_strict_receipt_rejects_resigned_semantic_mutations(
    tmp_path, monkeypatch, artifacts, mutation
):
    from scripts import validate_phase_4j87_receipt as validator

    root = tmp_path / "copy"
    shutil.copytree(artifacts, root)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, root)
    receipt = j87.read(root / "receipt.json")
    if mutation == "scope":
        receipt["scope"]["runtime_registration"] = True
    elif mutation == "decision":
        receipt["decision"] = READY
    else:
        name = {
            "missing": "graph-coverage-report.json",
            "locations": "source-locations.json",
            "conventions": "convention-comparison.json",
            "fixture": "qualification-request.json",
            "receipt": "qualification-receipt.json",
        }[mutation]
        obj = j87.read(root / name)
        if mutation == "missing":
            obj["reports"][0]["missing_terms"] = []
        elif mutation == "locations":
            obj["records"] = [{"location": "fabricated"}]
        elif mutation == "conventions":
            obj["performed"] = True
        elif mutation == "fixture":
            obj["fixture_files"][0]["sha256"] = "0" * 64
        else:
            obj["decision"] = READY
        j87.write(
            root / name, j87.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        manifest = j87.read(root / "manifest.json")
        manifest["artifacts"][name] = j87.pin(root / name)
        j87.write(
            root / "manifest.json",
            j87.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j87.pin(root / "manifest.json")
    j87.write(
        root / "receipt.json",
        j87.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, root)


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "lammps",
        "msi2lmp",
        "fit_parameters",
        "generate_parameters",
        "reserved_validation",
    ),
)
def test_prohibited_child_process_is_blocked(tmp_path, name):
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


def test_explicit_request_and_external_pin_required_before_ingestion(tmp_path):
    _, _, args = source_set(tmp_path / "source")
    with pytest.raises(ValueError, match="explicit_request_required"):
        j87.emit(
            REPO, tmp_path / "no-request", manifest=args[0], manifest_sha256=args[1]
        )
    assert not (tmp_path / "no-request").exists()
    with pytest.raises(ValueError, match="external_source_pin_required"):
        capture_request(REPO / j87.REQUEST_PATH, None)


def test_source_scope_cannot_be_extended_by_a_new_request(tmp_path):
    _, cap, _ = source_set(tmp_path / "source")
    specs = [{**SPECS[0], "degree_of_polymerization": 2}]
    path, digest = request_for(tmp_path / "new-graphs", cap, specs)
    result = qualify_opls(capture_request(path, digest), cap, j87.spec_from_data)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert "exact_graph_scope_not_evidenced" in result["blocking_reasons"]


def test_changed_audited_source_torsion_convention_still_needs_request_match(tmp_path):
    _, original, _ = source_set(tmp_path / "original")

    def change(_, rows):
        next(
            r
            for r in rows
            if r["kind"] == "convention" and r["value"]["name"] == "torsion_phase"
        )["value"]["value"] = "different_explicit_phase"

    _, changed, _ = source_set(tmp_path / "changed", change=change)
    path, digest = request_for(tmp_path / "request", original)
    result = qualify_opls(capture_request(path, digest), changed, j87.spec_from_data)
    assert result["decision"] == "parameter_bundle_incomplete"
    assert (
        "request_source_convention_mismatch:torsion_phase" in result["blocking_reasons"]
    )


@pytest.mark.parametrize(
    "mutation", ("fixture_bytes", "request_bytes", "duplicate", "symlink")
)
def test_fixture_capture_integrity_and_path_safety(tmp_path, mutation):
    _, cap, _ = source_set(tmp_path / "source")
    path, digest = request_for(tmp_path / "request", cap)
    captured = capture_request(path, digest)
    if mutation == "fixture_bytes":
        captured["files"][0]["raw_text"] += "\n"
        with pytest.raises(ValueError, match="captured_fixture_hash_mutation"):
            qualify_opls(captured, cap, j87.spec_from_data)
    elif mutation == "request_bytes":
        captured["request_raw"] += "\n"
        with pytest.raises(ValueError, match="captured_request_hash_mutation"):
            qualify_opls(captured, cap, j87.spec_from_data)
    elif mutation == "duplicate":
        captured["files"] *= 2
        with pytest.raises(ValueError, match="duplicate_captured_fixture"):
            qualify_opls(captured, cap, j87.spec_from_data)
    else:
        request = j87.read(path)
        file = path.parent / request["fixtures"][0]["path"]
        original = tmp_path / "original-fixture.json"
        file.rename(original)
        file.symlink_to(original)
        with pytest.raises(ValueError, match="fixture_path_symlink_forbidden"):
            capture_request(path, digest)
