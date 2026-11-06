"""J86 software controls: manufactured bundles/audits are not scientific sources."""

import json
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.forcefields.assignment_core.profiles import profile_contract
from island.forcefields.assignment_core.schema import FAMILIES
from island.forcefields.bundle_adapters import capture
from island.forcefields.bundle_adapters.scoped_qualification import (
    QUALIFIED,
    QualificationRequest,
    ScopedQualification,
    aggregate_qualifications,
    qualify_request,
    rank_qualifications,
)
from scripts import phase_4j86_qualification as j86
from tests.test_phase_4j82 import graph
from tests.test_phase_4j84 import WHEN
from tests.test_phase_4j85 import supplied

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j86.scope()


def request(g, family="OPLS-AA", profile_id=None):
    return QualificationRequest(
        family, profile_id or profile_contract(family).id, (g.identity,)
    )


def make_set(g, captures):
    results = []
    for family in FAMILIES:
        cap = captures.get(family, capture(capture_timestamp=WHEN))
        profile_id = cap["manifest"]["profile_id"] if cap["manifest"] else None
        r = request(g, family, profile_id)
        results.append(qualify_request(cap, (g,), request=r))
    return results


def aggregate(results):
    results = tuple(results)
    return aggregate_qualifications(
        results, qualification_set=[r.request for r in results]
    )


def inputs(caps):
    value = j86.read(REPO / j86.FIXTURE)
    value["bundle_captures"] = caps
    for declaration in value["required_profiles"]:
        for cap in caps:
            if declaration["family"] == cap["manifest"]["family"]:
                declaration["profile_id"] = cap["manifest"]["profile_id"]
    return value


def test_authorized_opls_and_four_missing_profiles_are_partial(tmp_path):
    g, cap, _ = supplied(tmp_path, authorized=True)
    results = make_set(g, {"OPLS-AA": cap})
    report = aggregate(results)
    assert report.decision == "partial_bundle_set"
    assert len(report.missing_or_unauthorized_profiles) == 4
    assert len(report.qualified_profiles) == 1 and not report.incomplete_profiles
    assert results[0].qualification.decision == QUALIFIED
    assert all(
        r.qualification.decision == "source_bundle_required" for r in results[1:]
    )
    assert aggregate(iter(results)) == report
    assert not report.runtime_registration and not report.automatic_assignment_enabled
    assert not report.production_validated and report.accepted_coefficients == 0
    assert report.simulation_readiness == "not_established"
    assert report.chemical_family_transfer == "not_established"


def test_incomplete_pcff_cannot_hide_behind_complete_opls(tmp_path):
    g, opls, _ = supplied(tmp_path / "opls", authorized=True)

    def remove_cross(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == "bond_angle")
        ]

    _, pcff, _ = supplied(
        tmp_path / "pcff", "PCFF", g=g, change=remove_cross, authorized=True
    )
    results = make_set(g, {"OPLS-AA": opls, "PCFF": pcff})
    report = aggregate(results)
    assert report.decision == "partial_bundle_set"
    assert report.incomplete_profiles == (("PCFF", pcff["manifest"]["profile_id"]),)
    subset = [r for r in results if r.request.family in ("OPLS-AA", "PCFF")]
    assert aggregate(subset).decision == "parameter_bundle_incomplete"
    pcff_result = next(r.qualification for r in results if r.request.family == "PCFF")
    assert pcff_result.decision == "parameter_bundle_incomplete"
    assert json.loads(pcff_result.graph_reports_json)[0]["missing_terms"]


@pytest.mark.parametrize("family", FAMILIES)
def test_one_explicit_profile_scope_can_qualify_without_family_transfer(
    tmp_path, family
):
    g, cap, _ = supplied(tmp_path, family, authorized=True)
    r = qualify_request(
        cap, (g,), request=request(g, family, cap["manifest"]["profile_id"])
    )
    report = aggregate((r,))
    assert report.decision == "family_scoped_qualification_ready_runtime_blocked"
    assert report.scope_of_qualification == "declared_exact_graphs_only"
    assert report.chemical_family_transfer == "not_established"
    assert report.accepted_coefficients == 0 and not report.runtime_registration


def test_all_required_profiles_must_be_complete(tmp_path):
    g = graph()
    caps = {f: supplied(tmp_path / f, f, g=g, authorized=True)[1] for f in FAMILIES}
    report = aggregate(make_set(g, caps))
    assert report.decision == "family_scoped_qualification_ready_runtime_blocked"
    assert len(report.qualified_profiles) == 5
    assert (
        not report.missing_or_unauthorized_profiles and not report.incomplete_profiles
    )
    assert not report.runtime_registration and report.accepted_coefficients == 0


def test_opls_and_gaff2_rank_separately(tmp_path):
    outcomes = {}
    for f in ("OPLS-AA", "GAFF2"):
        g, cap, _ = supplied(tmp_path / f, f, authorized=True)
        outcomes[f] = qualify_request(
            cap, (g,), request=request(g, f, cap["manifest"]["profile_id"])
        ).qualification
    assert rank_qualifications((outcomes["OPLS-AA"],), requested_family="opls") == (
        outcomes["OPLS-AA"],
    )
    assert rank_qualifications((outcomes["GAFF2"],), requested_family="gaff2") == (
        outcomes["GAFF2"],
    )
    for family in ("opls", "GAFF2"):
        with pytest.raises(ValueError, match="cross_family_ranking_forbidden"):
            rank_qualifications(outcomes.values(), requested_family=family)
    with pytest.raises(TypeError):
        rank_qualifications(outcomes.values())


@pytest.mark.parametrize("family", ("GAFF", "GAFF2", "PCFF", "COMPASS"))
def test_opls_cannot_select_other_family(family):
    r = make_set(graph(), {})[FAMILIES.index(family)].qualification
    with pytest.raises(ValueError, match="cross_family_ranking_forbidden"):
        rank_qualifications((r,), requested_family="opls")


def test_deterministic_within_family_tie_breaking(tmp_path):
    results = []
    for name in ("z", "a"):
        g, cap, _ = supplied(
            tmp_path / name,
            change=lambda m, _, n=name: m.update(profile_id=n),
            authorized=True,
        )
        results.append(
            qualify_request(
                cap, (g,), request=request(g, profile_id=name)
            ).qualification
        )
    ranked = rank_qualifications(iter(results), requested_family="OPLS-AA")
    assert ranked == rank_qualifications(reversed(results), requested_family="opls")
    assert [r.profile_id for r in ranked] == ["a", "z"]
    with pytest.raises(ValueError, match="duplicate_candidate_bundle"):
        rank_qualifications([results[0]] * 2, requested_family="opls")


@pytest.mark.parametrize("mode", ("missing", "extra", "duplicate", "empty"))
def test_declared_set_cannot_lose_required_results(mode):
    rows = make_set(graph(), {})
    declared = [r.request for r in rows]
    if mode == "missing":
        rows.pop()
    elif mode == "extra":
        declared.pop()
    elif mode == "duplicate":
        rows.append(rows[0])
    else:
        rows, declared = [], []
    with pytest.raises(
        ValueError,
        match="coverage_mismatch|duplicate_profile_result|qualification_set_required",
    ):
        aggregate_qualifications(rows, qualification_set=declared)


def test_duplicate_required_profile_rejected():
    rows = make_set(graph(), {})
    with pytest.raises(ValueError, match="duplicate_required_profile"):
        aggregate_qualifications(
            rows, qualification_set=[r.request for r in rows] + [rows[0].request]
        )


def test_wrong_graph_or_profile_or_family_cannot_qualify(tmp_path):
    g, cap, _ = supplied(tmp_path, authorized=True)
    r = request(g, profile_id=cap["manifest"]["profile_id"])
    result = qualify_request(cap, (g,), request=r)
    with pytest.raises(ValueError, match="declared_profile_identity_mismatch"):
        ScopedQualification(replace(r, profile_id="another"), result.qualification)
    with pytest.raises(ValueError, match="cross_family_qualification_forbidden"):
        ScopedQualification(replace(r, family="GAFF"), result.qualification)
    with pytest.raises(ValueError, match="declared_graph_scope_mismatch"):
        qualify_request(cap, (graph(dp=2),), request=r)
    with pytest.raises(ValueError, match="candidate_graph_scopes_differ"):
        rank_qualifications(
            (result.qualification, make_set(graph(dp=2), {})[0].qualification),
            requested_family="opls",
        )


@pytest.mark.parametrize("value", ("", "opls/gaff", "PCFF/COMPASS", None))
def test_implicit_family_switching_is_not_an_alias(value):
    with pytest.raises(ValueError, match="unknown_requested_family"):
        rank_qualifications((), requested_family=value)


@pytest.mark.parametrize(
    "field,value",
    (
        ("runtime_registration", True),
        ("automatic_assignment_enabled", True),
        ("production_validated", True),
        ("simulation_readiness", "ready"),
        ("accepted_coefficients", 1),
    ),
)
def test_even_partial_report_gates_cannot_be_enabled(field, value):
    report = aggregate(make_set(graph(), {}))
    with pytest.raises(ValueError, match="init=False"):
        replace(report, **{field: value})


def test_unauthorized_synthetic_does_not_count_as_qualified(tmp_path):
    g, cap, _ = supplied(tmp_path)
    report = aggregate(make_set(g, {"OPLS-AA": cap}))
    assert report.decision == "source_bundle_required"
    assert (
        not report.qualified_profiles
        and len(report.missing_or_unauthorized_profiles) == 5
    )


def test_partial_report_deterministic_and_relocated(tmp_path):
    g, cap, args = supplied(tmp_path / "sources", authorized=True)
    value = inputs([cap])
    outputs = j86.reconstruct(value)
    assert outputs["result.json"]["decision"] == "partial_bundle_set"
    value["required_profiles"].reverse()
    assert outputs == j86.reconstruct(value)
    declaration = tmp_path / "set.json"
    without_capture = inputs([])
    without_capture["required_profiles"][0]["profile_id"] = cap["manifest"][
        "profile_id"
    ]
    j86.write(declaration, without_capture)
    first = tmp_path / "intake"
    emitted = j86.emit(REPO, first, cap, declaration, j86.pin(declaration))
    assert (
        j86.replay(first, emitted["snapshot_sha256"])["decision"]
        == "partial_bundle_set"
    )
    relocated = tmp_path / "relocated"
    shutil.copytree(first, relocated)
    # Separate source re-ingestion hashes every original file after relocation.
    moved = tmp_path / "moved-source"
    shutil.copytree(args[0].parent, moved)
    audit = tmp_path / "review.json"
    shutil.copyfile(args[2], audit)
    fresh = capture(
        moved / "manifest.json", args[1], audit, args[3], capture_timestamp=WHEN
    )
    assert qualify_request(
        fresh, (g,), request=request(g, profile_id=cap["manifest"]["profile_id"])
    ) == qualify_request(
        cap, (g,), request=request(g, profile_id=cap["manifest"]["profile_id"])
    )
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j86_qualification",
            "replay",
            "--directory",
            str(relocated),
            "--snapshot-sha256",
            emitted["snapshot_sha256"],
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr + child.stdout
    assert json.loads(child.stdout) == emitted["result"]


def test_custom_declaration_requires_external_pin(tmp_path):
    declaration = tmp_path / "set.json"
    shutil.copyfile(REPO / j86.FIXTURE, declaration)
    with pytest.raises(ValueError, match="declaration_pin_mismatch"):
        j86.emit(REPO, tmp_path / "output", declaration=declaration)
    assert not (tmp_path / "output").exists()


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    target = tmp_path_factory.mktemp("j86") / "evidence"
    j86.emit(REPO, target)
    return target


def test_no_source_default_exposes_all_five_profiles(artifacts):
    result = j86.replay(artifacts)
    assert result["decision"] == "source_bundle_required"
    assert len(result["missing_or_unauthorized_profiles"]) == 5
    assert result["required_profiles"] == 5 and result["graph_fixtures"] == 1
    assert result["source_hashes"] == [] and result["accepted_coefficients"] == 0
    rankings = j86.read(artifacts / "source-qualification.json")[
        "rankings_by_requested_family_and_graph_scope"
    ]
    assert set(rankings) == set(FAMILIES)
    assert len(j86.read(artifacts / "class-i-completeness.json")["reports"]) == 3
    assert len(j86.read(artifacts / "class-ii-completeness.json")["reports"]) == 2


@pytest.mark.parametrize(
    "mutation", ("scope", "decision", "missing", "ranking", "class", "receipt")
)
def test_strict_validator_rejects_resigned_semantic_changes(
    tmp_path, monkeypatch, artifacts, mutation
):
    from scripts import validate_phase_4j86_receipt as validator

    root = tmp_path / "copy"
    shutil.copytree(artifacts, root)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, root)
    receipt = j86.read(root / "receipt.json")
    if mutation == "scope":
        receipt["scope"]["automatic_assignment_enabled"] = True
    elif mutation == "decision":
        receipt["decision"] = "family_scoped_qualification_ready_runtime_blocked"
    else:
        name = (
            "qualification-receipt.json"
            if mutation == "receipt"
            else "source-qualification.json"
            if mutation in ("missing", "ranking")
            else "class-ii-completeness.json"
        )
        obj = j86.read(root / name)
        if mutation == "receipt":
            obj["decision"] = "partial_bundle_set"
        elif mutation == "missing":
            obj["aggregate"]["missing_or_unauthorized_profiles"] = []
        elif mutation == "ranking":
            obj["rankings_by_requested_family_and_graph_scope"]["OPLS-AA"] = [
                {
                    "graph_identities": [],
                    "candidates": [{"family": "GAFF2", "profile_id": "forged"}],
                }
            ]
        else:
            obj["reports"] = []
        j86.write(
            root / name, j86.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        m = j86.read(root / "manifest.json")
        m["artifacts"][name] = j86.pin(root / name)
        j86.write(
            root / "manifest.json",
            j86.signed({k: v for k, v in m.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j86.pin(root / "manifest.json")
    j86.write(
        root / "receipt.json",
        j86.signed({k: v for k, v in receipt.items() if k != "identity"}),
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
def test_prohibited_children_are_rejected(tmp_path, name):
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


def test_multiple_missing_profiles_within_one_family_remain_visible():
    value = inputs([])
    value["required_profiles"] = [
        {"family": "OPLS-AA", "profile_id": name, "graphs": ["hydrocarbon"]}
        for name in ("source_a", "source_b")
    ]
    output = j86.reconstruct(value)
    result = output["result.json"]
    assert result["decision"] == "source_bundle_required"
    assert result["missing_or_unauthorized_profiles"] == [
        ["OPLS-AA", "source_a"],
        ["OPLS-AA", "source_b"],
    ]
    assert output["source-qualification.json"][
        "rankings_by_requested_family_and_graph_scope"
    ] == {"OPLS-AA": []}


def test_forged_individual_decision_and_changed_required_scope_rejected():
    rows = make_set(graph(), {})
    r = rows[0]
    with pytest.raises(ValueError, match="inconsistent_per_profile_decision"):
        ScopedQualification(r.request, replace(r.qualification, decision=QUALIFIED))
    declared = [replace(r.request, graph_identities=("0" * 64,))] + [
        r.request for r in rows[1:]
    ]
    with pytest.raises(ValueError, match="qualification_set_graph_scope_mismatch"):
        aggregate_qualifications(rows, qualification_set=declared)


@pytest.mark.parametrize("graphs", ((), ("x" * 64,), ("0" * 64, "0" * 64)))
def test_malformed_graph_declaration_rejected(graphs):
    with pytest.raises(
        ValueError,
        match="graphs_required|invalid_graph_identity|duplicate_declared_graph",
    ):
        QualificationRequest("OPLS-AA", "explicit", graphs)


def test_aggregate_rechecks_complete_counts_and_authorized_identity(tmp_path):
    g, cap, _ = supplied(tmp_path, authorized=True)
    r = request(g, profile_id=cap["manifest"]["profile_id"])
    q = qualify_request(cap, (g,), request=r).qualification
    reports = json.loads(q.graph_reports_json)
    reports[0]["coverage"]["bond"]["covered"] -= 1
    with pytest.raises(ValueError, match="inconsistent_per_profile_decision"):
        ScopedQualification(r, replace(q, graph_reports_json=json.dumps(reports)))
    with pytest.raises(ValueError, match="authorized_source_identity_required"):
        ScopedQualification(r, replace(q, source_hashes=()))
