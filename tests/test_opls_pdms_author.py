"""Author-source policies must not silently become the distinct paper Qeq model."""

import base64
import io
import json
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from hashlib import sha1, sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.forcefields.external_pdms import digest
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_author import (
    COMMIT,
    PINS,
    AuthorOPLSAudit,
    audit_author_opls,
    load_author_opls_audit,
)
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.forcefields.preparation import ForceFieldRequest, ForceFieldRequestError
from scripts.opls_pdms_author import fetch

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("mutation", [None, "commit_tree", "redirect"])
def test_explicit_fetch_binds_commit_tree_and_same_repository_only(
    inputs, monkeypatch, tmp_path, mutation
):
    options, _ = inputs
    tree_sha = "5b3bcc406226a9f5af9fcaee4abf5a465e9abc26"
    blobs, entries = {}, []
    for name in PINS:
        raw = (options["sources"] / name).read_bytes()
        blob_sha = sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        blobs[blob_sha] = {"content": base64.b64encode(raw).decode()}
        entries.append({"path": name, "sha": blob_sha, "size": len(raw)})
    calls = []

    def read(request, timeout):
        path = request.full_url.split("PDMS_benchmark", 1)[1]
        calls.append(path)
        if path == "/git/commits/" + COMMIT:
            payload = {
                "sha": COMMIT,
                "tree": {"sha": "0" * 40 if mutation == "commit_tree" else tree_sha},
            }
        elif path == "/git/trees/" + tree_sha + "?recursive=1":
            payload = {"sha": tree_sha, "truncated": False, "tree": entries}
        else:
            payload = blobs[path.rsplit("/", 1)[1]]
        response = io.BytesIO(json.dumps(payload).encode())
        response.url = (
            "https://api.github.com/repos/unrelated/parameters" + path
            if mutation == "redirect"
            else "https://api.github.com/repositories/902160948" + path
        )
        return response

    monkeypatch.setattr("urllib.request.urlopen", read)
    dest = tmp_path / "explicit_fetch"
    if mutation:
        reason = (
            "wrong_author_commit_tree"
            if mutation == "commit_tree"
            else "unexpected_source_redirect"
        )
        with pytest.raises(ValueError, match=reason):
            fetch(dest)
        assert not dest.exists() and len(calls) == 1
    else:
        result = fetch(dest)
        assert result["commit"] == COMMIT and result["tree"] == tree_sha
        assert calls[:2] == [
            "/git/commits/" + COMMIT,
            "/git/trees/" + tree_sha + "?recursive=1",
        ]
        for name, pin in PINS.items():
            assert sha256((dest / name).read_bytes()).hexdigest() == pin


@pytest.fixture(scope="module")
def inputs():
    sources = os.environ.get("ISLAND_PDMS_AUTHOR_SOURCES")
    docs = os.environ.get("ISLAND_LAMMPS_DOCS")
    if not sources or not docs:
        pytest.skip(
            "Requires explicit, already downloaded author sources and engine docs"
        )
    supplement = build_opls_pdms_literature(
        DomainSelection(
            (ROOT / "docs/evidence/phase_4j45/selection.json").read_text(),
            (ROOT / "docs/evidence/phase_4j44/source-audit.json").read_text(),
        )
    )
    options = {
        "sources": Path(sources),
        "supplement": supplement,
        "references": ROOT / "island_ffref",
        "engine_docs": Path(docs),
    }
    return options, audit_author_opls(**options)


def test_author_settings_resolve_policies_without_conflating_the_models(inputs):
    options, result = inputs
    p = result.payload
    assert p["commit"] == COMMIT
    assert p["author_variant"]["units"] == "real"
    assert p["author_variant"]["LJ_mixing"] == "geometric_sigma_and_epsilon"
    assert p["author_variant"]["LJ_scaling_by_bond_distance_1_2_3"] == [
        "0.0",
        "0.0",
        "0.5",
    ]
    assert p["author_variant"]["Coulomb_scaling_by_bond_distance_1_2_3"] == [
        "0.0",
        "0.0",
        "0.5",
    ]
    assert p["author_variant"]["pair_style"] == ["lj/charmm/coul/long", "9.0", "11.0"]
    assert p["active_commands"]["pair_style"]["line"] == 8
    assert p["paper_comparison"]["bonded_and_self_LJ_rows_match"]
    assert not p["paper_comparison"]["pair_style_matches"]
    assert not p["paper_comparison"]["charge_branch_matches_table_s4"]
    assert not p["paper_comparison"]["paper_Qeq_model_reproduced"]
    assert options["supplement"].payload["exact_charge_vector"] is None


def test_source_graph_partition_and_charge_precision_are_preserved(inputs):
    _, result = inputs
    g = result.payload["source_graph"]
    assert g["atom_count"] == 32448 and g["molecule_count"] == 64
    assert g["interaction_counts"] == {
        "Bonds": 32384,
        "Angles": 61504,
        "Dihedrals": 77568,
        "Impropers": 0,
    }
    assert all(
        len(m["terminal_silicon"]) == 2 and m["atoms"] == 507 for m in g["molecules"]
    )
    assert all(m["exact_serialized_charge_sum"] == "-4.50E-15" for m in g["molecules"])
    assert Decimal(g["exact_serialized_total_charge"]) == Decimal("-2.8800E-13")
    assert {q for c in g["type_charge_histograms"].values() for q in c} == {
        format(float(q), ".16e")
        for q in (
            "-1.261",
            "-1.255",
            "2.121",
            "1.874",
            "-1.102",
            "-1.083",
            "0.224",
            "0.2225",
        )
    }


def test_equation_preserving_conversions_and_evidenced_zero_rows(inputs):
    _, result = inputs
    rows = result.payload["selected_coefficient_rows"]
    assert len(rows) == 17
    bond = next(r for r in rows if r["family"] == "bond" and r["labels"] == ["o", "si"])
    assert Decimal(bond["converted_values"]["openmm_k_kJ_per_mol_nm2"]) == Decimal(
        "312963.2"
    )
    assert Decimal(bond["converted_values"]["r0_nm"]) == Decimal("0.164")
    torsion = next(
        r
        for r in rows
        if r["family"] == "proper_torsion" and Decimal(r["original_values"][2]) != 0
    )
    assert Decimal(torsion["converted_values"]["opls_K_kJ_per_mol"][2]) == Decimal(
        "0.75312"
    )
    assert "/2" in torsion["converted_values"]["equation"]
    for row in rows:
        assert (
            row["data_evidence"]["raw_text"]
            and row["selected_settings_evidence"]["raw_text"]
        )
        assert not row["execution_authorized"]
        assert row["evidence_checksum"] == digest(
            {k: v for k, v in row.items() if k != "evidence_checksum"}
        )


@pytest.mark.parametrize(
    "change",
    [
        "pair_style",
        "charge_branch",
        "conversion",
        "source",
        "missing_row",
        "zero",
        "parity",
        "authority",
        "profile",
        "status",
    ],
)
def test_resigning_cannot_forge_author_policies_or_authority(change, inputs):
    options, result = inputs
    p = result.payload
    if change == "pair_style":
        p["author_variant"]["pair_style"] = ["lj/cut/coul/long", "11", "11"]
    elif change == "charge_branch":
        p["author_variant"]["charge_method"] = "Qeq"
    elif change == "conversion":
        p["selected_coefficient_rows"][0]["converted_values"]["epsilon_kJ_per_mol"] = (
            "0"
        )
    elif change == "source":
        p["source_hashes"]["param/OPLS/opls_system.in.init"] = "f" * 64
    elif change == "missing_row":
        p["selected_coefficient_rows"].pop()
    elif change == "zero":
        next(r for r in p["selected_coefficient_rows"] if r["family"] == "bond")[
            "original_values"
        ] = ["0", "0"]
    elif change == "parity":
        p["paper_comparison"]["paper_Qeq_model_reproduced"] = True
    elif change == "authority":
        p["execution_authorized"] = True
    elif change == "profile":
        p["profile_identity"] = "native_oplsaa"
    else:
        p["status"] = "prepared"
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        AuthorOPLSAudit(
            pack(p), options["sources"], options["supplement"], options["engine_docs"]
        )
    assert error.value.status == "source_incompatible"


def test_offline_audit_does_not_call_network_or_numerical_preparation(
    inputs, monkeypatch
):
    options, _ = inputs

    def forbidden(*args, **kwargs):
        pytest.fail("offline source audit reached network or numerical preparation")

    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    monkeypatch.setattr("island.forcefields.preparation.prepare_forcefield", forbidden)
    result = audit_author_opls(**options)
    with pytest.raises(ForceFieldRequestError):
        ForceFieldRequest("oplsaa", result)
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        result.require_executable()
    assert (
        error.value.status == "scientifically_unqualified"
        and error.value.result is result
    )
    assert result.payload["profile_identity"] is None
    assert result.payload["backend_executions"] == []


def test_stale_author_source_rejected_even_if_only_a_comment_changes(inputs, tmp_path):
    options, result = inputs
    relocated = tmp_path / "source"
    for name in PINS:
        dest = relocated / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(options["sources"] / name, dest)
    altered = relocated / "param/OPLS/opls_system.in.init"
    altered.write_bytes(altered.read_bytes() + b"\n# altered\n")
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        audit_author_opls(**{**options, "sources": relocated})
    assert error.value.status == "source_incompatible"
    path = tmp_path / "record.json"
    path.write_text(result.json_text)
    with pytest.raises(OPLSPDMSLiteratureError):
        load_author_opls_audit(path, None, **options)


def test_fresh_process_relocation_requires_source_and_independent_record_pin(
    inputs, tmp_path
):
    options, result = inputs
    moved = tmp_path / "relocated"
    moved.mkdir()
    for label, original, names in (
        ("source", options["sources"], PINS),
        ("docs", options["engine_docs"], result.payload["engine_documentation"]),
        (
            "references",
            options["references"],
            result.payload["paper_comparison"]["source_evidence"]["source_hashes"],
        ),
    ):
        for name in names:
            dest = moved / label / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original / name, dest)
    record = moved / "record.json"
    record.write_text(result.json_text)
    selection, catalog = moved / "selection.json", moved / "catalog.json"
    selection.write_text(options["supplement"].selection.json_text)
    catalog.write_text(options["supplement"].selection.catalog_text)
    run = subprocess.run(
        [
            sys.executable,
            "scripts/opls_pdms_author.py",
            "--sources",
            str(moved / "source"),
            "--engine-docs",
            str(moved / "docs"),
            "--references",
            str(moved / "references"),
            "--selection",
            str(selection),
            "--catalog",
            str(catalog),
            "--reconstruct",
            str(record),
            "--sha256",
            sha256(record.read_bytes()).hexdigest(),
            "--require-executable",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 1 and not run.stderr
    report = json.loads(run.stdout)
    assert report["identity"] == result.payload["identity"]
    assert report["status"] == "scientifically_unqualified"
    assert report["execution_profiles"] == report["backend_executions"] == []
    assert not any(p.suffix in (".chk", ".dcd", ".xml") for p in moved.rglob("*"))
