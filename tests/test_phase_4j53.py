"""New-source contracts are not source authority or executable profiles."""

import copy
import json
import os
import shutil
import subprocess
import sys
import zipfile
from dataclasses import FrozenInstanceError, replace
from hashlib import sha256

import pytest

from island.charge_references.records import pack
from island.forcefields.opls_registry import OPLSAssets, OPLSResolutionError
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import (
    canonical,
    digest,
    final_graph,
    requirements,
)
from island.opls_final_resolution import (
    NINE_CONVENTIONS,
    ConventionSources,
    corpus_audit,
    use_convention_sources,
)
from island.opls_source_execution import reconstruct, run_md, save_preparation
from island.opls_source_intake import (
    STATUSES,
    SourcePackage,
    blocker_reopen_report,
    inspect_source,
    inventory,
    prepare_polymer_forcefield,
    use_source_package,
)
from island.opls_source_intake.schema import (
    REQUIRED_CONVENTIONS,
    SCHEMA,
    VERSION,
    source_intake_schema,
)
from island.polymer import use_opls_assets
from island.workflows.bundle import system_data
from tests.test_phase_4j50 import ARCHIVE, ARCHIVE_SHA256, fixture
from tests.test_phase_4j51 import PAPER, SUPPLEMENT


def pin(path):
    return sha256(path.read_bytes()).hexdigest()


def write_manifest(root, document):
    (root / "source-intake.json").write_text(canonical(document))
    return SourcePackage(root, pin(root / "source-intake.json"))


def refresh(root, document):
    roles = {r["path"]: r["role"] for r in document["files"]}
    document["files"] = [
        {**r, "role": roles.get(r["path"], "author_parameters")}
        for r in inventory(root)
        if r["path"] != "source-intake.json"
    ]
    document["authorization"]["covered_files"] = [
        {"file": r["path"], "sha256": r["sha256"]}
        for r in document["files"]
        if r["role"] != "permission"
    ]
    return write_manifest(root, document)


@pytest.fixture(scope="module")
def contract(tmp_path_factory):
    """Synthetic contract data, never scientific PDMS qualification or registration."""
    pytest.importorskip("rdkit")
    root = tmp_path_factory.mktemp("synthetic-contract")
    system = fixture(0)
    (root / "system.json").write_text(pack(system_data(system)))
    (root / "permission.txt").write_text(
        "TEST ONLY permission for synthetic contract\n"
    )
    (root / "version.txt").write_text("synthetic_v1\n")
    (root / "code.py").write_text(
        "# Fixture evidence; not executed or scientifically verified\n"
    )
    (root / "nested").mkdir()
    (root / "nested" / "empty.txt").write_bytes(b"")
    graph = final_graph(system)
    order = graph["atom_order"]
    values = []

    def record(key, value, units, labels):
        raw = canonical({"key": key, "value": value, "units": units})
        values.append(raw)
        return {
            "key": key,
            "family": "opls",
            "labels": labels,
            "functional_form": "synthetic_non_executable_contract",
            "resolution": "source_resolved",
            "value": value,
            "original_value": value,
            "units": units,
            "conversion": {
                "formula": "identity",
                "implementation": "test",
                "version": "1",
            },
            "line": len(values),
        }

    conventions = [
        record(
            key,
            "synthetic_v1" if key == "charge_method" else "TEST ONLY",
            "dimensionless",
            [],
        )
        for key in REQUIRED_CONVENTIONS
    ]
    rows = []
    for request in requirements(system, "opls"):
        kind, sites = request["family"], request["sites"]
        value, units = "TEST ONLY", "dimensionless"
        if kind == "atom_type":
            value = "synthetic_type"
        elif kind == "charge":
            value, units = "0.125", "elementary_charge"
        elif kind == "mass":
            value, units = str(system.topology.sites[sites[0]].mass), "dalton"
        rows.append(
            record(request["id"], value, units, ["synthetic_type"] * len(sites))
        )
    (root / "evidence.txt").write_text("\n".join(values) + "\n")
    doi, citation = "10.0000/contract-test", "Synthetic contract fixture, not a paper"

    def evidence(name, line):
        return {
            "file": name,
            "sha256": pin(root / name),
            "location": f"line {line}",
            "first_line": line,
            "raw_text": (root / name).read_text().splitlines()[line - 1],
            "doi": doi,
            "citation": citation,
        }

    for row in conventions + rows:
        line = row.pop("line")
        row["evidence"] = [evidence("evidence.txt", line)]
        row["reconstruction"] = {
            "code": evidence("code.py", 1),
            "result": evidence("evidence.txt", line),
            "implementation": "test",
            "version": "1",
            "value": row["value"],
            "units": row["units"],
        }
        row["evidence_checksum"] = digest(row)
    document = {
        "schema": VERSION,
        "family": "opls",
        "source_kind": "external_literature",
        "proposed_profile": "island_external_opls_synthetic_test_v1",
        "source_version": "synthetic_v1",
        "doi": doi,
        "citation": citation,
        "authorization": {
            "status": "permitted",
            "scope": "local_research",
            "rights_holder": "synthetic fixture",
            "redistribution": "not_established",
            "evidence": evidence("permission.txt", 1),
        },
        "version_evidence": evidence("version.txt", 1),
        "conventions": conventions,
        "domains": [
            {
                "name": "synthetic_graph_contract_only",
                "system_file": "system.json",
                "graph_identity": digest(graph),
                "typed_graph_identity": digest(
                    {
                        "graph_identity": digest(graph),
                        "atom_types": [(i, "synthetic_type") for i in order],
                    }
                ),
                "atom_types": ["synthetic_type"] * len(order),
                "charges": [0.125] * len(order),
                "charge_units": "elementary_charge",
                "charge_origin": "synthetic, never a PDMS model",
                "charge_method": "synthetic_v1",
                "rows": rows,
            }
        ],
        "files": [],
    }
    roles = {
        "system.json": "system",
        "permission.txt": "permission",
        "version.txt": "version",
        "code.py": "independent_reconstruction",
        "evidence.txt": "extracted_evidence",
    }
    document["files"] = [
        {**r, "role": roles.get(r["path"], "author_parameters")}
        for r in inventory(root)
    ]
    document["authorization"]["covered_files"] = [
        {"file": r["path"], "sha256": r["sha256"]}
        for r in document["files"]
        if r["role"] != "permission"
    ]
    package = write_manifest(root, document)
    assert inspect_source(package)["status"] == "source_complete_not_executed"
    return root, document


@pytest.fixture
def package(contract, tmp_path):
    root, document = contract
    path = tmp_path / "package"
    shutil.copytree(root, path)
    doc = copy.deepcopy(document)
    return path, doc, write_manifest(path, doc)


@pytest.fixture(scope="module")
def prepared():
    if not all(p.is_file() for p in (ARCHIVE, PAPER, SUPPLEMENT)):
        pytest.skip("Explicit pinned local source assets unavailable")
    system = fixture(2)
    with (
        use_opls_assets(OPLSAssets(ARCHIVE, ARCHIVE_SHA256)),
        use_convention_sources(ConventionSources(PAPER, SUPPLEMENT)),
    ):
        preparation = prepare_polymer_forcefield(system, forcefield="opls")
    return system, preparation


def test_exact_statuses_and_schema():
    assert STATUSES == (
        "not_supplied",
        "supplied_unreadable",
        "supplied_incomplete",
        "supplied_conflicting",
        "source_complete_not_executed",
        "source_complete_executed",
        "source_complete_independently_verified",
    )
    assert SCHEMA["additionalProperties"] is False
    assert set(NINE_CONVENTIONS) <= set(REQUIRED_CONVENTIONS)


def test_schema_copies_cannot_change_intake_contract(package):
    root, doc, _ = package
    schema = source_intake_schema()
    schema["properties"]["family"] = {"const": "pcff"}
    schema["additionalProperties"] = True
    doc["family"] = "pcff"
    assert source_intake_schema()["properties"]["family"] == {"const": "opls"}
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_incomplete"


def test_result_copies_cannot_reopen_or_promote_a_source(package):
    _, _, spec = package
    result = inspect_source(spec)
    result["status"] = "source_complete_independently_verified"
    result["MD_authorized"] = True
    original = inspect_source(spec)
    assert original["status"] == "source_complete_not_executed"
    assert original["MD_authorized"] is False


def test_not_supplied_and_missing_directory(tmp_path):
    assert inspect_source()["status"] == "not_supplied"
    assert (
        inspect_source(SourcePackage(tmp_path / "absent", "0" * 64))["status"]
        == "supplied_unreadable"
    )


def test_structurally_complete_is_not_executed_or_registered(package):
    _, _, spec = package
    result = inspect_source(spec)
    assert result["status"] == "source_complete_not_executed"
    assert result["new_profiles"] == [] and not result["automatically_selectable"]
    assert (
        result["execution_evidence"]
        is result["independent_verification_evidence"]
        is None
    )
    assert not result["MD_authorized"] and result["candidate_domains"]


@pytest.mark.parametrize("setting", NINE_CONVENTIONS)
def test_exact_permanent_blocker_preserved(setting):
    report = blocker_reopen_report(inspect_source())
    row = next(r for r in report["rows"] if r["setting"] == setting)
    original = next(
        r for r in corpus_audit()["final_resolutions"] if r["setting"] == setting
    )
    assert row["permanent_blocker"] == original
    assert row["required_new_artifact"] == original["exact_missing_artifact"]
    assert not row["reopened"] and row["candidate_evidence"] is None


@pytest.mark.parametrize("key", REQUIRED_CONVENTIONS)
def test_missing_each_execution_convention(package, key):
    root, doc, _ = package
    doc["conventions"] = [r for r in doc["conventions"] if r["key"] != key]
    result = inspect_source(write_manifest(root, doc))
    assert result["status"] == "supplied_incomplete"
    assert key in result["diagnostics"][0]["missing"]


@pytest.mark.parametrize(
    "key",
    [
        "distance_units",
        "dielectric",
        "special_bonds",
        "torsion_convention",
        "pair_cutoff_angstrom",
        "reciprocal_tolerance",
        "lj_mixing",
    ],
)
def test_resigned_convention_mutation_has_no_reconstruction(package, key):
    root, doc, _ = package
    row = next(r for r in doc["conventions"] if r["key"] == key)
    row["value"] = "changed"
    row["evidence_checksum"] = digest(
        {k: v for k, v in row.items() if k != "evidence_checksum"}
    )
    result = inspect_source(write_manifest(root, doc))
    assert result["status"] == "supplied_conflicting"
    assert result["diagnostics"][0]["code"] == "independent_reconstruction_conflict"


@pytest.mark.parametrize(
    "mutation",
    [
        "source_replace",
        "source_hash",
        "permission_file",
        "unit",
        "charge",
        "graph",
        "wrong_family",
        "manifest",
        "atom_type",
    ],
)
def test_negative_mutations(package, mutation):
    root, doc, spec = package
    if mutation == "source_replace":
        (root / "evidence.txt").write_text("replacement\n")
    elif mutation == "source_hash":
        doc["files"][0]["sha256"] = "0" * 64
    elif mutation == "permission_file":
        (root / "permission.txt").write_text("permission revoked\n")
    elif mutation == "unit":
        doc["conventions"][0]["units"] = "changed"
    elif mutation == "charge":
        doc["domains"][0]["charges"][0] += 0.01
    elif mutation == "graph":
        doc["domains"][0]["graph_identity"] = "0" * 64
    elif mutation == "wrong_family":
        doc["family"] = "pcff"
    elif mutation == "atom_type":
        doc["domains"][0]["atom_types"][0] = "other"
    else:
        doc["source_version"] = "replaced"
    if mutation not in ("source_replace", "permission_file"):
        new_spec = write_manifest(root, doc)
        if mutation != "manifest":
            spec = new_spec
    result = inspect_source(spec)
    assert result["status"] in ("supplied_conflicting", "supplied_incomplete")
    assert not result["automatically_selectable"]


@pytest.mark.parametrize("status", ["unknown", "prohibited"])
def test_permission_required_even_with_repin(package, status):
    root, doc, _ = package
    doc["authorization"]["status"] = status
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_incomplete"


def test_article_permission_cannot_implicitly_cover_an_archive(package):
    root, doc, _ = package
    doc["authorization"]["covered_files"].pop()
    result = inspect_source(write_manifest(root, doc))
    assert result["status"] == "supplied_conflicting"
    assert result["diagnostics"][0]["code"] == "permission_scope_incomplete"


@pytest.mark.parametrize("field", ["claimed_status", "qualified", "backend_defaults"])
def test_claimed_verification_or_defaults_not_accepted(package, field):
    root, doc, _ = package
    doc[field] = "source_complete_independently_verified"
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_incomplete"


def test_no_backend_resolution_promotion(package):
    root, doc, _ = package
    doc["conventions"][0]["resolution"] = "backend_resolved"
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_incomplete"


def test_invalid_graph_envelope_returns_structured_diagnostic(package):
    root, doc, _ = package
    (root / "system.json").write_text("{}")
    result = inspect_source(refresh(root, doc))
    assert result["status"] == "supplied_incomplete"
    assert result["diagnostics"][0]["code"] == "invalid_source_structure"


def test_duplicate_JSON_fields_are_conflicting(package):
    root, doc, _ = package
    raw = canonical(doc)
    (root / "source-intake.json").write_text('{"family":"opls",' + raw[1:])
    result = inspect_source(SourcePackage(root, pin(root / "source-intake.json")))
    assert result["status"] == "supplied_conflicting"
    assert result["diagnostics"][0]["code"] == "duplicate_JSON_key"


def test_incomplete_or_conflicting_interactions(package):
    root, doc, _ = package
    removed = doc["domains"][0]["rows"].pop()
    result = inspect_source(write_manifest(root, doc))
    assert removed["key"] in result["diagnostics"][0]["missing"]
    doc["domains"][0]["rows"].append(doc["domains"][0]["rows"][0])
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_conflicting"


def test_recursive_inventory_omission_symlink_and_path(package):
    root, doc, _ = package
    (root / "nested" / "hidden").write_text("uninventoried")
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_conflicting"
    (root / "nested" / "hidden").unlink()
    (root / "link").symlink_to(root / "evidence.txt")
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_conflicting"
    (root / "link").unlink()
    doc["files"][0]["path"] = "../outside"
    assert inspect_source(write_manifest(root, doc))["status"] == "supplied_conflicting"


@pytest.mark.parametrize("member", ["../escape", "/absolute"])
def test_archive_path_rejection(package, member):
    root, doc, _ = package
    with zipfile.ZipFile(root / "source.zip", "w") as z:
        z.writestr(member, "no extraction")
    assert inspect_source(refresh(root, doc))["status"] == "supplied_conflicting"


def test_archive_member_hashes_and_corruption(package):
    root, doc, _ = package
    with zipfile.ZipFile(root / "source.zip", "w") as z:
        z.writestr("nested/source.in", "test source")
    report = inspect_source(refresh(root, doc))
    assert report["status"] == "source_complete_not_executed"
    assert report["archive_members"][0]["sha256"] == sha256(b"test source").hexdigest()
    (root / "source.zip").write_bytes(b"not a ZIP")
    assert inspect_source(refresh(root, doc))["status"] == "supplied_unreadable"


def test_canonical_family_only_and_immutable_report(prepared):
    system, p = prepared
    p.validate_integrity(system)
    assert p.payload["family"] == "opls" and p.payload["selected_profile"] == PROFILE
    assert p.payload["source_intake"]["status"] == "not_supplied"
    assert p.payload["charge_origin"] == p.prior.payload["charge_origin"]
    assert p.payload["assignment_binding"] == {
        "canonical_assignment_identity": p.assignment.identity,
        "graph_identity": p.assignment.payload["graph_identity"],
        "typed_graph_identity": digest(p.assignment.payload["typed_graph"]),
        "charge_identity": digest(p.assignment.payload["graph_charges"]),
        "source_identity": p.assignment.payload["source_identity"],
    }
    assert not p.payload["new_source_selected"] and not p.payload["default_MD_enabled"]
    with pytest.raises(FrozenInstanceError):
        p.json_text = "{}"
    with pytest.raises(OPLSResolutionError):
        replace(p, json_text="{}").validate_integrity(system)


@pytest.mark.parametrize("action", ["minimize", "nvt", "npt"])
def test_MD_rejected_before_output_or_engine(prepared, tmp_path, action):
    system, p = prepared
    target = tmp_path / "must_not_exist"
    with pytest.raises(OPLSResolutionError, match="new_authoritative_source_required"):
        run_md(system, p, target, action=action, research_opt_in=True)
    assert not target.exists()


def test_installed_unverified_package_never_selected(package, prepared):
    _, _, spec = package
    system, prior = prepared
    with use_source_package(spec):
        p = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            assets=prior.assignment.assets,
            sources=prior.prior.sources,
        )
    assert not p.payload["new_source_selected"]
    assert p.payload["source_intake"]["status"] == "source_complete_not_executed"
    with pytest.raises(OPLSResolutionError, match="no_qualified_execution_adapter"):
        p.require_execution(system)


def test_charge_or_graph_mapping_mutation(prepared):
    system, p = prepared
    mutated = copy.deepcopy(system)
    mutated.topology.bonds.pop(next(iter(mutated.topology.bonds)))
    with pytest.raises(OPLSResolutionError):
        p.validate_integrity(mutated)
    assignment = p.assignment
    value = assignment.payload
    value["graph_charges"]["vector"][0][1] = "9"
    changed = replace(assignment, json_text=canonical(value))
    old = p.prior
    altered = replace(
        p,
        prior=replace(
            old,
            prior=replace(
                old.prior, prior=replace(old.prior.prior, assignment=changed)
            ),
        ),
    )
    with pytest.raises(OPLSResolutionError):
        altered.validate_integrity(system)


def test_source_replacement_after_preparation(package, prepared):
    root, _, spec = package
    system, old = prepared
    with use_source_package(spec):
        p = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            assets=old.assignment.assets,
            sources=old.prior.sources,
        )
    (root / "permission.txt").write_text("revoked")
    with pytest.raises(OPLSResolutionError, match="intake_source_or_report_mutation"):
        p.validate_integrity(system)


@pytest.mark.parametrize("seed", [1, 99397])
def test_contract_relocation_fresh_process_and_order(package, tmp_path, seed):
    root, doc, spec = package
    moved = tmp_path / "moved"
    # Physical creation order never participates in resolution.
    moved.mkdir()
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_file():
            target = moved / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(path.read_bytes())
    command = [
        sys.executable,
        "-m",
        "island.intake_opls_source",
        "inspect",
        "--package",
        str(moved),
        "--manifest-sha256",
        spec.manifest_sha256,
    ]
    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONHASHSEED": str(seed)},
    )
    assert process.returncode == 0, process.stderr
    assert json.loads(process.stdout) == inspect_source(spec)
    # Declared row/file order may change only the independent package pin; the
    # domain and qualification decisions must remain the same.
    doc["files"].reverse()
    doc["conventions"].reverse()
    assert (
        inspect_source(write_manifest(moved, doc))["status"]
        == "source_complete_not_executed"
    )


def test_preparation_relocation(prepared, tmp_path):
    system, p = prepared
    expected = save_preparation(system, p, tmp_path / "original")
    shutil.copytree(tmp_path / "original", tmp_path / "moved")
    _, rebuilt = reconstruct(
        tmp_path / "moved",
        manifest_sha256=expected,
        assets=p.assignment.assets,
        sources=p.prior.sources,
    )
    assert rebuilt.identity == p.identity
    (tmp_path / "moved" / "system.json").write_text("{}")
    with pytest.raises(OPLSResolutionError):
        reconstruct(tmp_path / "moved", manifest_sha256=expected)


def test_CLI_schema_and_no_source(tmp_path):
    schema = subprocess.run(
        [sys.executable, "-m", "island.intake_opls_source", "schema"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert schema.returncode == 0 and json.loads(schema.stdout) == SCHEMA
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.intake_opls_source",
            "inspect",
            "--json",
            str(tmp_path / "report.json"),
            "--blocker-reopen-json",
            str(tmp_path / "blockers.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert (
        json.loads((tmp_path / "report.json").read_text())["status"] == "not_supplied"
    )
    assert len(json.loads((tmp_path / "blockers.json").read_text())["rows"]) == 9
