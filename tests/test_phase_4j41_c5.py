"""External schema controls use real raw rows, never synthetic PDMS parameters."""

import json
import os
import subprocess
import sys
import tarfile
from dataclasses import FrozenInstanceError
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields.external_pdms import (
    AUTHORITY,
    FAMILIES,
    LIMITS,
    ROWS_SCHEMA,
    ExternalParameterBundle,
    bundle_template,
    load_external_pdms_bundle,
    pdms_final_graph,
    requirements,
    source_row,
)
from scripts.pcff_j5_reference import read_source
from scripts.phase_4j41_c5 import CASES, MUTATIONS, fixture, mutated_record


@pytest.fixture(scope="module")
def system():
    return fixture(CASES[1])


@pytest.fixture(scope="module")
def raw_pcff():
    path = (
        Path(
            os.environ.get(
                "ISLAND_J41_DISTRIBUTION", "../lammps/lammps/tools/msi2lmp/frc_files"
            )
        )
        / "pcff.frc"
    )
    if path.is_file():
        return path.read_bytes()
    with tarfile.open("docs/evidence/phase_4j41_c4/execution.tar.xz") as tar:
        return tar.extractfile("study/pcff/pcff.frc").read()


def stages(report):
    return {r["stage"] for r in report.payload["diagnostics"]}


@pytest.mark.parametrize("case", CASES, ids=[c[0] for c in CASES])
def test_pdms_graphs_are_explicit_h_validated_and_source_bound(case):
    s = fixture(case)
    r = requirements(s)
    assert len(r["atom_order"]) == s.number_of_sites
    assert r["graph"]["molecule_membership"]
    assert r["graph"]["components"] == [r["atom_order"]]
    assert {"quartic_bond", "quartic_angle", "angle-angle", "nonbond(9-6)"} <= {
        v["family"] for v in r["interactions"].values()
    }
    assert all(
        d in r["interactions"]
        for v in r["interactions"].values()
        for d in v["dependencies"]
    )
    if case[2]:
        assert s.metadata["polymer"]["source_psmiles"] == "[*:1]O[Si](C)(C)[*:2]"
        assert s.metadata["polymer"]["end_groups"]["head"] == "[*]C"
        assert s.metadata["polymer"]["end_groups"]["tail"] == "[*]O"


def test_missing_bundle_has_precise_asset_result(system, tmp_path):
    bundle, outcome = load_external_pdms_bundle(
        tmp_path / "absent.json", system, expected_sha256=None
    )
    assert bundle is None and outcome.payload["status"] == "asset_unavailable"
    assert outcome.payload["limitations"] == LIMITS and not outcome.complete


def test_incomplete_template_preserves_every_missing_interaction(system, tmp_path):
    b = ExternalParameterBundle(pack(bundle_template(system)))
    report = b.validate(system, assets=tmp_path, expected_sha256=b.identity)
    missing = [
        e["context"]["key"]
        for e in report.payload["diagnostics"]
        if e["stage"] == "missing_interaction"
    ]
    assert set(missing) == set(requirements(system)["interactions"])
    assert not report.complete and report.payload["execution_profile"] is None


@pytest.mark.parametrize("mutation", MUTATIONS)
def test_resigned_mutations_hit_their_own_gate(system, tmp_path, mutation):
    raw, expected = mutated_record(bundle_template(system), mutation)
    b = ExternalParameterBundle(raw)
    report = b.validate(system, assets=tmp_path, expected_sha256=b.identity)
    assert expected in stages(report)
    assert not report.complete and not report.payload["execution_authorized"]


def test_bundle_pin_is_external_and_cannot_be_self_asserted(system, tmp_path):
    b = ExternalParameterBundle(pack(bundle_template(system)))
    assert "bundle_pin_mismatch" in stages(
        b.validate(system, assets=tmp_path, expected_sha256=None)
    )
    assert "bundle_pin_mismatch" in stages(
        b.validate(system, assets=tmp_path, expected_sha256="0" * 64)
    )


def test_record_is_owned_and_frozen(system):
    b = ExternalParameterBundle(pack(bundle_template(system)))
    p = b.payload
    p["authority"]["execution_authorized"] = True
    assert b.payload["authority"] == AUTHORITY
    assert bundle_template(system)["authority"] == AUTHORITY
    with pytest.raises(FrozenInstanceError):
        b.json_text = "changed"


@pytest.mark.parametrize(
    "bad", ["[SiH2]O", "C[SiH]O", "C[Si](C)(C)[O-]", "COP(=O)(OC)OC"]
)
def test_wrong_valence_or_chemical_state_rejected(bad):
    with pytest.raises(ValueError):
        pdms_final_graph(from_smiles(bad, random_seed=2026))


def test_element_atomic_number_and_component_mutations(system):
    from island.workflows.bundle import system_data, system_from

    s = system_from(system_data(system))
    first = next(i for i, a in s.topology.sites.items() if a.element == "Si")
    s.topology.sites[first].atomic_number = 15
    with pytest.raises(ValueError, match="wrong_element"):
        pdms_final_graph(s)


def real_row(raw_pcff, family, *, zero=False):
    rows = read_source(raw_pcff)[family, "cff91"]
    row = next(
        r
        for r in rows
        if not any("*" in t for t in r["types"])
        and len(r["values"]) == len(FAMILIES[family][2])
        and (not zero or all(v == 0 for v in r["values"]))
    )
    values = list(map(str, row["values"]))
    # Contract control reconstructs literal FRC values. It is not an external
    # bundle, does not map PDMS, and is never passed to prepare_forcefield.
    count = len(row["types"])
    r = {
        "family": family,
        "types": row["types"],
        "values": values,
        "units": FAMILIES[family][2][: len(values)],
        "provenance": "literal pinned-PCFF row contract control only",
        "citations": [
            {
                "file": "pcff.frc",
                "line": row["line"],
                "raw": row["raw"],
                "columns": list(range(2 + count, 2 + count + len(values))),
            }
        ],
    }
    doc = {
        "schema": ROWS_SCHEMA,
        "source_family": "external:contract_control_only",
        "rows": {"r": r},
        "provenance": "schema unit control",
    }
    return r, doc


@pytest.mark.parametrize(
    "family", ["quartic_bond", "bond-angle", "torsion_3", "nonbond(9-6)"]
)
def test_independent_literal_source_reconstruction(raw_pcff, family):
    row, doc = real_row(raw_pcff, family)
    values = source_row(row, doc, {"pcff.frc": raw_pcff})
    assert list(map(str, values)) == row["values"]


def test_explicit_source_zero_allowed_but_inserted_zero_rejected(raw_pcff):
    row, doc = real_row(raw_pcff, "bond-angle", zero=True)
    assert all(v == 0 for v in source_row(row, doc, {"pcff.frc": raw_pcff}))
    row, doc = real_row(raw_pcff, "quartic_bond")
    row["values"][0] = "0"
    with pytest.raises(ValueError, match="source_coefficient_changed_or_inserted_zero"):
        source_row(row, doc, {"pcff.frc": raw_pcff})


def test_present_asset_tampering_is_not_hidden_by_other_missing_parameters(
    system, tmp_path
):
    p = bundle_template(system)
    p["source_family"] = "external:contract_control_only"
    (tmp_path / "evidence.txt").write_text("original\n")
    p["files"] = {
        "evidence.txt": {
            "sha256": sha256(b"original\n").hexdigest(),
            "family": p["source_family"],
            "role": "evidence",
            "provenance": "unit control",
        }
    }
    (tmp_path / "evidence.txt").write_text("tampered\n")
    b = ExternalParameterBundle(pack(p))
    assert "source_hash_mismatch" in stages(
        b.validate(system, assets=tmp_path, expected_sha256=b.identity)
    )


def test_wrong_primary_units_or_source_rows_rejected(system, tmp_path):
    p = bundle_template(system)
    family = "external:contract_control_only"
    doc = {
        "schema": ROWS_SCHEMA,
        "source_family": "COMPASS",
        "rows": {"unexpected": {}},
        "provenance": "negative control",
    }
    raw = json.dumps(doc).encode()
    (tmp_path / "rows.json").write_bytes(raw)
    p.update(
        source_family=family,
        primary_source="rows.json",
        files={
            "rows.json": {
                "sha256": sha256(raw).hexdigest(),
                "family": family,
                "role": "parameters",
                "provenance": "negative control",
            }
        },
    )
    b = ExternalParameterBundle(pack(p))
    report = b.validate(system, assets=tmp_path, expected_sha256=b.identity)
    assert {"primary_source_identity_mismatch", "source_rows_mismatch"} <= stages(
        report
    )


def test_path_escape_rejected(system, tmp_path):
    p = bundle_template(system)
    p["files"]["../escape"] = {
        "sha256": "0" * 64,
        "family": None,
        "role": "evidence",
        "provenance": "negative",
    }
    b = ExternalParameterBundle(pack(p))
    assert "malformed_bundle" in stages(
        b.validate(system, assets=tmp_path, expected_sha256=b.identity)
    )


def test_relocation_and_fresh_process_preserve_blocked_bundle(system, tmp_path):
    from island.workflows.bundle import system_data
    from scripts.phase_4j41_c2 import write

    a = tmp_path / "relocated"
    a.mkdir()
    b = ExternalParameterBundle(pack(bundle_template(system)))
    (a / "bundle.json").write_text(b.json_text)
    write(a / "system.json", system_data(system))
    command = """import sys
from pathlib import Path
from island.forcefields.external_pdms import load_external_pdms_bundle
from island.workflows.bundle import system_from
from scripts.phase_4j41_c2 import read
p=Path(sys.argv[1])
b,r=load_external_pdms_bundle(p/'bundle.json',system_from(read(p/'system.json')),expected_sha256=sys.argv[2])
print(r.json_text)
"""
    r = subprocess.run(
        [sys.executable, "-c", command, str(a), b.identity],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, r.stderr
    expected = b.validate(system, assets=a, expected_sha256=b.identity)
    assert r.stdout.strip() == expected.json_text.strip()


def test_public_cli_has_no_user_car_mdf_requirement():
    r = subprocess.run(
        [sys.executable, "scripts/phase_4j41_c5.py", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0
    assert "--bundle" in r.stdout and "--bundle-sha256" in r.stdout
    assert "--car" not in r.stdout and "--mdf" not in r.stdout


def test_sico_couplings_cannot_be_dropped(tmp_path):
    s = from_smiles("COC[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    template = bundle_template(s)
    atoms = {a["id"]: a["element"] for a in template["graph"]["sites"]}
    keys = {
        k
        for k, r in template["coverage"].items()
        if r["family"] in ("bond-angle", "bond-bond")
        and [atoms[i] for i in r["sites"]] in (["Si", "C", "O"], ["O", "C", "Si"])
    }
    assert keys
    b = ExternalParameterBundle(pack(template))
    report = b.validate(s, assets=tmp_path, expected_sha256=b.identity)
    missing = {
        e["context"]["key"]
        for e in report.payload["diagnostics"]
        if e["stage"] == "missing_interaction"
    }
    assert keys <= missing


def test_source_row_wrong_units_rejected(raw_pcff):
    row, doc = real_row(raw_pcff, "quartic_bond")
    row["units"][0] = "nm"
    with pytest.raises(ValueError, match="source_row_units_or_shape_invalid"):
        source_row(row, doc, {"pcff.frc": raw_pcff})


@pytest.mark.parametrize("target", ["receipt", "manifest", "archive"])
def test_tampered_publication_evidence_rejected(tmp_path, target):
    from scripts.phase_4j41_c5 import BASE, BRANCH
    from scripts.validate_phase_4j41_c5_receipt import validate
    from scripts.validate_phase_4j41_receipt import sha, signed

    m = signed(
        {"archives": {"execution.tar.xz": "0" * 64, "validation.tar.xz": "0" * 64}}
    )
    if target == "manifest":
        m["forged"] = True
    (tmp_path / "manifest.json").write_text(json.dumps(m))
    r = signed(
        {
            "schema": "island.j41-c5-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": sha(tmp_path / "manifest.json"),
        }
    )
    if target == "receipt":
        r["forged"] = True
    (tmp_path / "receipt.json").write_text(json.dumps(r))
    (tmp_path / "execution.tar.xz").write_bytes(b"tampered")
    with pytest.raises(
        ValueError,
        match={
            "receipt": "Receipt checksum",
            "manifest": "Manifest checksum",
            "archive": "Archive changed",
        }[target],
    ):
        validate(tmp_path)
