"""Real single-source PDMS rows, explicit missing couplings and replay guards."""

import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields.external_pdms import digest
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pdms_source_completion import (
    PDMSCompletionBlockedError,
    PDMSCompletionError,
    PDMSCompletionRecord,
    assess_pdms_source_completion,
    file_hash,
    load_pdms_source_completion,
)
from island.forcefields.siloxane_literature import load_siloxane_parameters
from scripts.audit_pdms_source_completion import audit
from scripts.complete_pdms_source_parameters import PAPER_PIN, corroborate, run

ROOT = Path(__file__).resolve().parents[1]
DISTRIBUTION = ROOT.parent / "lammps/lammps/tools/msi2lmp/frc_files"
CATALOG = ROOT / "docs/evidence/sun_siloxane_parameters/parameters.json"
REFERENCES = ROOT / "island_ffref"


@pytest.fixture(scope="module")
def distribution():
    if not (DISTRIBUTION / "pcff.frc").is_file():
        pytest.skip("authorized local FRC distribution unavailable")
    return load_pcff_sip_distribution(DISTRIBUTION)


@pytest.fixture(scope="module")
def hmds(distribution):
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    return system, assess_pdms_source_completion(system, distribution)


def test_all_printed_angle_fields_and_missing_couplings_preserved(hmds):
    system, record = hmds
    p = record.payload
    assert p["coverage"]["quartic_angle"] == {"assigned": 49}
    assert p["coverage"]["bond-bond"] == {"assigned": 49}
    assert p["coverage"]["bond-angle"] == {"assigned": 49}
    assert len(p["missing_dependencies"]) == 378
    assert p["native_charge_record"]["complete"]
    assert len(p["exact_native_charge_vector"]) == system.number_of_sites
    assert all(e["type"] != "si" for e in p["typing"]["entries"].values())
    row = p["source_rows"]["quartic_angle:cff91:2197"]
    assert row["types"] == ["h", "c", "sio"]
    assert row["original_values"] == ["111.5360", "30.2481", "-15.5255", "0.0000"]
    assert "quartic_angle:cff91:2196" not in p["source_rows"]
    for family in (
        "bond-bond_1_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
    ):
        assert p["coverage"][family] == {"missing": 60}
    assert p["coverage"]["angle-angle"] == {"assigned": 18, "missing": 78}
    assert all(
        q["path"].startswith(("direct", "ordinary_"))
        for q in p["queries"].values()
        if q["status"] == "assigned"
    )
    with pytest.raises(PDMSCompletionBlockedError) as exc:
        record.require_executable()
    assert exc.value.result is record
    assert exc.value.diagnostics == p["missing_dependencies"]
    with pytest.raises(FrozenInstanceError):
        record.json_text = "changed"
    p["source_rows"].clear()
    assert record.payload["source_rows"]


@pytest.mark.parametrize(
    "mutation",
    [
        "zero",
        "units",
        "row",
        "charges",
        "inventory",
        "order",
        "graph",
        "component",
        "source",
        "authority",
        "status",
        "profile",
        "evidence",
        "resolution",
    ],
)
def test_resigned_mutations_rejected(hmds, distribution, mutation):
    system, record = hmds
    p = record.payload
    first = next(iter(p["source_rows"].values()))
    if mutation == "zero":
        first["original_values"][0] = "0.0000"
    elif mutation == "units":
        first["units"][0] = "nm"
    elif mutation == "row":
        p["source_rows"].pop(first["id"])
    elif mutation == "charges":
        p["exact_native_charge_vector"][1] += 0.001
    elif mutation == "inventory":
        p["missing_dependencies"].clear()
    elif mutation == "order":
        p["requirements"]["atom_order"].reverse()
    elif mutation == "graph":
        p["requirements"]["graph"]["bonds"].pop()
    elif mutation == "component":
        p["requirements"]["graph"]["molecule_membership"]["1"] = "forged"
    elif mutation == "source":
        p["source_hashes"]["pcff.frc"] = "0" * 64
    elif mutation == "authority":
        p["execution_authorized"] = True
    elif mutation == "status":
        p["status"] = "prepared"
    elif mutation == "profile":
        p["profile_identity"] = "native_pcff"
    elif mutation == "evidence":
        first["raw"] += " forged"
    else:
        next(q for q in p["queries"].values() if q["status"] == "missing")["status"] = (
            "assigned"
        )
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(PDMSCompletionError) as exc:
        altered = PDMSCompletionRecord(pack(p))
        altered.validate_integrity(system, distribution)
    assert exc.value.status == "source_incompatible"


@pytest.mark.parametrize(
    "smiles",
    ["[SiH2]O[SiH3]", "[O-][Si](C)(C)C", "CO[Si](C)(C)C", "COP(=O)(O)OC", "C[SiH3]"],
)
def test_invalid_or_unsupported_graphs_are_classified(distribution, smiles):
    with pytest.raises(PDMSCompletionError) as exc:
        assess_pdms_source_completion(from_smiles(smiles), distribution)
    assert exc.value.status == "typing_blocked"


def test_stale_graph_and_molecule_rejected(hmds, distribution):
    system, record = hmds
    changed = system.copy()
    first = next(iter(changed.topology.sites))
    changed.topology.sites[first].metadata["chain_id"] = "changed"
    with pytest.raises(PDMSCompletionError):
        record.validate_integrity(changed, distribution)
    changed = system.copy()
    atom = changed.topology.sites[first]
    changed.topology.sites[first] = replace(atom, formal_charge=1)
    with pytest.raises(PDMSCompletionError):
        record.validate_integrity(changed, distribution)


def test_loader_requires_external_pin_and_present_record(hmds, distribution, tmp_path):
    system, record = hmds
    path = tmp_path / "record.json"
    with pytest.raises(PDMSCompletionError) as exc:
        load_pdms_source_completion(
            path, "0" * 64, system=system, distribution=distribution
        )
    assert exc.value.status == "asset_unavailable"
    path.write_text(record.json_text)
    with pytest.raises(PDMSCompletionError):
        load_pdms_source_completion(
            path, "0" * 64, system=system, distribution=distribution
        )
    assert (
        load_pdms_source_completion(
            path, file_hash(path), system=system, distribution=distribution
        )
        == record
    )


@pytest.mark.parametrize("asset", ["frc", "templates", "rlb"])
def test_tampered_distribution_cannot_supply_rows(hmds, distribution, asset):
    if asset == "frc":
        changed = replace(
            distribution,
            source=replace(distribution.source, raw=distribution.source.raw + b"\n"),
        )
    else:
        changed = replace(distribution, **{asset: b"changed source"})
    with pytest.raises(PDMSCompletionError) as exc:
        assess_pdms_source_completion(hmds[0], changed)
    assert exc.value.status == "source_incompatible"


def test_eight_paper_blanks_corroborated_without_rewriting(distribution):
    if not (REFERENCES / "pcff/ma00107a006.pdf").exists():
        pytest.skip("local PDFs unavailable")
    supplement = load_siloxane_parameters(CATALOG, PAPER_PIN, references=REFERENCES)
    result = corroborate(supplement, distribution)
    assert len(result["rows"]) == 8
    assert all(r["paper_original_values_unchanged"][-1] is None for r in result["rows"])
    assert all(
        r["complete_frc_row"]["original_values"][-1] == "0.0000" for r in result["rows"]
    )


def test_relocated_source_and_fresh_process_reconstruction(distribution, tmp_path):
    if not (REFERENCES / "pcff/ma00107a006.pdf").exists():
        pytest.skip("local PDFs unavailable")
    relocated = tmp_path / "sources"
    relocated.mkdir()
    for name in distribution.hashes:
        shutil.copyfile(DISTRIBUTION / name, relocated / name)
    output = tmp_path / "study"
    report = run(relocated, CATALOG, REFERENCES, output)
    assert [f["missing_interaction_requests"] for f in report["fixtures"]] == [
        54,
        378,
        378,
        528,
        828,
    ]
    raw = audit(output, relocated)
    assert raw["numerical_rows"] == 47
    cmd = [
        sys.executable,
        str(ROOT / "scripts/complete_pdms_source_parameters.py"),
        "--distribution",
        str(relocated),
        "--catalog",
        str(CATALOG),
        "--references",
        str(REFERENCES),
        "--output",
        str(output),
        "--reconstruct",
        "--require-executable",
    ]
    proc = subprocess.run(
        cmd, cwd=ROOT, capture_output=True, text=True, timeout=180, check=False
    )
    assert proc.returncode == 1, proc.stderr
    assert json.loads(proc.stdout) == report
    assert report["execution_profiles"] == report["backend_executions"] == []
    assert not list(output.glob("*.chk"))
    altered = output / "missing-interactions.json"
    altered.write_text("[]")
    proc = subprocess.run(
        cmd, cwd=ROOT, capture_output=True, text=True, timeout=180, check=False
    )
    assert proc.returncode == 1
    assert json.loads(proc.stdout)["status"] == "source_incompatible"
