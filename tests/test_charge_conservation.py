"""Synthetic software evidence, never live QM acceptance."""

import os
import subprocess
import sys
from copy import deepcopy
from math import fsum

import numpy as np
import pytest
from test_oligomer_charge_references import molecule

from island.charge_references import (
    ChargeProjection,
    RawChargeObservation,
    audit_charge_projections,
    load_raw_observation,
    project_charge_observation,
    project_molecular_charges,
    save_record,
)
from island.charge_references.observations import OBSERVATION_SCHEMA, _observation_data
from island.charge_references.projection import POLICY
from island.charge_references.records import pack
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.lineage import generated_atom_name
from island.workflows import storage
from island.workflows.bundle import system_data
from scripts import audit_charge_conservation as cli


def synthetic_observation(dp=3, seed=2026, oxygen=False, shift=0.0):
    s = molecule(dp, oxygen=oxygen, seed=seed)
    s.coordinates.translate((shift, 0, 0))
    s.metadata["test_evidence"] = "SYNTHETIC software fixture; no QM executed"
    ids = sorted(s.topology.sites)
    names = {
        i: generated_atom_name(s.topology.sites[i].element, k)
        for k, i in enumerate(ids)
    }
    index = {s: i + 1 for i, s in enumerate(ids)}
    charges = {i: (0.01 if k % 2 == 0 else -0.009) for k, i in enumerate(ids)}
    atoms = []
    for i in ids:
        xyz = " ".join(map(str, s.coordinates.get(i)))
        atoms.append(
            f"{index[i]} {names[i]} {xyz} {s.topology.sites[i].element} 1 MOL {charges[i]}"
        )
    bonds = [
        f"{k + 1} {index[a]} {index[b]} 1" for k, (a, b) in enumerate(s.topology.bonds)
    ]
    text = (
        "@<TRIPOS>MOLECULE\nTEST\n"
        + f"{len(ids)} {len(bonds)} 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n"
        + "\n".join(atoms)
        + "\n@<TRIPOS>BOND\n"
        + "\n".join(bonds)
        + "\n"
    )
    sqm_input = "&qmmm qm_theory='AM1', qmcharge=0, /\n" + "\n".join(
        f"{s.topology.sites[i].atomic_number} {names[i]} "
        + " ".join(f"{v:.3f}" for v in s.coordinates.get(i))
        for i in ids
    )
    sqm_output = (
        "Final Structure\n"
        + "\n".join(
            f"QMMM: {index[i]} {index[i]} {s.topology.sites[i].element} 0 0 0"
            for i in ids
        )
        + "\nCalculation Completed\n"
    )
    top = (
        "%FLAG ATOM_NAME\n%FORMAT(20a4)\n"
        + "".join(names[i] for i in ids)
        + "\n%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n"
        + "".join(f"{s.topology.sites[i].atomic_number:8d}" for i in ids)
        + "\n%FLAG CHARGE\n%FORMAT(5E16.8)\n"
        + "".join(f"{charges[i] * 18.2223:16.8E}" for i in ids)
        + "\n"
    )
    values = {
        "ANTECHAMBER_AC.AC": "\n".join(
            f"ATOM {index[i]} {names[i]}MOL 1 "
            + " ".join(f"{v:.3f}" for v in s.coordinates.get(i))
            + f" {charges[i]} {s.topology.sites[i].element}"
            for i in ids
        ),
        "ANTECHAMBER_AM1BCC.AC": "\n".join(
            f"ATOM {index[i]} {names[i]}MOL 1 0 0 0 {charges[i]} {s.topology.sites[i].element}"
            for i in ids
        ),
        "input.mol2": text,
        "typed.mol2": text,
        "result.prmtop": top,
        "lineage.json": storage.json_bytes(
            {
                "schema": "island_ambertools_lineage_v1",
                "input_name_to_site_id": {n: i for i, n in names.items()},
                "input_index_to_site_id": {str(k): i for i, k in index.items()},
            }
        ).decode(),
        "sqm.in": sqm_input,
        "sqm.out": sqm_output,
        "leap.in": 'source "/synthetic/leaprc.gaff2"\n',
        "leap.log": "Exiting LEaP: Errors = 0; Warnings = 0; Notes = 0.",
        "antechamber.stdout.log": "Running: /synthetic/atomtype -p gaff2\nRunning: /synthetic/sqm -O -i sqm.in -o sqm.out\nRunning: /synthetic/am1bcc -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC",
        "antechamber.stderr.log": "",
        "parmchk2.stdout.log": "",
        "parmchk2.stderr.log": "",
        "tleap.stdout.log": "",
        "tleap.stderr.log": "",
    }
    files = {"island-ambertools-test/" + k: v for k, v in values.items()}
    files["input-system.json"] = storage.json_bytes(
        storage.encode(system_data(s))
    ).decode()
    name = f"{'peo' if oxygen else 'pe'}-dp{dp}-seed{seed}"
    row = {
        "case": name,
        "chemistry": "peo" if oxygen else "pe",
        "dp": dp,
        "seed": seed,
        "psmiles": s.metadata["polymer"]["source_psmiles"],
        "sites": len(ids),
        "status": "failed",
        "failure": "SYNTHETIC import: Charge assignment is incomplete; not real QM",
        "source_checksums": {k: storage.checksum(v.encode()) for k, v in files.items()},
    }
    evidence = {
        "case": name,
        "manifest_json": storage.json_bytes(
            {"cases": [row], "label": "SYNTHETIC SOFTWARE TEST"}
        ).decode(),
        "files": files,
        "historical_reference": None,
    }
    return RawChargeObservation(
        pack(
            {
                "schema": OBSERVATION_SCHEMA,
                "evidence": evidence,
                "data": _observation_data(evidence),
            }
        )
    )


def test_analytical_projection_ownership_order_and_idempotence():
    s = molecule()
    ids = sorted(s.topology.sites)
    q = {i: float(k) / 100 for k, i in enumerate(ids)}
    before = deepcopy(q)
    system_before = system_data(s)
    p = project_molecular_charges(s, q, policy=POLICY)
    mean = sum(range(len(ids))) / 100 / len(ids)
    assert p["uniform_offset_e"] == pytest.approx(-mean, abs=1e-16)
    assert p["projected_charges"] == pytest.approx(
        {i: v - mean for i, v in q.items()}, abs=1e-16
    )
    assert abs(fsum(p["projected_charges"].values())) < 1e-12
    assert q == before and system_data(s) == system_before
    again = project_molecular_charges(s, p["projected_charges"], policy=POLICY)
    assert again["projected_charges"] == pytest.approx(
        p["projected_charges"], abs=1e-15
    )
    reverse = project_molecular_charges(
        molecule(reverse=True), dict(reversed(list(q.items()))), policy=POLICY
    )
    assert p == reverse
    constant = project_molecular_charges(s, {i: 0.04 for i in ids}, policy=POLICY)
    assert len(set(constant["projected_charges"].values())) == 1
    neutral = {i: 0.0 for i in ids}
    neutral[ids[0]] = 1.0
    neutral[ids[1]] = -1.0
    assert (
        project_molecular_charges(s, neutral, policy=POLICY)["projected_charges"]
        == neutral
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda q: q.pop(next(iter(q))),
        lambda q: q.update({-999: 0.1}),
        lambda q: q.update({next(iter(q)): True}),
        lambda q: q.update({next(iter(q)): np.bool_(False)}),
        lambda q: q.update({next(iter(q)): float("nan")}),
        lambda q: q.update({next(iter(q)): float("inf")}),
        lambda q: q.update({next(iter(q)): "0"}),
        lambda q: q.update({False: 0.0}),
    ],
)
def test_bad_charges(mutation):
    s = molecule()
    q = {i: 0.1 for i in s.topology.sites}
    mutation(q)
    with pytest.raises(ChargeReferenceError):
        project_molecular_charges(s, q, policy=POLICY)


@pytest.mark.parametrize("which", ["policy", "scope", "limit"])
def test_unsupported(which):
    s = molecule(dp=51 if which == "limit" else 3)
    if which == "scope":
        s.metadata["polymer"]["source_psmiles"] = "[*]CC[*]"
    with pytest.raises(ChargeReferenceError):
        project_molecular_charges(
            s,
            {i: 0.0 for i in s.topology.sites},
            policy="automatic" if which == "policy" else POLICY,
        )


def test_observation_projection_failed_status_and_repeat_groups(tmp_path):
    o = synthetic_observation()
    initial = o.json_text
    p = project_charge_observation(o, policy=POLICY)
    d = p.payload["projection"]
    assert p.payload["observation"]["data"]["historical_status"] == "failed"
    assert p.payload["observation"]["data"]["historical_reference_identity"] is None
    for row in d["repeat_changes"]:
        assert row["change_e"] == pytest.approx(
            row["atom_count"] * d["uniform_offset_e"], abs=1e-15
        )
    a = audit_charge_projections([p]).payload["audit"]
    raw = a["raw"]["summaries"][0]
    projected = a["projected"]["summaries"][0]
    for r, t in zip(raw["repeats"], projected["repeats"], strict=True):
        for g, h in zip(r["groups"], t["groups"], strict=True):
            assert h["hydrogen_count"] == g["hydrogen_count"]
            assert h["hydrogen_sum"] - g["hydrogen_sum"] == pytest.approx(
                g["hydrogen_count"] * d["uniform_offset_e"], abs=1e-15
            )
    save_record(o, tmp_path / "o.json")
    loaded = load_raw_observation(tmp_path / "o.json")
    assert loaded.identity == o.identity
    with pytest.raises(ChargeReferenceError):
        save_record(p, tmp_path / "o.json")
    assert o.json_text == initial and (tmp_path / "o.json").read_text() == initial
    x = p.payload
    x["projection"]["raw_charges"].clear()
    assert p.payload["projection"]["raw_charges"]


@pytest.mark.parametrize(
    "key,value",
    [
        ("uniform_offset_e", 1.0),
        ("policy", "other"),
        ("projected_residual_e", True),
        ("numerical_tolerance_e", 0.002),
        ("production_validated", True),
    ],
)
def test_rechecksummed_projection_tampering(key, value, tmp_path):
    p = project_charge_observation(synthetic_observation(), policy=POLICY).payload
    p["projection"][key] = value
    with pytest.raises(ChargeReferenceError):
        save_record(ChargeProjection(pack(p)), tmp_path / "invalid.json")
    assert not (tmp_path / "invalid.json").exists()


@pytest.mark.parametrize(
    "field", ["raw_charges", "source_manifest_sha256", "historical_status"]
)
def test_observation_recomputed_checksum_is_not_sufficient(field):
    p = synthetic_observation().payload
    p["data"][field] = "changed"
    with pytest.raises(ChargeReferenceError):
        RawChargeObservation(pack(p)).validate_integrity()


def test_neutral_oligomer_is_not_neutral_naive_transfer():
    from island.charge_references import repeat_correspondence
    from island.charge_references.audit import audit_summaries, summarize_charge_mapping

    s = molecule()
    c = repeat_correspondence(s)
    q = {i: 0.0 for i in s.topology.sites}
    for r, value in zip(c["repeats"], [-0.1, 0.2, -0.1], strict=True):
        q[r["groups"][0]["site_id"]] = value
    p = project_molecular_charges(s, q, policy=POLICY)
    assert p["projected_residual_e"] == 0
    summary = summarize_charge_mapping(c, p["projected_charges"], "synthetic", {}, 0.0)
    result = audit_summaries([summary])["diagnostic_extrapolations"][0]
    assert result["predicted_total_e"][20] == pytest.approx(3.4)


def test_offline_without_optional_dependencies(tmp_path):
    o = synthetic_observation()
    save_record(o, tmp_path / "observation.json")
    code = """
import sys
class Block:
 def find_spec(self,fullname,*args):
  if fullname.split('.')[0] in {'rdkit','parmed','openmm','scipy'}: raise AssertionError(fullname)
sys.meta_path.insert(0,Block())
from island.charge_references import load_raw_observation,project_charge_observation,audit_charge_projections
p=project_charge_observation(load_raw_observation(sys.argv[1]),policy='uniform_molecular_l2_v1')
audit_charge_projections([p]).validate_integrity()
"""
    subprocess.run(
        [sys.executable, "-c", code, str(tmp_path / "observation.json")],
        check=True,
        env=os.environ.copy(),
    )


def test_exact_manifest_paths_and_malformed_lineage(tmp_path):
    from island.charge_references import observe_retained_case

    o = synthetic_observation()
    e = o.payload["evidence"]
    root = tmp_path / "source"
    directory = root / e["case"]
    directory.mkdir(parents=True)
    for relative, text in e["files"].items():
        path = directory / relative
        path.parent.mkdir(exist_ok=True)
        path.write_text(text)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(e["manifest_json"])
    assert (
        observe_retained_case(root, manifest, e["case"]).payload["data"]
        == o.payload["data"]
    )
    (directory / "island-ambertools-unlisted").mkdir()
    with pytest.raises(ChargeReferenceError, match="Ambiguous"):
        observe_retained_case(root, manifest, e["case"])
    (directory / "island-ambertools-unlisted").rmdir()
    (directory / "island-ambertools-test/typed.mol2").write_text("substitution")
    with pytest.raises(ChargeReferenceError, match="checksum"):
        observe_retained_case(root, manifest, e["case"])


def test_cli_gate_and_existing_output_unchanged(tmp_path, monkeypatch):
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}")
    output = tmp_path / "output"
    args = [
        "--source",
        str(tmp_path / "source"),
        "--manifest",
        str(manifest),
        "--output",
        str(output),
        "--policy",
        POLICY,
    ]
    monkeypatch.setattr(cli, "execute", lambda *_: False)
    assert cli.main(args) == 1
    assert {p.name for p in output.iterdir()} == {"declaration.json"}
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    assert cli.main(args) == 1
    assert {p.name: p.read_bytes() for p in output.iterdir()} == before
    monkeypatch.setattr(cli, "execute", lambda *_: True)
    assert cli.main(args + ["--execute-declared"]) == 0


@pytest.mark.parametrize(
    "change",
    ["duplicate_atom", "bond_order", "sqm_failure", "element", "missing_artifact"],
)
def test_rechecksummed_scientific_evidence_rejected(change):
    p = synthetic_observation().payload
    e = p["evidence"]
    key = "island-ambertools-test/typed.mol2"
    if change == "duplicate_atom":
        lines = e["files"][key].splitlines()
        i = lines.index("@<TRIPOS>ATOM") + 1
        lines.insert(i, lines[i])
        e["files"][key] = "\n".join(lines) + "\n"
    elif change == "bond_order":
        e["files"][key] = e["files"][key].replace(
            "@<TRIPOS>BOND\n", "@<TRIPOS>BOND\n999 1 2 2\n"
        )
    elif change == "element":
        lines = e["files"][key].splitlines()
        i = lines.index("@<TRIPOS>ATOM") + 1
        r = lines[i].split()
        r[5] = "O"
        lines[i] = " ".join(r)
        e["files"][key] = "\n".join(lines) + "\n"
    elif change == "sqm_failure":
        key = "island-ambertools-test/sqm.out"
        e["files"][key] += "\nSCF convergence failure\n"
    else:
        e["files"].pop(key)
    from island.dynamics._checkpoint_data import strict_load

    manifest = strict_load(e["manifest_json"])
    if key in e["files"]:
        manifest["cases"][0]["source_checksums"][key] = storage.checksum(
            e["files"][key].encode()
        )
    e["manifest_json"] = storage.json_bytes(manifest).decode()
    with pytest.raises(ChargeReferenceError):
        RawChargeObservation(pack(p)).validate_integrity()


@pytest.mark.parametrize("key", ["raw_charges", "projected_charges"])
def test_projection_charge_values_recomputed(key):
    p = project_charge_observation(synthetic_observation(), policy=POLICY).payload
    p["projection"][key][next(iter(p["projection"][key]))] += 0.001
    with pytest.raises(ChargeReferenceError):
        ChargeProjection(pack(p)).validate_integrity()


def test_partial_matrix_execution_gate(tmp_path, monkeypatch):
    """Injected software outcomes remain failures, not executed reference evidence."""
    import json

    root = tmp_path / "source"
    root.mkdir()
    rows = []
    for c in cli.MATRIX:
        o = synthetic_observation(c["dp"], c["seed"], c["chemistry"] == "peo")
        e = o.payload["evidence"]
        row = json.loads(e["manifest_json"])["cases"][0]
        rows.append(row)
        directory = root / row["case"]
        directory.mkdir()
        for rel, text in e["files"].items():
            p = directory / rel
            p.parent.mkdir(exist_ok=True)
            p.write_text(text)
    manifest = tmp_path / "manifest.json"
    manifest.write_bytes(storage.json_bytes({"cases": rows}))
    evidence = tmp_path / "evidence.json"
    evidence.write_bytes(
        storage.json_bytes(
            {
                "original_manifest_sha256": storage.checksum(manifest.read_bytes()),
                "cases": rows,
            }
        )
    )
    output = tmp_path / "out"
    output.mkdir()
    original = cli.observe_retained_case

    def fail_one(source, path, name):
        if name == rows[0]["case"]:
            raise ChargeReferenceError("Injected failure")
        return original(source, path, name)

    monkeypatch.setattr(cli, "observe_retained_case", fail_one)
    assert cli.execute(root, manifest, evidence, output) is False
    outcomes = storage.read_json(output / "outcomes.json")
    assert outcomes["complete_diagnostic_matrix"] is False
    assert len([r for r in outcomes["cases"] if r["status"] == "passed"]) == 7
    assert outcomes["source_hashes_unchanged"] is True


def test_retained_pre_qm_ac_rounding_is_the_sqm_input():
    synthetic_observation(shift=0.0003).validate_integrity()
