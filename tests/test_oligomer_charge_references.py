"""Synthetic data-only charge evidence; these tests never claim a QM calculation."""

import subprocess
import sys
from copy import deepcopy
from dataclasses import replace

import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references import (
    ChargeReference,
    audit_charge_references,
    load_charge_reference,
    repeat_correspondence,
    save_record,
)
from island.charge_references.records import pack
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.models import digest
from island.forcefields.charges import ProvidedChargeEngine
from island.forcefields.typing.signatures import graph_signature
from island.workflows import storage
from island.workflows.bundle import record, system_data


def molecule(dp=3, oxygen=False, seed=2026, reverse=False):
    elements = ["C", "C", "O"] if oxygen else ["C", "C"]
    definition = "[*:1]CCO[*:2]" if oxygen else "[*:1]CC[*:2]"
    graph = Topology()
    rows = []
    bonds = []
    coords = {}
    h = 10000
    for r in range(dp):
        for j, e in enumerate(elements, 1):
            i = 101 + 17 * (r * len(elements) + j)
            meta = {
                "chain_id": "A",
                "repeat_unit_index": r,
                "repeat_unit_type": "A",
                "source_repeat_atom_index": j,
            }
            rows.append(
                AtomSite(
                    i,
                    e,
                    {"C": 12.011, "O": 15.999}[e],
                    metadata=meta,
                    element=e,
                    atomic_number={"C": 6, "O": 8}[e],
                )
            )
            coords[i] = (float(r), float(j), 0.0)
            index = r * len(elements) + j - 1
            neighbors = 1 if index in (0, dp * len(elements) - 1) else 2
            for n in range((4 if e == "C" else 2) - neighbors):
                h += 19
                rows.append(
                    AtomSite(
                        h,
                        "H",
                        1.008,
                        metadata={
                            k: v
                            for k, v in meta.items()
                            if k != "source_repeat_atom_index"
                        }
                        | {"generated_hydrogen": True},
                        element="H",
                        atomic_number=1,
                    )
                )
                coords[h] = (float(r), float(j), float(n + 1))
                bonds.append((i, h))
            if index:
                bonds.append((i - 17, i))
    for a in reversed(rows) if reverse else rows:
        graph.add_site(a)
    for a, b in reversed(bonds) if reverse else bonds:
        graph.add_bond(a, b, order=1)
    polymer = {
        "source_psmiles": definition,
        "architecture": "linear",
        "polymer_type": "homopolymer",
        "sequence": ["A"] * dp,
        "repeat_unit_definitions": {"A": definition},
        "composition_counts": {"A": dp},
        "degree_of_polymerization": dp,
        "number_of_repeat_units": dp,
        "number_of_inter_repeat_unit_bonds": dp - 1,
        "chain_id": "A",
        "head_site_id": 118,
        "tail_site_id": 101 + 17 * dp * len(elements),
        "coordinate_generation": {
            "template_seed": seed,
            "assembly_seed": seed,
            "method": "local_templates_self_avoiding_random_walk",
            "success": True,
        },
    }
    return MolecularSystem(graph, Coordinates(coords), metadata={"polymer": polymer})


def synthetic_reference(dp=3, seed=2026, interior=0.001):
    """Self-consistent synthetic signed envelope, not computational authenticity."""
    system = molecule(dp, seed=seed)
    corr = repeat_correspondence(system)
    assert corr["compatible"]
    charges = {}
    for r in corr["repeats"]:
        net = (
            0.05
            if r["role"] == "head"
            else -0.05 - (dp - 2) * interior
            if r["role"] == "tail"
            else interior
        )
        for k, g in enumerate(r["groups"]):
            charges.update({h: 0.01 for h in g["hydrogen_ids"]})
            charges[g["site_id"]] = -0.01 * len(g["hydrogen_ids"]) + (
                net if k == 0 else 0
            )
    cr = ProvidedChargeEngine().assign(
        system,
        charges,
        source="SYNTHETIC SOFTWARE FIXTURE; no QM executed",
        tolerance=0.002,
    )
    xyz = {str(s): system.coordinates.get(s).tolist() for s in system.topology.sites}
    mapping = dict(enumerate(sorted(charges)))
    r = {
        "schema": "island_ambertools_preparation_v2",
        "engine_version": "2",
        "requested_force_field": "gaff2",
        "charge_method": "am1bcc",
        "input_coordinates_angstrom": xyz,
        "input_coordinate_signature": digest(xyz),
        "charge_validation_tolerance_e": 0.002,
        "lineage": {"prmtop_index_to_site_id": mapping},
        "artifact_sha256": {"result.prmtop": "a" * 64},
        "charge_outcome": {
            "mode": "am1bcc",
            "qm_run": True,
            "convergence_marker": "Calculation Completed",
            "sqm_out_sha256": "b" * 64,
        },
        "stages": [
            {"stage": s, "returncode": 0} for s in ("antechamber", "parmchk2", "tleap")
        ],
    }
    content = {
        "charge_signature": cr.result_signature,
        "graph": graph_signature(system.topology),
        "mapping": sorted(mapping.items()),
        "source_sha256": "a" * 64,
        "provenance": {
            "ambertools_preparation": deepcopy(r),
            "force_field": "gaff2",
            "charge_method": "AM1-BCC",
        },
    }
    r["imported_result_signature"] = digest(content)
    source = system_data(system)
    p = {
        "schema": "island_oligomer_charge_reference_v1",
        "system": source,
        "system_sha256": storage.checksum(storage.json_bytes(storage.encode(source))),
        "correspondence": corr,
        "charge_result": record(cr),
        "import_content": content,
        "preparation": {"record": r, "record_signature": digest(r)},
        "seeds": {"template_seed": seed, "assembly_seed": seed},
        "unit": "elementary_charge",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    ref = ChargeReference(pack(p))
    ref.validate_integrity()
    return ref


@pytest.mark.parametrize("oxygen", [False, True])
def test_mapping_and_ordering(oxygen):
    a, b = molecule(5, oxygen), molecule(5, oxygen, reverse=True)
    before = deepcopy(a.to_dict())
    ca = repeat_correspondence(a)
    assert ca == repeat_correspondence(b) and ca["compatible"]
    assert a.to_dict() == before
    assert len(ca["repeats"][0]["groups"][0]["hydrogen_ids"]) == 3
    assert len(ca["repeats"][2]["groups"][0]["hydrogen_ids"]) == 2
    if oxygen:
        assert len(ca["repeats"][-1]["groups"][-1]["hydrogen_ids"]) == 1
        assert not ca["repeats"][2]["groups"][-1]["hydrogen_ids"]


@pytest.mark.parametrize(
    "mutation",
    [
        lambda s: s.metadata["polymer"].update(head_site_id=999),
        lambda s: s.metadata["polymer"].update(source_psmiles="[*:1]OCC[*:2]"),
        lambda s: s.topology.sites[118].metadata.pop("source_repeat_atom_index"),
        lambda s: s.topology.sites[135].metadata.update(source_repeat_atom_index=1),
        lambda s: setattr(s.topology.sites[118], "formal_charge", 1),
        lambda s: setattr(s.topology.sites[118], "mass", 13.0),
        lambda s: s.topology.sites[118].metadata.update(cip_label="R"),
        lambda s: s.topology.sites[10019].metadata.update(repeat_unit_index=1),
        lambda s: s.metadata["polymer"]["sequence"].__setitem__(1, "B"),
    ],
)
def test_ambiguous_or_unsupported_provenance(mutation):
    s = molecule()
    mutation(s)
    result = repeat_correspondence(s)
    assert not result["compatible"] and result["diagnostics"]


def test_repeat_sums_extrapolation_and_conformation():
    refs = [synthetic_reference(dp) for dp in (3, 5, 7)] + [
        synthetic_reference(5, 80317, interior=0.002)
    ]
    audit = audit_charge_references(refs)
    audit.validate_integrity()
    p = audit.payload["audit"]
    for s in p["summaries"]:
        assert sum(r["net_charge"] for r in s["repeats"]) == pytest.approx(0, abs=1e-15)
    first = p["diagnostic_extrapolations"][0]
    assert first["predicted_total_e"][20] == pytest.approx(0.017)
    assert first["predicted_total_e"][100] == pytest.approx(0.097)
    pair = next(
        c for c in p["comparisons"] if c["kind"] == "conformation_seed_variation"
    )
    assert pair["metrics"]["heavy_difference"]["maximum_absolute"] == pytest.approx(
        0.003
    )
    assert pair["metrics"]["hydrogen_sum_difference"]["maximum_absolute"] == 0


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p: p.update(production_validated=True),
        lambda p: p["charge_result"]["assignments"][118].update(charge=0.7),
        lambda p: p["charge_result"]["assignments"][118].update(charge=float("nan")),
        lambda p: p["system"]["sites"][0].update(formal_charge=1),
        lambda p: p["preparation"]["record"].update(imported_result_signature="a" * 64),
        lambda p: p["import_content"].update(graph="b" * 64),
        lambda p: p["seeds"].update(assembly_seed=4),
        lambda p: p["preparation"]["record"]["stages"][0].update(returncode=1),
    ],
)
def test_tamper_rejected_before_publication(tmp_path, mutation):
    p = synthetic_reference().payload
    mutation(p)
    with pytest.raises((ChargeReferenceError, ValueError)):
        save_record(ChargeReference(pack(p)), tmp_path / "bad.json")
    assert not (tmp_path / "bad.json").exists()


def test_ownership_strict_json_exclusive_publication(tmp_path):
    ref = synthetic_reference()
    p = ref.payload
    p["seeds"]["assembly_seed"] = 0
    assert ref.payload["seeds"]["assembly_seed"] == 2026
    assert deepcopy(ref) == ref and replace(ref) == ref
    path = tmp_path / "ref.json"
    save_record(ref, path)
    raw = path.read_bytes()
    assert load_charge_reference(path).identity == ref.identity
    with pytest.raises(ChargeReferenceError):
        save_record(ref, path)
    assert path.read_bytes() == raw
    for raw in ("not json", "{}", '{"payload":{},"payload":{}}', '{"payload":NaN}'):
        with pytest.raises(ChargeReferenceError):
            ChargeReference(raw).validate_integrity()


def test_simple_comparisons_without_scientific_dependencies(tmp_path):
    path = tmp_path / "ref.json"
    save_record(synthetic_reference(), path)
    code = """
import sys
class Block:
 def find_spec(self, fullname, *args):
  if fullname.split('.')[0] in {'rdkit','parmed','openmm','scipy'}: raise ImportError(fullname)
sys.meta_path.insert(0,Block())
from island.charge_references import load_charge_reference,audit_charge_references
r=load_charge_reference(sys.argv[1]);audit_charge_references([r]).validate_integrity()
"""
    subprocess.run([sys.executable, "-c", code, str(path)], check=True)


def test_audit_tamper_is_recomputed():
    from island.charge_references import ChargeAudit

    p = audit_charge_references([synthetic_reference()]).payload
    p["audit"]["diagnostic_extrapolations"][0]["predicted_total_e"][20] = 0
    with pytest.raises(ChargeReferenceError):
        ChargeAudit(pack(p)).validate_integrity()


def test_scrambled_stable_ids_do_not_choose_local_identity():
    s = molecule(5, oxygen=True)
    old = repeat_correspondence(s)
    ids = list(s.topology.sites)
    remap = {i: 90000 - 37 * n for n, i in enumerate(ids)}
    graph = Topology()
    for i in reversed(ids):
        site = deepcopy(s.topology.sites[i])
        site.id = remap[i]
        graph.add_site(site)
    for b in s.topology.bonds.values():
        graph.add_bond(remap[b.site1], remap[b.site2], order=b.order)
    meta = deepcopy(s.metadata)
    for key in ("head_site_id", "tail_site_id"):
        meta["polymer"][key] = remap[meta["polymer"][key]]
    other = MolecularSystem(
        graph, Coordinates({remap[i]: s.coordinates.get(i) for i in ids}), metadata=meta
    )
    new = repeat_correspondence(other)
    assert new["compatible"]
    for a, b in zip(old["repeats"], new["repeats"], strict=True):
        for ga, gb in zip(a["groups"], b["groups"], strict=True):
            assert gb["site_id"] == remap[ga["site_id"]]
            assert set(gb["hydrogen_ids"]) == {remap[h] for h in ga["hydrogen_ids"]}
            assert gb["heavy_environment"] == ga["heavy_environment"]


def test_correctly_resigned_failed_calculation_still_rejected():
    p = synthetic_reference().payload
    r = p["preparation"]["record"]
    r["stages"][0]["returncode"] = 1
    p["import_content"]["provenance"]["ambertools_preparation"] = {
        k: v for k, v in r.items() if k != "imported_result_signature"
    }
    r["imported_result_signature"] = digest(p["import_content"])
    p["preparation"]["record_signature"] = digest(r)
    with pytest.raises(ChargeReferenceError, match="Failed or incomplete"):
        ChargeReference(pack(p)).validate_integrity()


def test_cli_failed_matrix_and_exclusive_output(tmp_path):
    from scripts.validate_oligomer_charges import MATRIX, SETTINGS, main

    source = tmp_path / "source"
    source.mkdir()
    rows = [
        dict(
            c,
            case=f"{c['chemistry']}-dp{c['dp']}-seed{c['seed']}",
            status="failed",
            failure="Injected unavailable tools",
        )
        for c in MATRIX
    ]
    storage.publish(source / "declared-settings.json", storage.json_bytes(SETTINGS))
    storage.publish(
        source / "matrix.json", storage.json_bytes({"cases": rows, "complete": False})
    )
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    output = tmp_path / "audit"
    assert main(["--audit-existing", str(source), "--output", str(output)]) == 1
    assert storage.read_json(output / "outcomes.json")["validated_reference_count"] == 0
    assert {p.name: p.read_bytes() for p in source.iterdir()} == before
    with pytest.raises(FileExistsError):
        main(["--audit-existing", str(source), "--output", str(output)])


def test_creation_requires_actual_validated_preparation_type():
    from types import SimpleNamespace

    from island.charge_references import create_charge_reference

    with pytest.raises(ChargeReferenceError, match="AmberToolsPreparationResult"):
        create_charge_reference(
            molecule(), SimpleNamespace(validate_integrity=lambda _: None)
        )


def test_cli_declared_success_exit_is_zero_software_injection(tmp_path, monkeypatch):
    from scripts import validate_oligomer_charges as cli

    # Only exit dispatch is injected; this is not evidence of live matrix success.
    monkeypatch.setattr(cli, "audit_existing", lambda *_: True)
    assert (
        cli.main(
            [
                "--audit-existing",
                str(tmp_path / "source"),
                "--output",
                str(tmp_path / "new"),
            ]
        )
        == 0
    )
