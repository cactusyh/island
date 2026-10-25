"""C4 full-converter semantics, chemical boundaries and stale evidence controls."""

import copy
import os
import subprocess
import sys
import tarfile
from argparse import Namespace
from decimal import Decimal
from pathlib import Path

import pytest

from island.charge_references.records import pack, unpack
from island.chemistry import from_smiles
from island.forcefields.pcff.sip_variants import SiPSourceVariant
from scripts.pcff_j41_c4_reference import (
    angle_angle,
    dependencies,
    reconstruct,
    resolve,
    select,
    table,
    trace_audit,
)
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c3 import sources
from scripts.phase_4j41_c4 import BLOCKS, CASES, execute, mapping
from scripts.validate_phase_4j41_c4_receipt import validate_case


@pytest.fixture(scope="module")
def variants(tmp_path_factory):
    folder = Path(
        os.environ.get(
            "ISLAND_J41_DISTRIBUTION", "../lammps/lammps/tools/msi2lmp/frc_files"
        )
    )
    if folder.is_dir():
        return {v.filename: v for v in sources(folder)}
    result = {}
    with tarfile.open("docs/evidence/phase_4j41_c3/execution.tar.xz") as tar:
        for name in ("pcff.frc", "compass_published.frc"):
            prefix = "study/" + name.removesuffix(".frc") + "/"
            p = unpack(tar.extractfile(prefix + "inventory.json").read().decode())
            result[name] = SiPSourceVariant(
                name,
                tar.extractfile(prefix + name).read(),
                p["sha256"],
                pack(p["provenance"]),
                tuple(
                    (n, tar.extractfile(prefix + n).read(), h)
                    for n, h in p["companion_hashes"].items()
                ),
            )
    return result


def materialize(variant, folder):
    folder.mkdir(exist_ok=True)
    (folder / variant.filename).write_bytes(variant.raw)
    for n, raw, _ in variant.companions:
        (folder / n).write_bytes(raw)


def row(key, value, version="1", rid="test"):
    # Synthetic matcher unit controls do not enter a force-field inventory.
    return {
        "types": key.split(),
        "values": [Decimal(value)],
        "version": Decimal(version),
        "id": rid,
    }


def test_converter_version_table_keeps_first_slot_and_equal_version():
    rows = [
        row("a b", "1", rid="old"),
        row("b a", "2", rid="reverse"),
        row("a b", "3", "2", "new"),
        row("a b", "4", "2", "tie"),
    ]
    assert [r["id"] for r in table(rows)] == ["new", "reverse"]
    assert select(rows, ["a", "b"])["row_id"] == "new"
    assert select(rows, ["b", "a"])["row_id"] == "new"


def test_exact_pass_precedes_wildcards_but_wildcards_use_file_order():
    rows = [
        row("* b *", "1", rid="broad"),
        row("a b *", "2", rid="specific"),
        row("a b c", "3", rid="exact"),
    ]
    assert select(rows, ["a", "b", "c"])["row_id"] == "exact"
    assert select(rows[:2], ["a", "b", "c"])["row_id"] == "broad"
    assert select(rows[:2][::-1], ["a", "b", "c"])["row_id"] == "specific"


def test_role_preserving_angle_angle_does_not_reverse_center():
    from collections import defaultdict

    rows = defaultdict(list)
    rows["angle-angle", "cff91"] = [
        row("d b c a", "5", rid="role"),
        row("d c b a", "9", rid="bad"),
    ]
    result = angle_angle(rows, "cff91", ["a", "b", "c", "d"])
    assert result["row_id"] == "role"
    assert [r["row_id"] for r in result["caller_slots"]] == ["role", None, None]
    assert result["caller_success_is_not_completeness"]


def test_conflicting_wildcard_file_order_is_not_physical_authority():
    from collections import defaultdict

    rows = defaultdict(list)
    rows["quartic_angle", "cff91"] = [
        row("* b c", "1", rid="first"),
        row("a b *", "2", rid="second"),
    ]
    result = resolve(rows, "cff91", "quartic_angle", ["a", "b", "c"])
    assert result["row_id"] == "first"
    assert result["status"] == "source_row_ambiguous"


def test_auto_table_cannot_complete_the_full_converter(variants):
    from scripts.pcff_j5_reference import read_source

    rows = read_source(variants["pcff.frc"].raw)
    result = resolve(rows, "cff91", "quartic_angle", ["h_", "c_", "h_"])
    assert result["status"] == "source_row_absent"
    assert all(s["namespace"] == "cff91" for s in result["searches"])


@pytest.mark.parametrize("filename", ["pcff.frc", "compass_published.frc"])
@pytest.mark.parametrize(
    "name", ["methylsilane", "disiloxane", "hexamethyldisiloxane", "phosphazene_probe"]
)
def test_priority_domains_cannot_be_promoted_by_converter_labels(
    variants, tmp_path, filename, name
):
    _, smiles, domain = next(c for c in CASES if c[0] == name)
    v = variants[filename]
    materialize(v, tmp_path)
    system = from_smiles(smiles, random_seed=2026)
    labels, typing = mapping(system, v, domain, tmp_path)
    audit = reconstruct(system, v, labels, typing)
    assert audit["status"] in ("source_incomplete", "typing_blocked")
    assert not audit["execution_authorized"] and audit["execution_profile"] is None
    if labels is not None:
        assert audit["missing_couplings"]
        assert any(
            q["resolution"]["status"] == "source_row_absent"
            for q in audit["missing_couplings"]
        )
        assert audit["charges"]["complete"]
        assert not audit["charges"]["normalized"]


@pytest.mark.parametrize(
    "name",
    [
        "phosphate_ester",
        "phosphonate",
        "phosphine_oxide",
        "protonated_p",
        "pf5",
        "mixed_sip",
        "unsupported_si",
    ],
)
def test_chemical_boundary_is_not_reinterpreted(variants, tmp_path, name):
    _, smiles, domain = next(c for c in CASES if c[0] == name)
    s = from_smiles(smiles, random_seed=2026)
    for v in variants.values():
        materialize(v, tmp_path)
        labels, typing = mapping(s, v, domain, tmp_path)
        audit = reconstruct(s, v, labels, typing)
        assert labels is None and audit["status"] == "typing_blocked"
        assert audit["charges"] is None and audit["interactions"] == []


def test_skipped_bond13_zero_is_not_a_source_coefficient(variants, tmp_path):
    v = variants["pcff.frc"]
    materialize(v, tmp_path)
    s = from_smiles("C[SiH3]", random_seed=2026)
    labels, typing = mapping(s, v, "mixed_silane", tmp_path)
    p = reconstruct(s, v, labels, typing)
    skipped = [q for q in p["missing_couplings"] if q["family"] == "bond-bond_1_3"]
    assert len(skipped) == 9
    for q in skipped:
        assert q["resolution"]["status"] == "converter_zero_initialization"
        assert (
            q["resolution"]["independent_physical_lookup"]["status"]
            == "source_row_absent"
        )
    # A skipped lookup and initialized zero cannot replace the absent source row.
    assert p["status"] == "source_incomplete"


def test_equilibrium_dependencies_preserve_roles_and_cannot_supply_couplings():
    assert dependencies("bond-angle", [1, 2, 3]) == [
        ("quartic_bond", [1, 2]),
        ("quartic_bond", [2, 3]),
    ]
    assert dependencies("end_bond-torsion_3", [1, 2, 3, 4]) == [
        ("quartic_bond", [1, 2]),
        ("quartic_bond", [3, 4]),
    ]
    assert dependencies("angle-angle", [1, 2, 3, 4]) == [
        ("quartic_angle", [1, 2, 3]),
        ("quartic_angle", [3, 2, 4]),
        ("quartic_angle", [1, 2, 4]),
    ]


@pytest.mark.parametrize(
    "mutation",
    [
        "source",
        "companion",
        "graph",
        "order",
        "charge",
        "dependency",
        "coupling",
        "zero",
        "profile",
        "p_label",
    ],
)
def test_reserialized_stale_or_incomplete_evidence_rejected(
    variants, tmp_path, mutation
):
    from island.workflows.bundle import system_data

    v = variants["pcff.frc"]
    materialize(v, tmp_path)
    s = from_smiles("C[SiH3]", random_seed=2026)
    case = tmp_path / "case"
    case.mkdir()
    write(case / "system.json", system_data(s))
    labels, typing = mapping(s, v, "mixed_silane", tmp_path)
    p = copy.deepcopy(reconstruct(s, v, labels, typing))
    if mutation == "source":
        p["source_sha256"] = variants["compass_published.frc"].expected_sha256
    elif mutation == "companion":
        p["companion_hashes"]["pcff.rlb"] = "0" * 64
    elif mutation == "graph":
        p["system_sha256"] = "0" * 64
    elif mutation == "order":
        p["atom_order"].reverse()
    elif mutation == "charge":
        p["charges"]["partial_charges"][1] = "0"
    elif mutation == "dependency":
        p["dependencies"].pop()
    elif mutation == "coupling":
        p["missing_couplings"] = []
    elif mutation == "zero":
        p["interactions"][0]["resolution"]["values"] = ["0", "0"]
    elif mutation == "p_label":
        p["labels"][2] = "p="
    else:
        p["execution_profile"] = "invented_c4_profile"
    write(case / "raw-reconstruction.json", p)
    with pytest.raises(ValueError, match="Contradictory"):
        validate_case(case, v, "mixed_silane")


def test_unavailable_converter_is_a_classified_nonzero_result(tmp_path):
    args = Namespace(
        output=tmp_path / "out",
        source=tmp_path / "absent",
        lammps=tmp_path / "absent",
        distribution=tmp_path / "absent",
    )
    assert execute(args) == 1
    p = read(args.output / "report.json")
    assert p["status"] == "converter_unavailable"
    assert p["periodic"] == BLOCKS and p["acceptance_gate"] == 1


def test_full_converter_raw_trace_rejects_wrong_selection(variants):
    with pytest.raises(ValueError, match="disagrees"):
        trace_audit("C4_MATCH #quartic_bond 2 999 0 0 c si", variants["pcff.frc"])


def test_public_cli_requires_no_user_car_mdf():
    r = subprocess.run(
        [sys.executable, "scripts/phase_4j41_c4.py", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0
    assert "--car" not in r.stdout and "--mdf" not in r.stdout
