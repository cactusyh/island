"""Analytical fragment transfer and synthetic provided-charge software tests."""

import subprocess
import sys
from copy import deepcopy

import pytest

pytest.importorskip("rdkit")

from island.builders import build_linear_polymer
from island.charge_references.records import pack
from island.exceptions import FragmentChargeError
from island.fragment_charges import (
    AM1BCCFragmentChargeBackend,
    FragmentChargeAssignment,
    FragmentChargeTemplate,
    ProvidedFragmentChargeBackend,
    assign_fragment_charges,
    create_fragment_template,
    load_fragment_record,
    prepare_capped_fragment,
    save_fragment_record,
)
from island.workflows import storage
from island.workflows.bundle import system_data, system_from

PE = "[*:1]CC[*:2]"


@pytest.fixture
def analytical():
    f = prepare_capped_fragment(PE)
    p = f.payload
    m = p["mapping"]
    q = {s["id"]: 0.05 for s in p["system"]["sites"]}
    q[m["heavy"][1]] = -0.3
    q[m["heavy"][2]] = -0.5
    q[m["caps"]["head"]["site_id"]] = 0.2
    q[m["caps"]["tail"]["site_id"]] = 0.4
    result = ProvidedFragmentChargeBackend(
        q, source="synthetic analytical test"
    ).calculate(f)
    return (
        f,
        q,
        create_fragment_template(
            result, conservation_policy="strict", charge_tolerance=1e-12
        ),
    )


@pytest.mark.parametrize("dp", [1, 2, 3, 4, 50])
def test_analytical_cap_transfer(analytical, dp):
    f, _q, t = analytical
    s = build_linear_polymer(PE, dp=dp, generate_3d=False)
    before = system_data(s)
    result = assign_fragment_charges(t, s)
    d = result.payload["data"]
    charges = d["assignments"]
    assert d["total_e"] == pytest.approx(0, abs=1e-12)
    assert len(d["removed_cap_transfers"]) == 2 * (dp - 1)
    for g in d["groups"]:
        parent = s.topology.sites[g["parent_site_id"]]
        r = parent.metadata["repeat_unit_index"]
        j = parent.metadata["source_repeat_atom_index"]
        expected = (
            (-0.3 if r == 0 else -0.1) if j == 1 else (-0.5 if r == dp - 1 else -0.1)
        )
        assert charges[parent.id] == pytest.approx(expected, abs=1e-15)
        assert len(g["target_hydrogen_ids"]) == g["native_count"] + len(
            g["retained_caps"]
        )
        assert len({charges[h] for h in g["target_hydrogen_ids"]}) == 1
    assert (
        before == system_data(s)
        and f.payload["mapping"] == analytical[0].payload["mapping"]
    )
    copied = result.payload
    copied["data"]["assignments"].clear()
    assert result.payload["data"]["assignments"]


def test_shared_parent():
    definition = "[*:1]C[*:2]"
    f = prepare_capped_fragment(definition)
    m = f.payload["mapping"]
    q = {s["id"]: 0.15 for s in f.payload["system"]["sites"]}
    q[m["heavy"][1]] = -0.6
    q[m["caps"]["head"]["site_id"]] = 0.1
    q[m["caps"]["tail"]["site_id"]] = 0.2
    t = create_fragment_template(
        ProvidedFragmentChargeBackend(q).calculate(f), conservation_policy="strict"
    )
    s = build_linear_polymer(definition, dp=3, generate_3d=False)
    d = assign_fragment_charges(t, s).payload["data"]
    middle = next(g for g in d["groups"] if g["repeat"] == 1)
    assert set(middle["removed_caps"]) == {"head", "tail"}
    assert d["assignments"][middle["parent_site_id"]] == pytest.approx(-0.3)
    assert d["total_e"] == pytest.approx(0, abs=1e-12)


@pytest.mark.parametrize(
    "definition",
    ["[*:1]CCO[*:2]", "[*:1]CC(c1ccccc1)[*:2]", "[*:1]CC(C)(C(=O)OC)[*:2]"],
)
def test_generic_ring_branch_hetero(definition):
    f = prepare_capped_fragment(definition)
    q = {s["id"]: i * 0.01 for i, s in enumerate(f.payload["system"]["sites"])}
    result = ProvidedFragmentChargeBackend(q).calculate(f)
    with pytest.raises(FragmentChargeError):
        create_fragment_template(
            result, conservation_policy="strict", charge_tolerance=1e-12
        )
    t = create_fragment_template(result, conservation_policy="uniform_fragment_l2_v1")
    for dp in (2, 3):
        s = build_linear_polymer(definition, dp=dp, generate_3d=False)
        assert abs(assign_fragment_charges(t, s).payload["data"]["total_e"]) < 1e-12
    assert result.payload["raw_charges"] == q


@pytest.mark.parametrize(
    "definition",
    [
        "[*]CC[*]",
        "[*:1]=CC[*:2]",
        "[*:1][C@H](F)C[*:2]",
        "[*:1][13CH2]C[*:2]",
        "[*:1][CH]C[*:2]",
        "[*:1]C.[*:2]C",
        "[*:1]C/C=C/C[*:2]",
        "[*:1]C[NH2+]C[*:2]",
    ],
)
def test_unsupported(definition):
    with pytest.raises(FragmentChargeError):
        prepare_capped_fragment(definition)


def test_order_ids_cache_and_offline(analytical, tmp_path):
    _, _, t = analytical
    s = build_linear_polymer(PE, dp=4, generate_3d=False)
    original = assign_fragment_charges(t, s).payload["data"]["assignments"]
    p = system_data(s)
    rename = {a["id"]: a["id"] * 17 + 51 for a in p["sites"]}
    for a in p["sites"]:
        a["id"] = rename[a["id"]]
    p["sites"].reverse()
    for b in p["bonds"]:
        b["site1"] = rename[b["site1"]]
        b["site2"] = rename[b["site2"]]
    p["bonds"].reverse()
    p["coordinates"] = {rename[k]: v for k, v in p["coordinates"].items()}
    for k in ("head_site_id", "tail_site_id"):
        p["metadata"]["polymer"][k] = rename[p["metadata"]["polymer"][k]]
    changed = system_from(p)
    new = assign_fragment_charges(t, changed).payload["data"]["assignments"]
    assert new == {rename[k]: v for k, v in original.items()}
    save_fragment_record(t, tmp_path / "t.json")
    assert (
        load_fragment_record(
            tmp_path / "t.json", expected_cache_key=t.cache_key
        ).identity
        == t.identity
    )
    with pytest.raises(FragmentChargeError):
        load_fragment_record(tmp_path / "t.json", expected_cache_key="wrong")
    with pytest.raises(FragmentChargeError):
        save_fragment_record(t, tmp_path / "t.json")
    storage.publish(tmp_path / "s.json", storage.json_bytes(storage.encode(p)))
    code = """
import sys
for name in ('rdkit','openmm','parmed','scipy'):sys.modules[name]=None
from island.fragment_charges import load_fragment_record,assign_fragment_charges
from island.workflows import storage
from island.workflows.bundle import system_from
from pathlib import Path
t=load_fragment_record(sys.argv[1]);s=system_from(storage.decode(storage.read_json(Path(sys.argv[2]))))
assign_fragment_charges(t,s).validate_integrity()
"""
    subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(tmp_path / "t.json"),
            str(tmp_path / "s.json"),
        ],
        check=True,
    )


@pytest.mark.parametrize(
    "defect", ["graph", "orientation", "hydrogen", "stereo", "repeat"]
)
def test_target_mismatch(analytical, defect):
    t = analytical[2]
    s = build_linear_polymer(PE, dp=3, generate_3d=False)
    if defect == "orientation":
        s.metadata["polymer"]["head_site_id"] = s.metadata["polymer"]["tail_site_id"]
    if defect == "repeat":
        s.metadata["polymer"]["source_psmiles"] = "[*:1]CO[*:2]"
    if defect == "stereo":
        next(iter(s.topology.sites.values())).metadata["cip_label"] = "R"
    if defect == "graph":
        from dataclasses import replace

        key = next(iter(s.topology.bonds))
        s.topology.bonds[key] = replace(s.topology.bonds[key], order=2)
    if defect == "hydrogen":
        h = next(i for i, a in s.topology.sites.items() if a.element == "H")
        del s.topology.sites[h]
    before = deepcopy(s)
    with pytest.raises(FragmentChargeError):
        assign_fragment_charges(t, s)
    assert s.metadata == before.metadata


def test_malformed_and_compatibility(analytical):
    f, q, t = analytical
    for bad in (
        {},
        dict(q, unknown=0),
        {k: True for k in q},
        {k: float("nan") for k in q},
    ):
        with pytest.raises(FragmentChargeError):
            ProvidedFragmentChargeBackend(bad).calculate(f)
    for family in ("opls", "pcff"):
        with pytest.raises(FragmentChargeError):
            AM1BCCFragmentChargeBackend(force_field=family, work_root="unused")
    with pytest.raises(FragmentChargeError, match="Provided charges"):
        assign_fragment_charges(
            t,
            build_linear_polymer(PE, dp=2, generate_3d=False),
            force_field={"family": "gaff2", "data_sha256": "a" * 64},
        )
    p = t.payload
    p["data"]["transformed_charges"][next(iter(q))] += 0.1
    with pytest.raises(FragmentChargeError):
        FragmentChargeTemplate(pack(p)).validate_integrity()
    p = assign_fragment_charges(
        t, build_linear_polymer(PE, dp=2, generate_3d=False)
    ).payload
    p["data"]["removed_cap_transfers"] = []
    with pytest.raises(FragmentChargeError):
        FragmentChargeAssignment(pack(p)).validate_integrity()


def test_charge_backend_runner_and_semantic_failures(monkeypatch, tmp_path):
    """Injected tool output, explicitly synthetic; real QM is a separate gate."""
    from types import SimpleNamespace

    from island.forcefields.ambertools import engine
    from island.fragment_charges.core import FragmentChargeResult

    executable = tmp_path / "synthetic-tool"
    executable.write_text("SYNTHETIC SOFTWARE TEST")
    data = tmp_path / "gaff2.dat"
    data.write_text("SYNTHETIC SOFTWARE TEST")
    chain = SimpleNamespace(
        executables={k: str(executable) for k in ("antechamber", "parmchk2", "tleap")},
        amberhome=tmp_path,
        data_file=data,
        versions={"antechamber": "SYNTHETIC"},
        package={},
    )
    monkeypatch.setattr(engine, "_discover", lambda options: chain)
    calls = []

    def run(stage, command, directory, timeout, **kwargs):
        calls.append(command)
        (directory / "typed.mol2").write_text((directory / "input.mol2").read_text())
        (directory / "sqm.out").write_text("SYNTHETIC Calculation Completed")
        for stream in ("stdout", "stderr"):
            (directory / f"antechamber.{stream}.log").write_text("")
        return {
            "stage": stage,
            "command": command,
            "returncode": 0,
            "elapsed_seconds": 0.0,
            "stdout_tail": "",
            "stderr_tail": "",
            "stdout_sha256": storage.checksum(b""),
            "stderr_sha256": storage.checksum(b""),
        }

    monkeypatch.setattr(engine, "_run_stage", run)
    result = AM1BCCFragmentChargeBackend(
        force_field="gaff2", work_root=tmp_path
    ).calculate(prepare_capped_fragment(PE))
    assert len(calls) == 1
    t = create_fragment_template(result, conservation_policy="uniform_fragment_l2_v1")
    for dp in (2, 4):
        s = build_linear_polymer(PE, dp=dp, generate_3d=False)
        assign_fragment_charges(
            t, s, force_field=result.payload["backend"]["force_field"]
        ).validate_integrity()
    assert len(calls) == 1
    for defect in ("command", "completion", "force_field"):
        p = result.payload
        if defect == "command":
            p["evidence"]["stage"]["command"][
                p["evidence"]["stage"]["command"].index("bcc")
            ] = "rc"
        if defect == "completion":
            p["evidence"]["files"]["sqm.out"] = "SYNTHETIC incomplete"
            p["evidence"]["sha256"]["sqm.out"] = storage.checksum(
                b"SYNTHETIC incomplete"
            )
        if defect == "force_field":
            p["backend"]["force_field"]["data_sha256"] = "a" * 64
        with pytest.raises(FragmentChargeError):
            FragmentChargeResult(pack(p)).validate_integrity()
    for family in ("opls", "pcff", "gaff"):
        with pytest.raises(FragmentChargeError):
            assign_fragment_charges(
                t,
                s,
                force_field={
                    "family": family,
                    "data_sha256": result.payload["backend"]["force_field"][
                        "data_sha256"
                    ],
                },
            )

    def failed(*args, **kwargs):
        raise RuntimeError("injected tool failure")

    monkeypatch.setattr(engine, "_run_stage", failed)
    with pytest.raises(FragmentChargeError, match="injected tool failure"):
        AM1BCCFragmentChargeBackend(force_field="gaff2", work_root=tmp_path).calculate(
            prepare_capped_fragment(PE)
        )
    assert list(tmp_path.glob("island-fragment-*/input.mol2"))


def test_aggregate_gate_and_no_repeat_execution(tmp_path):
    from scripts import validate_fragment_charges as cli

    rows = [{"status": "passed"} for _ in range(4)]
    rows[0]["integration"] = {"status": "failed"}
    assert cli.summarize(rows)["complete_fragment_mapping"]
    assert not cli.summarize(rows)["complete"]
    rows[0]["integration"]["status"] = "passed"
    assert cli.summarize(rows)["complete"]
    rows[2]["status"] = "failed"
    assert not cli.summarize(rows)["complete"]
    assert not cli.summarize([])["complete"]
    destination = tmp_path / "new"
    assert cli.main(["--declare-only", "--output", str(destination)]) == 0
    before = (destination / "declaration.json").read_bytes()
    assert cli.main(["--declare-only", "--output", str(destination)]) == 1
    assert (destination / "declaration.json").read_bytes() == before


def test_offline_aromatic_evidence_rejects_all_single_ring():
    from island.forcefields.ambertools.lineage import prepare_input
    from island.fragment_charges.amber import _mol2

    f = prepare_capped_fragment("[*:1]CC(c1ccccc1)[*:2]")
    system = system_from(f.payload["system"])
    prepared = prepare_input(system, None)
    _mol2(prepared.mol2_text, system, prepared.names)
    malformed = prepared.mol2_text.replace(" ar\n", " 1\n")
    assert malformed != prepared.mol2_text
    with pytest.raises(FragmentChargeError, match="pi-valence"):
        _mol2(malformed, system, prepared.names)


@pytest.mark.parametrize(
    "definition",
    [
        "[*:1]Cc1ccoc1C[*:2]",
        "[*:1]Cc1ccsc1C[*:2]",
        "[*:1]Cc1cc[nH]c1C[*:2]",
        "[*:1]Cc1ccccc1C[*:2]",
    ],
)
def test_aromatic_parent_environments(definition, tmp_path):
    fragment = prepare_capped_fragment(definition)
    q = {s["id"]: 0.0 for s in fragment.payload["system"]["sites"]}
    template = create_fragment_template(
        ProvidedFragmentChargeBackend(q, source="synthetic software test").calculate(
            fragment
        ),
        conservation_policy="strict",
    )
    save_fragment_record(template, tmp_path / "template.json")
    for dp in (1, 3):
        system = build_linear_polymer(definition, dp=dp, generate_3d=False)
        before = system_data(system)
        assignment = assign_fragment_charges(
            load_fragment_record(tmp_path / "template.json"), system
        )
        assert set(assignment.payload["data"]["assignments"]) == set(
            system.topology.sites
        )
        assert assignment.payload["data"]["total_e"] == 0
        assert system_data(system) == before


def test_aromatic_oxygen_rechecksummed_invalid_environment():
    from island.fragment_charges import CappedFragment

    f = prepare_capped_fragment("[*:1]Cc1ccoc1C[*:2]")
    p = f.payload
    oxygen = next(s["id"] for s in p["system"]["sites"] if s["element"] == "O")
    system = system_from(p["system"])
    bond = next(
        b for b in system.topology.bonds.values() if oxygen in (b.site1, b.site2)
    )
    from dataclasses import replace

    key = next(k for k, b in system.topology.bonds.items() if b is bond)
    system.topology.bonds[key] = replace(bond, order=1.0, aromatic=False)
    p["system"] = system_data(system)
    with pytest.raises(FragmentChargeError, match="aromatic"):
        CappedFragment(pack(p)).validate_integrity()


def test_interoperability_serialization_and_corruption():
    from scripts.validate_fragment_interoperability import compare_charges

    q = {1: -0.1, 9: 0.1}
    assert (
        compare_charges(q, {1: -0.1 + 3e-10, 9: 0.1})["max_per_site_difference_e"]
        < 1e-8
    )
    with pytest.raises(ValueError):
        compare_charges(q, {1: -0.09, 9: 0.09})  # Neutrality cannot mask corruption.
    with pytest.raises(ValueError):
        compare_charges({i: 0.0 for i in range(200)}, {i: 9e-9 for i in range(200)})


def test_furan_offline_mapping(tmp_path):
    f = prepare_capped_fragment("[*:1]Cc1ccoc1C[*:2]")
    q = {s["id"]: 0.01 * i for i, s in enumerate(f.payload["system"]["sites"])}
    t = create_fragment_template(
        ProvidedFragmentChargeBackend(q, source="synthetic software data").calculate(f),
        conservation_policy="uniform_fragment_l2_v1",
    )
    save_fragment_record(t, tmp_path / "t.json")
    system = build_linear_polymer(f.payload["definition"], dp=3, generate_3d=False)
    storage.publish(
        tmp_path / "s.json", storage.json_bytes(storage.encode(system_data(system)))
    )
    subprocess.run(
        [
            sys.executable,
            "-c",
            """
import sys
for name in ('rdkit', 'parmed', 'openmm', 'scipy'): sys.modules[name] = None
from pathlib import Path
from island.fragment_charges import load_fragment_record, assign_fragment_charges
from island.workflows import storage
from island.workflows.bundle import system_from
s = system_from(storage.decode(storage.read_json(Path(sys.argv[2]))))
a = assign_fragment_charges(load_fragment_record(sys.argv[1]), s)
assert abs(a.payload['data']['total_e']) < 1e-12
""",
            str(tmp_path / "t.json"),
            str(tmp_path / "s.json"),
        ],
        check=True,
    )
