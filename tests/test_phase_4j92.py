"""Synthetic software contracts only; never substitute for the real PE5 bundle."""

import os
import shutil
from pathlib import Path

import pytest

from island.forcefields.assignment_core.schema import identity
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.opls_providers.intake92 import (
    AUDIT_FIELDS,
    FILES,
    READY,
    SOURCE_DIR,
    intake,
)
from island.forcefields.opls_providers.lammps92 import expand_pairs
from island.forcefields.opls_providers.polymer import finite_graph
from island.forcefields.opls_providers.registry import obligations
from scripts.phase_4j73_source_intake import pin, read, signed, write
from scripts.phase_4j92_intake import emit, replay


def control(root, *, approved=True, change=None):
    root.mkdir()
    spec = {
        "psmiles": "*CC*",
        "repeat_count": 5,
        "end_groups": ["*[H]", "*[H]"],
        "protonation_state": "as_encoded",
        "charge_state": 0,
        "stereochemistry": "achiral",
        "tacticity": "not_applicable",
        "stereo_seed": 2026,
        "intended_export_engine": "lammps",
    }
    graph = finite_graph(spec)
    mol = molecule(graph)
    needed, _ = obligations(graph)
    lines = [
        "LAMMPS data file generated with LigParGen software-control - Jorgensen Lab",
        "",
    ]
    count = len(graph.atom_order)
    for n, v in [
        ("atoms", count),
        ("bonds", len(needed["bond"])),
        ("angles", len(needed["angle"])),
        ("dihedrals", len(needed["torsion"])),
        ("atom types", count),
        ("bond types", len(needed["bond"])),
        ("angle types", len(needed["angle"])),
        ("dihedral types", len(needed["torsion"])),
    ]:
        lines.append(f"{v} {n}")
    lines += ["0 100 xlo xhi", "0 100 ylo yhi", "0 100 zlo zhi", "", "Masses", ""]
    lines += [f"{a.GetIdx() + 1} {a.GetMass()}" for a in mol.GetAtoms()]
    lines += ["", "Pair Coeffs", ""] + [f"{i} 1 1" for i in graph.atom_order]
    for section, family, values in [
        ("Bond Coeffs", "bond", "1 1"),
        ("Angle Coeffs", "angle", "1 100"),
        ("Dihedral Coeffs", "torsion", "1 1 1 1"),
    ]:
        lines += ["", section, ""] + [
            f"{i} {values}" for i, _ in enumerate(needed[family], 1)
        ]
    lines += ["", "Atoms", ""] + [f"{i} 1 {i} 0 {i} 0 0" for i in graph.atom_order]
    for section, family in [
        ("Bonds", "bond"),
        ("Angles", "angle"),
        ("Dihedrals", "torsion"),
    ]:
        lines += ["", section, ""] + [
            f"{i} {i} " + " ".join(map(str, s)) for i, s in enumerate(needed[family], 1)
        ]
    data = "\n".join(lines) + "\n"
    settings = """units real
atom_style full
boundary m m m
pair_style lj/cut/coul/cut 40.0 40.0
pair_modify mix geometric tail no
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style cvff
special_bonds lj/coul 0.0 0.0 0.5
read_data Lig.lammps.lmp
neighbor 2.0 bin
run 0
variable e equal pe
print "LIG: ${e} kcal/mol"
"""
    q = {
        "schema": "island.j92.request.v1",
        "target": "finite_PE_five_mer",
        "family": "OPLS-AA",
        "settings": spec,
        "backends": ["lammps", "openmm"],
        "atom_count": count,
        "graph_identity": graph.identity,
        "atom_mapping": [[i, i] for i in graph.atom_order],
    }
    provider = {
        "schema": "island.j92.provider.v1",
        "provider": "LigParGen",
        "version": "software-control",
        "family": "OPLS-AA",
        "source_authority": "manufactured software fixture",
        "license": {
            "identifier": "synthetic-control",
            "authorization": "authorized_local_use",
        },
        "provenance": "synthetic parser control, not scientific evidence",
        "atom_typing_origin": "software control site IDs",
        "charge_origin": "manufactured software control",
        "synthetic": True,
        "improper_policy": {
            "center_position": 1,
            "omissions": [
                {"sites": list(s), "reason": "software-only explicit omission"}
                for s in needed["improper"]
            ],
        },
    }
    if change:
        data, settings, q, provider = change(data, settings, q, provider)
    (root / "Lig.lammps.lmp").write_text(data)
    (root / "Lig.lammps.in").write_text(settings)
    write(root / "request.json", q)
    write(root / "provider_metadata.json", provider)
    audit = {k: None for k in AUDIT_FIELDS}
    audit.update(schema="island.j92.audit.v1", status="pending", independent=False)
    write(root / "independent_audit.json", audit)
    args = repin(root)
    if approved:
        r = intake(**args).report
        if "record_claims" in r:
            audit.update(
                status="approved",
                independent=True,
                reviewer="software reviewer control",
                authority="synthetic authority only",
                source_binding=r["source_binding"],
                graph_identity=q["graph_identity"],
                settings_identity=identity(r["conventions"]),
                expansion_hash=r["pair_expansion"]["expansion_hash"],
                claim_hashes=sorted(c["claim_hash"] for c in r["record_claims"]),
                improper_policy_identity=identity(provider["improper_policy"]),
            )
            write(root / "independent_audit.json", audit)
            args = repin(root)
    return args


def repin(root):
    manifest = {
        "schema": "island.j92.lammps-bundle.v1",
        "bundle_id": "software-only",
        "family": "OPLS-AA",
        "files": [
            {"path": n, "sha256": pin(root / n), "size": (root / n).stat().st_size}
            for n in FILES
        ],
    }
    write(root / "manifest.json", manifest)
    return {"source_dir": root, "manifest_sha256": pin(root / "manifest.json")}


def test_complete_two_file_adapter_synthetic(tmp_path):
    r = intake(**control(tmp_path / "source")).report
    assert r["decision"] == "audited_exact_graph", r["blockers"]
    assert (
        r["synthetic"]
        and r["scientific_graphs_qualified"] == 0
        and r["qualified_graph"] is None
    )
    assert r["pair14_count"] == 81 and r["exclusion_count"] == 91
    assert r["source_rows_inspected"] == 455 and r["source_files_inspected"] == 5
    assert all(r[k] == v for k, v in READY.items())
    assert not r["missing_terms"]
    assert r["conventions"]["conventions"]["units"] == "real"
    assert all(c["executed"] is False for c in r["conventions"]["inert_commands"])
    assert r["pair_expansion"]["source_command"]["source"]["line"] == 10


def test_pending_audit_is_provisional_not_qualified(tmp_path):
    r = intake(**control(tmp_path / "source", approved=False)).report
    assert r["decision"] == "provisional_parameterization"
    assert r["qualified_graph"] is None and r["accepted_coefficients"] == 0


@pytest.mark.parametrize(
    "old,new",
    [
        ("units real", "units metal"),
        ("mix geometric", "mix arithmetic"),
        ("0.0 0.0 0.5", "0.0 0.0 1.0"),
        ("dihedral_style opls", "dihedral_style harmonic"),
        ("40.0 40.0", "40.0 30.0"),
        ("improper_style cvff", "improper_style harmonic"),
        ("bond_style harmonic", "bond_style class2"),
        ("angle_style harmonic", "angle_style class2"),
        ("boundary m m m", "boundary p p p"),
        ("read_data Lig.lammps.lmp", "read_data Other.lmp"),
        ("tail no", "tail yes"),
    ],
)
def test_declared_conventions_are_verified(tmp_path, old, new):
    args = control(
        tmp_path / "source", change=lambda d, s, q, p: (d, s.replace(old, new), q, p)
    )
    r = intake(**args).report
    assert r["decision"] == "parameter_bundle_incomplete"
    assert any("convention_blocked" in b for b in r["blockers"])
    assert r["pair_expansion"] is None
    assert r["conventions"]["conventions"] == {"validated": False}


@pytest.mark.parametrize(
    "lj,coul,ok", [("0.5", "0.5", True), ("0.25", "0.5", False), ("0.5", "0.25", False)]
)
def test_separate_14_policies(tmp_path, lj, coul, ok):
    args = control(
        tmp_path / "source",
        change=lambda d, s, q, p: (
            d,
            s.replace("lj/coul 0.0 0.0 0.5", f"lj 0 0 {lj} coul 0 0 {coul}"),
            q,
            p,
        ),
    )
    r = intake(**args).report
    assert (r["decision"] == "audited_exact_graph") == ok


@pytest.mark.parametrize(
    "command",
    [
        "include other.in",
        "shell touch /tmp/no",
        "pair_modify mix arithmetic",
        "kspace_style pppm 1e-4",
        "run 1",
        "minimize 1e-4 1e-6 10 10",
        "special_bonds lj/coul 0 0 0.5 angle yes",
        "fix 1 all nve",
        "units ${u}",
    ],
)
def test_unknown_override_dynamic_commands_block_without_execution(tmp_path, command):
    args = control(
        tmp_path / "source",
        change=lambda d, s, q, p: (d, s + "\n" + command + "\n", q, p),
    )
    r = intake(**args).report
    assert r["decision"] == "parameter_bundle_incomplete" and r["blockers"]


@pytest.mark.parametrize("which", FILES + ("manifest.json",))
def test_all_source_hashes_rechecked(tmp_path, which):
    args = control(tmp_path / "source")
    out = tmp_path / "reports"
    e = emit(Path.cwd(), out, **args)
    with (args["source_dir"] / which).open("a") as f:
        f.write(" ")
    assert intake(**args).report["decision"] == "parameter_bundle_incomplete"
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(out, e["snapshot_sha256"], **args)


@pytest.mark.parametrize(
    "key", ["end_groups", "charge_state", "stereochemistry", "tacticity"]
)
def test_missing_exact_construction_metadata(tmp_path, key):
    def change(d, s, q, p):
        q["settings"].pop(key)
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert r["decision"] == "parameter_bundle_incomplete"


@pytest.mark.parametrize(
    "kind", ["identity", "atom_count", "mapping", "caps", "repeat"]
)
def test_graph_mutation(tmp_path, kind):
    def change(d, s, q, p):
        if kind == "identity":
            q["graph_identity"] = "0" * 64
        elif kind == "atom_count":
            q["atom_count"] += 1
        elif kind == "mapping":
            q["atom_mapping"][0] = q["atom_mapping"][1]
        elif kind == "caps":
            q["settings"]["end_groups"] = ["*C", "*C"]
        else:
            q["settings"]["repeat_count"] = 6
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert r["decision"] == "parameter_bundle_incomplete"


def test_improper_omissions_must_be_explicit(tmp_path):
    def change(d, s, q, p):
        p["improper_policy"]["omissions"].pop()
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert r["decision"] == "parameter_bundle_incomplete"
    assert any(t["family"] == "improper" for t in r["missing_terms"])


def test_parse_cvff_improper_coefficients_and_records(tmp_path):
    args = control(tmp_path / "source", approved=False)
    root = args["source_dir"]
    data = (root / "Lig.lammps.lmp").read_text()
    p = read(root / "provider_metadata.json")
    s = p["improper_policy"]["omissions"].pop()["sites"]
    data = data.replace("32 atoms", "32 atoms\n1 impropers\n1 improper types")
    data += (
        "\nImproper Coeffs\n\n1 1 -1 2\n\nImpropers\n\n1 1 "
        + " ".join(map(str, s))
        + "\n"
    )
    (root / "Lig.lammps.lmp").write_text(data)
    write(root / "provider_metadata.json", p)
    r = intake(**repin(root)).report
    assert r["decision"] == "provisional_parameterization", r["blockers"]
    assert not r["missing_terms"]


def test_missing_angle_remains_explicit(tmp_path):
    args = control(tmp_path / "source", approved=False)
    root = args["source_dir"]
    lines = (root / "Lig.lammps.lmp").read_text().splitlines()
    index = lines.index("Angles") + 2
    lines.pop(index)
    lines = [x.replace("60 angles", "59 angles") for x in lines]
    (root / "Lig.lammps.lmp").write_text("\n".join(lines) + "\n")
    r = intake(**repin(root)).report
    assert any(t["family"] == "angle" for t in r["missing_terms"])
    assert r["decision"] == "parameter_bundle_incomplete"


def test_missing_coefficient_references_rejected(tmp_path):
    args = control(tmp_path / "source", approved=False)
    root = args["source_dir"]
    data = (
        (root / "Lig.lammps.lmp")
        .read_text()
        .replace("Bonds\n\n1 1 ", "Bonds\n\n1 999 ")
    )
    (root / "Lig.lammps.lmp").write_text(data)
    assert "missing_coefficient_reference" in str(
        intake(**repin(root)).report["blockers"]
    )


def test_pair_expansion_shortest_distance_and_disconnected_guard(tmp_path):
    args = control(tmp_path / "source")
    r = intake(**args).report
    settings = r["conventions"]
    graph = "a" * 64
    e = expand_pairs(
        [1, 2, 3, 4, 5, 6],
        [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)],
        graph,
        settings,
    )
    assert len(e["pair14"]) == 3 and len(e["exclusions"]) == 12
    assert len({tuple(p["sites"]) for p in e["exclusions"]}) == 12
    with pytest.raises(ValueError, match="disconnected"):
        expand_pairs([1, 2, 3], [(1, 2)], graph, settings)
    with pytest.raises(ValueError, match="duplicate"):
        expand_pairs([1, 2], [(1, 2), (2, 1)], graph, settings)


@pytest.mark.parametrize("kind", ["extra", "missing", "symlink"])
def test_file_inventory(tmp_path, kind):
    args = control(tmp_path / "source")
    root = args["source_dir"]
    if kind == "extra":
        (root / "unlisted.txt").write_text("must not read")
    elif kind == "missing":
        (root / "independent_audit.json").unlink()
    else:
        (root / "link").symlink_to(root / "Lig.lammps.lmp")
    r = intake(**args).report
    assert (
        r["decision"] == "parameter_bundle_incomplete"
        and r["source_files_inspected"] == 0
    )


def test_missing_directory_no_fallback_or_inspection(tmp_path, monkeypatch):
    def forbidden(*a):
        raise AssertionError("source inspection")

    monkeypatch.setattr("island.forcefields.opls_providers.intake92.pinned", forbidden)
    r = intake(tmp_path / "absent").report
    assert (
        r["decision"] == "source_bundle_required" and r["source_files_inspected"] == 0
    )


def test_inside_worktree_rejected():
    r = intake(Path.cwd() / "bundle", manifest_sha256="a" * 64).report
    assert "licensed_files_must_remain_outside_Git_worktree" in r["blockers"]


def test_backend_scope_only_openmm_and_lammps(tmp_path):
    def change(d, s, q, p):
        q["backends"] = ["unsupported_engine"]
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert r["decision"] == "parameter_bundle_incomplete"


def test_deterministic_and_relocated_replay_requires_source_bytes(tmp_path):
    args = control(tmp_path / "source")
    root = tmp_path / "result"
    e = emit(Path.cwd(), root, **args)
    assert replay(root, e["snapshot_sha256"], **args) == e["result"]
    target = tmp_path / "relocated"
    shutil.copytree(args["source_dir"], target)
    shutil.rmtree(args["source_dir"])
    assert (
        replay(
            root,
            e["snapshot_sha256"],
            source_dir=target,
            manifest_sha256=args["manifest_sha256"],
        )
        == e["result"]
    )
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(root, e["snapshot_sha256"], **args)
    # Settings declarations are provenance; licensed coefficient arrays stay ephemeral.
    assert '"data"' not in "".join(p.read_text() for p in root.iterdir())
    assert all("tokens" not in c for c in e["result"]["record_claims"])


def test_strict_receipt_rejects_resigned_gate_mutation(tmp_path, monkeypatch):
    from scripts import validate_phase_4j92_receipt as validator

    root = tmp_path / "result"
    args = {"source_dir": tmp_path / "missing"}
    e = emit(Path.cwd(), root, **args)
    monkeypatch.setattr(
        validator, "verification", lambda repo, root: {"software_only": True}
    )
    monkeypatch.setattr(validator, "historical", lambda repo: {"control": "unchanged"})
    validator.seal(Path.cwd(), root, snapshot_sha=e["snapshot_sha256"], **args)
    assert (
        validator.validate(Path.cwd(), root, snapshot_sha=e["snapshot_sha256"], **args)[
            "status"
        ]
        == "passed"
    )
    r = read(root / "receipt.json")
    r["scope"]["runtime_registration"] = True
    write(
        root / "receipt.json", signed({k: v for k, v in r.items() if k != "identity"})
    )
    with pytest.raises(ValueError, match="scope_or_readiness"):
        validator.validate(Path.cwd(), root, snapshot_sha=e["snapshot_sha256"], **args)


def test_real_source_when_explicitly_available():
    root = Path(os.environ.get("ISLAND_OPLS_SOURCE_DIR", str(SOURCE_DIR)))
    if not root.is_dir():
        pytest.skip(
            "Explicit external J92 source directory absent; no real scientific coverage."
        )
    sha = os.environ.get("ISLAND_OPLS_MANIFEST_SHA256")
    if not sha:
        pytest.skip(
            "Independent manifest pin not supplied; no unpinned source inspection."
        )
    r = intake(root, manifest_sha256=sha).report
    assert not r["synthetic"]
    assert all(r[k] == v for k, v in READY.items())


@pytest.mark.parametrize(
    "field",
    [
        "source_binding",
        "graph_identity",
        "settings_identity",
        "expansion_hash",
        "claim_hashes",
        "improper_policy_identity",
        "independent",
        "reviewer",
    ],
)
def test_independent_audit_exact_binding(tmp_path, field):
    args = control(tmp_path / "source")
    root = args["source_dir"]
    audit = read(root / "independent_audit.json")
    audit[field] = (
        [] if field == "claim_hashes" else False if field == "independent" else ""
    )
    write(root / "independent_audit.json", audit)
    r = intake(**repin(root)).report
    assert r["decision"] == "parameter_bundle_incomplete"
    assert r["qualified_graph"] is None and r["accepted_coefficients"] == 0


@pytest.mark.parametrize(
    "kind",
    [
        "charge",
        "mass",
        "version",
        "authorization",
        "typing_origin",
        "charge_origin",
        "bond_graph",
        "improper_order",
        "non_PE_repeat",
    ],
)
def test_source_graph_authority_and_charge_guards(tmp_path, kind):
    def change(d, s, q, p):
        if kind == "charge":
            d = d.replace("Atoms\n\n1 1 1 0 ", "Atoms\n\n1 1 1 1 ")
        elif kind == "mass":
            d = d.replace("Masses\n\n1 12.011", "Masses\n\n1 1.008")
        elif kind == "version":
            p["version"] = "unmatched-version"
        elif kind == "authorization":
            p["license"]["authorization"] = "not_authorized"
        elif kind.endswith("origin"):
            p["atom_typing_origin" if kind == "typing_origin" else kind] = "unknown"
        elif kind == "bond_graph":
            start = d.index("Bonds\n\n")
            first = d[start:].splitlines()[2].split()
            first[-1] = "32"
            d = d[:start] + d[start:].replace(
                d[start:].splitlines()[2], " ".join(first), 1
            )
        elif kind == "improper_order":
            p["improper_policy"]["center_position"] = 0
        else:
            q["settings"]["psmiles"] = "*CCC*"
            graph = finite_graph(q["settings"])
            q["graph_identity"] = graph.identity
            q["atom_count"] = len(graph.atom_order)
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert r["decision"] == "parameter_bundle_incomplete", r
    assert r["blockers"] and r["accepted_coefficients"] == 0


@pytest.mark.parametrize(
    "section",
    ["Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Atoms"],
)
def test_incomplete_source_sections_block_without_fill(tmp_path, section):
    args = control(tmp_path / "source", approved=False)
    root = args["source_dir"]
    lines = (root / "Lig.lammps.lmp").read_text().splitlines()
    lines.pop(lines.index(section) + 2)
    (root / "Lig.lammps.lmp").write_text("\n".join(lines) + "\n")
    r = intake(**repin(root)).report
    assert r["decision"] == "parameter_bundle_incomplete"
    assert r["missing_terms"] and all(v["covered"] == 0 for v in r["coverage"].values())


@pytest.mark.parametrize(
    "name",
    [
        "request.json",
        "provider_metadata.json",
        "independent_audit.json",
        "manifest.json",
    ],
)
def test_malformed_metadata_structured_failure(tmp_path, name):
    args = control(tmp_path / "source")
    root = args["source_dir"]
    write(root / name, {"schema": []})
    if name == "manifest.json":
        args["manifest_sha256"] = pin(root / name)
    else:
        args = repin(root)
    r = intake(**args).report
    assert r["decision"] == "parameter_bundle_incomplete" and r["blockers"]


def test_source_style_must_precede_read_data(tmp_path):
    def change(d, s, q, p):
        s = s.replace("dihedral_style opls\n", "") + "dihedral_style opls\n"
        return d, s, q, p

    r = intake(**control(tmp_path / "source", change=change)).report
    assert "setting_after_read_data:dihedral_style" in str(r["blockers"])
