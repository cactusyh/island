"""Fast, network-free AmberTools boundary tests; no fake GAFF claim."""

import copy
import json
import subprocess
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import pytest

from island.builders import build_linear_polymer
from island.exceptions import (
    AmberToolsInputError,
    AmberToolsStageError,
    AmberToolsUnavailableError,
)
from island.forcefields.ambertools import (
    AmberToolsOptions,
    AmberToolsParameterizationEngine,
)
from island.forcefields.ambertools import engine as backend
from island.forcefields.ambertools.lineage import (
    parse_and_validate_mol2,
    prepare_input,
)
from island.forcefields.ambertools.models import digest


@pytest.fixture(scope="module")
def polymer():
    pytest.importorskip("rdkit")
    return build_linear_polymer(
        "[*]CC[*]", dp=2, coordinate_method="local_templates",
        template_seed=2026, assembly_seed=2026,
    )


def test_optional_import_is_lazy():
    script = (
        "import sys; import island.forcefields.ambertools; "
        "assert 'parmed' not in sys.modules; assert 'rdkit' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", script], check=True)


def test_options_preflight_and_missing_tools(polymer, monkeypatch):
    with pytest.raises(AmberToolsInputError, match="force_field"):
        AmberToolsOptions("not_gaff", "provided", {})
    with pytest.raises(AmberToolsInputError, match="exact stable-site"):
        AmberToolsOptions("gaff2", "provided")
    with pytest.raises(AmberToolsInputError, match="AM1-BCC"):
        AmberToolsOptions("gaff", "am1bcc", {})
    charges = {site_id: 0.0 for site_id in polymer.topology.sites}
    source = polymer.to_dict()
    options = AmberToolsOptions("gaff", "provided", charges)
    charges[next(iter(charges))] = 2.0
    assert next(iter(options.provided_charges.values())) == 0.0
    copied_options = copy.deepcopy(options)
    assert dict(copied_options.provided_charges) == dict(options.provided_charges)
    with pytest.raises(AmberToolsInputError, match="coverage"):
        AmberToolsParameterizationEngine().parameterize(
            polymer, AmberToolsOptions("gaff", "provided", {1: 0.0})
        )
    monkeypatch.setattr(backend.shutil, "which", lambda _: None)
    with pytest.raises(AmberToolsUnavailableError, match="antechamber"):
        AmberToolsParameterizationEngine().parameterize(polymer, options)
    assert polymer.to_dict() == source


def test_input_scope_and_mol2_lineage(polymer, tmp_path):
    prepared = prepare_input(polymer, None)
    path = tmp_path / "input.mol2"
    path.write_text(prepared.mol2_text)
    output_index, _, _ = parse_and_validate_mol2(
        path, polymer, prepared.names, prepared.expected_cip, stage="input"
    )
    assert set(output_index.values()) == set(polymer.topology.sites)
    text = path.read_text()
    before_atoms, rest = text.split("@<TRIPOS>ATOM\n", 1)
    atoms_text, rest = rest.split("@<TRIPOS>BOND\n", 1)
    bonds_text, after_bonds = rest.split("@<TRIPOS>SUBSTRUCTURE", 1)
    atom_lines = atoms_text.splitlines()
    permutation = {int(line.split()[0]): len(atom_lines) - int(line.split()[0]) + 1
                   for line in atom_lines}
    new_atoms = []
    for line in reversed(atom_lines):
        fields = line.split()
        fields[0] = str(permutation[int(fields[0])])
        new_atoms.append(" ".join(fields))
    new_bonds = []
    for line in bonds_text.splitlines():
        fields = line.split()
        fields[1] = str(permutation[int(fields[1])])
        fields[2] = str(permutation[int(fields[2])])
        new_bonds.append(" ".join(fields))
    path.write_text(before_atoms + "@<TRIPOS>ATOM\n" + "\n".join(new_atoms)
                    + "\n@<TRIPOS>BOND\n" + "\n".join(new_bonds)
                    + "\n@<TRIPOS>SUBSTRUCTURE" + after_bonds)
    reordered, _, _ = parse_and_validate_mol2(
        path, polymer, prepared.names, prepared.expected_cip, stage="reordered"
    )
    assert reordered != output_index
    assert set(reordered.values()) == set(output_index.values())
    corrupted = path.read_text().replace(next(iter(prepared.names)), "BAD0", 1)
    path.write_text(corrupted)
    with pytest.raises(AmberToolsInputError, match="unknown atom name"):
        parse_and_validate_mol2(path, polymer, prepared.names, {}, stage="corrupt")
    no_h = polymer.copy()
    hydrogens = [i for i, site in no_h.topology.sites.items() if site.element == "H"]
    for hydrogen in hydrogens[:2]:
        no_h.topology.remove_site(hydrogen)
    from island.core import Coordinates

    no_h.coordinates = Coordinates({
        site_id: no_h.coordinates.get(site_id)
        for site_id in no_h.topology.sites
    })
    with pytest.raises(AmberToolsInputError, match="hydrogen"):
        prepare_input(no_h, None)
    planar = polymer.copy()
    planar.metadata["coordinate_source"] = "rdkit_2d"
    with pytest.raises(AmberToolsInputError, match="2D"):
        prepare_input(planar, None)


def test_output_bond_semantics_and_assigned_stereochemistry(tmp_path):
    from island.chemistry import from_smiles

    aromatic = from_smiles("Oc1ccccc1", random_seed=2026)
    input_record = prepare_input(aromatic, None)
    path = tmp_path / "phenol.mol2"
    path.write_text(input_record.mol2_text)
    assert parse_and_validate_mol2(
        path, aromatic, input_record.names, input_record.expected_cip,
        stage="input",
    )
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    converted = system_to_rdkit_graph(aromatic)
    kekule = Chem.Mol(converted.mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    by_index = {index + 1: site_id
                for index, site_id in enumerate(sorted(aromatic.topology.sites))}
    lines = path.read_text().splitlines()
    in_bonds = False
    for index, line in enumerate(lines):
        if line == "@<TRIPOS>BOND":
            in_bonds = True
            continue
        if line.startswith("@<TRIPOS>") and in_bonds:
            break
        if in_bonds:
            fields = line.split()
            if len(fields) < 4 or fields[3] != "ar":
                continue
            left, right = by_index[int(fields[1])], by_index[int(fields[2])]
            bond = kekule.GetBondBetweenAtoms(
                converted.site_id_to_rdkit_index[left],
                converted.site_id_to_rdkit_index[right],
            )
            fields[3] = str(int(bond.GetBondTypeAsDouble()))
            lines[index] = " ".join(fields)
    path.write_text("\n".join(lines) + "\n")
    assert parse_and_validate_mol2(
        path, aromatic, input_record.names, input_record.expected_cip,
        stage="equivalent Kekule output",
    )
    path.write_text(input_record.mol2_text)
    changed = path.read_text().replace(" 1\n", " 2\n", 1)
    path.write_text(changed)
    with pytest.raises(AmberToolsInputError, match="bond-order/aromaticity|bond semantics"):
        parse_and_validate_mol2(
            path, aromatic, input_record.names, input_record.expected_cip,
            stage="changed hydroxyl bond",
        )


def test_closed_shell_and_formal_charge_preflight():
    from island.chemistry import from_smiles

    ammonium = from_smiles("[NH4+]", random_seed=2026)
    assert prepare_input(ammonium, None).formal_charge == 1
    radical = from_smiles("[CH3]", random_seed=2026)
    with pytest.raises(AmberToolsInputError, match="even-electron|Open-shell"):
        prepare_input(radical, None)


def test_assigned_stereochemistry_survives_exchange(tmp_path):
    from island.chemistry import from_smiles

    chiral = from_smiles("F[C@H](Cl)Br", random_seed=2026)
    chiral_input = prepare_input(chiral, None)
    assert chiral_input.expected_cip
    chiral_site = next(iter(chiral_input.expected_cip))
    first, second = sorted(chiral.topology.neighbors(chiral_site))[:2]
    site_to_name = {site_id: name for name, site_id in chiral_input.names.items()}
    path = tmp_path / "chiral.mol2"
    path.write_text(chiral_input.mol2_text)
    lines = path.read_text().splitlines()
    entries = {}
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) > 5 and fields[1] in {site_to_name[first], site_to_name[second]}:
            entries[fields[1]] = index
    left, right = entries[site_to_name[first]], entries[site_to_name[second]]
    left_fields, right_fields = lines[left].split(), lines[right].split()
    left_fields[2:5], right_fields[2:5] = right_fields[2:5], left_fields[2:5]
    lines[left], lines[right] = " ".join(left_fields), " ".join(right_fields)
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(AmberToolsInputError, match="stereochemistry expected"):
        parse_and_validate_mol2(
            path, chiral, chiral_input.names, chiral_input.expected_cip,
            stage="mirrored stereocenter",
        )


def test_stage_exit_timeout_and_artifact_failures(tmp_path):
    with pytest.raises(AmberToolsStageError, match="exit code 3") as failed:
        backend._run_stage(
            "bad", [sys.executable, "-c", "import sys;sys.exit(3)"], tmp_path, 3,
        )
    assert failed.value.stage == "bad" and failed.value.returncode == 3
    with pytest.raises(AmberToolsStageError, match="timed out") as timed:
        backend._run_stage(
            "slow", [sys.executable, "-c", "import time;time.sleep(10)"],
            tmp_path, 0.1,
        )
    assert timed.value.artifact_dir == str(tmp_path)
    with pytest.raises(AmberToolsStageError, match="missing or empty"):
        backend._require_artifact("tleap", tmp_path / "missing.prmtop", tmp_path)
    with pytest.raises(AmberToolsStageError, match="unresolved"):
        backend._run_stage(
            "parmchk2", [sys.executable, "-c", "print('ATTN, need revision')"],
            tmp_path, 3,
        )


def test_force_field_data_discovery_is_family_specific(tmp_path, monkeypatch):
    cmd_dir = tmp_path / "dat/leap/cmd"
    parm_dir = tmp_path / "dat/leap/parm"
    cmd_dir.mkdir(parents=True)
    parm_dir.mkdir(parents=True)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name in ("antechamber", "parmchk2", "tleap"):
        (bin_dir / name).write_text("mock executable\n")
    for family in ("gaff", "gaff2"):
        (parm_dir / f"{family}.dat").write_text(f"{family} version-test header\n")
        (cmd_dir / f"leaprc.{family}").write_text(
            f"loadamberparams {family}.dat\n"
        )
    monkeypatch.setattr(backend.shutil, "which", lambda name: str(bin_dir / name))
    monkeypatch.setattr(backend, "_probe_version", lambda _: "unavailable (mock)")
    for family in ("gaff", "gaff2"):
        chain = backend._discover(AmberToolsOptions(
            family, "am1bcc", amberhome=tmp_path,
        ))
        assert chain.data_file.name == f"{family}.dat"
        assert chain.leaprc.name == f"leaprc.{family}"
        assert chain.versions["antechamber"] == "unavailable (mock)"
    (cmd_dir / "leaprc.gaff2").write_text("loadamberparams gaff.dat\n")
    with pytest.raises(AmberToolsUnavailableError, match="cannot be verified"):
        backend._discover(AmberToolsOptions(
            "gaff2", "am1bcc", amberhome=tmp_path,
        ))


@pytest.mark.parametrize("failure,stage", [
    ("sqm_missing", "antechamber"),
    ("sqm_failed", "antechamber"),
    ("sqm_incomplete", "antechamber"),
    ("frcmod_unresolved", "parmchk2"),
    ("typed_missing", "antechamber"),
    ("prmtop_missing", "tleap"),
    ("leap_error", "tleap"),
])
def test_mocked_qm_and_parameter_failures_retain_artifacts(
    polymer, tmp_path, monkeypatch, failure, stage,
):
    data = tmp_path / "gaff2.dat"
    leaprc = tmp_path / "leaprc.gaff2"
    data.write_text("GAFF2 mock\n")
    leaprc.write_text("loadamberparams gaff2.dat\n")
    monkeypatch.setattr(backend, "_discover", lambda _: backend._Toolchain(
        {name: name for name in ("antechamber", "parmchk2", "tleap")},
        tmp_path, leaprc, data, {},
    ))

    def fake_stage(current, command, directory, _timeout, *, environment=None):
        assert environment["AMBERHOME"] == str(tmp_path)
        if current == "antechamber" and failure != "typed_missing":
            (directory / "typed.mol2").write_text((directory / "input.mol2").read_text())
            if failure != "sqm_missing":
                (directory / "sqm.out").write_text(
                    "SCF failed to converge" if failure == "sqm_failed"
                    else "SCF converged but final result absent"
                    if failure == "sqm_incomplete"
                    else "Calculation Completed"
                )
        if current == "parmchk2":
            (directory / "typed.frcmod").write_text(
                "ATTN, need revision" if failure == "frcmod_unresolved" else "MASS\n"
            )
        if current == "tleap":
            if failure != "prmtop_missing":
                (directory / "result.prmtop").write_text("synthetic placeholder\n")
            (directory / "result.rst7").write_text("synthetic placeholder\n")
            (directory / "leap.log").write_text(
                "Errors = 1\n" if failure == "leap_error" else "Errors = 0\n"
            )
        return {"stage": current, "command": command}

    monkeypatch.setattr(backend, "_run_stage", fake_stage)
    original = polymer.to_dict()
    with pytest.raises(AmberToolsStageError) as failure_result:
        AmberToolsParameterizationEngine().parameterize(
            polymer, AmberToolsOptions(
                "gaff2", "am1bcc", work_root=tmp_path,
            ),
        )
    assert failure_result.value.stage == stage
    assert Path(failure_result.value.artifact_dir).is_dir()
    assert (Path(failure_result.value.artifact_dir) / "lineage.json").is_file()
    assert polymer.to_dict() == original


def test_provided_mode_rejects_unexpected_qm_run(polymer, tmp_path, monkeypatch):
    data, leaprc = tmp_path / "gaff.dat", tmp_path / "leaprc.gaff"
    data.write_text("GAFF mock\n")
    leaprc.write_text("loadamberparams gaff.dat\n")
    monkeypatch.setattr(backend, "_discover", lambda _: backend._Toolchain(
        {name: name for name in ("antechamber", "parmchk2", "tleap")},
        tmp_path, leaprc, data, {},
    ))

    def unexpected_qm(stage, _command, directory, _timeout, *, environment=None):
        assert environment["AMBERHOME"] == str(tmp_path)
        (directory / "typed.mol2").write_text((directory / "input.mol2").read_text())
        (directory / "sqm.out").write_text("Calculation Completed\n")
        return {"stage": stage}

    monkeypatch.setattr(backend, "_run_stage", unexpected_qm)
    charges = {site_id: 0.0 for site_id in polymer.topology.sites}
    with pytest.raises(AmberToolsStageError, match="unexpectedly ran SQM"):
        AmberToolsParameterizationEngine().parameterize(
            polymer, AmberToolsOptions(
                "gaff", "provided", charges, work_root=tmp_path,
            ),
        )


@dataclass(frozen=True)
class _FakeImported:
    provenance: dict
    result_signature: str = ""

    def content_signature(self):
        return digest(self.provenance)

    def validate_integrity(self, _system):
        assert self.result_signature == self.content_signature()


@pytest.mark.parametrize("family,mode", [
    ("gaff", "provided"), ("gaff2", "provided"), ("gaff2", "am1bcc"),
])
def test_mocked_workflow_commands_and_record(polymer, tmp_path, monkeypatch, family, mode):
    commands = []
    supplied = {site_id: 0.0 for site_id in polymer.topology.sites}
    first, second = sorted(supplied)[:2]
    if mode == "provided":
        supplied[first], supplied[second] = 0.2, -0.2
    data = tmp_path / f"{family}.dat"
    data.write_text(f"{family} mock version 1\n")
    leaprc = tmp_path / f"leaprc.{family}"
    leaprc.write_text(f"loadamberparams {family}.dat\n")
    executables = {}
    for name in ("antechamber", "parmchk2", "tleap"):
        path = tmp_path / name
        path.write_text("mock executable\n")
        executables[name] = str(path)
    chain = backend._Toolchain(
        executables,
        tmp_path, leaprc, data, {name: "mock tool" for name in (
            "antechamber", "parmchk2", "tleap",
        )},
    )
    monkeypatch.setattr(backend, "_discover", lambda _: chain)

    def fake_stage(stage, command, directory, _timeout, *, environment=None):
        assert environment["AMBERHOME"] == str(tmp_path)
        commands.append((stage, command))
        if stage == "antechamber":
            (directory / "typed.mol2").write_text((directory / "input.mol2").read_text())
            if mode == "am1bcc":
                (directory / "sqm.out").write_text("SCF converged\nCalculation Completed\n")
        elif stage == "parmchk2":
            (directory / "typed.frcmod").write_text("BOND\nca-ca 0.000 1.500\n")
        else:
            for name in ("result.prmtop", "result.rst7", "leap.log"):
                (directory / name).write_text("mock, not real Amber output\n")
        return {"stage": stage, "command": command}

    monkeypatch.setattr(backend, "_run_stage", fake_stage)
    monkeypatch.setattr(backend, "_prmtop_lineage", lambda _p, _r, system, names, _c, _d: (
        {index: site_id for index, site_id in enumerate(names.values())},
        supplied if mode == "provided" else
        {site_id: 0.0 for site_id in system.topology.sites},
    ))
    monkeypatch.setattr(backend, "import_amber_prmtop", lambda *_a, **_kw:
                        _FakeImported({"mocked": True}))
    opts = AmberToolsOptions(
        family, mode, supplied if mode == "provided" else None,
        work_root=tmp_path,
    )
    original = polymer.to_dict()
    result = AmberToolsParameterizationEngine().parameterize(polymer, opts)
    assert result.record["requested_force_field"] == family
    assert result.record["charge_method"] == mode
    if mode == "provided":
        assert result.record["charge_outcome"]["qm_run"] is False
        assert result.record["input_coordinates_angstrom"]
    assert result.record["tool_versions"]["antechamber"] == "mock tool"
    assert result.record["imported_result_signature"] == result.imported_result.result_signature
    assert commands[0][1][commands[0][1].index("-at") + 1] == family
    assert commands[1][1][commands[1][1].index("-s") + 1] == (
        "1" if family == "gaff" else "2"
    )
    assert ("-c", "rc" if mode == "provided" else "bcc") in list(zip(
        commands[0][1], commands[0][1][1:], strict=False
    ))
    assert polymer.to_dict() == original
    result.validate_integrity(polymer)
    assert result.record["artifact_dir"] is None
    assert not list(tmp_path.glob("island-ambertools-*"))


def test_mocked_mapping_corruption_and_failure_artifacts(polymer, tmp_path, monkeypatch):
    prepared = prepare_input(polymer, None)
    path = tmp_path / "typed.mol2"
    path.write_text(prepared.mol2_text.replace(next(iter(prepared.names)), "BAD0", 1))
    with pytest.raises(AmberToolsInputError):
        parse_and_validate_mol2(path, polymer, prepared.names, {}, stage="antechamber")
    stage = tmp_path / "job"
    stage.mkdir()
    with pytest.raises(AmberToolsStageError, match="per-site charges changed"):
        backend._check_charge_values({1: 0.2}, {1: 0.0}, "tleap", stage)


def test_preparation_record_detects_coordinate_changes(polymer, tmp_path):
    # The imported parameter signature is graph-based; the preparation record
    # independently binds the exact conformer sent to the external tools.
    from island.forcefields.ambertools.models import AmberToolsPreparationResult

    coordinates = {str(i): polymer.coordinates.get(i).tolist()
                   for i in sorted(polymer.topology.sites)}
    record = {"schema": "island_ambertools_preparation_v1",
              "input_coordinate_signature": digest(coordinates),
              "imported_result_signature": digest({"mocked": True})}
    imported = _FakeImported({"mocked": True}, digest({"mocked": True}))
    result = AmberToolsPreparationResult(imported, record, digest(record))
    result.validate_integrity(polymer)
    copy.deepcopy(result).validate_integrity(polymer)
    replace(result).validate_integrity(polymer)
    changed = polymer.copy()
    first = next(iter(changed.topology.sites))
    changed.coordinates.set(first, (5, 6, 7))
    from island.exceptions import InvalidAmberImportResultError

    with pytest.raises(InvalidAmberImportResultError, match="coordinates changed"):
        result.validate_integrity(changed)
    assert json.loads(json.dumps(dict(result.record))) == record
