"""Real provider comparisons for the PolyJarvis SI trial-build chemistries."""

import os
import re
import subprocess
from importlib.util import find_spec
from math import isfinite

import pytest

from island.exceptions import ForceFieldError
from island.forcefields.assignment import (
    assign_emc,
    assign_polymer_parameters,
    assign_radonpy,
)

# The repeat units used in PolyJarvis's archived *.esh inputs.
REPEATS = {
    "PE": "*CC*",
    "PEG": "*CCO*",
    "PLLA": "*[C@@H](C)C(=O)O*",
    "aPS": "*CC(*)c1ccccc1",
    "sPVC": "*C[C@H](Cl)C[C@@H](Cl)*",
    "PEEK": "*Oc1ccc(C(=O)c2ccc(Oc3ccc(*)cc3)cc2)cc1",
    "PSU": "*Oc1ccc(C(C)(C)c2ccc(Oc3ccc(S(=O)(=O)c4ccc(*)cc4)cc3)cc2)cc1",
    "PTFE": "*C(F)(F)C(*)(F)F",
}

PAPER_DP = {
    "PE": 179,
    "PEG": 114,
    "PLLA": 70,
    "aPS": 49,
    "sPVC": 41,
    "PEEK": 18,
    "PSU": 12,
    "PTFE": 50,
}

requires_emc = pytest.mark.skipif(
    find_spec("rdkit") is None
    or (not os.environ.get("EMC_ROOT") and find_spec("pyemc") is None),
    reason="requires an EMC installation",
)


@pytest.mark.parametrize("family", ("gaff", "gaff2"))
@pytest.mark.parametrize("name", REPEATS)
def test_radonpy_direct_comparison(family, name):
    pytest.importorskip("radonpy")
    from radonpy.ff.gaff import GAFF
    from radonpy.ff.gaff2 import GAFF2

    from island.builders import build_linear_polymer
    from island.chemistry import to_rdkit

    psmiles = REPEATS[name]
    dp = PAPER_DP[name] if os.environ.get("ISLAND_ASSIGNMENT_PAPER_DP") == "1" else 3
    system = build_linear_polymer(
        psmiles, dp=dp, generate_3d=False, head_end_group="*[H]", tail_end_group="*[H]"
    )
    graph = to_rdkit(system)
    ref = graph.mol
    engine = GAFF() if family == "gaff" else GAFF2()
    assert engine.ff_assign(ref, charge="gasteiger")
    result = assign_polymer_parameters(psmiles, family, dp=dp, charge="gasteiger")
    assert result.site_ids == tuple(graph.rdkit_index_to_site_id.values())
    assert (result.lj14_scale, result.coulomb14_scale) == (0.5, 5 / 6)
    for observed, reference in zip(result.mol.GetAtoms(), ref.GetAtoms(), strict=True):
        assert observed.GetProp("ff_type") == reference.GetProp("ff_type")
        for prop in ("ff_sigma", "ff_epsilon", "AtomicCharge"):
            assert observed.GetDoubleProp(prop) == reference.GetDoubleProp(prop)
    for observed, reference in zip(result.mol.GetBonds(), ref.GetBonds(), strict=True):
        for prop in ("ff_type", "ff_r0", "ff_k"):
            assert observed.GetProp(prop) == reference.GetProp(prop)
    for kind in ("angles", "dihedrals", "impropers"):
        assert [term.to_dict() for term in getattr(result.mol, kind)] == [
            term.to_dict() for term in getattr(ref, kind)
        ]


@requires_emc
@pytest.mark.parametrize(
    ("family", "name"),
    [("pcff", name) for name in REPEATS]
    + [("oplsaa", name) for name in ("PE", "aPS", "PTFE", "PSU")],
)
def test_emc_assignment(family, name, tmp_path):
    pytest.importorskip("rdkit")
    from island.builders import build_linear_polymer

    dp = PAPER_DP[name] if os.environ.get("ISLAND_ASSIGNMENT_PAPER_DP") == "1" else 3
    system = build_linear_polymer(
        REPEATS[name],
        dp=dp,
        generate_3d=False,
        head_end_group="*[H]",
        tail_end_group="*[H]",
    )
    result = assign_polymer_parameters(REPEATS[name], family, dp=dp)
    assert re.search(rf"(?m)^\s*{system.number_of_sites}\s+atoms$", result.data)
    assert re.search(rf"(?m)^\s*{system.number_of_bonds}\s+bonds$", result.data)
    atoms = result.data.split("\nAtoms\n")[1].lstrip().split("\n\n", 1)[0].splitlines()
    charges = [float(row.split()[3]) for row in atoms]
    assert len(charges) == system.number_of_sites
    assert sum(charges) == pytest.approx(0, abs=5e-4)
    assert "pair_coeff" in result.parameters
    if family == "pcff":
        for term in ("bb", "ba", "mbt", "ebt", "at", "aat", "bb13", "aa"):
            assert re.search(rf"_coeff\s+\d+ {term}\s", result.parameters)
    if lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS"):
        single_point(lmp, tmp_path, result.data, result.parameters)


def single_point(lmp, tmp_path, data, params):
    (tmp_path / "chain.data").write_text(data)
    (tmp_path / "chain.params").write_text(params)
    long_range = (
        "kspace_style pppm 0.001\n"
        if re.search(r"(?m)^pair_style\s+\S*/long\b", params)
        else ""
    )
    run = subprocess.run(
        [lmp, "-log", "none"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        input="units real\natom_style full\nread_data chain.data\n"
        f"include chain.params\n{long_range}"
        'run 0\nprint "assignment_energy=$(pe)"\n',
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    energy = float(re.search(r"assignment_energy=(\S+)", run.stdout)[1])
    assert isfinite(energy)
    return energy


def test_emc_input_semantics():
    pytest.importorskip("rdkit")
    from rdkit import Chem

    from island.forcefields.emc import _repeat

    def field_input(s):
        u, _ = _repeat(s, "strict")
        return Chem.MolToSmiles(
            u.mol, rootedAtAtom=u.head.dummy_atom_index, canonical=False
        )

    assert field_input("[*:2]OC[*:1]") == field_input("*CO*")
    with pytest.raises(ValueError, match="stereochemistry"):
        field_input("*/C=C/*")
    with pytest.raises(ValueError, match="whitespace"):
        field_input("*CC*\nITEM END")
    with pytest.raises(ValueError, match="connected"):
        field_input("*CC*.C")
    with pytest.raises(ValueError, match="positive integer"):
        assign_emc("*CC*", "pcff", dp=True)


@requires_emc
def test_emc_missing_parameters():
    pytest.importorskip("rdkit")
    # Stock OPLS/2024 lacks the alkyl/ether increment in this release.
    with pytest.raises(ForceFieldError, match=r"(?is)increment.*missing"):
        assign_emc(REPEATS["PEG"], "oplsaa", dp=3)


def test_radonpy_noncontiguous_site_ids():
    pytest.importorskip("radonpy")
    from island.chemistry import from_smiles
    from island.core import Coordinates

    system = from_smiles("C=O", site_ids=[90, 2, 541, 7])
    before = system.to_dict()
    result = assign_radonpy(system, "gaff2")
    assert [system.topology.sites[i].element for i in result.site_ids] == [
        a.GetSymbol() for a in result.mol.GetAtoms()
    ]
    assert all(not a.HasProp("AtomicCharge") for a in result.mol.GetAtoms())
    assert system.to_dict() == before
    system.coordinates = Coordinates()
    graph_only = assign_radonpy(system, "gaff2", charge="gasteiger")
    assert graph_only.mol.GetNumConformers() == 0
    assert graph_only.site_ids == result.site_ids
    with pytest.raises(ForceFieldError):
        assign_radonpy(from_smiles("CC", add_hydrogens=False), "gaff")


@pytest.mark.parametrize(
    ("repeat", "end_groups", "expected"),
    [
        ("*CC/C=C/CC*", {}, "CC/C=C/CCCC/C=C/CCCC/C=C/CC"),
        (r"*CC/C=C\CC*", {}, r"CC/C=C\CCCC/C=C\CCCC/C=C\CC"),
        ("*/C=C/*", {}, "C=C/C=C/C=C"),
        (r"*/C=C\*", {}, r"C=C/C=C\C=C"),
        ("[3H]CO[3H]", {}, "COCOCO"),
        ("[3H]CO*", {}, "COCOCO"),
        ("*CC*", {"head_end_group": "*/C=C/C"}, "C/C=C/CCCCCC"),
        (
            "*[C@H](*)Cl",
            {"head_end_group": "*F", "tail_end_group": "*Br"},
            "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)Br",
        ),
        ("[*:2]OC[*:1]", {"head_end_group": "*F", "tail_end_group": "*C"}, "FCOCOCOC"),
    ],
)
def test_radonpy_polymer_structure(repeat, end_groups, expected):
    pytest.importorskip("radonpy")
    from rdkit import Chem

    result = assign_polymer_parameters(repeat, "gaff2", dp=3, **end_groups)
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == Chem.MolToSmiles(
        Chem.MolFromSmiles(expected)
    )


def test_native_radonpy_input(tmp_path):
    pytest.importorskip("radonpy")
    from hashlib import sha256
    from pathlib import Path

    from radonpy.ff import gaff2_mod
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(r"F/C=C\C(F)(F)F"))
    assert AllChem.EmbedMolecule(mol, randomSeed=14) == 0
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetDoubleProp("AtomicCharge", i / 10)
    database = tmp_path / "custom.json"
    database.write_bytes(
        (Path(gaff2_mod.__file__).parent / "ff_dat/gaff2_mod.json").read_bytes() + b"\n"
    )
    ref = Chem.Mol(mol)
    assert gaff2_mod.GAFF2_mod(db_file=str(database)).ff_assign(ref)
    result = assign_radonpy(mol, "gaff2_mod", db_file=database)
    assert result.database_sha256 == sha256(database.read_bytes()).hexdigest()
    assert result.site_ids == tuple(range(mol.GetNumAtoms()))
    assert Chem.MolToSmiles(result.mol) == Chem.MolToSmiles(ref)
    assert (
        result.mol.GetConformer().GetPositions().tolist()
        == mol.GetConformer().GetPositions().tolist()
    )
    for i, atom in enumerate(result.mol.GetAtoms()):
        assert atom.GetProp("ff_type") == ref.GetAtomWithIdx(i).GetProp("ff_type")
        assert atom.GetDoubleProp("AtomicCharge") == i / 10
        assert not mol.GetAtomWithIdx(i).HasProp("ff_type")


@requires_emc
def test_emc_stereo_and_chain_settings():
    from rdkit import Chem

    from island.forcefields.emc import _repeat

    unit, omitted = _repeat("*C[C@H](Cl)CC/C=C/CC*", "polyjarvis")
    assert omitted and "@" in Chem.MolToSmiles(unit.mol)
    r = assign_emc(r"*/C=C\CC*", "pcff", dp=3, nchains=2, stereo_policy="polyjarvis")
    assert r.bond_stereo_omitted and r.nchains == 2
    assert not any(c in r.emc_smiles for c in "/\\")
    assert re.search(r"(?m)^\s*64\s+atoms$", r.data)  # two C12H20 chains
    ua = assign_emc("*CC*", "opls/2012/opls-ua", dp=2, nchains=2)
    assert ua.field == "opls/2012/opls-ua"
    assert re.search(r"(?m)^\s*12\s+atoms$", ua.data)  # two methyl-capped C6 chains


@requires_emc
@pytest.mark.parametrize(
    "repeat",
    [
        "*C(C(*)(c1c(c(c(c(c1[2H])[2H])[2H])[2H])[2H])[2H])([2H])[2H]",
        "*C([2H])([2H])C(*)(C(=O)OC([2H])([2H])[2H])C([2H])([2H])[2H]",
    ],
)
def test_emc_isotope_mass_and_energy(repeat, tmp_path, monkeypatch):
    from collections import Counter

    from island.forcefields import emc
    from island.forcefields._emc_isotopes import _section

    original = emc._build
    native = []

    def capture(*args, **kwargs):
        raw = original(*args, **kwargs)
        native.append(raw)
        return raw

    monkeypatch.setattr(emc, "_build", capture)
    result = assign_emc(repeat, "pcff", dp=2, nchains=2)
    types = {
        v.split()[0]
        for v in _section(result.data, "Masses").splitlines()
        if float(v.split()[1]) == pytest.approx(2.014101778, abs=1e-7)
    }
    rows = [r.split() for r in _section(result.data, "Atoms").splitlines()]
    assert Counter(r[1] for r in rows if r[2] in types) == {"1": 16, "2": 16}
    if lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS"):
        initial = single_point(lmp, tmp_path, *native[-1])
        assert single_point(
            lmp, tmp_path, result.data, result.parameters
        ) == pytest.approx(initial, rel=1e-10, abs=1e-7)


@requires_emc
def test_pcff_fused_aromatic_repeat():
    from rdkit import Chem

    repeat = (
        "*c1ccc(-c2nc3cc(-c4ccc5nc(*)c(-c6ccc7ccccc7c6)nc5c4)ccc3nc2-c2ccc3ccccc3c2)cc1"
    )
    result = assign_emc(repeat, "pcff", dp=3)
    assert Chem.MolToSmiles(Chem.MolFromSmiles(result.emc_smiles)) == Chem.MolToSmiles(
        Chem.MolFromSmiles(repeat)
    )
    assert re.search(r"(?m)^\s*212\s+atoms$", result.data)


@requires_emc
def test_polyjarvis_opls_installation(tmp_path):
    import shutil

    from island.forcefields.polyjarvis_fields import PATCHED, install_polyjarvis_opls

    if not (checkout := os.environ.get("ISLAND_POLYJARVIS_CHECKOUT")):
        pytest.skip("requires published PolyJarvis stock files and patches")
    field_dir = install_polyjarvis_opls(tmp_path / "field", checkout=checkout)
    assert install_polyjarvis_opls(field_dir, checkout=checkout) == field_dir
    result = assign_emc("*[Si](C)(C)O*", "opls", dp=3, field_dir=field_dir)
    assert result.field_sha256 == {"opls/2024/" + k: v for k, v in PATCHED.items()}
    assert re.search(r"(?m)^\s*32\s+atoms$", result.data)
    if lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS"):
        single_point(lmp, tmp_path, result.data, result.parameters)
    from pathlib import Path

    fixture = tmp_path / "checkout" / "emc_fields"
    shutil.copytree(Path(checkout) / "emc_fields", fixture)
    patch = fixture / "opls-2024-opls-aa.prm.patch"
    patch.write_bytes(patch.read_bytes() + b"\n")
    with pytest.raises(ForceFieldError, match="patch hash"):
        install_polyjarvis_opls(tmp_path / "broken", checkout=fixture.parent)
    assert not (tmp_path / "broken/opls/2024").exists()


@pytest.mark.parametrize("smiles", ["*c1ccc(-c2nc3cc4nc(*)oc4cc3o2)cc1", "*/C=C/*"])
def test_native_atom_order(smiles):
    pytest.importorskip("radonpy")
    from radonpy.core.poly import polymerize_MolFromSmiles
    from radonpy.ff.gaff2 import GAFF2
    from rdkit import Chem

    native = polymerize_MolFromSmiles(smiles, n=3, terminal="[H]")
    assert GAFF2().ff_assign(native, charge="gasteiger")
    result = assign_polymer_parameters(
        smiles, "gaff2", dp=3, atom_order="radonpy", charge="gasteiger"
    )
    standard = assign_polymer_parameters(smiles, "gaff2", dp=3)
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == Chem.MolToSmiles(
        Chem.RemoveHs(standard.mol)
    )
    by_id = {
        site: (
            atom.GetSymbol(),
            standard.mol.GetConformer().GetPositions()[atom.GetIdx()].tolist(),
        )
        for site, atom in zip(standard.site_ids, standard.mol.GetAtoms(), strict=True)
    }
    for i, atom in enumerate(result.mol.GetAtoms()):
        assert atom.GetProp("ff_type") == native.GetAtomWithIdx(i).GetProp("ff_type")
        assert atom.GetDoubleProp("AtomicCharge") == native.GetAtomWithIdx(
            i
        ).GetDoubleProp("AtomicCharge")
        assert by_id[result.site_ids[i]] == (
            atom.GetSymbol(),
            result.mol.GetConformer().GetPositions()[i].tolist(),
        )
    for kind in ("angles", "dihedrals", "impropers"):
        assert [term.to_dict() for term in getattr(result.mol, kind)] == [
            term.to_dict() for term in getattr(native, kind)
        ]


@requires_emc
@pytest.mark.parametrize(
    ("field", "smiles"),
    [
        ("opls/2012/opls-aa", "*C(C*)Cl"),
        ("opls/2012/opls-ua", "*C(=CCC*)CCCCCCC"),
        ("opls/2017/opls-aa", "*CC*"),
        ("opls/2017/opls-ua", "*CC*"),
        ("opls/2024/opls-ua", "*C1C(*)C2CC1C1C=CCC12"),
        ("pcff_ore", "*C(=CCC*)CCCCCCC"),
    ],
)
def test_emc_field_versions(field, smiles, tmp_path):
    r = assign_polymer_parameters(smiles, field, dp=4, nchains=2, density=0.9)
    assert r.field == field
    if lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS"):
        single_point(lmp, tmp_path, r.data, r.parameters)


@requires_emc
@pytest.mark.parametrize(
    ("field", "repeat"), [("pcff", "*CCO*"), ("opls/2012/opls-aa", "*C(C*)Cl")]
)
def test_polyjarvis_direct_lammps_parameters(field, repeat, tmp_path):
    import sys
    from pathlib import Path

    from island.forcefields.emc import _installation

    if not (checkout := os.environ.get("ISLAND_POLYJARVIS_CHECKOUT")):
        pytest.skip("requires a PolyJarvis checkout")
    builder = Path(checkout) / "mcp-servers/mcp-emc-server/smiles_to_emc.py"
    if not builder.is_file():
        pytest.skip("requires the PolyJarvis SMILES builder")
    root, _ = _installation(None)
    options = [
        "--field",
        field,
        "--dp",
        "3",
        "--nchains",
        "1",
        "--density",
        "0.1",
        "--seed",
        "2026",
    ]
    subprocess.run(
        [sys.executable, str(builder.resolve()), repeat, str(tmp_path), *options],
        check=True,
        env={**os.environ, "EMC_ROOT": str(root)},
    )
    result = assign_emc(repeat, field, dp=3, seed=2026)

    def coeffs(s):
        return [
            line.split("#")[0].split()
            for line in s.splitlines()
            if re.match(r"^\w+_coeff\s", line)
        ]

    assert coeffs(result.parameters) == coeffs(
        (tmp_path / "emc_build.params").read_text()
    )

    def entries(content, section):
        rows = (
            content.partition("\n" + section + "\n\n")[2]
            .split("\n\n", 1)[0]
            .splitlines()
        )
        values = [line.split("#")[0].split() for line in rows]
        return [v[:4] for v in values] if section == "Atoms" else values

    reference = (tmp_path / "emc_build.data").read_text()
    for label in ("Atoms", "Masses", "Bonds", "Angles", "Dihedrals", "Impropers"):
        assert entries(result.data, label) == entries(reference, label)


@pytest.mark.parametrize(
    ("repeat", "dp", "expected"), [("*CCO*", 1, "CCO"), ("*[CH2:3]O*", 3, "COCOCO")]
)
def test_native_order_capped_monomer_and_mapped_atoms(repeat, dp, expected):
    pytest.importorskip("radonpy")
    from rdkit import Chem

    result = assign_polymer_parameters(repeat, "gaff2", dp=dp, atom_order="radonpy")
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == expected


@requires_emc
def test_emc_macrocycle_depth(tmp_path):
    result = assign_emc(
        "*CC(*)C(=O)OC1CCCCCCCCCCC1", "pcff", dp=4, nchains=2, density=0.9
    )
    assert result.automatic_ring_depth
    if lmp := os.environ.get("ISLAND_ASSIGNMENT_LAMMPS"):
        single_point(lmp, tmp_path, result.data, result.parameters)
