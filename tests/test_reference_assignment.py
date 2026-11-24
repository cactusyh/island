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
    not os.environ.get("EMC_ROOT") and find_spec("pyemc") is None,
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
        (tmp_path / "chain.data").write_text(result.data)
        (tmp_path / "chain.params").write_text(result.parameters)
        run = subprocess.run(
            [lmp, "-log", "none"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            input="units real\natom_style full\nread_data chain.data\n"
            "include chain.params\nkspace_style pppm 0.001\n"
            'run 0\nprint "assignment_energy=$(pe)"\n',
            check=False,
        )
        assert run.returncode == 0, run.stdout + run.stderr
        assert isfinite(float(re.search(r"assignment_energy=(\S+)", run.stdout)[1]))


def test_emc_input_semantics():
    pytest.importorskip("rdkit")
    from island.forcefields.emc import _field_input

    assert _field_input("[*:2]OC[*:1]") == _field_input("*CO*")
    with pytest.raises(ValueError, match="stereochemistry"):
        _field_input("*/C=C/*")
    with pytest.raises(ValueError, match="whitespace"):
        _field_input("*CC*\nITEM END")
    with pytest.raises(ValueError, match="connected"):
        _field_input("*CC*.C")
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
