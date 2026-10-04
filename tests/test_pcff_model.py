"""Synthetic software records and analytical kernels, not source accuracy claims."""

from itertools import product
from math import cos, sin

import numpy as np
import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_class2 import assigned

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFModelSpecification,
    automatic,
    define_pcff_model,
    load_pcff_model,
    model,
    save_pcff_model,
    source,
    special_pair_policy,
)
from island.forcefields.pcff.class2 import FAMILIES
from island.forcefields.pcff.source import PCFFSource, digest
from island.forcefields.pcff.terms import (
    Class2Term,
    finite_difference_forces,
    pair_energy,
)


@pytest.fixture
def parameters(synthetic, monkeypatch):  # noqa: F811
    rows = []
    for family, (arity, _, units) in FAMILIES.items():
        if family == "bond-bond_1_3":
            continue
        rows.append("#" + family + " cff91")
        values = {
            "quartic_bond": [1.2, 5, -2, 1],
            "quartic_angle": [109, 5, -2, 1],
            "torsion_3": [1, 0, 2, 180, 3, 0],
            "nonbond(9-6)": [3.5, 0.1],
            "wilson_out_of_plane": [1, 0],
            "bond-bond": [0],
        }.get(family, [0.2] * len(units))
        for labels in product(("c", "h"), repeat=arity):
            rows.append("1 1 " + " ".join(labels) + " " + " ".join(map(str, values)))
    raw = synthetic.raw.replace(b"#end", ("\n".join(rows) + "\n#end").encode())
    sha = digest(raw)
    monkeypatch.setitem(source.PIN, "sha256", sha)
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", sha)
    monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
    return assigned(explicit([(0, 1)], [3, 3]), PCFFSource(raw, sha))


def policy(lj=(0, 0, 1), coulomb=(0, 0, 1)):
    return special_pair_policy(lj=lj, coulomb=coulomb)


def test_model_policy_zeros_and_historical_ownership(parameters, tmp_path):
    original = parameters.json_text
    result = define_pcff_model(parameters, special_pairs=policy())
    p = result.payload
    assert p["model_definition_complete"] and not p["raw_parameter_coverage_complete"]
    assert p["term_origins"]["policy_derived_zero"] == 9
    assert all(
        t["source_rows"] == [] and t["coefficients"] == [0]
        for t in p["terms"]
        if t["origin"] == "policy_derived_zero"
    )
    assert any(
        t["origin"] == "source_row" and t["coefficients"] == [0] for t in p["terms"]
    )
    assert len(result.numerical_terms()) == len(p["terms"])
    assert parameters.json_text == original
    assert p["numerical_verification"] == "not_attested_by_definition"
    assert not p["simulation_snapshot_available"]
    p["terms"].clear()
    assert result.payload["terms"]
    save_pcff_model(result, tmp_path / "model.json")
    assert (
        load_pcff_model(tmp_path / "model.json", parameters).identity == result.identity
    )
    with pytest.raises(PCFFError):
        save_pcff_model(result, tmp_path / "model.json")


def test_policy_identity_and_tampering(parameters):
    a = define_pcff_model(parameters, special_pairs=policy())
    b = define_pcff_model(parameters, special_pairs=policy(coulomb=(0, 0, 0.5)))
    assert a.identity != b.identity
    for field in (
        "compatibility_profile",
        "special_pair_inventory",
        "terms",
        "nonbonded",
    ):
        p = unpack(a.json_text)
        p[field] = []
        with pytest.raises(PCFFError):
            PCFFModelSpecification(pack(p), parameters).validate_integrity()
    p = unpack(a.json_text)
    p["model_definition_complete"] = False
    with pytest.raises(PCFFError):
        PCFFModelSpecification(pack(p), parameters).validate_integrity()
    p = unpack(a.json_text)
    p["terms"][0]["coefficients"][0] += 0.1
    with pytest.raises(PCFFError):
        PCFFModelSpecification(pack(p), parameters).validate_integrity()


def test_unrelated_missing_rows_block_definition(parameters, monkeypatch):
    raw = parameters.source.raw.replace(b"1 1 c h 1.2 5 -2 1\n", b"").replace(
        b"1 1 h c 1.2 5 -2 1\n", b""
    )
    sha = digest(raw)
    monkeypatch.setitem(source.PIN, "sha256", sha)
    monkeypatch.setitem(automatic.PROFILE, "source_sha256", sha)
    monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
    broken = assigned(explicit([(0, 1)], [3, 3]), PCFFSource(raw, sha))
    result = define_pcff_model(broken, special_pairs=policy())
    assert not result.payload["model_definition_complete"]
    assert result.payload["diagnostics"]
    with pytest.raises(PCFFError):
        result.numerical_terms()


@pytest.mark.parametrize(
    "values", [(True, 0, 1), (0, -1, 1), (0, float("nan"), 1), (0, 1), (0, 0, 2)]
)
def test_invalid_policy(values):
    with pytest.raises(PCFFError):
        special_pair_policy(lj=values, coulomb=(0, 0, 1))


def test_analytic_bond_force_and_torsion_convention():
    term = Class2Term("quartic_bond", (17, 93), (1.0, 2.0, -3.0, 4.0))
    xyz = {17: (0.0, 0.0, 0.0), 93: (1.4, 0.0, 0.0)}
    assert term.energy(xyz) == pytest.approx(2 * 0.4**2 - 3 * 0.4**3 + 4 * 0.4**4)
    f = finite_difference_forces(term, xyz, displacement=1e-5)
    assert f[17][0] == pytest.approx(4 * 0.4 - 9 * 0.4**2 + 16 * 0.4**3, abs=1e-8)
    phi = 0.7
    x = {2: (0, 1, 0), 5: (0, 0, 0), 19: (1, 0, 0), 38: (1, cos(phi), sin(phi))}
    p = (1.0, 0.2, 2.0, -0.3, 3.0, 0.5)
    tor = Class2Term("torsion_3", (2, 5, 19, 38), p)
    expected = sum(p[2 * n - 2] * (1 - cos(n * phi - p[2 * n - 1])) for n in (1, 2, 3))
    assert tor.energy(x) == pytest.approx(expected)
    assert Class2Term("torsion_3", (38, 19, 5, 2), p).energy(x) == pytest.approx(
        expected
    )


def test_asymmetric_reversal_id_remapping_and_pair_algebra():
    x = {
        11: (0.1, 1.3, 0.2),
        22: (0.0, 0.0, 0.0),
        44: (1.5, 0.2, 0.1),
        99: (2.0, 1.1, 1.4),
    }
    p = (1.0, -2.0, 0.3, 4.0, 0.5, -0.6)
    term = Class2Term("end_bond-torsion_3", (11, 22, 44, 99), p, (1.1, 1.4))
    rev = Class2Term("end_bond-torsion_3", (99, 44, 22, 11), p[3:] + p[:3], (1.4, 1.1))
    assert term.energy(x) == pytest.approx(rev.energy(x))
    moved = {-i: v for i, v in reversed(list(x.items()))}
    assert Class2Term(
        term.family, tuple(-i for i in term.sites), p, term.equilibria
    ).energy(moved) == term.energy(x)
    r = ((2**6 + 4**6) / 2) ** (1 / 6)
    eps = 2 * np.sqrt(0.2 * 0.8) * 2**3 * 4**3 / (2**6 + 4**6)
    result = pair_energy(
        r,
        rmin_i=2,
        epsilon_i=0.2,
        rmin_j=4,
        epsilon_j=0.8,
        charge_i=0.1,
        charge_j=-0.2,
        lj_weight=0.4,
        coulomb_weight=0.7,
    )
    assert result["lj_9_6"] == pytest.approx(-0.4 * eps)
    assert result["coulomb"] == pytest.approx(-0.7 * 332.06371 * 4.184 * 0.02 / r)
    zero = pair_energy(
        r,
        rmin_i=2,
        epsilon_i=0,
        rmin_j=4,
        epsilon_j=0.8,
        charge_i=0.1,
        charge_j=-0.2,
        lj_weight=0.4,
        coulomb_weight=0.7,
    )
    assert zero["lj_9_6"] == 0 and zero["coulomb"] == result["coulomb"]


def test_bad_kernel_inputs():
    with pytest.raises(PCFFError):
        Class2Term("wilson_out_of_plane", (1, 2, 3, 4), (1.0, 0.0))
    term = Class2Term("quartic_bond", (1, 2), (1.0, 2.0, 3.0, 4.0))
    for coords in (
        {True: (0, 0, 0), 2: (1, 0, 0)},
        {1: (True, 0, 0), 2: (1, 0, 0)},
        {1: (0, 0, 0), 2: (0, 0, 0)},
        {1: (1e308, 0, 0), 2: (-1e308, 0, 0)},
    ):
        with pytest.raises(PCFFError):
            term.energy(coords)
    with pytest.raises(PCFFError):
        term.energy({1: (0, 0, 0), 2: (1, 0, 0)}, coordinate_unit="nm")


def test_retained_actual_lammps_term_oracle():
    import json
    from pathlib import Path

    records = json.loads(
        (Path(__file__).parent / "fixtures/pcff_class2_terms.json").read_text()
    )
    for row in records["fixtures"]:
        term = Class2Term(
            row["family"],
            tuple(row["sites"]),
            tuple(row["coefficients"]),
            tuple(row["equilibria"]),
        )
        coords = dict(enumerate(row["xyz_angstrom"]))
        assert term.energy(coords) == pytest.approx(
            row["energy_kcal_mol"] * 4.184, abs=1e-10, rel=1e-12
        )
        forces = finite_difference_forces(term, coords, displacement=1e-4)
        assert np.allclose(
            [forces[i] for i in sorted(forces)],
            np.array(row["forces_kcal_mol_angstrom"]) * 4.184,
            atol=3e-5,
            rtol=1e-7,
        )


def test_offline_model_validation(parameters, tmp_path):
    import os
    import subprocess
    import sys

    (tmp_path / "source.frc").write_bytes(parameters.source.raw)
    (tmp_path / "assignment.json").write_text(parameters.json_text)
    m = define_pcff_model(parameters, special_pairs=policy())
    (tmp_path / "model.json").write_text(m.json_text)
    code = """
import sys
class Block:
    def find_spec(self,name,*args):
        if name.split('.')[0] in {'rdkit','openmm','foyer','parmed','scipy'}:
            raise RuntimeError('unexpected optional dependency: '+name)
sys.meta_path.insert(0,Block())
from pathlib import Path
from island.forcefields.pcff import source,automatic,model,load_pcff_parameters,load_pcff_model
root=Path(sys.argv[1]);raw=(root/'source.frc').read_bytes();sha=source.digest(raw)
source.PIN['sha256']=sha;automatic.PROFILE['source_sha256']=sha;model.PROFILE['frc_sha256']=sha
a=load_pcff_parameters(root/'assignment.json',source.PCFFSource(raw,sha))
m=load_pcff_model(root/'model.json',a)
assert m.payload['model_definition_complete']
assert m.numerical_terms()
"""
    completed = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path)],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
