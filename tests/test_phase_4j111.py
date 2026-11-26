"""Source-bound profile changes; no fabricated increments or missing-term bypass."""

import pytest
from rdkit import Chem

from island.exceptions import ForceFieldError
from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.emc_graph import _field_elements
from island.forcefields.emc_profiles import install
from scripts.check_opls_polymer_profile import SOURCE_MODELS, source_atoms, typed_atoms
from scripts.validate_j107b_coverage import SOURCE_MODELS as V2_MODELS

ADDITIONAL_SOURCE_MODELS = [
    ("MECA", "CCOC(=O)NC"),
    ("PMCA", "COC(=O)Nc1ccccc1"),
    ("MEAL", "CNc1ccccc1"),
    ("DMAL", "CN(C)c1ccccc1"),
    ("BNAM", "NCc1ccccc1"),
    ("DMES", "CSC"),
    ("EMES", "CCSC"),
    ("MESH", "CS"),
    ("MSME", "CS(=O)(=O)C"),
]


@pytest.fixture(scope="module")
def fields(tmp_path_factory):
    root = tmp_path_factory.mktemp("j111-fields")
    for family, profile in [("pcff", "pcff-polymer-v2"), ("opls", "opls-polymer-v3")]:
        manifest = install(profile, root / family)
        numerical = "cff/pcff/pcff.frc" if family == "pcff" else "opls/2024/opls-aa.prm"
        assert manifest["base"][numerical] == manifest["installed"][numerical]
    return root


@pytest.mark.parametrize(
    "residue,smiles", SOURCE_MODELS + V2_MODELS + ADDITIONAL_SOURCE_MODELS
)
def test_source_types_and_charges(fields, residue, smiles):
    field = fields / "opls"
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    expected = source_atoms(
        (field / "opls/2024/opls-aa.define").read_text(),
        residue,
        mol,
        _field_elements(field, "opls/2024/opls-aa"),
    )
    actual = typed_atoms(assign_polymer_parameters(mol, "opls", field_dir=field))
    for atom, (kind, q) in expected.items():
        assert actual[atom][0] == kind
        assert actual[atom][1] == pytest.approx(q, abs=1e-6)


@pytest.mark.parametrize(
    "smiles,kind,charge",
    [
        ("O=CNc1ccccc1", "n", -0.4395),
        ("CCOC(=O)Nc1ccccc1", "n_2", -0.726),
        ("COC(=O)NCC", "n_2", -0.699),
    ],
)
def test_pcff_native_increment_sums(fields, smiles, kind, charge):
    mol = Chem.MolFromSmiles(smiles)
    atoms = typed_atoms(
        assign_polymer_parameters(mol, "pcff", field_dir=fields / "pcff")
    )
    n = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    assert atoms[n][0] == kind
    assert atoms[n][1] == pytest.approx(charge, abs=1e-6)


@pytest.mark.parametrize(
    "smiles", ["CC(=O)N(CC)CC", "CCSCC", "CCS(=O)(=O)CC", "CCOC1CCCCC1", "NC1CCC(N)CC1"]
)
def test_development_junctions(fields, smiles):
    assign_polymer_parameters(
        Chem.MolFromSmiles(smiles), "opls", field_dir=fields / "opls"
    )


@pytest.mark.parametrize("smiles", ["CCC#N", "CS(=O)C", "CP(=O)(C)C"])
def test_unsupported_chemistry_remains_failure(fields, smiles):
    with pytest.raises(ForceFieldError):
        assign_polymer_parameters(
            Chem.MolFromSmiles(smiles),
            "opls",
            field_dir=fields / "opls",
            backend_timeout=10,
        )
