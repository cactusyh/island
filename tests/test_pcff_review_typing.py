"""Synthetic source software regressions; real pinned-source checks are separate."""

from copy import deepcopy

import pytest
from test_pcff_automatic import synthetic  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield
from island.forcefields.pcff import (
    assign_pcff_source_types,
    expanded,
    source,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import POLICY
from island.forcefields.pcff.source import PCFFSource, digest


@pytest.fixture
def chemical_source(synthetic, monkeypatch, tmp_path):  # noqa: F811
    raw = synthetic.raw.replace(
        b"#equivalence cff91",
        b"1 1 o_1 15.999 O 1 synthetic\n1 1 nt 14.007 N 1 synthetic\n"
        b"1 1 ct 12.011 C 2 synthetic\n1 1 c_0 12.011 C 3 synthetic\n"
        b"1 1 cl 35.453 Cl 1 synthetic\n1 1 o* 15.999 O 2 synthetic\n"
        b"1 1 hw 1.008 H 1 synthetic\n#equivalence cff91",
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    path = tmp_path / "synthetic.frc"
    path.write_bytes(raw)
    return PCFFSource(raw, digest(raw)), path


@pytest.mark.parametrize("smiles,label", [("O=O", "o_1"), ("N#N", "nt")])
def test_elemental_on_rejected_by_typing_explicit_and_preparation(
    chemical_source, smiles, label
):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    src, path = chemical_source
    system = from_smiles(smiles, random_seed=2026)
    before = deepcopy(system)
    for version in (1, 2):
        profile = f"island_pcff_source_graph_v{version}"
        result = type_pcff_atoms(system, src, profile=profile)
        assert not result.complete
        assert not result.payload["assignments"]
        with pytest.raises(PCFFError, match="unresolved"):
            assign_pcff_source_types(
                system,
                src,
                dict.fromkeys(system.topology.sites, label),
                provenance="synthetic review bypass test",
                profile=profile,
            )
        with pytest.raises(PCFFError):
            prepare_forcefield(
                system,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        path,
                        (0, 0, 1),
                        (0, 0, 1),
                        typing_profile=profile,
                        resolution_policy=POLICY,
                    ),
                ),
            )
        # Rechecksummed contradictory evidence must not become valid by loading.
        data = unpack(result.json_text)
        data["assignments"] = dict.fromkeys(system.topology.sites, label)
        with pytest.raises(PCFFError):
            type(result)(pack(data), src).validate_integrity(system)
    assert system.topology.sites == before.topology.sites
    assert system.topology.bonds == before.topology.bonds


@pytest.mark.parametrize(
    "smiles,expected", [("C=O", "o_1"), ("CC#N", "nt"), ("O", "o*"), ("ClCl", "cl")]
)
def test_declared_positive_environments_unchanged(chemical_source, smiles, expected):
    pytest.importorskip("rdkit")
    from island.chemistry import from_smiles

    src, _ = chemical_source
    system = from_smiles(smiles, random_seed=2026)
    v2 = type_pcff_atoms(system, src, profile="island_pcff_source_graph_v2")
    assert v2.complete and expected in v2.assignments.values()
    explicit = assign_pcff_source_types(
        system,
        src,
        v2.assignments,
        provenance="synthetic software control",
        profile="island_pcff_source_graph_v2",
    )
    assert explicit.assignments == v2.assignments
    v1 = type_pcff_atoms(system, src, profile="island_pcff_source_graph_v1")
    if smiles == "ClCl":
        assert not v1.complete
    else:
        assert v1.assignments == v2.assignments
