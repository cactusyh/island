"""C1 source completeness, phosphorus state evidence and immutable graph binding."""

import copy
import os
import tarfile
from dataclasses import replace
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff import assess_pcff_sip, prepare_pcff_sip
from island.forcefields.pcff.mixed_sip import (
    DISTRIBUTION_PINS,
    PCFFSiPDistribution,
    assess_pcff_mixed_sip,
    load_pcff_mixed_sip_assessment,
    load_pcff_sip_distribution,
    prepare_pcff_mixed_sip,
)
from island.forcefields.pcff.source import PCFFSource
from scripts.pcff_j41_c1_reference import mixed_inventory
from scripts.pcff_j41_reference import audit_queries


@pytest.fixture(scope="module")
def distribution():
    path = Path(
        os.environ.get(
            "ISLAND_J41_DISTRIBUTION", "../lammps/lammps/tools/msi2lmp/frc_files"
        )
    )
    if path.is_dir():
        return load_pcff_sip_distribution(path)
    archive = Path("docs/evidence/phase_4j41/execution.tar.xz")
    if not archive.is_file():
        pytest.skip("Pinned reviewed source distribution unavailable")
    with tarfile.open(archive) as tar:
        raw = {n: tar.extractfile("study/" + n).read() for n in DISTRIBUTION_PINS}
    return PCFFSiPDistribution(
        PCFFSource(raw["pcff.frc"], DISTRIBUTION_PINS["pcff.frc"]),
        raw["pcff.rlb"],
        raw["pcff_templates.dat"],
    )


@pytest.fixture(scope="module")
def cases(distribution):
    values = {}
    for name, smiles, profile in (
        ("silane", "[SiH4]", "mixed_silane"),
        ("methyl", "C[SiH3]", "mixed_silane"),
        ("disiloxane", "[SiH3]O[SiH3]", "siloxane"),
        ("oligomer", "C[SiH2]C[SiH2]C", "organosilicon_polymer"),
    ):
        system = from_smiles(smiles, random_seed=2026)
        values[name] = (
            system,
            assess_pcff_mixed_sip(
                system, distribution, profile="island_pcff_" + profile + "_c1_v1"
            ),
        )
    return values


def test_historical_silane_identity_unchanged(cases, distribution):
    system, a = cases["silane"]
    old = assess_pcff_sip(
        system, distribution.source, profile="island_pcff_silane_source_v1"
    )
    assert a.payload["complete"] and a.payload["execution_authorized"]
    assert (
        prepare_pcff_mixed_sip(system, a).identity
        == prepare_pcff_sip(system, old).identity
    )


@pytest.mark.parametrize(
    "name,count", [("methyl", 18), ("disiloxane", 54), ("oligomer", 99)]
)
def test_required_mixed_couplings_never_zero(cases, name, count):
    system, a = cases[name]
    p = a.payload
    assert not p["execution_authorized"]
    missing = [
        d
        for d in p["diagnostics"]
        if d["reason"] == "required_source_row_or_dependency_unresolved"
    ]
    assert len(missing) == count
    assert all(
        d["interaction"]["sites"] and d["interaction"]["supplied_types"]
        for d in missing
    )
    assert audit_queries(system, a, a.distribution.source.raw)["inventory_equal"]
    with pytest.raises(PCFFError, match="Incomplete or unqualified"):
        prepare_pcff_mixed_sip(system, a)


def test_wildcard_ambiguity_and_exact_row_policy(cases):
    _, a = cases["oligomer"]
    queries = a.payload["inspection"]["interaction_queries"]
    assert any(q["resolution"]["status"] == "ambiguous" for q in queries)
    assert any(
        d["reason"] == "exact_si_p_source_row_required"
        for d in a.payload["diagnostics"]
    )
    assert a.payload["policy"]["missing_cross_term"].startswith("reject")


@pytest.mark.parametrize(
    "smiles,kind,reason",
    [
        (
            "COP(=O)(OC)OC",
            "phosphate_ester_pentavalent",
            "p_oxo_template_connections_conflict",
        ),
        (
            "CP(=O)(OC)OC",
            "phosphonate_pentavalent",
            "p_oxo_template_connections_conflict",
        ),
        ("CP(C)(C)=O", "other_oxo_pentavalent", "p_oxo_template_connections_conflict"),
        ("P", "phosphine_like_trivalent", "phosphine_no_degree_three_source_type"),
        ("CP(C)C", "phosphine_like_trivalent", "phosphine_no_degree_three_source_type"),
        (
            "FP(F)(F)(F)F",
            "unsupported_connectivity_or_valence",
            "unsupported_phosphorus_oxidation_or_connectivity",
        ),
        (
            "[PH4+]",
            "unsupported_connectivity_or_valence",
            "charged_phosphorus_unsupported",
        ),
    ],
)
def test_p_chemistry_disambiguated_without_guessing_type(
    distribution, smiles, kind, reason
):
    s = from_smiles(smiles, random_seed=2026)
    a = assess_pcff_mixed_sip(
        s, distribution, profile="island_pcff_phosphorus_audit_c1_v1"
    )
    p = a.payload["phosphorus_environments"][0]
    assert p["classification"] == kind and p["reason"] == reason
    assert p["selected_type"] is None
    assert not p["execution_authorized"] and a.payload["typing"] is None
    assert "connections=5" in p["evidence"]["p="]
    assert "neighbor count" in p["evidence"]["interpretation"]
    with pytest.raises(PCFFError):
        prepare_pcff_mixed_sip(s, a)


@pytest.mark.parametrize(
    "smiles,reason",
    [
        ("[SiH2]", "unsupported_silicon_valence_or_identity"),
        ("N[SiH3]", "mixed_si_neighbor_type_not_independently_qualified"),
        ("S[SiH3]", "mixed_si_neighbor_type_not_independently_qualified"),
        ("F[SiH3]", "mixed_si_neighbor_type_not_independently_qualified"),
        ("Cl[SiH3]", "mixed_si_neighbor_type_not_independently_qualified"),
        ("[SiH3-]", "charged_isotopic_or_radical_unsupported"),
        ("C", "outside_selected_profile"),
    ],
)
def test_unsupported_si_and_neighbors(distribution, smiles, reason):
    s = from_smiles(smiles, random_seed=2026)
    a = assess_pcff_mixed_sip(s, distribution, profile="island_pcff_mixed_silane_c1_v1")
    assert any(d["reason"] == reason for d in a.payload["diagnostics"])
    assert not a.payload["complete"]


def test_all_source_environment_families_are_visible(distribution):
    inv = mixed_inventory(distribution.source.raw)
    assert set(inv["environments"]) == {
        "Si-" + n for n in ("C", "H", "O", "N", "S", "F", "Cl", "Br", "I", "Si", "P")
    }
    for env in inv["environments"].values():
        assert {
            "bond_increments",
            "angle-angle",
            "bond-bond_1_3",
            "angle-torsion_3",
            "end_bond-torsion_3",
            "middle_bond-torsion_3",
            "nonbond(9-6)",
            "auto_equivalence",
            "torsion-torsion",
        } <= set(env["families"])
        assert env["families"]["torsion-torsion"] == []
    assert inv["environments"]["Si-C"]["counts"]["bond_increments"] > 0
    assert inv["environments"]["Si-N"]["counts"]["bond_increments"] > 0


@pytest.mark.parametrize(
    "field",
    [
        "graph",
        "atom_order",
        "coordinates",
        "distribution_hashes",
        "policy",
        "charge_origin",
        "inspection",
        "charge_vector",
        "component_totals",
        "execution_authorized",
        "production_validated",
        "simulation_readiness",
    ],
)
def test_persistence_mutations_rejected(cases, field, tmp_path):
    system, a = cases["methyl"]
    p = copy.deepcopy(a.payload)
    if field == "atom_order":
        p[field].reverse()
    elif field == "inspection":
        p[field]["interaction_queries"].pop()
    elif field == "coordinates":
        p[field][0][1][0] += 0.1
    elif field == "charge_vector":
        p["inspection"]["native_charge_record"]["native_charge_record"][
            "partial_charges"
        ][1] += 0.01
    elif field == "component_totals":
        p["inspection"]["native_charge_record"]["native_charge_record"]["components"][
            0
        ]["charge"] = 1
    else:
        p[field] = "mutated"
    path = tmp_path / "assessment.json"
    path.write_text(pack(p))
    with pytest.raises(PCFFError, match="Contradictory"):
        load_pcff_mixed_sip_assessment(path, a.distribution, system=system)


@pytest.mark.parametrize("which", ["source", "rlb", "templates"])
def test_wrong_distribution_bytes_reject(distribution, which):
    if which == "source":
        from hashlib import sha256

        raw = distribution.source.raw + b"\n"
        changed = replace(distribution, source=PCFFSource(raw, sha256(raw).hexdigest()))
    else:
        changed = replace(distribution, **{which: getattr(distribution, which) + b"\n"})
    with pytest.raises(PCFFError):
        changed.validate_integrity()


def test_graph_connectivity_mutation_and_roundtrip(cases, tmp_path):
    s, a = cases["silane"]
    path = tmp_path / "assessment.json"
    path.write_text(a.json_text)
    assert (
        load_pcff_mixed_sip_assessment(path, a.distribution, system=s).identity
        == a.identity
    )
    changed = s.copy()
    changed.topology.remove_site(5)
    with pytest.raises(PCFFError):
        a.validate_integrity(changed)
    with pytest.raises(PCFFError, match="Unregistered"):
        assess_pcff_mixed_sip(s, a.distribution, profile="universal_sip")


def test_p_connectivity_and_evidence_mutation(distribution, tmp_path):
    s = from_smiles("COP(=O)(OC)OC", random_seed=2026)
    a = assess_pcff_mixed_sip(
        s, distribution, profile="island_pcff_phosphorus_audit_c1_v1"
    )
    assert a.payload["phosphorus_environments"][0]["formal_oxidation_number"] == 5
    p = a.payload
    p["phosphorus_environments"][0]["selected_type"] = "p="
    path = tmp_path / "record.json"
    path.write_text(pack(p))
    with pytest.raises(PCFFError, match="Contradictory"):
        load_pcff_mixed_sip_assessment(path, distribution, system=s)
