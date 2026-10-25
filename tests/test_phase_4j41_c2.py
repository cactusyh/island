"""C2 exact source-completion boundaries and independent negative controls."""

import copy
import os
import tarfile
from dataclasses import replace
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff.mixed_sip import (
    DISTRIBUTION_PINS,
    PCFFSiPDistribution,
    assess_pcff_mixed_sip,
    load_pcff_sip_distribution,
    prepare_pcff_mixed_sip,
)
from island.forcefields.pcff.sip_completion import (
    EXECUTABLE_PROFILES,
    SILANE_PROFILE,
    assess_pcff_sip_completion,
    load_pcff_sip_completion,
    prepare_pcff_sip_completion,
)
from island.forcefields.pcff.source import PCFFSource
from scripts.pcff_j41_c2_reference import source_audit


@pytest.fixture(scope="module")
def distribution():
    folder = Path(
        os.environ.get(
            "ISLAND_J41_DISTRIBUTION", "../lammps/lammps/tools/msi2lmp/frc_files"
        )
    )
    if folder.is_dir():
        return load_pcff_sip_distribution(folder)
    with tarfile.open("docs/evidence/phase_4j41/execution.tar.xz") as archive:
        raw = {n: archive.extractfile("study/" + n).read() for n in DISTRIBUTION_PINS}
    return PCFFSiPDistribution(
        PCFFSource(raw["pcff.frc"], DISTRIBUTION_PINS["pcff.frc"]),
        raw["pcff.rlb"],
        raw["pcff_templates.dat"],
    )


@pytest.fixture(scope="module")
def cases(distribution):
    out = {}
    for name, smiles, domain in (
        ("silane", "[SiH4]", "mixed_silane"),
        ("methyl", "C[SiH3]", "mixed_silane"),
        ("disiloxane", "[SiH3]O[SiH3]", "siloxane"),
        ("oligomer", "C[SiH2]C[SiH2]C", "organosilicon_polymer"),
        ("phosphate", "COP(=O)(OC)OC", "phosphorus_audit"),
    ):
        s = from_smiles(smiles, random_seed=2026)
        out[name] = (s, assess_pcff_sip_completion(s, distribution, domain=domain))
    return out


def test_only_verified_silane_has_execution_profile(cases, distribution):
    s, decision = cases["silane"]
    old = assess_pcff_mixed_sip(
        s, distribution, profile="island_pcff_mixed_silane_c1_v1"
    )
    assert decision.payload["c1_identity"] == old.identity
    assert EXECUTABLE_PROFILES == (SILANE_PROFILE,)
    assert (
        prepare_pcff_sip_completion(s, decision, profile=SILANE_PROFILE).identity
        == prepare_pcff_mixed_sip(s, old).identity
    )
    for name, (s, decision) in cases.items():
        if name != "silane":
            assert decision.payload["execution_profile"] is None
            with pytest.raises(PCFFError, match="Source-incomplete"):
                prepare_pcff_sip_completion(s, decision, profile=SILANE_PROFILE)


@pytest.mark.parametrize("name", ["methyl", "disiloxane", "oligomer"])
def test_exhaustive_source_tables_do_not_complete_absent_couplings(cases, name):
    system, decision = cases[name]
    audit = source_audit(system, decision)
    missing = [
        q
        for q in audit["queries"]
        if q["independent_resolution"]["status"] == "missing"
    ]
    assert missing and decision.payload["blocking_interactions"]
    for q in missing:
        searches = q["independent_resolution"]["searches"]
        assert searches and all(not s["candidate_ids"] for s in searches)
        for search in searches:
            table = audit["search_tables"][search["family"] + ":" + search["namespace"]]
            assert table["row_count"] == len(table["row_ids"])
            assert len(table["raw_table_sha256"]) == 64
    assert audit["charges"]["complete"]  # Neutrality does not supply parameters.


@pytest.mark.parametrize(
    "smiles,degree,oxo,h",
    [
        ("COP(=O)(OC)OC", 4, True, 0),
        ("CP(=O)(OC)OC", 4, True, 0),
        ("CP(C)(C)=O", 4, True, 0),
        ("P", 3, False, 3),
        ("CP(C)C", 3, False, 0),
        ("FP(F)(F)(F)F", 5, False, 0),
        ("[PH4+]", 4, False, 4),
        ("OP(=O)(O)O", 4, True, 0),
    ],
)
def test_p_candidates_require_chemical_evidence(distribution, smiles, degree, oxo, h):
    system = from_smiles(smiles, random_seed=2026)
    decision = assess_pcff_sip_completion(
        system, distribution, domain="phosphorus_audit"
    )
    audit = source_audit(system, decision)
    state = audit["phosphorus"]["sites"][0]
    assert state["degree"] == degree and state["oxo_template_matches"] == oxo
    assert state["selected_type"] is None
    p = decision.payload["phosphorus_decisions"][0]
    assert p["hydrogen_neighbors"] == h
    assert all(not c["authorized"] and c["rejection_reasons"] for c in p["candidates"])
    if oxo:
        assert state["oxo_template_connections_conflict"]
        assert p["oxygen_substitution"]


@pytest.mark.parametrize(
    "smiles",
    ["[SiH2]", "N[SiH3]", "S[SiH3]", "F[SiH3]", "Cl[SiH3]", "Br[SiH3]", "I[SiH3]"],
)
def test_unsupported_si_has_no_c2_execution(distribution, smiles):
    s = from_smiles(smiles, random_seed=2026)
    d = assess_pcff_sip_completion(s, distribution, domain="mixed_silane")
    assert d.payload["execution_profile"] is None and not d.payload["source_complete"]
    assert d.payload["c1_assessment"]["diagnostics"]


@pytest.mark.parametrize("file", list(DISTRIBUTION_PINS))
def test_wrong_source_hash_rejected(distribution, file):
    if file == "pcff.frc":
        d = replace(
            distribution,
            source=PCFFSource(distribution.source.raw + b"\n", DISTRIBUTION_PINS[file]),
        )
    else:
        attr = "rlb" if file == "pcff.rlb" else "templates"
        d = replace(distribution, **{attr: getattr(distribution, attr) + b"\n"})
    with pytest.raises(PCFFError, match="hash mismatch"):
        assess_pcff_sip_completion(
            from_smiles("[SiH4]", random_seed=2026), d, domain="mixed_silane"
        )


@pytest.mark.parametrize(
    "mutation",
    ["graph", "order", "charge", "coupling", "p_type", "source", "pressure", "policy"],
)
def test_resigned_mutations_do_not_bypass_reconstruction(cases, mutation, tmp_path):
    s, decision = cases["phosphate" if mutation == "p_type" else "methyl"]
    p = copy.deepcopy(decision.payload)
    if mutation == "graph":
        p["c1_assessment"]["graph"]["bonds"].pop()
    elif mutation == "order":
        p["c1_assessment"]["atom_order"].reverse()
    elif mutation == "charge":
        p["c1_assessment"]["inspection"]["native_charge_record"][
            "native_charge_record"
        ]["partial_charges"][1] += 0.1
    elif mutation == "coupling":
        p["blocking_interactions"].clear()
        p["execution_authorized"] = True
    elif mutation == "p_type":
        p["phosphorus_decisions"][0]["selected_type"] = "p="
    elif mutation == "source":
        p["source_evidence"]["atom_rows"][0]["connections"] = 3
    elif mutation == "pressure":
        p["pressure_status"] = "validated"
    else:
        p["policy"]["missing_source_row"] = "assume_zero"
    f = tmp_path / "decision.json"
    f.write_text(pack(p))
    with pytest.raises(PCFFError, match="Contradictory C2"):
        load_pcff_sip_completion(f, decision.distribution, system=s)


def test_wrong_final_graph_and_unknown_profile_rejected(cases):
    s, d = cases["silane"]
    with pytest.raises(PCFFError, match="Contradictory C2"):
        d.validate_integrity(cases["methyl"][0])
    with pytest.raises(PCFFError, match="Unqualified C2"):
        prepare_pcff_sip_completion(s, d, profile="island_pcff_phosphate_c2_v1")
