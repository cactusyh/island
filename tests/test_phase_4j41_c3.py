"""C3 single-source audit, missing-coupling and phosphorus-scope controls."""

import copy
import os
import tarfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from island.charge_references.records import pack, unpack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff.sip_variants import (
    SiPSourceVariant,
    validate_sip_variant_inventory,
)
from scripts.pcff_j41_c2_reference import compiled_searches
from scripts.pcff_j41_c3_reference import (
    audit_inventory,
    compass_trial,
    trial_searches,
    validate_compass_reconstruction,
)
from scripts.pcff_j41_c3_silane import raw_coefficients, raw_energy
from scripts.phase_4j41_c3 import sources


@pytest.fixture(scope="module")
def variants():
    folder = Path(
        os.environ.get(
            "ISLAND_J41_DISTRIBUTION", "../lammps/lammps/tools/msi2lmp/frc_files"
        )
    )
    if folder.is_dir():
        return {v.filename: v for v in sources(folder)}
    # Published raw bytes remain portable; never fetch a force field implicitly.
    with tarfile.open("docs/evidence/phase_4j41_c3/execution.tar.xz") as archive:
        out = {}
        for filename in ("pcff.frc", "compass_published.frc"):
            root = "study/" + filename.removesuffix(".frc") + "/"
            p = unpack(archive.extractfile(root + "inventory.json").read().decode())
            out[filename] = SiPSourceVariant(
                filename,
                archive.extractfile(root + filename).read(),
                p["sha256"],
                pack(p["provenance"]),
                tuple(
                    (n, archive.extractfile(root + n).read(), h)
                    for n, h in p["companion_hashes"].items()
                ),
            )
        return out


@pytest.mark.parametrize("filename", ["pcff.frc", "compass_published.frc"])
def test_every_raw_row_and_coefficient_reconstructs(variants, filename):
    v = variants[filename]
    audit = audit_inventory(v)
    p = v.payload
    assert audit["raw_rows_checked"] == len(p["rows"])
    assert audit["numeric_rows_checked"] > 500
    assert audit["source_sha256"] == v.expected_sha256
    assert not p["execution_authorized"]
    assert p["production_validated"] is False
    assert p["simulation_readiness"] == "not_established"
    assert p["sip_families"]["nonbond(9-6)"]
    assert p["sip_families"]["bond_increments"]


def test_compass_schema_and_p_scope_remain_distinct(variants):
    pcff, compass = (
        variants["pcff.frc"].payload,
        variants["compass_published.frc"].payload,
    )
    assert pcff["sip_labels"] == ["p", "p=", "si", "sio", "sz"]
    assert compass["sip_labels"] == ["p4=", "si4", "si4c"]
    assert pcff["end_marker_present"] and not compass["end_marker_present"]
    assert pcff["companion_hashes"] and not compass["companion_hashes"]
    assert all(r["connections"] is None for r in compass["atom_types"])
    assert any(r["connections"] == 4 for r in pcff["atom_types"])
    assert "auto_equivalence" not in compass["family_counts"]
    assert compass["sip_families"]["bond-bond_1_3"] == []
    assert any(
        "phosphazenes" in b["raw"]
        for b in compass["templates"]["compass_published.frc"]
        if b["type"] == "p4="
    )


@pytest.mark.parametrize(
    "smiles", ["C[SiH3]", "[SiH3]O[SiH3]", "C[Si](C)(C)O[Si](C)(C)C"]
)
def test_other_source_does_not_supply_missing_couplings(variants, smiles):
    v = variants["compass_published.frc"]
    s = from_smiles(smiles, random_seed=2026)
    trial = compass_trial(s, v)
    assert not trial["source_complete"] and not trial["execution_authorized"]
    assert trial["execution_profile"] is None
    missing = [
        q for q in trial["interactions"] if q["resolution"]["status"] == "missing"
    ]
    assert missing and any(q["family"] == "bond-bond_1_3" for q in missing)
    for q in missing:
        assert q["sites"] and q["types"]
        assert all(
            not search["result"]["candidate_ids"]
            for search in q["resolution"]["searches"]
        )
    assert trial["charges"]["complete"] and not trial["charges"]["normalized"]


@pytest.mark.parametrize(
    "smiles", ["COP(=O)(OC)OC", "CP(=O)(OC)OC", "CP(C)C", "[PH4+]", "FP(F)(F)(F)F"]
)
def test_broad_p_template_cannot_authorize_other_chemical_states(variants, smiles):
    trial = compass_trial(
        from_smiles(smiles, random_seed=2026), variants["compass_published.frc"]
    )
    assert trial["mapping"]["status"] == "typing_blocked"
    assert trial["charges"] is None and not trial["interactions"]
    assert not trial["source_complete"]
    assert any(r["element"] == "P" for r in trial["mapping"]["blockers"])


def test_phosphazene_template_is_a_candidate_not_a_complete_profile(variants):
    s = from_smiles("[PH2]1=N[PH2]=N[PH2]=N1", random_seed=2026)
    trial = compass_trial(s, variants["compass_published.frc"])
    assert trial["mapping"]["status"] == "diagnostic_candidate_only"
    assert trial["charges"]["complete"] and trial["charges"]["total_e"] == "0.0000"
    assert not trial["source_complete"] and not trial["execution_authorized"]
    assert any(
        q["family"] == "bond-bond_1_3" and q["resolution"]["status"] == "missing"
        for q in trial["interactions"]
    )


@pytest.mark.parametrize(
    "mutation", ["hash", "bytes", "family", "authorization", "companion", "duplicate"]
)
def test_unauthorized_or_inconsistent_variant_rejected(variants, mutation):
    v = variants["pcff.frc"]
    if mutation == "hash":
        v = replace(v, expected_sha256="0" * 64)
    elif mutation == "bytes":
        v = replace(v, raw=v.raw + b"\n")
    elif mutation in ("family", "authorization"):
        p = unpack(v.provenance_json)
        p[mutation] = "COMPASS" if mutation == "family" else "unknown"
        v = replace(v, provenance_json=pack(p))
    elif mutation == "companion":
        n, raw, h = v.companions[0]
        v = replace(v, companions=((n, raw + b"\n", h),))
    else:
        v = replace(v, companions=v.companions * 2)
    with pytest.raises(PCFFError):
        v.validate_integrity()


@pytest.mark.parametrize(
    "mutation", ["row", "source", "family", "template", "zero", "execution"]
)
def test_reserialized_inventory_cannot_mix_sources(variants, tmp_path, mutation):
    v = variants["compass_published.frc"]
    p = copy.deepcopy(v.payload)
    if mutation == "row":
        p["rows"].append(variants["pcff.frc"].payload["rows"][0])
    elif mutation == "source":
        p["sha256"] = variants["pcff.frc"].expected_sha256
    elif mutation == "family":
        p["provenance"]["family"] = "PCFF"
    elif mutation == "template":
        p["templates"] = variants["pcff.frc"].payload["templates"]
    elif mutation == "zero":
        p["sip_families"]["bond-bond_1_3"].append("invented-zero")
    else:
        p["execution_authorized"] = True
    path = tmp_path / "inventory.json"
    path.write_text(pack(p))
    with pytest.raises(PCFFError, match="Contradictory"):
        validate_sip_variant_inventory(path, v)


@pytest.mark.parametrize(
    "mutation",
    [
        "source",
        "graph",
        "order",
        "coordinates",
        "charge",
        "coupling",
        "p_type",
        "authority",
    ],
)
def test_independent_reconstruction_rejects_mutated_evidence(variants, mutation):
    v = variants["compass_published.frc"]
    s = from_smiles("C[SiH3]", random_seed=2026)
    p = copy.deepcopy(compass_trial(s, v))
    if mutation == "source":
        p["source_sha256"] = variants["pcff.frc"].expected_sha256
    elif mutation in ("graph", "order"):
        p["system_sha256"] = "0" * 64
    elif mutation == "coordinates":
        s.coordinates.set(1, np.asarray(s.coordinates.get(1)) + [0.1, 0, 0])
    elif mutation == "charge":
        p["charges"]["vector_e"][1] = "0"
    elif mutation == "coupling":
        p["interactions"] = [
            q for q in p["interactions"] if q["resolution"]["status"] != "missing"
        ]
    elif mutation == "p_type":
        p["mapping"]["labels"][1] = "p4="
    else:
        p["execution_authorized"] = True
    with pytest.raises(ValueError, match="Contradictory"):
        validate_compass_reconstruction(s, v, p)


def test_silane_raw_energies_are_translation_invariant_and_source_bound(variants):
    s = from_smiles("[SiH4]", random_seed=2026)
    xyz = np.array([s.coordinates.get(i) for i in sorted(s.topology.sites)])
    xyz[1] += [0.02, -0.01, 0.03]
    for v in variants.values():
        coeff, evidence = raw_coefficients(v)
        energy = raw_energy(xyz, 0, coeff)
        np.testing.assert_allclose(
            raw_energy(xyz + [7, 3, -2], 0, coeff), energy, atol=1e-10
        )
        assert energy[0] > 0 and energy[-2:].tolist() == [0, 0]
        assert all(ids for ids in evidence.values())


def test_finite_silane_source_completeness_is_not_profile_authority(variants):
    trial = compass_trial(
        from_smiles("[SiH4]", random_seed=2026), variants["compass_published.frc"]
    )
    assert trial["source_complete"]
    assert not trial["execution_authorized"] and trial["execution_profile"] is None
    assert trial["charges"]["complete"]
    assert all(c["complete"] for c in trial["charges"]["components"])


def test_compiled_matcher_observes_compass_absence_and_equivalences(variants):
    matcher = Path(
        os.environ.get(
            "ISLAND_J41_C3_MATCHER", "/tmp/j41-c3-validation/matcher/leaf.so"
        )
    )
    if not matcher.is_file():
        pytest.skip("Optional compiled msi2lmp leaf matcher not supplied")
    v = variants["compass_published.frc"]
    t = compass_trial(from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026), v)
    result = compiled_searches(v.raw, matcher, trial_searches(t))
    assert result and any(q["compiled_row"] is None for q in result)
    assert any(q["compiled_row"] and q["family"] == "quartic_bond" for q in result)


def test_inventory_positive_round_trip_has_no_model_authority(variants, tmp_path):
    for v in variants.values():
        path = tmp_path / (v.filename + ".json")
        path.write_text(pack(v.payload))
        assert validate_sip_variant_inventory(path, v) == v.identity
