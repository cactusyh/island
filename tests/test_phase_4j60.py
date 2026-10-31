"""J60 executable Hamiltonian, native data and identifiable joint design."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, signed
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.convergence_importer import native
from island.pcff_qm_pilot.hamiltonian import (
    accounting,
    cos_phi,
    evaluate,
    geometry,
    harmonic,
    load_model,
    parts,
)
from island.pcff_qm_pilot.hamiltonian_analysis import numerical
from island.pcff_qm_pilot.hamiltonian_fit import joint_fit
from scripts.independent_phase_4j60 import check

REPO = Path(__file__).resolve().parents[1]
E = REPO / "docs/evidence/phase_4j60"


def read(p):
    return json.loads(Path(p).read_text())


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def plan():
    return read(E / "plan/reference-plan.json")


@pytest.mark.parametrize("name", ["center", "angle-angle-0000", "bond-bond-0024"])
def test_independent_energy_and_all_force_columns(model, plan, name):
    result = check(model, plan["geometries"][name])
    assert result["finite_differences"][1]["fixed_force_max_error"] < 2e-5
    assert result["finite_differences"][1]["all_column_force_max_error"] < 2e-6


def test_rigid_reversal_and_periodicity_reconstruct():
    assert numerical(REPO, E) == read(E / "numerical-checks.json")


@pytest.mark.parametrize(
    "phi", [-np.pi - 1e-7, -np.pi + 1e-7, -1.1, 0, 0.7, np.pi, 2 * np.pi + 0.7]
)
def test_main_phase_and_torsion_derivative(phi):
    xyz = np.array(
        [[1, 0, 0], [0, 0, 0], [0, 0, 1], [np.cos(phi), np.sin(phi), 1]]
        + [[5 + i, 6, 7] for i in range(5)]
    )
    c = cos_phi(geometry(xyz), [1, 2, 3, 4])
    rev = cos_phi(geometry(xyz), [4, 3, 2, 1])
    assert c.v == pytest.approx(np.cos(phi), abs=1e-12)
    assert c.d == pytest.approx(rev.d, abs=1e-12)
    for n in (1, 2, 3):
        v = harmonic(c, n)
        tangent = np.array([-np.sin(phi), np.cos(phi), 0])
        assert v.d.reshape(9, 3)[3] @ tangent == pytest.approx(
            -n * np.sin(n * phi), abs=1e-11
        )
        delta = 0.23
        assert 1 + np.cos(n * phi - delta) == pytest.approx(
            1 - np.cos(n * phi - (delta + np.pi))
        )
        for alpha in (0, np.pi, 0.37):
            # Conditional identity only: retained source does not establish alpha.
            phase = delta - n * alpha + np.pi
            assert 1 + np.cos(n * (phi + alpha) - delta) == pytest.approx(
                1 - np.cos(n * phi - phase)
            )


def test_AAT_source_and_engine_never_equated():
    r = read(E / "hamiltonian-definition.json")
    aat = next(c for c in r["conventions"] if c["name"] == "angle_angle_torsion")
    assert "Phi - Phi1(0)" in aat["evidence"][0]["raw_text"]
    assert "\\cos (\\phi)" in aat["evidence"][1]["raw_text"]
    assert aat["authority"] == "declared_research_variant_conflicts_with_literal_FRC"
    assert 0.03 * (0.7 + 2 * np.pi) - 0.03 * 0.7 == pytest.approx(0.03 * 2 * np.pi)
    assert np.cos(0.7 + 2 * np.pi) == pytest.approx(np.cos(0.7))


def test_complete_accounting_does_not_supply_values(model, plan):
    table = accounting(model)
    assert table == read(E / "interaction-accounting.json")
    assert len(table["table"]) == len(model["table"])
    assert all(
        r["J60_status"]
        in (
            "verified_fixed",
            "explicit_fitted_unknown",
            "model_inapplicable",
            "structurally_inapplicable",
        )
        for r in table["table"]
    )
    assert table["numerical_parameter_completeness"] == "incomplete"
    assert model["author_confirmations"]["source_artifact_status"] == "not_yet_supplied"
    with pytest.raises(
        CompletionError, match="explicit_all_unknown_coefficients_required"
    ):
        evaluate(model, plan["geometries"]["center"])
    with pytest.raises(CompletionError):
        evaluate(model, plan["geometries"]["center"], {"one_partial_coefficient": 1})


def test_exclusions_scaling_and_no_pair_fit_weights(model, plan):
    r = parts(model, plan["geometries"]["center"])
    assert len(r["pairs"]) == 36
    assert sum(p["weight"] == 0 for p in r["pairs"]) == 21
    assert sum(p["weight"] == 1 for p in r["pairs"]) == 15
    assert all(p["weight"] == (p["graph_distance"] >= 3) for p in r["pairs"])
    assert len(r["columns"]) == 20
    assert all("Coulomb" not in c and "LJ" not in c for c in r["columns"])


@pytest.mark.parametrize(
    "field",
    [
        "family",
        "source_hashes",
        "charge_vector",
        "graph_identity",
        "atom_order",
        "energy_unit",
        "typing",
        "interaction_inventory",
    ],
)
def test_resigned_source_mapping_unit_mutation_rejected(model, plan, field):
    changed = copy.deepcopy(model)
    changed[field] = "changed"
    changed = signed({k: v for k, v in changed.items() if k != "identity"})
    with pytest.raises(CompletionError):
        parts(changed, plan["geometries"]["center"])


def test_design_is_complete_and_reduces_real_correlation():
    d = read(E / "joint-target-design.json")
    assert d["independent_coefficient_count"] == 20
    assert d["selected_complete_design"]["rank"] == 20
    assert d["selected_complete_design"]["condition"] < 6
    assert abs(d["selected_complete_design"]["end_middle_n3_cosine"]) < 0.01
    assert d["old_single_conformer_grid"]["end_middle_n3_cosine"] > 0.9999
    assert len(d["held_out"]) == 2
    assert not {x["geometry_identity"] for x in d["selected_training"]} & {
        x["geometry_identity"] for x in d["held_out"]
    }


@pytest.fixture(scope="module")
def manufactured(model, plan):
    """Independent nonzero software-control coefficients, never QM data."""
    d = read(E / "joint-target-design.json")
    names = parts(model, plan["geometries"]["center"])["columns"]
    values = {n: (i + 1) * 0.31 - 2.7 for i, n in enumerate(names)}
    out = []
    for sample in d["selected_training"] + d["held_out"]:
        r = evaluate(model, sample["coordinates"], values)
        out.append(
            {
                **sample,
                "graph_identity": model["graph_identity"],
                "atom_order": model["atom_order"],
                "conformer_identity": "training_manifold"
                if sample["split"] == "train"
                else sample["id"],
                "energy": r["energy"] + 105.0,
                "forces": r["forces"].tolist(),
                "reference_kind": "synthetic_software_control_only",
            }
        )
    return (
        values,
        [x for x in out if x["split"] == "train"],
        [x for x in out if x["split"] == "validation_only"],
    )


def test_identifiable_joint_solver_recovers_all_nonzero_control_values(
    model, plan, manufactured
):
    values, tr, va = manufactured
    r = joint_fit(
        model, tr, va, reference_gate_passed=True, criteria=plan["fit_criteria"]
    )
    assert r["coefficients"] == pytest.approx(values, abs=1e-9)
    assert r["untouched_validation"]["energy_RMSE"] < 1e-10
    assert r["untouched_validation"]["force_RMS"] < 1e-10
    assert r["scope"].startswith("experimental research Si-H")


@pytest.mark.parametrize(
    "failure", ["reference", "geometry_leak", "conformer_leak", "rigid_leak", "rank"]
)
def test_fit_rejects_unsafe_inference(model, plan, manufactured, failure):
    _, tr, va = copy.deepcopy(manufactured)
    passed = True
    if failure == "reference":
        passed = False
    if failure == "geometry_leak":
        va[0]["coordinates"] = tr[0]["coordinates"]
    if failure == "conformer_leak":
        va[0]["conformer_identity"] = tr[0]["conformer_identity"]
    if failure == "rigid_leak":
        va[0]["coordinates"] = (
            np.array(tr[0]["coordinates"])
            @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            + [1, 2, 3]
        ).tolist()
    if failure == "rank":
        tr = tr[:1]
    with pytest.raises(CompletionError):
        joint_fit(
            model, tr, va, reference_gate_passed=passed, criteria=plan["fit_criteria"]
        )


@pytest.fixture(scope="module")
def raw(tmp_path_factory):
    p = tmp_path_factory.mktemp("j60") / "native"
    unpack(E / "native-outputs.tar.xz", p)
    return p


@pytest.mark.parametrize(
    "key",
    [
        "ccq/center",
        "cc5/angle-angle-0000",
        "def2/bond-bond-0024",
        "tight/center",
        "method/center",
    ],
)
def test_real_native_QM_import_and_independent_units(raw, key):
    from island.pcff_qm.reference import CONVERSION
    from scripts.independent_phase_4j59 import native_values

    request = read(raw / "requests" / (key + ".json"))
    r, g = native(raw / "outputs" / key, request)
    energy, force = native_values(raw / "outputs" / key)
    assert energy == pytest.approx(
        r["energy_hartree"] * CONVERSION["energy_factor"], abs=2e-7
    )
    assert force == pytest.approx(
        -g * CONVERSION["energy_factor"] / CONVERSION["distance_factor"], abs=1e-10
    )


@pytest.mark.parametrize(
    "field", ["units", "elements", "charge", "scf_converged", "basis", "energy_hartree"]
)
def test_actual_native_mutations_rejected(raw, tmp_path, field):
    key = "ccq/center"
    dest = tmp_path / "output"
    shutil.copytree(raw / "outputs" / key, dest)
    r = read(dest / "result.json")
    r[field] = {
        "units": {},
        "elements": ["H"] * 9,
        "charge": 1,
        "scf_converged": False,
        "basis": "mutated",
        "energy_hartree": r["energy_hartree"] + 1,
    }[field]
    (dest / "result.json").write_text(json.dumps(r))
    with pytest.raises(CompletionError):
        native(dest, read(raw / "requests" / (key + ".json")))


def test_relocated_fresh_process_receipt(tmp_path):
    dest = tmp_path / "relocated"
    shutil.copytree(E, dest)
    p = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/validate_phase_4j60_receipt.py"),
            str(dest),
            "--repo",
            str(REPO),
        ],
        env={
            **os.environ,
            "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
            "OPENBLAS_NUM_THREADS": "1",
            "PYTHONHASHSEED": "460",
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=240,
    )
    assert p.returncode == 0, p.stdout + p.stderr


def test_actual_lammps_static_dihedral_controls():
    r = read(E / "backend-kernel-controls.json")
    assert r["static_kernel_check"] == "passed" and len(r["controls"]) == 4
    assert not r["physical_parameter_values_assigned"] and not r["PCFF_MD_executed"]
    assert max(c["force_max_error"] for c in r["controls"]) < 1e-9


@pytest.mark.parametrize(
    "mutation", ["source_replacement", "source_hash", "unit", "phase_rule"]
)
def test_source_and_contract_mutations_rejected(tmp_path, mutation):
    import lzma

    from scripts.validate_phase_4j60_receipt import source_checks

    root = tmp_path / "sources"
    root.mkdir()
    shutil.copytree(E / "sources", root / "sources")
    definition = read(E / "hamiltonian-definition.json")
    if mutation == "source_replacement":
        first = next(iter(definition["sources"].values()))
        (root / "sources" / (first["file"].replace("/", "__") + ".xz")).write_bytes(
            lzma.compress(b"replaced source")
        )
    elif mutation == "source_hash":
        next(iter(definition["sources"].values()))["sha256"] = "0" * 64
        definition = signed({k: v for k, v in definition.items() if k != "identity"})
    elif mutation == "unit":
        definition["conventions"][0]["units"] = "degree"
    else:
        next(c for c in definition["conventions"] if c["name"] == "main_torsion")[
            "interpreted_value"
        ] = "copy unchanged phase"
    (root / "hamiltonian-definition.json").write_text(json.dumps(definition))
    with pytest.raises(CompletionError):
        source_checks(root)


def test_held_out_conformers_are_not_rotated_or_permuted_training_points():
    r = read(E / "conformer-separation.json")
    assert r["automorphisms"] == 72
    assert (
        min(
            x["minimum_proper_rotation_and_graph_automorphism_RMS_A"]
            for x in r["held_out"]
        )
        > 0.02
    )
    assert r["held_out_mutual_minimum_RMS_A"] > 0.02


def test_real_special_pair_and_coulomb_unit_control():
    r = read(E / "backend-pair-control.json")
    assert r["special_bonds"] == [0, 0, 1]
    assert r["included_pairs"] == [[1, 4]]
    assert len(r["excluded_pairs"]) == 5
    assert r["energy_error"] < 1e-11 and r["force_max_error"] < 1e-10


def test_bootstrap_uses_only_training_and_has_finite_uncertainty(model, manufactured):
    from island.pcff_qm_pilot.hamiltonian_fit import bootstrap_training

    _, training, _ = manufactured
    r = bootstrap_training(model, training)
    assert r["valid_replicates"] >= 90 and r["held_out_used"] is False
    assert np.isfinite(np.array(r["correlation"])).all()
    assert max(r["sigma"]) < 1e-9


def test_partial_reference_gate_cannot_authorize_fit():
    from scripts.refine_phase_4j60 import reference_gates_pass

    records = [
        {
            "effect": k,
            "predeclared_gate_passed": True,
            "mixed_derivative_differences": {"angle-angle": {}, "bond-bond": {}},
        }
        for k in ("cardinal_QZ_5Z", "numerical_grid", "SCF_threshold")
    ]
    assert reference_gates_pass({"comparisons": records})
    assert not reference_gates_pass({"comparisons": records[:2]})
    del records[1]["mixed_derivative_differences"]["bond-bond"]
    assert not reference_gates_pass({"comparisons": records})


def test_complete_graph_synthetic_export_matches_native_energy_and_forces():
    from scripts.full_graph_phase_4j60 import reconstruct

    r = reconstruct(REPO, E)
    assert r["energy_error"] < 1e-10 and r["force_max_error"] < 1e-10
    report = read(E / "backend-full-graph-control.json")
    assert report["AA_only_carrier_count"] == 8
    assert len(report["synthetic_coefficients"]) == 20
    assert not report["source_physical_coefficients_fitted"]


def test_actual_full_rank_fit_is_rejected_without_promoting_source_rows():
    fit = read(E / "experimental-joint-fit.json")
    diagnostic = read(E / "model-discrepancy-diagnostics.json")
    decision = read(E / "fitting-decision.json")
    assert fit["rank"] == 20 and fit["condition"] < 6
    assert fit["training"]["force_RMS"] > 2.5
    assert fit["untouched_validation"]["force_RMS"] > 5.4
    assert not fit["held_out_thresholds_passed"]
    assert not fit["bootstrap_uncertainty_gate_passed"]
    assert diagnostic["active_bound_coefficients"] == []
    assert diagnostic["fraction_residual_SS_outside_all20_columns"] > 0.999999
    assert not diagnostic["held_out_used_for_probe_selection"]
    assert diagnostic["additional_coefficients_created"] == []
    assert decision["accepted_coefficients"] == []
    assert decision["parameter_completeness"] == "incomplete"
