"""Synthetic software controls cannot qualify a PCFF-compatible PDMS extension."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_pdms.extension import FLAGS, PROFILE, digest, signed
from island.pcff_pdms.extension.fitting import (
    LIMITS,
    STAGES,
    UNCERTAINTY,
    check_bounds,
    split_check,
    synthetic_stage,
)
from island.pcff_pdms.extension.inventory import (
    chemical_identity,
    scope_check,
    validate_inventory,
)
from island.pcff_pdms.extension.qualification import (
    GATES,
    evaluate_gates,
    prepare_polymer_forcefield,
)
from island.pcff_pdms.extension.targets import validate_observation, validate_target
from scripts import phase_4j72_c1_fence as fence
from scripts import phase_4j78_extension as j78
from scripts import validate_phase_4j78_receipt as validator

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / j78.ROOT


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j78.scope()


@pytest.fixture(scope="module")
def graphs():
    return j78.inputs(ROOT)[0]


@pytest.fixture(scope="module")
def documents():
    return j78.reconstruct(ROOT)


@pytest.mark.parametrize("name", ("hmds", "dp2", "dp4", "dp6"))
def test_exact_declared_scope_accepted(graphs, name):
    s = scope_check(graphs[name])
    assert s.box is None
    assert not any(a.formal_charge for a in s.topology.sites.values())


@pytest.mark.parametrize(
    "mutation",
    (
        "SiH",
        "SiOH",
        "SiOC",
        "cyclic",
        "crosslinked",
        "charge",
        "radical",
        "triplet",
        "Br",
        "P",
        "Se",
        "periodic",
        "unsupported_length",
    ),
)
def test_unsupported_chemical_graphs_rejected(graphs, mutation):
    p = copy.deepcopy(graphs["hmds"])
    if mutation in ("SiH", "SiOH", "SiOC", "Br", "P", "Se"):
        index = next(
            i
            for i, a in enumerate(p["sites"])
            if a["element"] == ("O" if mutation in ("SiOH", "SiOC") else "C")
        )
        p["sites"][index]["element"] = {"SiH": "H", "SiOH": "H", "SiOC": "C"}.get(
            mutation, mutation
        )
    elif mutation == "charge":
        p["sites"][0]["formal_charge"] = 1
    elif mutation == "radical":
        p["sites"][0]["metadata"]["radical_electrons"] = 1
    elif mutation == "triplet":
        p["metadata"]["multiplicity"] = 3
    elif mutation == "periodic":
        p["box"] = {"invalid": "periodic_not_allowed"}
    elif mutation == "unsupported_length":
        p["sites"][0]["element"] = "Si"
    else:
        p["bonds"].append({"site1": 1, "site2": 6, "order": 1.0, "aromatic": False})
    with pytest.raises(ValueError):
        scope_check(p)


def test_graph_identity_ignores_coordinates_but_not_mapping(graphs):
    p = copy.deepcopy(graphs["hmds"])
    first = min(p["coordinates"])
    p["coordinates"][first][0] += 1
    assert chemical_identity(p) == chemical_identity(graphs["hmds"])
    p["sites"][0]["formal_charge"] = 1
    with pytest.raises(ValueError):
        chemical_identity(p)


def test_actual_cyclic_siloxane_with_valid_local_valences_is_rejected(graphs):
    p = copy.deepcopy(graphs["dp2"])
    s = scope_check(p)
    topo = s.topology
    end_si = [
        i
        for i, a in topo.sites.items()
        if a.element == "Si"
        and sum(topo.sites[n].element == "O" for n in topo.neighbors(i)) == 1
    ]
    removed = set()
    for si in end_si:
        c = next(n for n in topo.neighbors(si) if topo.sites[n].element == "C")
        removed.update({c, *(topo.neighbors(c) - {si})})
    new_id = max(topo.sites) + 1
    oxygen = copy.deepcopy(next(a for a in p["sites"] if a["element"] == "O"))
    oxygen.update(id=new_id, name="cycle_O")
    p["sites"] = [a for a in p["sites"] if a["id"] not in removed] + [oxygen]
    p["coordinates"] = {i: x for i, x in p["coordinates"].items() if i not in removed}
    p["coordinates"][new_id] = [0, 0, 0]
    p["bonds"] = [b for b in p["bonds"] if not {b["site1"], b["site2"]} & removed]
    p["bonds"].extend(
        {"site1": si, "site2": new_id, "order": 1.0, "aromatic": False} for si in end_si
    )
    with pytest.raises(ValueError, match="linear_tms_capped"):
        scope_check(p)


def test_unique_terms_instances_and_reversal_aliases(documents):
    inv = documents["unique-term-inventory.json"]
    assert len(inv["terms"]) == inv["unique_terms"]
    assert len(inv["instances"]) == inv["site_instances"] > inv["unique_terms"]
    assert len({t["id"] for t in inv["terms"]}) == len(inv["terms"])
    for t in inv["terms"]:
        assert not set(t["field_aliases"]) & set(t["fields"])
        assert t["origin"] == "unqualified_extension_target" and t["value"] is None
    assert any(t["field_aliases"] for t in inv["terms"])
    assert all(i["local_environment"] for i in inv["instances"])
    assert "Class-II" not in inv["structural_nonapplicability"]


def test_duplicate_term_and_missing_cross_term_block(documents):
    inv = documents["unique-term-inventory.json"]
    with pytest.raises(ValueError, match="duplicate_parameter_alias"):
        validate_inventory([*inv["terms"], inv["terms"][0]], inv["instances"])
    cross = next(t for t in inv["terms"] if t["family"] == "angle-angle")
    with pytest.raises(ValueError, match="missing_Class_II"):
        validate_inventory([t for t in inv["terms"] if t != cross], inv["instances"])


def test_unresolved_or_missing_equilibrium_dependency_blocks(documents):
    inv = documents["unique-term-inventory.json"]
    instances = copy.deepcopy(inv["instances"])
    i = next(i for i in instances if i["equilibrium_parameter_dependencies"])
    i["equilibrium_parameter_dependencies"][0] = "unidentified_cross_term_reference/r0"
    with pytest.raises(ValueError, match="missing_Class_II_equilibrium_dependency"):
        validate_inventory(inv["terms"], instances)


@pytest.mark.parametrize(
    "family",
    (
        "COMPASS",
        "OPLS-AA",
        "GAFF",
        "GAFF2",
        "external_PDMS",
        "CFF91",
        "analog",
        "native_PCFF",
        "zero_fill",
        "wildcard",
    ),
)
def test_foreign_or_default_values_never_become_extension_targets(documents, family):
    inv = documents["unique-term-inventory.json"]
    term = copy.deepcopy(inv["terms"][0])
    term["origin"] = family
    term["value"] = 1.0
    with pytest.raises(ValueError, match="foreign_family_or_analog"):
        validate_inventory([term, *inv["terms"][1:]], inv["instances"])


@pytest.mark.parametrize(
    "mutation",
    (
        "atom_order",
        "coordinate_order",
        "atom_mapping",
        "method",
        "basis_set",
        "dispersion_treatment",
        "charge",
        "multiplicity",
        "input_hash",
        "source_location",
    ),
)
def test_resigned_target_protocol_order_and_mapping_mutations_rejected(
    documents, mutation
):
    expected = documents["qm-target-manifest.json"]["targets"][0]
    t = copy.deepcopy(expected)
    if mutation in ("atom_order", "coordinate_order"):
        t[mutation][0], t[mutation][1] = t[mutation][1], t[mutation][0]
    elif mutation == "atom_mapping":
        t[mutation][0][1] = 999
    elif mutation in (
        "method",
        "basis_set",
        "dispersion_treatment",
        "charge",
        "multiplicity",
    ):
        t["qm_protocol"][mutation] = "mutated"
    else:
        t[mutation] = "fabricated"
    t = signed({k: v for k, v in t.items() if k != "identity"})
    with pytest.raises(ValueError, match="target_protocol_mapping"):
        validate_target(t, expected)


def test_target_templates_are_not_executed_or_independent_observations(documents):
    manifest = documents["qm-target-manifest.json"]
    assert (
        manifest["completed_QM_targets"] == 0 and manifest["execution_profile"] is None
    )
    assert all(
        t["output_hash"] is None and not t["physical_observation"]
        for t in manifest["targets"]
    )
    kinds = {t["specification"]["kind"] for t in manifest["targets"]}
    assert {
        "equilibrium_geometry",
        "stationary_point_gradient",
        "Hessian_or_projected_curvature",
        "coupled_two_dimensional_deformations",
        "constrained_bond_scan",
        "constrained_angle_scan",
        "torsion_or_methyl_rotation_scan",
    } <= kinds
    t = manifest["targets"][0]
    validate_target(t, t)
    with pytest.raises(ValueError, match="QM_software_version_missing"):
        validate_observation(t, t, expected_protocol=t["qm_protocol"])
    assert not documents["training-validation-split.json"]["independence_verified"]


def mathematical_fixture(*, a=None, validation_error=0, scale=1):
    a = (
        np.asarray(a if a is not None else [[1, 0], [0, 1], [1, 1], [1, -1]], float)
        * scale
    )
    values = np.array([0.3, -0.2])
    va = np.array([[1, 2], [2, -1]], float) * scale
    y = a @ values
    vy = va @ values + validation_error
    parameters = ["synthetic_a", "synthetic_b"]
    bounds = [{"parameter": p, "lower": -1, "upper": 1} for p in parameters]
    metadata = {
        "evidence_kind": "synthetic_contract_control",
        "parameter_order": parameters,
        "bounds_identity": digest(bounds),
        "design_identity": digest([a.tolist(), y.tolist(), va.tolist(), vy.tolist()]),
        "training": [
            {
                "id": "training",
                "geometry_identity": "synthetic_training",
                "conformer_group": "synthetic_T",
                "role": "training",
            }
        ],
        "validation": [
            {
                "id": "validation",
                "geometry_identity": "synthetic_validation",
                "conformer_group": "synthetic_V",
                "role": "validation",
            }
        ],
    }
    return (
        a,
        y,
        va,
        vy,
        {"parameters": parameters, "bounds": bounds, "metadata": metadata},
    )


def test_bounded_synthetic_fit_reuses_solver_and_never_qualifies():
    a, y, va, vy, kw = mathematical_fixture()
    first = synthetic_stage(a, y, va, vy, **kw)
    assert first["software_gates_passed"] and first["coefficients"] == pytest.approx(
        [0.3, -0.2]
    )
    assert first["solver_identity"]["backend"] == "scipy.optimize.lsq_linear"
    assert first["diagnostics"]["rank"] == 2 and first["raw_solver_state"]["success"]
    assert len(first["uncertainty"]["attempts"]) == UNCERTAINTY["replicates"]
    assert not first["scientific_qualification"] and first["accepted_coefficients"] == 0
    assert first == synthetic_stage(a, y, va, vy, **kw)


@pytest.mark.parametrize(
    "mutation",
    ("bound", "bound_order", "parameter_order", "design", "real_data", "uncertainty"),
)
def test_fitting_contract_mutation_rejected(mutation):
    a, y, va, vy, kw = mathematical_fixture()
    if mutation == "bound":
        kw["bounds"][0]["upper"] = 10
    elif mutation == "bound_order":
        kw["bounds"].reverse()
    elif mutation == "parameter_order":
        kw["parameters"].reverse()
    elif mutation == "design":
        a[0, 0] += 1
    elif mutation == "real_data":
        kw["metadata"]["evidence_kind"] = "QM"
    elif mutation == "uncertainty":
        kw["uncertainty"] = {**UNCERTAINTY, "noise_scale": 0}
    with pytest.raises(ValueError):
        synthetic_stage(a, y, va, vy, **kw)


@pytest.mark.parametrize(
    "bounds",
    (
        [{"parameter": "a", "lower": None, "upper": None}],
        [{"parameter": "a", "lower": -np.inf, "upper": 1}],
        [{"parameter": "a", "lower": 1, "upper": 1}],
    ),
)
def test_unconstrained_or_missing_bounds_rejected(bounds):
    with pytest.raises(ValueError):
        check_bounds(["a"], bounds)


@pytest.mark.parametrize("key", ("id", "geometry_identity", "conformer_group"))
def test_training_validation_leakage_rejected(key):
    *_, kw = mathematical_fixture()
    m = kw["metadata"]
    m["validation"][0][key] = m["training"][0][key]
    with pytest.raises(ValueError, match="leakage"):
        split_check(m["training"], m["validation"])


def test_rank_deficiency_and_ill_conditioning_block_before_solve():
    for a, blocker in (
        ([[1, 1], [2, 2], [3, 3], [4, 4]], "rank_deficient_design"),
        ([[1, 1], [2, 2 + 1e-7], [3, 3], [4, 4]], "ill_conditioned_design"),
    ):
        a, y, va, vy, kw = mathematical_fixture(a=a)
        r = synthetic_stage(a, y, va, vy, **kw)
        assert blocker in r["blockers"] and r["coefficients"] is None


def test_heldout_and_uncertainty_failures_do_not_accept_coefficients():
    for cfg, blocker in (
        ({"validation_error": 1}, "heldout_error_gate_failed"),
        ({"scale": 1e-6}, "uncertainty_gate_failed"),
    ):
        a, y, va, vy, kw = mathematical_fixture(**cfg)
        r = synthetic_stage(a, y, va, vy, **kw)
        assert blocker in r["blockers"] and not r["software_gates_passed"]
        assert r["accepted_coefficients"] == 0


def test_all_stages_require_frozen_predecessors():
    a, y, va, vy, kw = mathematical_fixture()
    predecessors = []
    for stage in STAGES:
        if stage != "primary":
            with pytest.raises(ValueError, match="stage_order"):
                synthetic_stage(a, y, va, vy, **kw, stage=stage)
        r = synthetic_stage(a, y, va, vy, **kw, stage=stage, predecessors=predecessors)
        assert r["frozen_predecessor_identities"] == [
            p["identity"] for p in predecessors
        ]
        predecessors.append(r)


def complete_synthetic_gates():
    ev = {k: True for k in GATES}
    ev.update(
        {
            "limits": copy.deepcopy(LIMITS),
            "evidence_kind": "synthetic_contract_control",
            "metrics": {k: 0.0 for k in GATES},
            "geometry_errors": {
                "equilibrium_bond_error_A": 0.0,
                "equilibrium_angle_error_rad": 0.0,
            },
            "rank": 2,
            "columns": 2,
        }
    )
    ev["metrics"]["conditioning"] = 1.0
    return ev


def test_passing_synthetic_gates_cannot_claim_scientific_qualification():
    r = evaluate_gates(complete_synthetic_gates())
    assert r["software_contract_passed"] and not r["scientific_qualification"]


@pytest.mark.parametrize("gate", GATES)
def test_every_acceptance_gate_fail_closed(gate):
    ev = complete_synthetic_gates()
    if gate == "equilibrium_geometry":
        ev["geometry_errors"] = {}
    elif gate in (
        "training_energy",
        "heldout_energy",
        "training_force",
        "heldout_force",
        "maximum_energy_error",
        "maximum_local_force_error",
        "stationarity",
        "torsion_scan",
        "curvature",
        "conditioning",
        "parameter_uncertainty",
    ):
        ev["metrics"][gate] = None
    else:
        ev[gate] = False
    r = evaluate_gates(ev)
    assert gate in r["failed_gates"] and not r["scientific_qualification"]


def test_threshold_widening_rejected():
    ev = complete_synthetic_gates()
    ev["limits"]["force_RMS_kcal_mol_A"] = 50
    with pytest.raises(ValueError, match="threshold_mutation"):
        evaluate_gates(ev)


def test_impossible_condition_number_rejected():
    ev = complete_synthetic_gates()
    ev["metrics"]["conditioning"] = 0.0
    assert "conditioning" in evaluate_gates(ev)["failed_gates"]


def test_real_candidate_remains_blocked_and_native_registry_unchanged(
    graphs, documents
):
    from island.forcefields.pcff.registry import list_pcff_profiles
    from island.pcff_pdms import REGISTERED_PROFILES

    assert REGISTERED_PROFILES == ()
    assert PROFILE not in {p["name"] for p in list_pcff_profiles()}
    q = documents["qualification-report.json"]
    b = documents["candidate-parameter-bundle.json"]
    assert not q["qualified"] and q["decision"] == "derived_profile_blocked"
    assert not q["gate_report"]["gates"]["Class_II_coverage"]
    assert all(b[k] == v for k, v in FLAGS.items())
    assert all(p["values"] is None for p in b["parameters"])
    with pytest.raises(ValueError, match="derived_profile_not_qualified"):
        prepare_polymer_forcefield(
            graphs["hmds"], forcefield="pcff", bundle=b, qualification=q
        )
    fake = {**q, "qualified": True}
    with pytest.raises(ValueError, match="derived_profile_not_qualified"):
        prepare_polymer_forcefield(
            graphs["hmds"], forcefield="pcff", bundle=b, qualification=fake
        )


def test_historical_native_bytes_preserved():
    assert validator.history(REPO, ROOT) > 9952


def test_replay_deterministic(documents):
    assert documents == j78.reconstruct(ROOT)
    assert j78.replay(ROOT)["decision"] == "derived_profile_blocked"


def test_resigned_qualification_report_rejected_by_reconstruction(tmp_path, documents):
    for name in (
        "graph-inputs.json",
        "historical-files.json",
        "implementation-reuse.json",
        "frozen-inputs.json",
        *documents,
    ):
        shutil.copyfile(ROOT / name, tmp_path / name)
    report = j78.read(tmp_path / "qualification-report.json")
    report["qualified"] = True
    report["decision"] = "derived_profile_qualified"
    report["accepted_coefficients"] = 365
    j78.write(
        tmp_path / "qualification-report.json",
        signed({k: v for k, v in report.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="resigned_extension_semantic_mutation"):
        j78.replay(tmp_path)


@pytest.mark.parametrize(
    "mutation", ("source_hash", "graph", "profile", "source_family")
)
def test_resigned_source_freeze_mutation_rejected(tmp_path, mutation):
    for n in (
        "graph-inputs.json",
        "historical-files.json",
        "implementation-reuse.json",
        "frozen-inputs.json",
    ):
        shutil.copyfile(ROOT / n, tmp_path / n)
    d = j78.read(tmp_path / "graph-inputs.json")
    if mutation == "source_hash":
        d["graphs"]["hmds"]["original_system_sha256"] = "0" * 64
    else:
        d[mutation] = "fabricated"
    j78.write(
        tmp_path / "graph-inputs.json",
        signed({k: v for k, v in d.items() if k != "identity"}),
    )
    freeze = j78.read(tmp_path / "frozen-inputs.json")
    freeze["capture_files"]["graph-inputs.json"] = j78.pin(
        tmp_path / "graph-inputs.json"
    )
    j78.write(
        tmp_path / "frozen-inputs.json",
        signed({k: v for k, v in freeze.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="fabricated_or_resigned_source_freeze"):
        j78.inputs(tmp_path)


def test_resigned_qualification_receipt_mutations_rejected(documents):
    receipt = signed(
        {
            "schema": "island.j78-receipt.v1",
            "base": j78.BASE,
            "branch": j78.BRANCH,
            "scope": j78.SCOPE,
            "decision": "derived_profile_blocked",
        }
    )
    manifest = signed(
        {
            "schema": "island.j78-manifest.v1",
            "base": j78.BASE,
            "branch": j78.BRANCH,
            "scope": j78.SCOPE,
            "freeze_sha256": j78.FREEZE_SHA,
        }
    )
    validator.semantic(receipt, manifest)
    for field in (
        "native_pcff_complete",
        "executable",
        "automatic_assignment_enabled",
        "accepted_coefficients",
        "production_validated",
    ):
        r = copy.deepcopy(receipt)
        r["scope"][field] = True
        r = signed({k: v for k, v in r.items() if k != "identity"})
        with pytest.raises(ValueError):
            validator.semantic(r, manifest)


def test_relocated_bundle_replay_without_sources(tmp_path, documents):
    from island.pcff_qm_pilot.bundle import pack, unpack

    reuse = j78.read(ROOT / "implementation-reuse.json")
    relocated = tmp_path / "offline"
    for n in {*reuse["replay_source_files"], *j78.IMPLEMENTATION}:
        dst = relocated / n
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / n, dst)
    for n in (
        "graph-inputs.json",
        "historical-files.json",
        "implementation-reuse.json",
        "frozen-inputs.json",
        *documents,
    ):
        dst = relocated / j78.ROOT / n
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / n, dst)
    packed = tmp_path / "relocatable.tar.xz"
    pack(relocated, packed)
    relocated_again = tmp_path / "reconstructed"
    unpack(packed, relocated_again)
    p = subprocess.run(
        [sys.executable, "-m", "scripts.phase_4j78_extension", "replay"],
        cwd=relocated_again,
        env={
            **os.environ,
            "PYTHONPATH": str(relocated_again)
            + os.pathsep
            + str(relocated_again / "src"),
        },
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(p.stdout)["decision"] == "derived_profile_blocked"
    assert not (relocated / "island_ffref").exists()


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "psi4",
        "lammps",
        "msi2lmp",
        "converter",
        "run_md",
        "evaluate_energy",
        "reserved_validation",
    ),
)
def test_forbidden_fake_child_cannot_execute(tmp_path, name):
    fake = tmp_path / name
    marker = tmp_path / "ran"
    fake.write_text(
        "#!/usr/bin/env python\nfrom pathlib import Path\nPath("
        + repr(str(marker))
        + ").write_text('synthetic')\n"
    )
    fake.chmod(0o755)
    code = (
        "from scripts.phase_4j78_extension import scope; scope(); import subprocess; from scripts.phase_4j72_c1_fence import ScopeDenied;\ntry: subprocess.run(["
        + repr(str(fake))
        + "])\nexcept ScopeDenied: print('rejected')\nelse: raise RuntimeError('fence bypass')"
    )
    p = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert "rejected" in p.stdout and not marker.exists()


def test_missing_fence_blocks_capture_before_read(monkeypatch, tmp_path):
    def missing():
        raise RuntimeError("fence unavailable")

    monkeypatch.setattr(fence, "assert_installed", missing)
    with pytest.raises(RuntimeError, match="fence unavailable"):
        j78.capture(tmp_path, tmp_path / "new")


@pytest.mark.parametrize(
    "control",
    (
        "real_QM",
        "MD",
        "LAMMPS_OpenMM",
        "converter",
        "scientific_fitting",
        "evaluator",
        "reserved_validation",
    ),
)
def test_real_scientific_controls_not_launched(control):
    pytest.skip(
        control + " not authorized by an execution profile; synthetic contracts only"
    )
