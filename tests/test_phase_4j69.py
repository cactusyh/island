"""Source audit controls; synthetic mutations never replace retained QM evidence."""

import copy
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, signed
from island.pcff_qm_pilot.hmds_local_baseline import derivative_check, fixed_families
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.hmds_source_audit import (
    source_graph,
    source_sum,
    validate_graph,
)
from scripts.audit_phase_4j69 import ROOT, inputs
from scripts.phase_4j69_test_policy import denied_command
from scripts.validate_phase_4j69_receipt import (
    DECLARATIONS,
    decision_semantics,
    semantics,
)

REPO = Path(__file__).resolve().parents[1]
E = REPO / ROOT


def resign(r):
    return signed({k: v for k, v in r.items() if k != "identity"})


@pytest.fixture(scope="module")
def retained():
    from island.pcff_qm_pilot.hmds_hamiltonian import load_model

    _, raw = inputs(REPO, E)
    model = load_model(REPO)
    frc = raw["tools/msi2lmp/frc_files/pcff.frc"]
    graph, selected = source_graph(model, frc)
    rows = read(REPO / "docs/evidence/phase_4j67/fitting-dataset.json")["observations"]
    return model, frc, graph, selected, rows


def test_source_graph_matches_retained_and_missing_are_explicit(retained):
    _, _, graph, _, _ = retained
    validate_graph(read(E / "source-dependency-graph.json"), graph)
    assert len(graph["missing_interactions"]) == 378
    assert graph["source_complete"] is False
    missing = [
        r for r in graph["nodes"] if r.get("kind") == "explicit_research_unknown"
    ]
    assert all(
        r["internal_values"] is None and r["LAMMPS_coefficients"] is None
        for r in missing
    )


@pytest.mark.parametrize("replacement", [b"306.1233", b"1.6563", b"157.0270"])
def test_source_row_SHA_mutation(retained, replacement):
    model, frc, *_ = retained
    before = (
        b"157.0260"
        if replacement.startswith(b"157")
        else b"1.6562"
        if replacement.startswith(b"1.")
        else b"306.1232"
    )
    assert before in frc
    with pytest.raises(CompletionError, match="SHA"):
        source_graph(model, frc.replace(before, replacement, 1))


@pytest.mark.parametrize(
    "field,value",
    [
        ("source_rows", ["quartic_bond:cff91:1776"]),
        ("source_units", ["nm", "kJ"]),
        ("execution_values", [1.6563, 306, -517, 673]),
        ("reference_order", "permute any sites"),
        ("equation", "0.5*K2*dr^2"),
        ("LAMMPS_coefficients", [1.6562, 1280, -2164, 2818]),
        ("dependencies", []),
    ],
)
def test_resigned_dependency_semantic_mutations(retained, field, value):
    graph = copy.deepcopy(retained[2])
    node = next(
        r
        for r in graph["nodes"]
        if r["id"]
        == ("bond-angle:2,5,6" if field == "dependencies" else "quartic_bond:2,5")
    )
    node[field] = value
    with pytest.raises(CompletionError, match="semantic"):
        validate_graph(resign(graph), retained[2])


@pytest.mark.parametrize("mutation", ["units", "phase"])
def test_actual_parser_conversion_mutation(retained, monkeypatch, mutation):
    from hashlib import sha256

    import island.forcefields.pcff.class2 as parser
    from island.forcefields.pcff.source import PCFFSource

    # Manufactured nonzero phases isolate the phase channel. They are never
    # assigned to HMDS or described as native QM/source evidence.
    catalog = parser.inspect_pcff_class2(
        PCFFSource(retained[1], sha256(retained[1]).hexdigest())
    )["records"]
    row = catalog["torsion_3:cff91:3133"]
    words = row["raw"].split()
    words[7], words[9], words[11] = "30.0", "180.0", "-90.0"
    manufactured = retained[1].replace(row["raw"].encode(), " ".join(words).encode(), 1)
    specimen = PCFFSource(manufactured, sha256(manufactured).hexdigest())
    phase_values = parser.inspect_pcff_class2(specimen)["records"][row["id"]][
        "normalized_values"
    ][1::2]
    np.testing.assert_allclose(phase_values, [np.pi / 6, np.pi, -np.pi / 2], atol=1e-15)

    original = parser.conversion

    def changed(unit):
        if mutation == "phase" and unit == "degree":
            return 1.0
        if mutation == "units" and unit.startswith("E"):
            return 1.0
        return original(unit)

    monkeypatch.setattr(parser, "conversion", changed)
    if mutation == "phase":
        wrong = parser.inspect_pcff_class2(specimen)["records"][row["id"]][
            "normalized_values"
        ][1::2]
        assert wrong == [30.0, 180.0, -90.0]
    with pytest.raises(CompletionError, match="values_units_order"):
        source_graph(retained[0], retained[1])


@pytest.mark.parametrize(
    "mutation", ["atom_order", "sites", "missing", "equilibrium", "prefactor"]
)
def test_stable_mapping_and_parameter_mutations(retained, mutation):
    model = copy.deepcopy(retained[0])
    if mutation == "atom_order":
        model["atom_order"][0:2] = [2, 1]
    else:
        node = next(
            r for r in model["interaction_inventory"] if r["id"] == "quartic_bond:2,5"
        )
        if mutation == "sites":
            node["sites"] = [2, 6]
        elif mutation == "missing":
            node["assignment"] = None
        elif mutation == "equilibrium":
            node["assignment"]["values"][0] += 0.001
        else:
            node["assignment"]["values"][1] *= 0.5
    model = resign(model)
    with pytest.raises(CompletionError):
        source_graph(model, retained[1])


@pytest.mark.parametrize("index", [0, 17, 18, 21, 22, 23, 24, 25])
def test_source_scalar_derivatives_and_family_closure(retained, index):
    model, _, _, selected, rows = retained
    xyz = np.array(rows[index]["coordinates"])
    energy, _ = source_sum(model, selected, xyz, jets=True)
    target = fixed_families(model, xyz)["families"]
    for family, e in energy.items():
        assert e.v == pytest.approx(target[family]["energy"], abs=1e-10)
        np.testing.assert_allclose(
            -e.d.reshape(27, 3), target[family]["forces"], atol=1e-10, rtol=0
        )
    # Independent directional difference of the scalar sum at all81 components
    # via a deterministic dense direction. Full component checks are in evidence.
    direction = np.cos(np.arange(81)).reshape(27, 3)
    h = 1e-5
    plus, _ = source_sum(model, selected, xyz + h * direction)
    minus, _ = source_sum(model, selected, xyz - h * direction)
    for k, e in energy.items():
        assert (plus[k] - minus[k]) / (2 * h) == pytest.approx(
            e.d @ direction.ravel(), abs=2e-5
        )


@pytest.mark.parametrize("index", [0, 17, 18, 21, 22, 23, 24, 25])
def test_equivalent_order_rotation_phase_and_prefactors(retained, index):
    model, _, _, selected, rows = retained
    xyz = np.array(rows[index]["coordinates"])
    original, _ = source_sum(model, selected, xyz, jets=True)
    reverse_model, reverse_values = copy.deepcopy(model), copy.deepcopy(selected)
    for row in reverse_model["interaction_inventory"]:
        f, s = row["family"], row["sites"]
        if f == "angle-angle":
            row["sites"] = [s[3], s[1], s[2], s[0]]
        elif f not in {"charge", "mass", "atom_type"}:
            row["sites"] = s[::-1]
            if f == "bond-angle":
                reverse_values[row["id"]] = reverse_values[row["id"]][::-1]
    reversed_sum, _ = source_sum(reverse_model, reverse_values, xyz, jets=True)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    rotated, _ = source_sum(
        model, selected, xyz @ rotation + [1.2, 0.2, -1.4], jets=True
    )
    for k, e in original.items():
        assert e.v == pytest.approx(reversed_sum[k].v, abs=1e-10)
        np.testing.assert_allclose(e.d, reversed_sum[k].d, atol=1e-10, rtol=0)
        assert rotated[k].v == pytest.approx(e.v, abs=1e-10)
        np.testing.assert_allclose(
            rotated[k].d.reshape(27, 3),
            e.d.reshape(27, 3) @ rotation,
            atol=1e-9,
            rtol=0,
        )
    # Plus/minus phase mapping is an algebraic identity including constants.
    phi = np.linspace(-np.pi, np.pi, 101)
    for n in (1, 2, 3):
        np.testing.assert_allclose(
            1 + np.cos(n * phi), 1 - np.cos(n * phi - np.pi), atol=2e-15
        )


@pytest.mark.parametrize("variant", ["SiO_stretch", "bridge_angle"])
def test_research_diagnostics_derivatives_not_global_equilibrium_change(
    retained, variant
):
    m, _, _, _, rows = retained
    result = derivative_check(m, rows[23]["coordinates"], variant)
    assert result["force_finite_difference_error"] < 2e-8
    assert (
        read(E / "equilibrium-phase-prefactor-audit.json")["J68_diagnostics"][variant][
            "propagates_global_source_equilibrium"
        ]
        is False
    )


@pytest.mark.parametrize("kind", ["duplicate", "null"])
def test_duplicate_and_null_diagnostic_columns_rejected(retained, kind):
    from island.pcff_qm_pilot.hmds_constrained import physical_problem
    from island.pcff_qm_pilot.hmds_local_baseline import span_sensitivity
    from scripts.analyze_phase_4j68 import reconstruct

    _, _, _, d, _, _ = reconstruct(REPO)
    p = physical_problem(d)
    a = copy.deepcopy(p)
    cols = np.zeros((len(p["y"]), 2)) if kind == "null" else p["w"][:, :2]
    a["w"] = np.column_stack((p["w"][:, :-1], cols, p["w"][:, -1]))
    assert span_sensitivity(p, a)["independent"] is False


def test_historical_bounded_and_free_failures_and_producers():
    r = read(E / "reconstruction.json")
    assert r["observations"] == 26 and r["protected"] == 22
    assert len(r["58_column_names"]) == 58 and r["one_energy_only_intercept"]
    models = r["J68_comparison"]
    assert models["bounded-Clarabel"]["rho"] == pytest.approx(1.196374166, abs=1e-8)
    assert models["bound-free-Clarabel"]["lower"] > 1
    assert all(m["lower"] > 1 for m in models.values())
    assert r["new_QM"] == r["MD_steps"] == 0 and r["reserved_validation_untouched"]


@pytest.mark.parametrize("field", list(DECLARATIONS))
def test_resigned_scientific_declarations_rejected(field):
    d = read(E / "decision.json")
    v = d[field]
    d[field] = (
        not v
        if type(v) is bool
        else v + 1
        if type(v) is int
        else ["forged"]
        if isinstance(v, list)
        else "forged"
    )
    with pytest.raises(CompletionError, match="declaration"):
        decision_semantics(resign(d))


def test_resigned_manifest_schema_rejected():
    # Synthetic envelope exercises schema semantics before the final receipt seal.
    n = {
        "offline": True,
        "observations": 26,
        "protected": 22,
        "physical_columns": 58,
        "energy_only_intercepts": 1,
        "source_available_subset_closed": True,
        "source_complete": False,
        "new_QM": 0,
        "MD_steps": 0,
        "validation_untouched": True,
        "hypotheses_executed": 0,
        "accepted_coefficients": 0,
    }
    from scripts.audit_phase_4j69 import BASE, REQUIRED

    receipt = signed(
        {
            **DECLARATIONS,
            "schema": "island.j69-receipt.v1",
            "base_commit": BASE,
            "required_J68": REQUIRED,
            "branch": "codex/phase-4j69-authoritative-sio-source-audit",
            "manifest_sha256": "a" * 64,
            "software_counts": {},
            "terminal_unique_tests": 1,
            "numerical_reconstruction": n,
            "decision_identity": "a" * 64,
            "historical_files_preserved": 1,
        }
    )
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {"a": "a"},
            "historical_files": {"a": "a"},
            "artifacts": {"a": "a"},
        }
    )
    semantics(receipt, manifest)
    forged = copy.deepcopy(receipt)
    forged["production_validated"] = True
    with pytest.raises(CompletionError, match="declaration"):
        semantics(resign(forged), manifest)
    forged = copy.deepcopy(receipt)
    forged["numerical_reconstruction"]["source_complete"] = True
    with pytest.raises(CompletionError, match="nested_scientific"):
        semantics(resign(forged), manifest)
    manifest["scientific_validation"] = "established"
    with pytest.raises(CompletionError, match="manifest_schema"):
        semantics(receipt, resign(manifest))


def test_native_nonactivation_and_execution_policy():
    d = read(E / "decision.json")
    assert d["validation_activated"] is False and d["validation_consumed"] is False
    assert d["hypotheses_executed"] == d["research_solver_calls"] == 0
    for command in (
        ["sqm", "-i", "x"],
        ["psi4", "input"],
        [sys.executable, "-m", "scripts.validate_phase_4j37"],
        [sys.executable, "-c", "integrator.step(1)"],
    ):
        assert denied_command(command)
    assert not denied_command([sys.executable, "scripts/replay_phase_4j69.py"])


def test_fresh_process_relocation_offline(tmp_path):
    # Copy bytes across filesystems into a physically independent tree.
    moved = tmp_path / "relocated"
    moved.mkdir()
    shutil.copytree(REPO / "docs/evidence", moved / "docs/evidence")
    for directory in ("src", "scripts"):
        shutil.copytree(
            REPO / directory,
            moved / directory,
            ignore=shutil.ignore_patterns("__pycache__"),
        )
    code = """
import importlib.abc,sys
from pathlib import Path
class NoNative(importlib.abc.MetaPathFinder):
 def find_spec(self, fullname, path=None, target=None):
  if fullname.split('.')[0] in {'psi4','clarabel','scs'}: raise AssertionError('forbidden execution dependency')
sys.meta_path.insert(0,NoNative())
from scripts.replay_phase_4j69 import replay
print(replay(Path(sys.argv[1]),Path(sys.argv[1])/'docs/evidence/phase_4j69'))
"""
    p = subprocess.run(
        [sys.executable, "-c", code, str(moved)],
        capture_output=True,
        text=True,
        check=False,
        env={
            **os.environ,
            "PYTHONPATH": str(moved) + ":" + str(moved / "src"),
            "OPENBLAS_NUM_THREADS": "1",
        },
        timeout=600,
        cwd=moved,
    )
    assert p.returncode == 0, p.stdout + p.stderr
    assert "'validation_untouched': True" in p.stdout
