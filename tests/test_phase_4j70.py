"""Synthetic intake controls are separate from actual source/QM evidence."""

import copy
import os
import shutil
import subprocess
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, canonical, signed
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.source_intake import (
    SourceRow,
    completeness,
    converter_gate,
    intake_row,
    missing_rows,
    orientations,
    require_complete,
)
from scripts.intake_phase_4j70 import BASE, REQUIRED, ROOT, inputs
from scripts.validate_phase_4j70_receipt import DECLARATIONS, declarations, semantics

REPO = Path(__file__).resolve().parents[1]
E = REPO / ROOT


def resign(r):
    return signed({k: v for k, v in r.items() if k != "identity"})


@pytest.fixture(scope="module")
def retained():
    model = load_model(REPO)
    return model, missing_rows(model)


def raw_frc(row, reverse=False, duplicate=False):
    # Manufactured nonzero source cells; not literature, native HMDS or QM evidence.
    t = orientations(row["family"], row["ordered_types"])[int(reverse)]
    n = len(row["contract"]["fields"])
    line = "1.0 0 " + " ".join(t) + " " + " ".join(str(i + 1) for i in range(n)) + "\n"
    return (
        "!BIOSYM forcefield 1\n"
        "#atom_types\n1.0 0 zz 12 C 4 synthetic\n"
        "#equivalence\n1.0 0 zz zz zz zz zz zz\n"
        "#auto_equivalence\n1.0 0 zz zz zz zz zz zz zz zz zz zz\n"
        "#bond_increments\n1.0 0 zz zz 0.1 -0.1\n#"
        + row["family"]
        + " cff91\n"
        + line
        + (line if duplicate else "")
        + "#end\n"
    ).encode()


def test_exact_378_matrix_no_values_or_registration(retained):
    model, rows = retained
    assert read(E / "completeness-matrix.json")["instances"] == rows
    assert len(rows) == 378 and len({r["id"] for r in rows}) == 378
    assert not completeness(model, rows)
    assert all(r["coefficients"] is None and r["source"] is None for r in rows)
    with pytest.raises(CompletionError, match="parameter_incomplete"):
        require_complete(
            model,
            rows,
            {"complete": True, "deterministic": True},
            full_inventory_closed=True,
        )


@pytest.mark.parametrize(
    "kind",
    [
        "partial",
        "duplicate",
        "site",
        "dependency",
        "identity",
        "zero",
        "wildcard",
        "units",
        "prefactor",
        "phase",
        "multiplicity",
        "equation",
        "status",
    ],
)
def test_resigned_row_and_bundle_mutations(retained, kind):
    model, original = retained
    rows = copy.deepcopy(original)
    r = rows[0]
    if kind == "partial":
        rows.pop()
    elif kind == "duplicate":
        rows[-1] = copy.deepcopy(r)
    elif kind == "site":
        r["stable_sites"][0] = 27
        r["id"] = r["family"] + ":" + ",".join(map(str, r["stable_sites"]))
    elif kind == "dependency":
        r["dependencies"][0]["equilibrium"] += 0.01
    elif kind == "identity":
        r["id"] += "forged"
    elif kind == "zero":
        r["coefficients"] = [0.0]
    elif kind == "wildcard":
        r["ordered_types"][0] = "*"
    elif kind == "units":
        r["contract"]["source_units"] = ["kJ"]
    elif kind == "prefactor":
        r["contract"]["prefactor"] = 0.5
    elif kind in ("phase", "multiplicity", "equation"):
        r["contract"]["scalar_energy" if kind == "equation" else kind] = "forged"
    elif kind == "status":
        r["status"] = "native_source"
    if kind not in ("partial", "duplicate"):
        rows[0] = resign(r)
    with pytest.raises(CompletionError):
        completeness(model, rows)


@pytest.mark.parametrize(
    "family",
    [
        "angle-angle",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-torsion_3",
        "bond-bond_1_3",
        "angle-angle-torsion_1",
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_manufactured_source_normalization_orientation(retained, family, reverse):
    row = next(r for r in retained[1] if r["family"] == family)
    raw = raw_frc(row, reverse=reverse)
    value = intake_row(
        row,
        raw,
        expected_sha256=sha256(raw).hexdigest(),
        category="literature_external",
    ).payload
    assert value["source"]["category"] == "literature_external"
    assert value["source"]["row_ids"] == [family + ":cff91:11"]
    expected = np.arange(1, len(row["contract"]["fields"]) + 1, dtype=float)
    # If ordered labels are palindromic, conflicting asymmetric blocks are rejected
    # by the existing resolver; selected fixtures here have unambiguous end order.
    if reverse and family in ("end_bond-torsion_3", "angle-torsion_3"):
        expected = np.r_[expected[3:], expected[:3]]
    np.testing.assert_allclose(
        value["normalized"]["values"], expected * 4.184, atol=1e-12
    )
    assert value["eligible_static"] == ("torsion" not in family)
    assert value["coefficients"][0] == list(range(1, len(expected) + 1))
    # Immutable storage and copied payload.
    value["coefficients"][0][0] = 999
    assert (
        SourceRow(
            canonical(
                intake_row(
                    row,
                    raw,
                    expected_sha256=sha256(raw).hexdigest(),
                    category="literature_external",
                ).payload
            )
        ).payload["coefficients"][0][0]
        == 1
    )


@pytest.mark.parametrize(
    "kind",
    [
        "hash",
        "row_id",
        "normalized",
        "raw_coefficient",
        "cross_family",
        "duplicate",
        "missing",
    ],
)
def test_source_proof_mutations(retained, kind):
    row = retained[1][0]
    raw = raw_frc(row)
    if kind == "hash":
        with pytest.raises(CompletionError, match="hash"):
            intake_row(
                row,
                raw + b"\n",
                expected_sha256=sha256(raw).hexdigest(),
                category="native_pcff",
            )
        return
    if kind == "cross_family":
        with pytest.raises(CompletionError, match="cross_family"):
            intake_row(
                row,
                raw,
                expected_sha256=sha256(raw).hexdigest(),
                category="native_pcff",
                source_family="COMPASS",
            )
        return
    if kind in ("duplicate", "missing"):
        raw = (
            raw_frc(row, duplicate=True)
            if kind == "duplicate"
            else b"!BIOSYM forcefield 1\n#end\n"
        )
        with pytest.raises(CompletionError):
            intake_row(
                row,
                raw,
                expected_sha256=sha256(raw).hexdigest(),
                category="native_pcff",
            )
        return
    p = intake_row(
        row, raw, expected_sha256=sha256(raw).hexdigest(), category="native_pcff"
    ).payload
    if kind == "row_id":
        p["source"]["row_ids"] = ["angle-angle:cff91:4"]
    elif kind == "normalized":
        p["normalized"]["conversion_factors"] = [1.0]
    else:
        p["coefficients"][0][0] = 0.0
    with pytest.raises(CompletionError):
        SourceRow(canonical(resign(p)))


def test_native_external_mixing_rejected(retained):
    model, rows = retained
    rows = copy.deepcopy(rows)
    for i, category in enumerate(("native_pcff", "literature_external")):
        raw = raw_frc(rows[i])
        rows[i] = intake_row(
            rows[i], raw, expected_sha256=sha256(raw).hexdigest(), category=category
        ).payload
    with pytest.raises(CompletionError, match="mixing"):
        completeness(model, rows)


def test_manufactured_complete_converter_table_gate():
    # Software-only complete-table transcript; not an executed native success.
    required = ["Bond", "Angle", "Dihedral", "AA", "EBT", "MBT", "AT", "AAT", "BB13"]
    assert converter_gate(0, "complete tables", required, required)
    assert not converter_gate(17, "Unable to find endbonddihedral", required, required)
    assert not converter_gate(0, "complete tables", required[:-1], required)
    assert not converter_gate(0, "missing row", required, required)


def test_actual_strict_converter_and_static_closure():
    c = read(E / "converter/report.json")
    assert [r["returncode"] for r in c["attempts"]] == [17, 17]
    assert c["deterministic_failure"] and not c["complete"]
    s = read(E / "static-closure.json")
    assert s["configurations"] == 26 and s["available_subset_closed"]
    assert s["max_force_component_difference_kcal_mol_A"] < 1e-9
    assert s["complete_source_static_closure"] == "blocked"


def test_canonical_pcff_fail_closed_exact_missing_diagnostics():
    p = read(E / "canonical-capability.json")
    assert p["ordinary_api"] == 'forcefield="pcff"' and p["selected_profile"] is None
    assert not p["automatic_assignment"] and len(p["missing_rows"]) == 378
    assert not p["manual_profile_required"] and not p["manual_charge_method_required"]
    assert p["ordinary_public_status"] != "prepared"


def test_actual_family_only_polymer_selection():
    from island.pcff_pdms import PCFFPDMSAssets, prepare_polymer

    assets = (
        REPO.parent
        / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
    )
    if not assets.is_dir():
        pytest.skip(
            "original declared source installation not available; offline replay remains supported"
        )
    system, prepared = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="pcff",
        dp=1,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        assets=PCFFPDMSAssets(assets),
    )
    p = prepared.validate_integrity(system)
    assert len(system.topology.sites) == 27 and len(p["missing_terms"]) == 378
    assert p["selected_profile"] is None and p["parameter_complete"] is False
    with pytest.raises(Exception, match="parameter_complete"):
        prepared.require()


def test_resigned_extra_row_schema_rejected(retained):
    p = copy.deepcopy(retained[1][0])
    p["source_complete"] = True
    with pytest.raises(CompletionError, match="schema_fields"):
        SourceRow(canonical(resign(p)))


@pytest.mark.parametrize("field", list(DECLARATIONS))
def test_resigned_receipt_declaration_mutations(field):
    p = signed(DECLARATIONS)
    v = p[field]
    p[field] = (
        not v
        if type(v) is bool
        else v + 1
        if type(v) is int
        else ["forged"]
        if isinstance(v, list)
        else "forged"
    )
    with pytest.raises(CompletionError, match="declaration"):
        declarations(resign(p))


@pytest.mark.parametrize(
    "mutation",
    [
        "none",
        "receipt_schema",
        "receipt_extra",
        "ancestry",
        "manifest_schema",
        "manifest_extra",
        "manifest_base",
        "nested_count",
        "nested_boolean_type",
    ],
)
def test_resigned_receipt_manifest_envelopes(mutation):
    # Manufactured envelopes exercise semantic validation independently of hashes.
    receipt = {
        **DECLARATIONS,
        "schema": "island.j70-receipt.v1",
        "base_commit": BASE,
        "required_J69": REQUIRED,
        "branch": "codex/phase-4j70-authoritative-pcff-missing-row-intake",
        "manifest_sha256": "0" * 64,
        "software_counts": {},
        "terminal_unique_tests": 0,
        "historical_files_preserved": 0,
        "numerical_reconstruction": {
            "offline": True,
            "instances": 378,
            "observations_retained": 26,
            "source_complete": False,
            "converter_complete": False,
            "parameter_assignment_capability": False,
            "new_QM": 0,
            "MD_steps": 0,
            "production_fits": 0,
            "accepted_coefficients": 0,
            "validation_untouched": True,
        },
    }
    manifest = {
        "schema": "island.numerical-evidence-manifest.v1",
        "base_commit": BASE,
        "implementation_files": {},
        "historical_files": {},
        "artifacts": {},
    }
    if mutation == "receipt_schema":
        receipt["schema"] = "forged"
    elif mutation == "receipt_extra":
        receipt["unvalidated_declaration"] = True
    elif mutation == "ancestry":
        receipt["required_J69"] = "0" * 40
    elif mutation == "manifest_schema":
        manifest["schema"] = "forged"
    elif mutation == "manifest_extra":
        manifest["unvalidated_declaration"] = True
    elif mutation == "manifest_base":
        manifest["base_commit"] = "0" * 40
    elif mutation == "nested_count":
        receipt["numerical_reconstruction"]["new_QM"] = 1
    elif mutation == "nested_boolean_type":
        receipt["numerical_reconstruction"]["validation_untouched"] = 1
    if mutation == "none":
        semantics(signed(receipt), signed(manifest))
    else:
        with pytest.raises(CompletionError):
            semantics(signed(receipt), signed(manifest))


def test_inventory_and_history_preserved():
    inputs(REPO, E)
    d = read(E / "decision.json")
    assert d["reserved_validation_activated"] is False
    assert d["new_QM"] == d["MD_steps"] == d["production_fits"] == 0
    assert d["J68_bound_free_incompatibility"].startswith("preserved")


def test_fresh_process_offline_relocation(tmp_path):
    moved = tmp_path / "relocated"
    moved.mkdir()
    for directory in ("src", "scripts", "docs/evidence"):
        shutil.copytree(
            REPO / directory,
            moved / directory,
            ignore=shutil.ignore_patterns("__pycache__"),
        )
    code = """
import importlib.abc,sys
from pathlib import Path
class BlockNative(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in {'psi4','clarabel','scs'}:raise AssertionError('native execution forbidden')
sys.meta_path.insert(0,BlockNative())
from scripts.replay_phase_4j70 import replay
print(replay(Path.cwd(),Path('docs/evidence/phase_4j70'),deep=True))
"""
    p = subprocess.run(
        [sys.executable, "-c", code],
        cwd=moved,
        env={
            **os.environ,
            "PYTHONPATH": str(moved) + ":" + str(moved / "src"),
            "OPENBLAS_NUM_THREADS": "1",
        },
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    assert p.returncode == 0, p.stdout + p.stderr
    assert "'validation_untouched': True" in p.stdout


@pytest.mark.parametrize(
    "mutation",
    [
        "static_metric",
        "static_count",
        "static_parent",
        "cap_status",
        "cap_profile",
        "cap_source_lookup",
        "cap_dependencies",
        "cap_diagnostics",
        "review_downloads",
        "review_parsed_namespace",
        "review_overlay",
        "converter_summary_emission",
        "converter_sidecar",
        "converter_ignore",
        "converter_source",
    ],
)
def test_resigned_nested_report_mutations(tmp_path, mutation):
    from island.pcff_completion import write
    from scripts.replay_phase_4j70 import replay

    # Copied metadata mutation controls; no native calculations or source changes.
    root = tmp_path / "tampered-evidence"
    shutil.copytree(E, root)
    if mutation.startswith("static_"):
        name = "static-closure.json"
        p = read(root / name)
        if mutation == "static_metric":
            p["max_energy_difference_kcal_mol"] = 1.0
        elif mutation == "static_count":
            p["configurations"] = 25
        else:
            p["parent_reconstruction"]["accepted_coefficients"] = 1
    elif mutation.startswith("cap_"):
        name = "canonical-capability.json"
        p = read(root / name)
        if mutation == "cap_status":
            p["ordinary_public_status"] = "prepared"
        elif mutation == "cap_profile":
            p["registered_bounded_profiles"] = ["forged"]
        elif mutation == "cap_source_lookup":
            p["missing_rows"][0]["source_lookup"]["status"] = "assigned"
        elif mutation == "cap_dependencies":
            p["missing_rows"][0]["dependencies"] = ["forged"]
        else:
            p["ordinary_public_diagnostics"] = []
    elif mutation.startswith("review_"):
        name = "source-review.json"
        p = read(root / name)
        if mutation == "review_downloads":
            p["network_downloads"] = 1
        elif mutation == "review_overlay":
            p["sources"][0]["eligible_native_overlay"] = True
        else:
            next(r for r in p["parsed_PCFF_files"] if r["category"] != "native_pcff")[
                "category"
            ] = "native_pcff"
    elif mutation == "converter_sidecar":
        name = "converter/attempt-1/attempt.json"
        p = read(root / name)
        p["returncode"] = 0
    else:
        name = "converter/report.json"
        p = read(root / name)
        if mutation == "converter_summary_emission":
            p["full_emission"] = "complete"
        else:
            a = p["attempts"][0]
            if mutation == "converter_ignore":
                a["argv"].append("-ignore")
            else:
                a["source_sha256"] = "0" * 64
            p["attempts"][0] = resign(a)
            write(root / "converter/attempt-1/attempt.json", p["attempts"][0])
    write(root / name, resign(p))
    with pytest.raises((CompletionError, ValueError)):
        replay(REPO, root)
