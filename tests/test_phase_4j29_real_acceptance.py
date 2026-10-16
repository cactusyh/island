import json
import os
import re
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.analysis import (
    AnalysisSelection,
    calculate_observables,
    load_analysis_record,
    load_verified_trajectory,
    make_analysis_record,
    relocate_analysis_record,
    write_analysis_record,
)
from island.exceptions import AnalysisError
from scripts.validate_phase_4j29_receipt import validate_receipt

DEFAULT_FIXTURE = Path("../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full")


def _graph(atomic_numbers, bonds):
    return {
        "graph_identity": "a" * 64,
        "sites": [{"id": index, "atomic_number": number} for index, number in enumerate(atomic_numbers, 1)],
        "bonds": [{"sites": list(bond)} for bond in bonds],
        "molecule_membership": {str(index): "molecule-A" for index in range(1, len(atomic_numbers) + 1)},
    }


def _mapping_identity():
    return "b" * 64


def _fixture():
    path = Path(os.environ.get("ISLAND_J29_WORKFLOW", DEFAULT_FIXTURE))
    if not (path / "manifest.json").exists():
        pytest.skip(f"retained J23 fixture is unavailable: {path}")
    return path


def test_real_j23_loader_and_all_observables():
    trajectory = load_verified_trajectory(_fixture())
    assert trajectory.selection is not None
    assert trajectory.topology_bonds
    assert trajectory.molecule_membership
    names = ("radius_of_gyration", "end_to_end_distance", "re_over_rg", "gyration_tensor", "chain_anisotropy", "chain_orientation", "nematic_order", "density", "network_connectivity")
    values = calculate_observables(trajectory, names, block_size=2)
    assert {name for name, value in values.items() if value["status"] == "ok"} == set(names)
    assert values["radius_of_gyration"]["value"] > 0
    assert values["end_to_end_distance"]["value"] > 0
    assert values["re_over_rg"]["value"] > 0


def test_real_record_relocation_and_separate_process_payload(tmp_path):
    trajectory = load_verified_trajectory(_fixture())
    values = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg", "density"), block_size=2)
    record = make_analysis_record(trajectory, values, provenance={"phase": "4J29"})
    source = write_analysis_record(record, tmp_path / "analysis.json")
    destination = relocate_analysis_record(source, tmp_path / "relocated.json")
    loaded = load_analysis_record(destination, expected_trajectory=trajectory)
    assert loaded.identity == record.identity
    assert loaded.payload["selection"]["identity"] == trajectory.selection.identity
    assert json.dumps(loaded.payload["observables"], sort_keys=True) == json.dumps(record.payload["observables"], sort_keys=True)
    code = """
import json
import sys
from island.analysis import load_analysis_record, load_verified_trajectory
trajectory = load_verified_trajectory(sys.argv[2])
record = load_analysis_record(sys.argv[1], expected_trajectory=trajectory, expected_selection_identity=sys.argv[3])
print(json.dumps(record.payload, sort_keys=True, separators=(\",\", \":\")))
"""
    completed = subprocess.run(
        [sys.executable, "-c", code, str(destination), str(_fixture()), trajectory.selection.identity],
        cwd=Path.cwd(), check=True, capture_output=True, text=True,
    )
    assert completed.stdout.strip() == json.dumps(record.payload, sort_keys=True, separators=(",", ":"))


def test_selection_mutations_reject_or_change_identity():
    trajectory = load_verified_trajectory(_fixture())
    selection = trajectory.selection
    chain = dict(selection.chains[0])
    duplicate = dict(chain, rg_atom_ids=list(chain["rg_atom_ids"]) + [chain["rg_atom_ids"][0]])
    with pytest.raises(AnalysisError, match="duplicate"):
        replace(selection, chains=(duplicate,))
    missing_endpoint = dict(chain, rg_atom_ids=list(chain["rg_atom_ids"])[1:])
    with pytest.raises(AnalysisError, match="endpoints"):
        replace(selection, chains=(missing_endpoint,))
    wrong_graph = replace(selection, final_graph_identity="f" * 64)
    with pytest.raises(AnalysisError, match="final graph identity"):
        wrong_graph.validate(trajectory)
    with pytest.raises(AnalysisError, match="Only atomic-mass"):
        replace(selection, mass_weighting="uniform")
    with pytest.raises(AnalysisError, match="Unsupported PBC"):
        replace(selection, pbc_geometry_policy="triclinic_v1")


def test_periodic_real_workflow_requires_graph_bound_selection():
    trajectory = load_verified_trajectory(_fixture())
    unbound = replace(trajectory, selection=None)
    with pytest.raises(AnalysisError, match="requires an AnalysisSelection"):
        calculate_observables(unbound, ("radius_of_gyration",), block_size=2)


@pytest.mark.parametrize(
    ("name", "graph", "message"),
    (
        ("ring", _graph((6, 6, 6), ((1, 2), (2, 3), (3, 1))), "cycle"),
        ("branched", _graph((6, 6, 6, 6), ((1, 2), (1, 3), (1, 4))), "branched"),
        ("disconnected", _graph((6, 6, 6, 6), ((1, 2), (3, 4))), "disconnected"),
    ),
)
def test_from_final_graph_rejects_non_linear_heavy_topology(name, graph, message):
    del name
    with pytest.raises(AnalysisError, match=message):
        AnalysisSelection.from_final_graph(graph, atom_mapping_identity=_mapping_identity())


def test_from_final_graph_accepts_valid_explicit_h_linear_chain():
    selection = AnalysisSelection.from_final_graph(
        _graph((6, 6, 1, 1), ((1, 2), (1, 3), (2, 4))),
        atom_mapping_identity=_mapping_identity(),
    )
    assert selection.chains[0]["backbone_path"] == [1, 2]
    assert selection.chains[0]["endpoint_ids"] == [1, 2]
    assert selection.chains[0]["rg_atom_ids"] == [1, 2, 3, 4]


def test_j29_receipt_and_manifest_have_consistent_exact_head_schema():
    receipt = json.loads(Path("docs/evidence/phase_4j29.json").read_text())
    manifest = json.loads(Path("docs/paper/reproducibility_manifest.json").read_text())
    phase = manifest["phase_4j29"]
    assert receipt["schema"] == "island.phase-4j29-evidence.v1"
    assert re.fullmatch(r"[0-9a-f]{40}", receipt["validated_commit"])
    assert receipt["validated_commit"] == phase["validated_commit"]
    assert receipt["branch"] == phase["branch"]
    assert receipt["base_commit"] == phase["base_commit"]
    mutated = dict(receipt, validated_commit="0" * 40)
    assert mutated["validated_commit"] != receipt["validated_commit"]


def test_receipt_validator_rejects_real_mutated_files(tmp_path):
    receipt_path = Path("docs/evidence/phase_4j29.json")
    manifest_path = Path("docs/paper/reproducibility_manifest.json")
    validate_receipt(receipt_path, manifest_path)
    receipt = json.loads(receipt_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    mutated_receipt = tmp_path / "receipt-mutated.json"
    changed_receipt = dict(receipt, branch="0" * 40)
    mutated_receipt.write_text(json.dumps(changed_receipt))
    with pytest.raises(AnalysisError, match="branch is invalid"):
        validate_receipt(mutated_receipt, manifest_path)
    mutated_manifest = tmp_path / "manifest-mutated.json"
    changed_manifest = dict(manifest, phase_4j29=dict(manifest["phase_4j29"], validated_commit="0" * 40))
    mutated_manifest.write_text(json.dumps(changed_manifest))
    with pytest.raises(AnalysisError, match="receipt/manifest mismatch"):
        validate_receipt(receipt_path, mutated_manifest)
    changed_manifest = dict(manifest, production_validated=True)
    mutated_manifest.write_text(json.dumps(changed_manifest))
    with pytest.raises(AnalysisError, match="readiness"):
        validate_receipt(receipt_path, mutated_manifest)
    changed_manifest = dict(manifest, environment=dict(manifest["environment"], focused_tests="16 passed"))
    mutated_manifest.write_text(json.dumps(changed_manifest))
    with pytest.raises(AnalysisError, match="focused test count"):
        validate_receipt(receipt_path, mutated_manifest)
    changed_manifest = dict(manifest, environment=dict(manifest["environment"], ordinary_tests="2091 passed, 10 skipped, 0 failed, 0 errors"))
    mutated_manifest.write_text(json.dumps(changed_manifest))
    with pytest.raises(AnalysisError, match="ordinary test count"):
        validate_receipt(receipt_path, mutated_manifest)


@pytest.mark.parametrize(
    ("graph", "message"),
    (
        (_graph((6, 6), ((1, 2),)), "duplicate site IDs"),
        (_graph((6, 6), ((1, 3),)), "unknown site ID"),
        (_graph((6, 6), ((1, 1),)), "self-bond"),
        (_graph((6, 6), ((1, 2), (2, 1))), "duplicate bonds"),
        (_graph((6, 6), ((1, 2),)), "malformed"),
        (_graph((6, 6), ()), "heavy-atom topology"),
    ),
)
def test_from_final_graph_rejects_malformed_records(graph, message):
    if message == "duplicate site IDs":
        graph["sites"].append(dict(graph["sites"][0]))
        graph["molecule_membership"]["1"] = "molecule-A"
    elif message == "unknown site ID":
        graph["bonds"] = [{"sites": [1, 9]}]
    elif message == "self-bond":
        graph["bonds"] = [{"sites": [1, 1]}]
    elif message == "duplicate bonds":
        graph["bonds"] = [{"sites": [1, 2]}, {"sites": [2, 1]}]
    elif message == "malformed":
        graph["bonds"] = [{}]
    with pytest.raises(AnalysisError, match=message):
        AnalysisSelection.from_final_graph(graph, atom_mapping_identity=_mapping_identity())
