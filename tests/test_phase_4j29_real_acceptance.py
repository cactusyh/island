import json
import os
from dataclasses import replace
from pathlib import Path

import pytest

from island.analysis import (
    calculate_observables,
    load_analysis_record,
    load_verified_trajectory,
    make_analysis_record,
    relocate_analysis_record,
    write_analysis_record,
)
from island.exceptions import AnalysisError

DEFAULT_FIXTURE = Path("../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full")


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
