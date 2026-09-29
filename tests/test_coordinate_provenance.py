"""Current coordinate provenance follows coordinate application across modules."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("rdkit")
from test_minimization import Harmonic

from island.builders import build_linear_polymer
from island.minimization import minimize_geometry
from island.minimization.models import coordinate_hash


class Target(Harmonic):
    def __init__(self, system):
        super().__init__()
        self.target = {
            site: system.coordinates.get(site) + 0.02 for site in system.topology.sites
        }

    def evaluate(self, coordinates=None, **kwargs):
        shifted = {
            site: np.array(xyz) - self.target[site] for site, xyz in coordinates.items()
        }
        return replace(
            super().evaluate(shifted, **kwargs),
            coordinate_fingerprint=coordinate_hash(coordinates),
        )


def test_local_template_minimization_current_provenance():
    system = build_linear_polymer(
        "[*:1]CC[*:2]",
        dp=3,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    original = deepcopy(system.to_dict())
    result = minimize_geometry(system, Target(system))
    assert result.converged
    output = result.to_system(system)
    output.validate()
    assert output.metadata["coordinate_source"] == "local_minimization"
    assert output.metadata["polymer"]["coordinates"] == "local_minimization"
    assert (
        output.metadata["polymer"]["coordinate_generation"]["coordinate_source"]
        == "local_minimization"
    )
    assert (
        output.metadata["local_minimization"]["original_coordinate_source"]
        == "local_templates_etkdg_incremental"
    )
    assert system.to_dict() == original


def test_history_repeated_minimization_and_conformation_application():
    from island.conformations import ConformationResult

    system = build_linear_polymer(
        "[*:1]CC[*:2]",
        dp=3,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    system.metadata["ambertools_preparation"] = {
        "signature": "historical",
        "charges": [0.1, -0.1],
    }
    system.metadata["parameter_provenance"] = {"signature": "parameters"}
    historical = deepcopy(system.metadata["polymer"]["coordinate_generation"])
    first_result = minimize_geometry(system, Target(system))
    first = first_result.to_system(system)
    previous = first.metadata["coordinate_history"][-1]
    assert previous["polymer_records"]["coordinate_generation"] == historical
    assert (
        "minimum_nonbonded_distance"
        not in first.metadata["polymer"]["coordinate_generation"]
    )
    assert (
        first.metadata["local_minimization"]["initial_coordinate_fingerprint"]
        == first_result.initial_coordinate_fingerprint
    )
    assert (
        first.metadata["local_minimization"]["coordinate_fingerprint"]
        == first_result.coordinate_fingerprint
    )
    second_result = minimize_geometry(first, Target(first))
    second = second_result.to_system(first)
    assert (
        second.metadata["local_minimization"]["original_coordinate_source"]
        == "local_minimization"
    )
    assert (
        len(second.metadata["coordinate_history"])
        == len(first.metadata["coordinate_history"]) + 1
    )
    assert (
        second.metadata["coordinate_history"][-1]["records"]["local_minimization"]
        == first.metadata["local_minimization"]
    )
    conformation = ConformationResult(
        second.coordinates.copy(), "test_conformation", 1, True, 1, 0, 0
    )
    third = conformation.apply_to(second)
    assert third.metadata["coordinate_source"] == "test_conformation"
    assert third.metadata["polymer"]["coordinates"] == "test_conformation"
    assert (
        third.metadata["polymer"]["coordinate_generation"]["coordinate_source"]
        == "test_conformation"
    )
    assert "local_minimization" not in third.metadata
    assert (
        third.metadata["coordinate_history"][-1]["records"]["local_minimization"]
        == second.metadata["local_minimization"]
    )
    for output in (first, second, third):
        output.validate()
        assert (
            output.metadata["ambertools_preparation"]
            == system.metadata["ambertools_preparation"]
        )
        assert (
            output.metadata["parameter_provenance"]
            == system.metadata["parameter_provenance"]
        )
        for site in system.topology.sites:
            assert output.topology.sites[site] == system.topology.sites[site]
        assert (
            output.metadata["polymer"]["sequence"]
            == system.metadata["polymer"]["sequence"]
        )
    assert (
        first.metadata["polymer"]["coordinate_generation"]["coordinate_source"]
        == "local_minimization"
    )


def test_nonpolymer_generation_namespace_and_diagnostic_source():
    from test_minimization import molecule

    from island.conformations import ConformationResult
    from island.minimization import MinimizationOptions

    system = molecule()
    system.metadata.pop("coordinate_source")
    system.metadata["coordinate_generation"] = {
        "coordinate_source": "external_coordinates",
        "old_diagnostic": 3.5,
    }
    result = minimize_geometry(
        system, Target(system), MinimizationOptions(max_evaluations=1)
    )
    output = result.to_system(system, allow_unconverged=True)
    assert (
        output.metadata["local_minimization"]["original_coordinate_source"]
        == "external_coordinates"
    )
    assert (
        output.metadata["coordinate_source"]
        == output.metadata["coordinate_generation"]["coordinate_source"]
        == "local_minimization_diagnostic"
    )
    assert "old_diagnostic" not in output.metadata["coordinate_generation"]
    assert (
        output.metadata["coordinate_history"][-1]["records"]["coordinate_generation"][
            "old_diagnostic"
        ]
        == 3.5
    )
    conformation = ConformationResult(
        output.coordinates.copy(), "next", 1, True, 1, 0, 0
    )
    conformation.apply_to(output, copy=False)
    assert (
        output.metadata["coordinate_source"]
        == output.metadata["coordinate_generation"]["coordinate_source"]
        == "next"
    )
    assert "local_minimization" not in output.metadata
    output.validate()
