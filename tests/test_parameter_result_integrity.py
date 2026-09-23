from copy import deepcopy
from dataclasses import replace

import pytest
from test_parameter_assignment_graph_only import graph_library, graph_typing

from island import Coordinates, MolecularSystem
from island.exceptions import InvalidParameterAssignmentResultError
from island.forcefields import ParameterizedSystem
from island.forcefields.parameters import (
    FamilyCoverage,
    ParameterAssignmentEngine,
    ParameterAssignmentResult,
)


def valid_result() -> tuple[object, object, object, object, ParameterAssignmentResult]:
    topology, typing, ruleset = graph_typing()
    library = graph_library(ruleset)
    result = ParameterAssignmentEngine().assign(topology, typing, ruleset, library)
    return topology, typing, ruleset, library, result


def molecular_system(topology: object) -> MolecularSystem:
    return MolecularSystem(
        topology=topology,
        coordinates=Coordinates(
            {
                site_id: [float(index), 0.0, 0.0]
                for index, site_id in enumerate(topology.sites)
            }
        ),
    )


def replace_bond_selection(
    result: ParameterAssignmentResult, **changes: object
) -> ParameterAssignmentResult:
    assignments = dict(result.bond_assignments)
    key = next(iter(assignments))
    assignments[key] = replace(assignments[key], **changes)
    return replace(
        result,
        bond_assignments=assignments,
        metadata=dict(result.metadata),
    )


def assert_malformed_result_rejected(
    malformed: ParameterAssignmentResult,
    topology: object,
    typing: object,
    library: object,
) -> None:
    system = molecular_system(topology)
    before = system.to_dict()
    assert malformed.is_input_compatible_with(topology, typing, library)
    assert not malformed.is_compatible_with(topology, typing, library)
    with pytest.raises(InvalidParameterAssignmentResultError):
        malformed.validate_integrity(topology, typing_result=typing, library=library)
    with pytest.raises(InvalidParameterAssignmentResultError):
        malformed.to_parameterized_system(system)
    with pytest.raises(InvalidParameterAssignmentResultError):
        ParameterizedSystem.from_assignment(system, malformed)
    assert system.to_dict() == before


def test_mapping_key_and_selection_site_ids_must_agree() -> None:
    topology, typing, _, library, result = valid_result()
    malformed = replace_bond_selection(result, site_ids=(999, 1000))
    assert_malformed_result_rejected(malformed, topology, typing, library)


@pytest.mark.parametrize(
    "change",
    [
        {"family": "angle"},
        {"atom_types": ("wrong", "types")},
        {"parameter_id": "wrong_parameter"},
        {"source": "wrong source"},
    ],
)
def test_selection_wrapper_inconsistencies_are_rejected(
    change: dict[str, object],
) -> None:
    topology, typing, _, library, result = valid_result()
    malformed = replace_bond_selection(result, **change)
    assert_malformed_result_rejected(malformed, topology, typing, library)


def test_selected_record_family_and_library_identity_are_validated() -> None:
    topology, typing, _, library, result = valid_result()
    key = next(iter(result.bond_assignments))
    selection = result.bond_assignments[key]
    wrong_family_record = next(iter(result.angle_assignments.values())).parameter
    malformed_family = replace_bond_selection(result, parameter=wrong_family_record)
    assert_malformed_result_rejected(malformed_family, topology, typing, library)

    wrong_library_record = replace(selection.parameter, library_name="other_library")
    malformed_library = replace_bond_selection(result, parameter=wrong_library_record)
    assert_malformed_result_rejected(malformed_library, topology, typing, library)


def test_coverage_and_completeness_must_match_result_content() -> None:
    topology, typing, _, library, result = valid_result()
    coverage = dict(result.coverage)
    original = coverage["bond"]
    coverage["bond"] = FamilyCoverage(
        required=original.required,
        assigned=original.assigned - 1,
        missing=1,
        ambiguous=0,
    )
    malformed = replace(result, coverage=coverage, metadata=dict(result.metadata))
    assert_malformed_result_rejected(malformed, topology, typing, library)

    false_complete = replace(
        result,
        complete_supported_scope=False,
        metadata=dict(result.metadata),
    )
    assert_malformed_result_rejected(false_complete, topology, typing, library)


def test_changed_selected_numerical_value_invalidates_stale_result_signature() -> None:
    topology, typing, _, library, result = valid_result()
    key = next(iter(result.bond_assignments))
    selection = result.bond_assignments[key]
    changed_record = replace(
        selection.parameter,
        force_constant=selection.parameter.force_constant + 1.0,
    )
    malformed = replace_bond_selection(result, parameter=changed_record)
    assert malformed.assignment_signature == result.assignment_signature
    assert_malformed_result_rejected(malformed, topology, typing, library)


def test_noop_replace_and_deepcopy_preserve_valid_protected_result() -> None:
    topology, typing, _, library, result = valid_result()
    reconstructed = replace(result)
    copied = deepcopy(result)

    for candidate in (reconstructed, copied):
        candidate.validate_integrity(topology, typing_result=typing, library=library)
        assert candidate.assignment_signature == result.assignment_signature
        assert candidate.is_compatible_with(topology, typing, library)
        with pytest.raises(TypeError):
            candidate.bond_assignments[(999, 1000)] = next(  # type: ignore[index]
                iter(candidate.bond_assignments.values())
            )


def test_metadata_isolated_from_callers_and_deep_copies() -> None:
    topology, typing, _, library, result = valid_result()
    caller_metadata = {"nested": {"labels": ["initial"]}}
    reconstructed = replace(result, metadata=caller_metadata)
    caller_metadata["nested"]["labels"].append("caller mutation")
    assert reconstructed.metadata["nested"]["labels"] == ["initial"]

    copied = deepcopy(reconstructed)
    copied.metadata["nested"]["labels"].append("copy mutation")
    assert reconstructed.metadata["nested"]["labels"] == ["initial"]
    copied.validate_integrity(topology, typing_result=typing, library=library)


def test_library_record_content_is_checked_even_with_matching_parameter_id() -> None:
    topology, typing, ruleset, library, result = valid_result()
    key = next(iter(result.bond_assignments))
    selection = result.bond_assignments[key]
    changed_record = replace(
        selection.parameter,
        equilibrium_length=selection.parameter.equilibrium_length + 0.001,
    )
    records = tuple(
        changed_record if record.parameter_id == changed_record.parameter_id else record
        for record in library.records
    )
    changed_library = replace(library, records=records)
    changed_result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, changed_library
    )
    assert changed_result.assignment_signature != result.assignment_signature
    assert changed_result.library_signature != result.library_signature
    assert not result.is_compatible_with(topology, typing, changed_library)
