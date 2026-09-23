from copy import deepcopy
from dataclasses import replace

import pytest

pytest.importorskip("rdkit")

from island import Coordinates, MolecularSystem
from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.exceptions import (
    IncompleteParameterAssignmentError,
    InvalidParameterDefinitionError,
    InvalidTypingResultError,
    UnsupportedParameterRequirementError,
)
from island.forcefields import (
    ParameterAssignmentEngine,
    ParameterizedSystem,
    RDKitSmartsAtomTypingEngine,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)


def typed_system(smiles: str) -> tuple[MolecularSystem, object, object]:
    system = from_smiles(smiles, random_seed=41)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    return system, typing, ruleset


def test_known_methane_assignments_have_independent_expected_counts() -> None:
    system, typing, ruleset = typed_system("C")
    result = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1()
    )
    assert result.complete_supported_scope
    assert {
        family: (coverage.required, coverage.assigned)
        for family, coverage in result.coverage.items()
    } == {
        "site": (5, 5),
        "bond": (4, 4),
        "angle": (6, 6),
        "proper_torsion": (0, 0),
    }
    assert {item.parameter_id for item in result.site_assignments.values()} == {
        "demo_lj_c",
        "demo_lj_hc",
    }
    assert {item.parameter_id for item in result.bond_assignments.values()} == {
        "demo_bond_ch"
    }
    assert {item.parameter_id for item in result.angle_assignments.values()} == {
        "demo_angle_hch"
    }
    assert result.charges_status == "unassigned"
    assert not result.production_validated
    assert result.simulation_readiness == "not_established"
    with pytest.raises(TypeError):
        result.site_assignments[999] = next(  # type: ignore[index]
            iter(result.site_assignments.values())
        )


def test_multiterm_torsion_is_one_selected_record() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    result = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1()
    )
    selections = [
        item
        for item in result.proper_torsion_assignments.values()
        if item.parameter_id == "demo_torsion_cccc"
    ]
    assert selections
    assert all(len(item.parameter.terms) == 2 for item in selections)
    assert result.coverage["proper_torsion"].ambiguous == 0


def test_missing_and_conflicting_records_have_structured_diagnostics() -> None:
    system, typing, ruleset = typed_system("CO")
    engine = ParameterAssignmentEngine()
    partial = engine.assign(
        system,
        typing,
        ruleset,
        island_demo_parameters_v1(),
        strict=False,
    )
    assert not partial.complete_supported_scope
    oxygen_id = next(
        site_id
        for site_id, assignment in typing.assignments.items()
        if assignment.atom_type == "demo_o_hydroxyl"
    )
    oxygen_diagnostic = next(
        item
        for item in partial.diagnostics
        if item.family == "site" and item.site_ids == (oxygen_id,)
    )
    assert oxygen_diagnostic.reason == "missing"
    assert oxygen_diagnostic.atom_types == ("demo_o_hydroxyl",)
    assert oxygen_diagnostic.candidate_parameter_ids == ()
    with pytest.raises(IncompleteParameterAssignmentError) as captured:
        engine.assign(system, typing, ruleset, island_demo_parameters_v1())
    assert captured.value.result.diagnostics == partial.diagnostics

    library = island_demo_parameters_v1()
    duplicate = replace(
        library.records[0], parameter_id="demo_lj_c_conflict", epsilon=0.31
    )
    conflicting = replace(library, records=(*library.records, duplicate))
    ambiguous = engine.assign(
        build_linear_polymer("[*]CC[*]", dp=1, generate_3d=False),
        RDKitSmartsAtomTypingEngine().type_system(
            build_linear_polymer("[*]CC[*]", dp=1, generate_3d=False), ruleset
        ),
        ruleset,
        conflicting,
        strict=False,
    )
    assert ambiguous.coverage["site"].ambiguous == 2
    assert any(
        item.candidate_parameter_ids == ("demo_lj_c", "demo_lj_c_conflict")
        for item in ambiguous.diagnostics
    )


def test_assignment_is_coordinate_independent_and_does_not_mutate_input() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    before = system.to_dict()
    engine = ParameterAssignmentEngine()
    first = engine.assign(system, typing, ruleset, island_demo_parameters_v1())
    assert system.to_dict() == before

    moved = system.copy()
    moved.coordinates.translate([50.0, -20.0, 3.0])
    moved.metadata["coordinate_generation"] = {"method": "replacement"}
    second = engine.assign(moved, typing, ruleset, island_demo_parameters_v1())
    assert first.assignment_signature == second.assignment_signature
    assert first.site_assignments == second.site_assignments
    assert first.bond_assignments == second.bond_assignments
    assert first.is_compatible_with(moved, typing, island_demo_parameters_v1())


@pytest.mark.parametrize("corruption", ["missing_key", "wrong_site", "diagnostic"])
def test_structurally_inconsistent_typing_result_is_rejected(corruption: str) -> None:
    system, typing, ruleset = typed_system("C")
    broken = deepcopy(typing)
    site_id = min(broken.assignments)
    if corruption == "missing_key":
        broken.assignments.pop(site_id)
    elif corruption == "wrong_site":
        broken.assignments[site_id] = replace(broken.assignments[site_id], site_id=9999)
    else:
        broken.diagnostics[site_id] = replace(
            broken.diagnostics[site_id], status="untyped"
        )
    with pytest.raises(InvalidTypingResultError):
        ParameterAssignmentEngine().assign(
            system, broken, ruleset, island_demo_parameters_v1()
        )


def test_misleading_complete_flag_does_not_hide_incomplete_typing() -> None:
    system, typing, ruleset = typed_system("C")
    broken = deepcopy(typing)
    site_id = min(broken.assignments)
    broken.assignments.pop(site_id)
    assert broken.complete is True
    assert broken.untyped_site_ids == ()
    with pytest.raises(InvalidTypingResultError, match="exactly cover"):
        ParameterAssignmentEngine().assign(
            system, broken, ruleset, island_demo_parameters_v1()
        )


def test_stale_graph_ruleset_and_mutated_typing_content_are_rejected() -> None:
    system, typing, ruleset = typed_system("C")
    stale_graph = system.copy()
    stale_graph.topology.sites[min(stale_graph.topology.sites)].formal_charge = 1
    with pytest.raises(InvalidTypingResultError, match="graph signature"):
        ParameterAssignmentEngine().assign(
            stale_graph, typing, ruleset, island_demo_parameters_v1()
        )

    stale_rules = replace(ruleset, description="changed content")
    with pytest.raises(InvalidTypingResultError, match="ruleset"):
        ParameterAssignmentEngine().assign(
            system, typing, stale_rules, island_demo_parameters_v1()
        )

    mutated = deepcopy(typing)
    site_id = min(mutated.assignments)
    mutated.assignments[site_id] = replace(
        mutated.assignments[site_id], atom_type="fabricated"
    )
    with pytest.raises(InvalidTypingResultError, match="disagree"):
        ParameterAssignmentEngine().assign(
            system, mutated, ruleset, island_demo_parameters_v1()
        )


def test_parameterized_system_is_an_owned_non_md_ready_snapshot() -> None:
    system, typing, ruleset = typed_system("C")
    before = system.to_dict()
    result = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1()
    )
    parameterized = ParameterizedSystem.from_assignment(system, result)
    assert system.to_dict() == before
    assert parameterized.system.to_dict() == before
    original_position = parameterized.system.coordinates.get(1)
    system.coordinates.translate([10.0, 0.0, 0.0])
    assert parameterized.system.coordinates.get(1) == pytest.approx(original_position)
    assert set(parameterized.site_assignments) == set(system.topology.sites)
    assert parameterized.metadata["parameter_assignment"]["charges_status"] == (
        "unassigned"
    )
    assert not parameterized.metadata["parameter_assignment"]["production_validated"]
    assert (
        parameterized.metadata["parameter_assignment"]["simulation_readiness"]
        == "not_established"
    )
    assert all(
        "atom_type" not in site.metadata
        for site in parameterized.system.topology.sites.values()
    )


def test_incomplete_result_cannot_create_parameterized_system() -> None:
    system, typing, ruleset = typed_system("CO")
    incomplete = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1(), strict=False
    )
    with pytest.raises(InvalidParameterDefinitionError, match="incomplete"):
        incomplete.to_parameterized_system(system)


def test_library_representation_requirement_is_enforced() -> None:
    system, typing, ruleset = typed_system("C")
    unsupported = system.copy()
    unsupported.representation = "coarse_grained"
    with pytest.raises(UnsupportedParameterRequirementError, match="atomistic"):
        ParameterAssignmentEngine().assign(
            unsupported, typing, ruleset, island_demo_parameters_v1()
        )


@pytest.mark.parametrize("dp", [50, 100])
def test_polyethylene_coverage_and_independent_inventory_counts(dp: int) -> None:
    system = build_linear_polymer("[*]CC[*]", dp=dp, generate_3d=False)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    result = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1()
    )
    assert result.complete_supported_scope
    assert {
        family: coverage.required for family, coverage in result.coverage.items()
    } == {
        "site": 6 * dp + 2,
        "bond": 6 * dp + 1,
        "angle": 12 * dp,
        "proper_torsion": 18 * dp - 9,
    }
    assert all(coverage.complete for coverage in result.coverage.values())
    assert not result.diagnostics


def test_failure_does_not_mutate_input() -> None:
    system, typing, ruleset = typed_system("CO")
    before = system.to_dict()
    with pytest.raises(IncompleteParameterAssignmentError):
        ParameterAssignmentEngine().assign(
            system, typing, ruleset, island_demo_parameters_v1()
        )
    assert system.to_dict() == before


def test_assignment_works_with_empty_coordinate_container() -> None:
    system, typing, ruleset = typed_system("C")
    coordinate_free = MolecularSystem(
        topology=system.topology,
        coordinates=Coordinates(),
        representation="atomistic",
    )
    result = ParameterAssignmentEngine().assign(
        coordinate_free, typing, ruleset, island_demo_parameters_v1()
    )
    assert result.complete_supported_scope
