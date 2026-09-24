"""Regressions for malformed charge results and shared typing validation."""

import os
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from math import inf, nan

import pytest
from test_atom_type_charges import charge_table
from test_parameter_assignment_graph_only import graph_library, graph_typing
from test_parameterized_composition import components, policy

from island.exceptions import (
    InvalidChargeAssignmentResultError,
    InvalidTypingResultError,
)
from island.forcefields import (
    AtomTypeChargeEngine,
    ParameterAssignmentEngine,
    ParameterizedSystem,
    ProvidedChargeEngine,
    compose_parameterized_system,
)


def _malformed_charge(result: object, field: str, value: object) -> object:
    if field == "assignment":
        assignments = dict(result.assignments)
        assignments[min(assignments)] = value
        return replace(result, assignments=assignments)
    return replace(result, **{field: value})


@pytest.mark.parametrize(
    ("field", "value", "context"),
    [
        ("assignment", None, "site"),
        ("tolerance", "bad", "tolerance"),
        ("diagnostics", (None,), "diagnostic"),
        ("coverage", None, "coverage"),
        ("target_charge", nan, "target_charge"),
        ("total_assigned_charge", inf, "total_assigned_charge"),
        ("total_within_tolerance", 1, "total_within_tolerance"),
    ],
)
def test_malformed_charge_result_has_focused_boundary_errors(
    field: str, value: object, context: str
) -> None:
    system, parameters, charges = components()
    before = system.to_dict()
    malformed = _malformed_charge(charges, field, value)

    assert malformed.is_input_compatible_with(system)
    assert not malformed.is_compatible_with(system)
    with pytest.raises(InvalidChargeAssignmentResultError, match=context):
        malformed.validate_integrity(system)
    with pytest.raises(InvalidChargeAssignmentResultError, match=context):
        compose_parameterized_system(system, parameters, malformed, policy())
    with pytest.raises(InvalidChargeAssignmentResultError, match=context):
        ParameterizedSystem.from_components(system, parameters, malformed, policy())
    assert system.to_dict() == before
    assert charges.is_compatible_with(system)


@pytest.mark.parametrize(
    ("field", "value", "context"),
    [
        ("diagnostics", (None,), "diagnostic"),
        ("coverage", None, "coverage"),
        ("tolerance", True, "tolerance"),
        ("total_formal_charge", nan, "total_formal_charge"),
        ("total_charge_residual", inf, "total_charge_residual"),
    ],
)
def test_malformed_numeric_and_diagnostic_fields_are_rejected(
    field: str, value: object, context: str
) -> None:
    system, _, charges = components()
    malformed = replace(charges, **{field: value})
    assert not malformed.is_compatible_with(system)
    with pytest.raises(InvalidChargeAssignmentResultError, match=context):
        malformed.validate_integrity(system)


def test_malformed_component_and_diagnostic_records_are_rejected() -> None:
    system, _, charges = components()
    changes = (
        ("component_diagnostics", (None,), "component"),
        (
            "component_diagnostics",
            (replace(charges.component_diagnostics[0], within_tolerance=False),),
            "component",
        ),
        (
            "diagnostics",
            (replace(charges.diagnostics[0], reason="bad"),)
            if charges.diagnostics
            else (None,),
            "diagnostic",
        ),
    )
    for field, value, context in changes:
        malformed = replace(charges, **{field: value})
        with pytest.raises(InvalidChargeAssignmentResultError, match=context):
            malformed.validate_integrity(system)
        assert not malformed.is_compatible_with(system)


def test_invalid_charge_mismatch_reason_and_component_values_are_rejected() -> None:
    system, _, charges = components()
    provided = ProvidedChargeEngine().assign(
        system,
        {
            site_id: (0.25 if site_id == min(system.topology.sites) else 0.0)
            for site_id in system.topology.sites
        },
        source="synthetic invalid-charge fixture",
        strict=False,
    )
    assert not provided.complete
    malformed_reason = replace(
        provided,
        diagnostics=(replace(provided.diagnostics[0], reason="invalid"),),
    )
    with pytest.raises(InvalidChargeAssignmentResultError, match="reason"):
        malformed_reason.validate_integrity(system)
    assert not malformed_reason.is_compatible_with(system)

    bad_component = replace(
        charges,
        component_diagnostics=(
            replace(charges.component_diagnostics[0], residual=nan),
        ),
    )
    with pytest.raises(InvalidChargeAssignmentResultError, match="residual"):
        bad_component.validate_integrity(system)


def test_typing_record_scalar_types_are_validated_by_both_engines() -> None:
    topology, typing, ruleset = graph_typing()
    malformed = deepcopy(typing)
    malformed.assignments[10] = replace(malformed.assignments[10], atom_type=["A"])
    for engine, args in (
        (ParameterAssignmentEngine(), (graph_library(ruleset),)),
        (AtomTypeChargeEngine(), (charge_table(ruleset),)),
    ):
        with pytest.raises(InvalidTypingResultError, match="atom_type"):
            engine.assign(topology, malformed, ruleset, *args)


@pytest.mark.parametrize(
    "mutation", ["unknown_rule", "diagnostic_id", "assignment_type", "diagnostic_type"]
)
def test_both_engines_reject_the_same_invalid_typing_result(mutation: str) -> None:
    topology, typing, ruleset = graph_typing()
    damaged = deepcopy(typing)
    if mutation == "unknown_rule":
        assignment = damaged.assignments[10]
        ids = (*assignment.matched_rule_ids, "zzz_unknown")
        damaged.assignments[10] = replace(
            assignment, selected_rule_ids=ids, matched_rule_ids=ids
        )
        damaged.diagnostics[10] = replace(
            damaged.diagnostics[10], surviving_rule_ids=ids, matched_rule_ids=ids
        )
    elif mutation == "diagnostic_id":
        damaged.diagnostics[10] = replace(damaged.diagnostics[10], site_id=999)
    elif mutation == "assignment_type":
        damaged.assignments[10] = None
    else:
        damaged.diagnostics[10] = None

    with pytest.raises(InvalidTypingResultError):
        ParameterAssignmentEngine().assign(
            topology, damaged, ruleset, graph_library(ruleset)
        )
    with pytest.raises(InvalidTypingResultError):
        AtomTypeChargeEngine().assign(topology, damaged, ruleset, charge_table(ruleset))
    assert typing.assignments[10].atom_type == "A"


def test_valid_incomplete_charge_result_remains_a_valid_diagnostic_object() -> None:
    topology, typing, ruleset = graph_typing()
    table = charge_table(ruleset)
    reduced = replace(
        table,
        entries=tuple(entry for entry in table.entries if entry.atom_type != "A"),
    )
    result = AtomTypeChargeEngine().assign(
        topology, typing, ruleset, reduced, strict=False
    )
    assert not result.complete
    assert result.coverage.missing == 1
    result.validate_integrity(topology, typing_result=typing, table=reduced)
    assert result.is_compatible_with(topology, typing_result=typing, table=reduced)


def test_valid_copying_and_coordinate_independence_remain_supported() -> None:
    system, _, charges = components()
    moved = system.copy()
    moved.coordinates.translate([13.0, -7.0, 2.0])
    assert charges.is_compatible_with(moved)
    assert deepcopy(charges).is_compatible_with(moved)
    assert replace(charges).is_compatible_with(moved)

    provided = ProvidedChargeEngine().assign(
        system,
        {site_id: 0.0 for site_id in system.topology.sites},
        source="synthetic validation fixture",
    )
    assert provided.is_compatible_with(moved)


def test_shared_typing_validation_imports_without_rdkit() -> None:
    code = """
import importlib.abc
import sys

class BlockRDKit(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'rdkit' or fullname.startswith('rdkit.'):
            raise RuntimeError('typing validation imported RDKit')
        return None

sys.meta_path.insert(0, BlockRDKit())
from island.forcefields.typing import validate_complete_typing_result
from island.forcefields.parameters import ParameterAssignmentEngine
from island.forcefields.charges import AtomTypeChargeEngine
assert callable(validate_complete_typing_result)
assert 'rdkit' not in sys.modules
"""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=False,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
