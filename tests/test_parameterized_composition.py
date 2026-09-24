from dataclasses import replace

import pytest

pytest.importorskip("rdkit")

from island.builders import build_linear_polymer
from island.exceptions import (
    InvalidChargeAssignmentResultError,
    InvalidParameterAssignmentResultError,
)
from island.forcefields import (
    NonbondedPolicy,
    ParameterAssignmentEngine,
    ParameterizedSystem,
    ProvidedChargeEngine,
    RDKitSmartsAtomTypingEngine,
    compose_parameterized_system,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)


def policy(*, coulomb_14: float = 0.833333333333) -> NonbondedPolicy:
    return NonbondedPolicy(
        "synthetic_composition_policy",
        "1",
        "lorentz_berthelot",
        0.0,
        0.0,
        0.0,
        0.0,
        0.5,
        coulomb_14,
        "SYNTHETIC SOFTWARE-TEST POLICY",
    )


def components() -> tuple[object, object, object]:
    system = build_linear_polymer("[*]CC[*]", dp=2, generate_3d=False)
    ruleset = island_demo_v1_ruleset()
    typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    parameters = ParameterAssignmentEngine().assign(
        system, typing, ruleset, island_demo_parameters_v1()
    )
    charges = ProvidedChargeEngine().assign(
        system,
        {site_id: 0.0 for site_id in system.topology.sites},
        source="SYNTHETIC SOFTWARE-TEST CHARGES",
    )
    return system, parameters, charges


def test_composition_creates_owned_snapshot_with_separate_statuses() -> None:
    system, parameters, charges = components()
    before = system.to_dict()
    composed = compose_parameterized_system(system, parameters, charges, policy())
    assert system.to_dict() == before
    assert composed.system.to_dict() == before
    assert set(composed.charge_assignments) == set(system.topology.sites)
    assert composed.nonbonded_policy == policy()
    assert len(composed.aggregate_signature) == 64
    assert composed.metadata["aggregate"] == {
        "aggregate_signature": composed.aggregate_signature,
        "parameter_coverage_complete": True,
        "charge_assignment_complete": True,
        "nonbonded_policy_available": True,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    assert parameters.charges_status == "unassigned"
    mixed = policy().mixed_parameters_for_types(
        parameters, "demo_c_aliphatic", "demo_h_on_carbon"
    )
    assert mixed.sigma == pytest.approx((0.34 + 0.25) / 2)
    assert mixed.epsilon == pytest.approx((0.30 * 0.10) ** 0.5)

    original = composed.system.coordinates.get(1)
    system.coordinates.translate([20.0, 0.0, 0.0])
    assert composed.system.coordinates.get(1) == pytest.approx(original)


def test_parameterized_system_classmethod_is_equivalent() -> None:
    system, parameters, charges = components()
    direct = compose_parameterized_system(system, parameters, charges, policy())
    classmethod_result = ParameterizedSystem.from_components(
        system, parameters, charges, policy()
    )
    assert classmethod_result.aggregate_signature == direct.aggregate_signature
    assert classmethod_result.site_assignments == direct.site_assignments
    assert classmethod_result.charge_assignments == direct.charge_assignments


def test_coordinate_changes_preserve_aggregate_signature() -> None:
    system, parameters, charges = components()
    first = compose_parameterized_system(system, parameters, charges, policy())
    moved = system.copy()
    moved.coordinates.translate([30.0, -4.0, 2.0])
    second = compose_parameterized_system(moved, parameters, charges, policy())
    assert first.aggregate_signature == second.aggregate_signature


def test_policy_changes_invalidate_aggregate_reuse_signature() -> None:
    system, parameters, charges = components()
    first = compose_parameterized_system(system, parameters, charges, policy())
    second = compose_parameterized_system(
        system, parameters, charges, policy(coulomb_14=0.75)
    )
    assert first.aggregate_signature != second.aggregate_signature


def test_malformed_parameter_or_charge_result_is_rejected_before_composition() -> None:
    system, parameters, charges = components()
    parameter_assignments = dict(parameters.bond_assignments)
    parameter_assignments.pop(next(iter(parameter_assignments)))
    malformed_parameters = replace(
        parameters,
        bond_assignments=parameter_assignments,
        metadata=dict(parameters.metadata),
    )
    before = system.to_dict()
    with pytest.raises(InvalidParameterAssignmentResultError):
        compose_parameterized_system(system, malformed_parameters, charges, policy())

    charge_assignments = dict(charges.assignments)
    first_site = min(charge_assignments)
    charge_assignments[first_site] = replace(charge_assignments[first_site], charge=0.1)
    malformed_charges = replace(
        charges,
        assignments=charge_assignments,
        metadata=dict(charges.metadata),
    )
    with pytest.raises(InvalidChargeAssignmentResultError):
        compose_parameterized_system(system, parameters, malformed_charges, policy())
    assert system.to_dict() == before
