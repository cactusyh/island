"""Validated composition of numerical parameters, charges, and policy."""

import hashlib
import json
from copy import deepcopy

from island.core import MolecularSystem
from island.exceptions import ParameterCompositionError
from island.forcefields.charges.models import ChargeAssignmentResult
from island.forcefields.nonbonded import NonbondedPolicy, nonbonded_policy_signature
from island.forcefields.parameterized import ParameterizedSystem
from island.forcefields.parameters import ParameterAssignmentResult


def compose_parameterized_system(
    system: MolecularSystem,
    parameter_result: ParameterAssignmentResult,
    charge_result: ChargeAssignmentResult,
    nonbonded_policy: NonbondedPolicy,
) -> ParameterizedSystem:
    """Create an owned validated snapshot from three independent components."""
    if not isinstance(system, MolecularSystem):
        raise TypeError("system must be a MolecularSystem")
    if not isinstance(parameter_result, ParameterAssignmentResult):
        raise TypeError("parameter_result must be a ParameterAssignmentResult")
    if not isinstance(charge_result, ChargeAssignmentResult):
        raise TypeError("charge_result must be a ChargeAssignmentResult")
    if not isinstance(nonbonded_policy, NonbondedPolicy):
        raise TypeError("nonbonded_policy must be a NonbondedPolicy")
    parameter_result.validate_integrity(system)
    charge_result.validate_integrity(system)
    if not parameter_result.complete_supported_scope:
        raise ParameterCompositionError(
            "Numerical parameter assignment is incomplete for supported scope"
        )
    if not charge_result.complete:
        raise ParameterCompositionError("Charge assignment is incomplete")
    if parameter_result.graph_signature != charge_result.graph_signature:
        raise ParameterCompositionError(
            "Parameter and charge results describe different chemical graphs"
        )
    policy_digest = nonbonded_policy_signature(nonbonded_policy)
    aggregate_signature = _aggregate_signature(
        parameter_result.assignment_signature,
        charge_result.result_signature,
        policy_digest,
        parameter_result.graph_signature,
    )
    return ParameterizedSystem(
        system=deepcopy(system),
        backend_name=(
            f"{parameter_result.library_name}:{parameter_result.library_version}"
        ),
        site_assignments=deepcopy(dict(parameter_result.site_assignments)),
        interaction_assignments={
            "bond": deepcopy(dict(parameter_result.bond_assignments)),
            "angle": deepcopy(dict(parameter_result.angle_assignments)),
            "proper_torsion": deepcopy(
                dict(parameter_result.proper_torsion_assignments)
            ),
        },
        charge_assignments=deepcopy(dict(charge_result.assignments)),
        nonbonded_policy=nonbonded_policy,
        aggregate_signature=aggregate_signature,
        metadata={
            "parameter_assignment": {
                "assignment_signature": parameter_result.assignment_signature,
                "complete_supported_scope": True,
            },
            "charge_assignment": {
                "result_signature": charge_result.result_signature,
                "method": charge_result.method,
                "method_version": charge_result.method_version,
                "complete": True,
                "unit": charge_result.unit,
            },
            "nonbonded_policy": {
                "name": nonbonded_policy.name,
                "version": nonbonded_policy.version,
                "policy_signature": policy_digest,
                "mixing_rule": nonbonded_policy.mixing_rule,
            },
            "aggregate": {
                "aggregate_signature": aggregate_signature,
                "parameter_coverage_complete": True,
                "charge_assignment_complete": True,
                "nonbonded_policy_available": True,
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
        },
    )


def _aggregate_signature(
    parameter_signature: str,
    charge_signature: str,
    policy_signature: str,
    graph_signature: str,
) -> str:
    payload = {
        "schema": "island_parameterized_components_v1",
        "parameter_assignment_signature": parameter_signature,
        "charge_assignment_signature": charge_signature,
        "nonbonded_policy_signature": policy_signature,
        "graph_signature": graph_signature,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
