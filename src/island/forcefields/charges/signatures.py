"""Deterministic signatures for charge tables, inputs, and results."""

import hashlib
import json
from dataclasses import asdict

from island.forcefields.charges.models import AtomTypeChargeTable


def charge_table_signature(table: AtomTypeChargeTable) -> str:
    return _digest(
        {
            "name": table.name,
            "version": table.version,
            "required_ruleset_name": table.required_ruleset_name,
            "required_ruleset_version": table.required_ruleset_version,
            "required_ruleset_signature": table.required_ruleset_signature,
            "source": table.source,
            "supported_representation": table.supported_representation,
            "unit": table.unit,
            "entries": [
                asdict(entry)
                for entry in sorted(table.entries, key=lambda item: item.entry_id)
            ],
        }
    )


def charge_input_signature(payload: object) -> str:
    return _digest(payload)


def charge_result_signature(result: object) -> str:
    return _digest(
        {
            "schema": "island_charge_assignment_result_v1",
            "method": result.method,
            "method_version": result.method_version,
            "source": result.source,
            "representation": result.representation,
            "assignments": [
                {"key": key, **asdict(assignment)}
                for key, assignment in sorted(result.assignments.items())
            ],
            "diagnostics": [asdict(item) for item in result.diagnostics],
            "coverage": asdict(result.coverage),
            "component_diagnostics": [
                asdict(item) for item in result.component_diagnostics
            ],
            "complete": result.complete,
            "unit": result.unit,
            "tolerance": result.tolerance,
            "total_assigned_charge": result.total_assigned_charge,
            "total_formal_charge": result.total_formal_charge,
            "total_charge_residual": result.total_charge_residual,
            "total_within_tolerance": result.total_within_tolerance,
            "target_charge": result.target_charge,
            "graph_signature": result.graph_signature,
            "typing_signature": result.typing_signature,
            "typing_assignment_signature": result.typing_assignment_signature,
            "table_signature": result.table_signature,
            "input_signature": result.input_signature,
        }
    )


def _digest(value: object) -> str:
    serialized = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()
