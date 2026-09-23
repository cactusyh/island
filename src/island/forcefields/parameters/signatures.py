"""Deterministic signatures for parameter libraries and assignments."""

import hashlib
import json
from dataclasses import asdict

from island.forcefields.parameters.models import ParameterLibrary


def parameter_library_signature(library: ParameterLibrary) -> str:
    """Hash normalized library identity, requirements, records, values, and units."""
    records = []
    for record in sorted(library.records, key=lambda item: item.parameter_id):
        content = asdict(record)
        content["family"] = record.family
        content["functional_form"] = record.functional_form
        records.append(content)
    return _digest(
        {
            "name": library.name,
            "version": library.version,
            "required_ruleset_name": library.required_ruleset_name,
            "required_ruleset_version": library.required_ruleset_version,
            "required_ruleset_signature": library.required_ruleset_signature,
            "supported_representation": library.supported_representation,
            "description": library.description,
            "source": library.source,
            "records": records,
        }
    )


def parameter_result_content_signature(result: object) -> str:
    """Hash complete assignment-result content except descriptive metadata.

    Version 2 covers selected numerical records, wrapper provenance, diagnostics,
    coverage, completeness semantics, and all authoritative input fingerprints.
    """

    def selections(mapping: object) -> list[dict[str, object]]:
        return [
            {
                "key": key,
                "family": selection.family,
                "site_ids": selection.site_ids,
                "atom_types": selection.atom_types,
                "parameter_id": selection.parameter_id,
                "source": selection.source,
                "parameter": parameter_record_content(selection.parameter),
            }
            for key, selection in sorted(mapping.items())
        ]

    diagnostics = sorted(
        result.diagnostics,
        key=lambda item: (item.family, item.site_ids, item.atom_types, item.reason),
    )
    coverage = result.coverage
    return _digest(
        {
            "signature_schema": "island_parameter_assignment_result_v2",
            "engine_name": result.engine_name,
            "engine_version": result.engine_version,
            "library_name": result.library_name,
            "library_version": result.library_version,
            "representation": result.representation,
            "site_assignments": selections(result.site_assignments),
            "bond_assignments": selections(result.bond_assignments),
            "angle_assignments": selections(result.angle_assignments),
            "proper_torsion_assignments": selections(result.proper_torsion_assignments),
            "diagnostics": [asdict(item) for item in diagnostics],
            "coverage": {
                family: asdict(coverage[family]) for family in sorted(coverage)
            },
            "complete_supported_scope": result.complete_supported_scope,
            "charges_status": result.charges_status,
            "production_validated": result.production_validated,
            "simulation_readiness": result.simulation_readiness,
            "graph_signature": result.graph_signature,
            "typing_signature": result.typing_signature,
            "typing_assignment_signature": result.typing_assignment_signature,
            "ruleset_signature": result.ruleset_signature,
            "library_signature": result.library_signature,
        }
    )


def parameter_record_content(record: object) -> dict[str, object]:
    """Return canonical dataclass content including fixed family/form fields."""
    content = asdict(record)
    content["family"] = record.family
    content["functional_form"] = record.functional_form
    return content


def typing_assignment_content_signature(result: object) -> str:
    """Hash mutable Phase 4A assignment and diagnostic content."""
    assignments = getattr(result, "assignments", {})
    diagnostics = getattr(result, "diagnostics", {})
    return _digest(
        {
            "complete": getattr(result, "complete", None),
            "untyped_site_ids": getattr(result, "untyped_site_ids", None),
            "ambiguous_site_ids": getattr(result, "ambiguous_site_ids", None),
            "assignments": [
                {
                    "key": site_id,
                    "site_id": assignment.site_id,
                    "atom_type": assignment.atom_type,
                    "selected_rule_ids": assignment.selected_rule_ids,
                    "matched_rule_ids": assignment.matched_rule_ids,
                }
                for site_id, assignment in sorted(assignments.items())
            ],
            "diagnostics": [
                {
                    "key": site_id,
                    "site_id": diagnostic.site_id,
                    "status": diagnostic.status,
                    "matched_rule_ids": diagnostic.matched_rule_ids,
                    "eliminated_rule_ids": diagnostic.eliminated_rule_ids,
                    "surviving_rule_ids": diagnostic.surviving_rule_ids,
                    "surviving_atom_types": diagnostic.surviving_atom_types,
                    "eliminated_by": diagnostic.eliminated_by,
                }
                for site_id, diagnostic in sorted(diagnostics.items())
            ],
        }
    )


def _digest(value: object) -> str:
    serialized = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()
