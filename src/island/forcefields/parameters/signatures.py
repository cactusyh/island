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


def parameter_assignment_signature(
    *,
    graph_digest: str,
    typing_digest: str,
    typing_assignment_digest: str,
    library_digest: str,
    selected_parameter_ids: tuple[tuple[str, tuple[int, ...], str], ...],
    engine_name: str,
    engine_version: str,
) -> str:
    """Hash all authoritative inputs and deterministic selected records."""
    return _digest(
        {
            "graph_signature": graph_digest,
            "typing_signature": typing_digest,
            "typing_assignment_signature": typing_assignment_digest,
            "library_signature": library_digest,
            "selected_parameters": selected_parameter_ids,
            "engine_name": engine_name,
            "engine_version": engine_version,
        }
    )


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
