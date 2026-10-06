"""Recomputed source-versus-model blocker assessments, never usable potentials."""

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .automatic import PCFFAutomaticTypingResult, graph_system
from .charges import identity
from .model import PROFILE, bb13_policy_applies
from .operational import inspect_pcff_operational_support
from .source import FLAGS, PCFFSource, boundary, require

SCHEMA = "island_pcff_resolution_assessment_v1"


def classify(query, model_diagnostics):
    """Separate source lookup results from the established executable policy."""
    result = query["resolution"]
    status = result["status"]
    policy = None
    if status == "not_applicable":
        category = "structurally_not_applicable"
        blocking = False
    elif bb13_policy_applies(
        query["family"],
        query["supplied_types"],
        query["equilibrium_dependencies"],
        resolution_policy=result.get("selection_policy"),
    ):
        category = "existing_policy_derived_zero"
        policy = PROFILE["bb13_policy"]
        blocking = False
    elif status == "ambiguous":
        category = "interpretation_unresolved"
        blocking = True
    elif status == "missing":
        category = "source_parameter_missing_under_declared_searches"
        blocking = True
    elif not query["dependencies_complete"]:
        category = "equilibrium_dependency_unresolved"
        blocking = True
    elif (
        query["family"] == "wilson_out_of_plane" and result["normalized_values"][1] != 0
    ):
        category = "nonzero_wilson_semantics_unresolved"
        blocking = True
    elif query["request_id"] in model_diagnostics:
        category = "native_model_rejection"
        blocking = True
    else:
        require(status == "assigned", "Unknown source-resolution status")
        category = "source_assigned"
        blocking = False
    return {
        "request_id": query["request_id"],
        "family": query["family"],
        "sites": query["sites"],
        "supplied_types": query["supplied_types"],
        "raw_source_status": status,
        "classification": category,
        "structural_model_blocker": blocking,
        "policy": policy,
        "dependencies": query["equilibrium_dependencies"],
        "selected_source_rows": [r["record_id"] for r in result.get("selected", [])],
        "candidate_source_rows": sorted(
            {r["record_id"] for r in result.get("candidates", [])}
        ),
        "native_diagnostic": model_diagnostics.get(query["request_id"]),
        "note": "Policy coverage is conditional on all native charge and model gates; never a source parameter row"
        if policy
        else result.get("reason"),
    }


def derive(typing, resolution_policy, special_pairs):
    p = typing.payload
    inspection = inspect_pcff_operational_support(
        graph_system(p["graph"]),
        typing,
        resolution_policy=resolution_policy,
        special_pairs=special_pairs,
    )
    diagnostics = {d["id"]: d for d in inspection["model_diagnostics"] if "id" in d}
    entries = [classify(q, diagnostics) for q in inspection["interaction_queries"]]
    charge = inspection["native_charge_record"]["native_charge_record"]
    return {
        "schema": SCHEMA,
        "diagnostic_only": True,
        "typing": p,
        "source": typing.source.identity,
        "typing_identity": typing.identity,
        "resolution_policy": resolution_policy,
        "special_pairs": deepcopy(special_pairs),
        "inspection_identity": identity(inspection),
        "native_identities": {
            k: inspection[k]
            for k in ("charge_identity", "assignment_identity", "model_identity")
        },
        "charges_complete": inspection["charges_complete"],
        "charge_diagnostics": charge["diagnostics"],
        "charge_components": charge["components"],
        "native_model_complete": inspection["model_complete"],
        "entries": entries,
        "counts": dict(Counter(e["classification"] for e in entries)),
        "structural_blocker_count": sum(e["structural_model_blocker"] for e in entries),
        "scope": "Current named native policy only; no new coefficient, equivalence, precedence, applicability or charge convention",
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFResolutionAssessment:
    """Owned, data-only diagnostic; cannot be adopted or evaluated as a model."""

    json_text: str
    source: PCFFSource

    @boundary
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        require(p["schema"] == SCHEMA, "Unsupported resolution assessment schema")
        typing = PCFFAutomaticTypingResult(pack(p["typing"]), self.source)
        typing.validate_integrity(system)
        expected = derive(typing, p["resolution_policy"], p["special_pairs"])
        require(pack(p) == pack(expected), "Contradictory PCFF resolution assessment")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def assess_pcff_resolution(system, typing, *, resolution_policy, special_pairs):
    require(
        type(typing) is PCFFAutomaticTypingResult, "Validated expanded typing required"
    )
    typing.validate_integrity(system)
    result = PCFFResolutionAssessment(
        pack(derive(typing, resolution_policy, special_pairs)), typing.source
    )
    result.validate_integrity(system)
    return result


@boundary
def save_pcff_resolution_assessment(result, path):
    require(
        type(result) is PCFFResolutionAssessment, "Expected PCFF resolution assessment"
    )
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_resolution_assessment(path, source, *, system=None):
    result = PCFFResolutionAssessment(Path(path).read_text(), source)
    result.validate_integrity(system)
    return result
