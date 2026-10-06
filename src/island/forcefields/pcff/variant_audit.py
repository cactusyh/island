"""Recomputed cross-variant diagnostics for the original J6 physical requests."""

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .charges import identity
from .model import PROFILE, bb13_policy_applies
from .resolution import PCFFResolutionAssessment
from .source import FLAGS, PIN, boundary, records, require, select
from .variants import query_compiled, query_context, selections

SCHEMA = "island_pcff_source_variant_audit_v1"


def component_sites(graph):
    adjacency = {s["id"]: set() for s in graph["sites"]}
    for b in graph["bonds"]:
        a, c = b["sites"]
        adjacency[a].add(c)
        adjacency[c].add(a)
    remaining = set(adjacency)
    components = []
    while remaining:
        todo = [min(remaining)]
        group = set()
        while todo:
            site = todo.pop()
            if site in group:
                continue
            group.add(site)
            todo.extend(adjacency[site] - group)
        remaining -= group
        components.append(sorted(group))
    return components


def charge_diagnostic(context, graph, labels, resolution_policy):
    charges = {s["id"]: Decimal(0) for s in graph["sites"]}
    contributions = []
    missing = []
    cache = {}
    atom_rows = records(context["inventory"], "atom_types")
    atom_evidence = []
    neighbors = Counter(i for b in graph["bonds"] for i in b["sites"])
    compatible = True
    for atom in graph["sites"]:
        row = select([r for r in atom_rows if r["data"]["type"] == labels[atom["id"]]])
        valid = (
            row is not None
            and row["record"]["data"]["element"] == atom["element"]
            and row["record"]["data"]["connections"] == neighbors[atom["id"]]
        )
        compatible = compatible and valid
        atom_evidence.append(
            {
                "site": atom["id"],
                "type": labels[atom["id"]],
                "source_atom_record": row,
                "element_connection_compatible": valid,
            }
        )
    for bond in graph["bonds"]:
        sites = bond["sites"]
        types = [labels[i] for i in sites]
        key = tuple(types)
        if key not in cache:
            cache[key] = query_compiled(
                context, "bond_increments", "cff91_auto", types, resolution_policy
            )
        q = cache[key]
        if q["classification"] != "source_row_present":
            missing.append({"sites": sites, "types": types, "query": q})
            continue
        vals = [Decimal(str(x)) for x in q["forward"]["normalized_values"]]
        for sid, value in zip(sites, vals, strict=True):
            charges[sid] += value
        contributions.append(
            {
                "sites": sites,
                "types": types,
                "endpoint_increments_e": [str(v) for v in vals],
                "source_selection": context["variant"]["selection"],
                "query": q,
            }
        )
    byid = {s["id"]: s for s in graph["sites"]}
    comps = []
    for group in component_sites(graph):
        formal = sum(byid[i]["formal_charge"] for i in group)
        total = sum((charges[i] for i in group), Decimal(0))
        coverage = not any(set(m["sites"]) & set(group) for m in missing)
        comps.append(
            {
                "sites": group,
                "formal_charge_e": formal,
                "known_contribution_total_e": str(total),
                "residual_e": str(total - Decimal(formal)) if coverage else None,
                "coverage_complete": coverage,
                "charge_check_passed": coverage
                and abs(total - Decimal(formal)) <= Decimal(str(PIN["tolerance_e"])),
            }
        )
    complete = (
        compatible and not missing and all(c["charge_check_passed"] for c in comps)
    )
    return {
        "base_charge_e": "0",
        "source_selection": context["variant"]["selection"],
        "original_pinned_type_labels": labels,
        "atom_record_evidence": atom_evidence,
        "atom_element_connection_compatible": compatible,
        "chemical_typing_equivalence_established": context["variant"]["selection"][
            "sha256"
        ]
        == PIN["sha256"],
        "contributions": contributions,
        "missing_bonds": missing,
        "components": comps,
        "charges_e": {i: str(v) for i, v in charges.items()},
        "source_charge_calculation_complete": complete,
        "native_charge_validation_authorized": context["variant"]["selection"]["sha256"]
        == PIN["sha256"],
        "classification": "source_row_present"
        if complete
        else "native_charge_incomplete",
        "scope": "Explicit comparison using original PCFF labels; no alternate-source automatic typing or charge assignment is certified",
    }


def derive(assessment, variants):
    assessment.validate_integrity()
    original = unpack(assessment.json_text)
    selections(variants)
    contexts = []
    seen = set()
    for variant in variants:
        if variant.expected_sha256 in seen:
            continue
        seen.add(variant.expected_sha256)
        contexts.append(query_context(variant))
    require(
        PIN["sha256"] in seen, "Original pinned source required for cross-source audit"
    )
    contexts.sort(key=lambda c: c["variant"]["selection"]["sha256"])
    graph = original["typing"]["graph"]
    labels = original["typing"]["assignments"]
    case_results = []
    for context in contexts:
        result_by_id = {}
        cache = {}
        for e in original["entries"]:
            key = (e["family"], tuple(e["supplied_types"]))
            if key not in cache:
                cache[key] = query_compiled(
                    context,
                    e["family"],
                    "cff91",
                    e["supplied_types"],
                    original["resolution_policy"],
                )
            result_by_id[e["request_id"]] = cache[key]
        entries = []
        for e in original["entries"]:
            q = result_by_id[e["request_id"]]
            dependencies = []
            for d in e["dependencies"]:
                dq = result_by_id[d["request_id"]]
                dependencies.append(
                    {
                        "request_id": d["request_id"],
                        "sites": d["sites"],
                        "status": "assigned"
                        if dq["classification"] == "source_row_present"
                        else dq["forward"]["status"],
                        "source_selection": context["variant"]["selection"],
                        "equilibrium_value": dq["forward"].get(
                            "normalized_values", [None]
                        )[0],
                        "query": dq,
                    }
                )
            status = q["classification"]
            policy = None
            pinned = context["variant"]["selection"]["sha256"] == PIN["sha256"]
            if e["raw_source_status"] == "not_applicable":
                status = "source_row_absent"
                reason = "Original structurally inapplicable interaction; not a missing required parameter"
            elif pinned and bb13_policy_applies(
                e["family"], e["supplied_types"], dependencies
            ):
                status = "policy_derived_zero"
                policy = PROFILE["bb13_policy"]
                reason = "Existing pinned operational profile only; never a source row"
            elif status == "source_row_present" and any(
                d["status"] != "assigned" for d in dependencies
            ):
                status = "equilibrium_dependency_missing"
                reason = "Required source-bound equilibrium dependency unresolved"
            else:
                reason = q.get("interpretation_diagnostic")
            entries.append(
                {
                    "request_id": e["request_id"],
                    "family": e["family"],
                    "sites": e["sites"],
                    "types": e["supplied_types"],
                    "classification": status,
                    "query": q,
                    "dependencies": dependencies,
                    "policy": policy,
                    "reason": reason,
                    "structurally_not_applicable": e["raw_source_status"]
                    == "not_applicable",
                }
            )
        charge = charge_diagnostic(
            context, graph, labels, original["resolution_policy"]
        )
        missing = [
            e
            for e in entries
            if e["classification"] not in ("source_row_present", "policy_derived_zero")
            and not e["structurally_not_applicable"]
        ]
        model_complete = (
            pinned
            and original["native_model_complete"]
            and not missing
            and charge["source_charge_calculation_complete"]
        )
        case_results.append(
            {
                "source": context["variant"]["selection"],
                "variant_identity": context["identity"],
                "entries": entries,
                "classification_counts": dict(
                    Counter(e["classification"] for e in entries)
                ),
                "required_interaction_blockers": len(missing),
                "charge_audit": charge,
                "native_model_complete": model_complete,
                "model_classification": "source_row_present"
                if model_complete
                else "model_incomplete",
                "candidate_model_semantics_authorized": pinned,
                "policy_derived_zeros": "Only applied to original pinned profile; never inherited by an alternate hash",
            }
        )
    primary = next(c for c in case_results if c["source"]["sha256"] == PIN["sha256"])
    primary_entries = {e["request_id"]: e for e in primary["entries"]}
    additions = []
    for result in case_results:
        if result is primary:
            continue
        for e in result["entries"]:
            pe = primary_entries[e["request_id"]]
            if (
                pe["query"]["classification"] == "source_row_absent"
                and e["query"]["classification"] == "source_row_present"
            ):
                additions.append(
                    {
                        "classification": "source_variant_only",
                        "request_id": e["request_id"],
                        "sites": e["sites"],
                        "types": e["types"],
                        "primary_source": primary["source"],
                        "alternate_source": result["source"],
                        "alternate_query": e["query"],
                        "dependencies_complete": all(
                            d["status"] == "assigned" for d in e["dependencies"]
                        ),
                        "model_assembly_authorized": False,
                    }
                )
    return {
        "schema": SCHEMA,
        "original_assessment": original,
        "original_assessment_identity": identity(original),
        "variant_identities": sorted(v.identity for v in variants),
        "variant_results": case_results,
        "source_variant_only": additions,
        "source_substitution_performed": False,
        "operational_gate": False,
        "full_source_complete": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFSourceVariantAudit:
    json_text: str
    assessment: PCFFResolutionAssessment
    variants: tuple

    @boundary
    def validate_integrity(self):
        require(
            type(self.assessment) is PCFFResolutionAssessment,
            "Validated original resolution assessment required",
        )
        p = unpack(self.json_text)
        require(p["schema"] == SCHEMA, "Unsupported source-variant audit schema")
        require(
            pack(p) == pack(derive(self.assessment, self.variants)),
            "Contradictory cross-source audit",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def audit_pcff_source_variants(assessment, variants):
    require(
        type(assessment) is PCFFResolutionAssessment,
        "Validated resolution assessment required",
    )
    variants = tuple(variants)
    result = PCFFSourceVariantAudit(
        pack(derive(assessment, variants)), assessment, variants
    )
    result.validate_integrity()
    return result


@boundary
def save_pcff_source_variant_audit(result, path):
    require(
        type(result) is PCFFSourceVariantAudit, "Expected cross-source diagnostic audit"
    )
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_source_variant_audit(path, *, assessment, variants):
    result = PCFFSourceVariantAudit(Path(path).read_text(), assessment, tuple(variants))
    result.validate_integrity()
    return result
