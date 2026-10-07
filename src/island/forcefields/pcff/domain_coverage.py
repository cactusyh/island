"""All-source J3 ledger: implementations and bounded observations stay distinct."""

from copy import deepcopy

from .catalog import inspect_pcff_full_source
from .domains import PROFILE_NAME, RULES
from .expanded import RULE_TYPES
from .source import FLAGS, boundary, records, require


@boundary
def pcff_domain_coverage(source, *, profile=PROFILE_NAME):
    """Owned offline source/rule ledger, without inferred global chemical coverage.

    No molecule or supplied receipt can mark a type universally validated. Every
    family retains its own required evidence and source/equivalence row identities.
    """
    require(
        profile
        in (PROFILE_NAME, "island_pcff_source_graph_v4", "island_pcff_source_graph_v5"),
        "Unsupported ledger profile",
    )
    rules = dict(RULES)
    if profile in ("island_pcff_source_graph_v4", "island_pcff_source_graph_v5"):
        from .organic_domains import RULES as ORGANIC_RULES

        rules.update(ORGANIC_RULES)
    if profile == "island_pcff_source_graph_v5":
        from .amine_domains import EVIDENCE

        rules.update(
            {
                "na": EVIDENCE["neutral"],
                "hn": EVIDENCE["neutral"],
                "n4": EVIDENCE["cation"],
                "h+": EVIDENCE["cation"],
            }
        )
    inventory = source.inventory
    catalog = inspect_pcff_full_source(source)
    families = [
        {k: deepcopy(v) for k, v in s.items() if k != "raw_lines"}
        for s in catalog["sections"]
        if s["name"]
        not in (
            "define",
            "reference",
            "version",
            "end",
            "atom_types",
            "equivalence",
            "auto_equivalence",
        )
    ]
    types = []
    for row in records(inventory, "atom_types"):
        label = row["data"]["type"]
        implemented = label in RULE_TYPES or label in rules
        eq = [
            {k: deepcopy(v) for k, v in r.items() if k != "raw"}
            for kind in ("equivalence", "auto_equivalence")
            for r in records(inventory, kind)
            if r["data"]["type"] == label
        ]
        aliases = sorted(
            {v for r in eq for v in r["data"]["families"].values()} | {label}
        )
        family_evidence = []
        for section in families:
            rows = [
                r
                for r in catalog["records"]
                if r["section"] == section["name"]
                and r["namespace"] == section["namespace"]
            ]
            candidates = [
                r["id"] for r in rows if any(t in aliases for t in r["types"])
            ]
            family_evidence.append(
                {
                    "family": section["name"],
                    "namespace": section["namespace"],
                    "status": (
                        "intentionally_out_of_operational_scope"
                        if section["name"] == "torsion-torsion_1"
                        else "interpretation_unresolved"
                        if section["name"] == "wilson_out_of_plane"
                        else "implemented_but_unverified"
                        if implemented
                        else "algorithm_not_implemented"
                    ),
                    "applicability_note": "No torsion-torsion source rows; operational scope exclusion is not a universal zero"
                    if section["name"] == "torsion-torsion_1"
                    else "Only zero Wilson equilibrium currently executable; nonzero signed/permutation semantics unresolved"
                    if section["name"] == "wilson_out_of_plane"
                    else "Requires complete graph-local interaction evidence",
                    "source_candidate_rows": candidates,
                    "candidate_meaning": "label/equivalence occurrence only; not an interaction assignment or coverage proof",
                    "next_action": f"Verify oriented {section['name']} interactions for {label} ({row['data']['comment']}) using the listed equivalences; retain missing rows after exact declared searches.",
                }
            )
        ambiguous = label in {
            "c+",
            "cr",
            "ci",
            "c_a",
            "cg",
            "nb",
            "n1",
            "n2",
            "nr",
            "nh+",
            "nho",
            "ni",
            "nz",
            "az",
            "oah",
            "oas",
            "ob",
            "osh",
            "oss",
            "sz",
            "hb",
            "hoa",
            "hos",
            "ca+",
            "sf",
            "s-",
        }
        unresolved = (
            "chemically_ambiguous" if ambiguous else "algorithm_not_implemented"
        )
        types.append(
            {
                "type": label,
                "source_row": row["id"],
                "source_line": row["line"],
                "source_description": row["data"]["comment"],
                "source_identity": source.identity,
                "graph_predicate_status": "implemented_but_unverified"
                if implemented
                else unresolved,
                "declared_predicate": rules.get(
                    label, "frozen J1/J2 predicate" if implemented else None
                ),
                "explicit_label_validation": "same chemical predicate; no bypass"
                if implemented
                else "unresolved/rejected",
                "native_charges": "interaction-dependent; zero base only; formal/component checks required",
                "atom_nonbonded": "source/equivalence candidate inventory; not globally verified",
                "families": family_evidence,
                "equivalence_rows": eq,
                "executable_evaluator_status": "required interactions must form a complete validated model",
                "independent_numerical_verification": "reference_validation_unavailable_for_full_declared_domain",
                "classification": "implemented_but_unverified"
                if implemented
                else unresolved,
                "next_action": f"{'Extend bounded evidence for' if implemented else 'Establish a deterministic graph/charge-state predicate for'} {label}: {row['data']['comment']}. Validate the ordinary/automatic paths listed here and every applicable coupling before claiming operational coverage.",
                "fixture_observations": [],
            }
        )
    return {
        "schema": "island_pcff_domain_coverage_v1",
        "source": source.identity,
        "profile": profile,
        "counts": catalog["counts"],
        "type_rows": types,
        "sections": families,
        "rule_labels": sum(t["declared_predicate"] is not None for t in types),
        "denominator": 133,
        "full_source_complete": False,
        "global_verified_labels": 0,
        "interpretation": "bounded fixture verification does not certify every environment of a source label",
        **FLAGS,
    }
