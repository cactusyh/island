"""Source-row accounting; presence and bounded observations are not full coverage."""

from collections import Counter
from hashlib import sha256

from .catalog import EXTRA, inspect_pcff_full_source
from .class2 import FAMILIES
from .domain_coverage import pcff_domain_coverage
from .source import FLAGS, SEMANTIC, boundary, require


@boundary
def pcff_source_row_ledger(source, *, verified_rows=None):
    source.require_assignment()
    verified_rows = {} if verified_rows is None else verified_rows
    require(
        type(verified_rows) is dict, "Explicit row/case verification mapping required"
    )
    catalog = inspect_pcff_full_source(source)
    numerical = {r["id"]: r for r in catalog["records"]}
    domains = pcff_domain_coverage(source, profile="island_pcff_source_graph_v4")
    labels = {r["type"]: r for r in domains["type_rows"]}
    rows, sections = [], []
    for section in source.inventory["sections"]:
        name, namespace = section["name"], section["namespace"]
        semantic = {r["line"]: r for r in section["records"]}
        count = 0
        for line in section["lines"]:
            fields = line["raw"].split("!", 1)[0].split()
            if (
                not fields
                or not fields[0][0].isdigit()
                or len(fields) < 2
                or not fields[1].isdigit()
            ):
                continue
            if name in ("version", "reference", "end"):
                continue
            rid = f"{name}:{namespace}:{line['line']}"
            parsed = semantic.get(line["line"], numerical.get(rid))
            rule = None
            if name == "atom_types":
                rule = labels[parsed["data"]["type"]]
            implemented = name in FAMILIES or name in EXTRA
            blocking = (
                rule["next_action"]
                if rule and rule["declared_predicate"] is None
                else "Charge/equivalence applicability requires a complete chemical graph and formal-charge contract"
                if name in SEMANTIC
                else "No source-backed executable interpretation"
                if not implemented
                else "Automatic namespace may supplement only explicitly authorized families; missing cross terms are not zero"
                if namespace == "cff91_auto"
                else "Bounded row observations do not verify every applicable graph/interaction"
            )
            if name == "wilson_out_of_plane":
                if parsed and parsed.get("original_values", [0, 0])[1] != 0:
                    blocking = "Nonzero signed Wilson equilibrium and legal peripheral permutation semantics unresolved; strict model disabled"
                elif namespace != "cff91":
                    blocking = "Automatic out-of-plane resolution not authorized by current native interpretation"
            rows.append(
                {
                    "id": rid,
                    "source_sha256": source.identity["sha256"],
                    "line": line["line"],
                    "raw_line_sha256": sha256(line["raw"].encode()).hexdigest(),
                    "section": name,
                    "namespace": namespace,
                    "version": fields[0],
                    "reference": fields[1],
                    "parser": "semantic_record"
                    if line["line"] in semantic
                    else "numerical_record"
                    if rid in numerical
                    else "raw_record_only",
                    "typing": rule["graph_predicate_status"]
                    if rule
                    else "interaction_or_lookup_record_not_atom_rule",
                    "typing_rule": rule["declared_predicate"] if rule else None,
                    "charge_rule": "direct/oriented ordinary/declared automatic increments; component checks; no normalization"
                    if name == "bond_increments"
                    else "not_a_charge_increment",
                    "parameter_family": name,
                    "executable_implementation": "Class2/source kernel; source-dependent applicability"
                    if name in FAMILIES
                    else "explicit quadratic/periodic fallback kernel"
                    if name in EXTRA
                    else "not_an_energy_term"
                    if name in SEMANTIC
                    else "implementation_unresolved",
                    "bounded_verified_cases": sorted(verified_rows.get(rid, [])),
                    "globally_verified": False,
                    "blocking_reason": blocking,
                }
            )
            count += 1
        sections.append(
            {
                "family": name,
                "namespace": namespace,
                "header_line": section["line"],
                "rows": count,
                "empty": count == 0,
                "empty_meaning": "No source coefficient; does not establish universal physical zero"
                if count == 0
                else None,
            }
        )
    known = {r["id"] for r in rows}
    require(
        set(verified_rows).issubset(known), "Verification names nonexistent source rows"
    )
    return {
        "schema": "island_pcff_source_row_coverage_v1",
        "source": source.identity,
        "rows": rows,
        "sections": sections,
        "counts": dict(Counter(r["section"] for r in rows)),
        "type_ledger": domains,
        "full_source_complete": False,
        "coverage_meaning": "All numerical rows accounted; no promotion from parsing/queryability or one bounded fixture",
        **FLAGS,
    }
