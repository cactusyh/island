"""Executable full-source ledger. Unresolved inventory is never counted complete."""

from collections import Counter
from copy import deepcopy

from .catalog import inspect_pcff_full_source
from .expanded import RULE_TYPES
from .source import FLAGS, boundary, records

STAGES = (
    "ingestion",
    "chemical_recognition",
    "charge_assignment",
    "parameter_resolution",
    "energy_force",
    "independent_verification",
)
STATUSES = (
    "implemented_and_verified",
    "implementation_missing",
    "source_parameter_missing",
    "chemical_rule_ambiguous",
    "interpretation_unresolved",
    "reference_validation_unavailable",
)


@boundary
def pcff_coverage_ledger(source):
    """Inventory-only baseline. Actual fixture receipts are a separate evidence file.

    This deliberately does not label code availability as numerical verification.
    It is deterministic, contains every source type/namespace, and can gate CI.
    """
    catalog = inspect_pcff_full_source(source)
    inv = source.inventory
    types = []
    for row in records(inv, "atom_types"):
        label = row["data"]["type"]
        entry = {
            "type": label,
            "source_row": row["id"],
            "element": row["data"]["element"],
            "description": row["data"]["comment"],
            "stages": {stage: "reference_validation_unavailable" for stage in STAGES},
            "rule_implemented": label in RULE_TYPES,
            "fixtures": [],
            "issues": [],
        }
        entry["stages"]["ingestion"] = "implemented_and_verified"
        if label not in RULE_TYPES:
            ambiguous = {
                "c+",
                "cr",
                "ci",
                "c_a",
                "cg",
                "nb",
                "nn",
                "n1",
                "n2",
                "nr",
                "nh+",
                "nho",
                "ni",
                "npc",
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
                "dw",
            }
            entry["stages"]["chemical_recognition"] = (
                "chemical_rule_ambiguous"
                if label in ambiguous
                else "implementation_missing"
            )
            entry["issues"].append(
                "No audited graph predicate/independent typed fixture implemented; see source row. Ambiguous specialized environments require chemical evidence; other missing rules remain implementation work."
            )
        else:
            entry["issues"].append(
                "Graph rule implemented; case-specific evidence and remaining interaction gaps must be read from acceptance receipts."
            )
        types.append(entry)
    families = []
    for section in catalog["sections"]:
        name, namespace = section["name"], section["namespace"]
        if name in ("define", "reference", "version", "end"):
            continue
        families.append(
            {
                "family": name,
                "namespace": namespace,
                "source_line": section["line"],
                "numeric_record_count": section["numerical_records"],
                "semantic_record_count": section["semantic_records"],
                "record_ids": section["record_ids"],
                "duplicate_ordered_keys": section["repeated_ordered_keys"],
                "stages": {
                    stage: (
                        "implemented_and_verified"
                        if stage == "ingestion"
                        else "reference_validation_unavailable"
                    )
                    for stage in STAGES
                },
                "issues": (
                    [
                        "Empty source section; no source parameter and no universal physical-zero inference."
                    ]
                    if name == "torsion-torsion_1"
                    else [
                        "Automatic namespace interpreted and explicitly queryable; no implicit quadratic/torsion/charge fallback in models."
                    ]
                    if namespace == "cff91_auto"
                    else []
                ),
            }
        )
    for row in families:
        if row["family"] in ("torsion_1", "quadratic_bond", "quadratic_angle"):
            row["stages"]["energy_force"] = "implementation_missing"
        if row["family"] == "torsion-torsion_1":
            row["stages"]["parameter_resolution"] = "source_parameter_missing"
            row["stages"]["energy_force"] = "interpretation_unresolved"
    return {
        "schema": "island_pcff_full_source_coverage_v1",
        "source": source.identity,
        "target_counts": catalog["counts"],
        "distinct_atom_types": len({r["type"] for r in types}),
        "type_rows": types,
        "families": families,
        "status_vocabulary": list(STATUSES),
        "completion": {
            stage: dict(Counter(r["stages"][stage] for r in types)) for stage in STAGES
        },
        "full_source_complete": False,
        "meaning": "Baseline ledger tracks implementation and evidence gaps. Separate acceptance receipts contain executed fixture outcomes; parsing does not imply model support.",
        **FLAGS,
    }


def attach_fixture_evidence(ledger, outcomes):
    """Attach traceable experimental observations; do not promote entire types.

    A passing molecule verifies only its exercised environments, not every use of
    its atom types. The conservative baseline completion measures stay unchanged.
    """
    result = deepcopy(ledger)
    for row in result["type_rows"]:
        row["fixtures"] = [
            {
                "case": c["name"],
                "status": c["status"],
                "evidence_directory": c["name"],
                "stages": c.get("stages", {}),
            }
            for c in outcomes
            if row["type"] in c.get("types", [])
        ]
    result["executed_case_counts"] = dict(Counter(c["status"] for c in outcomes))
    result["bounded_fixture_coverage"] = {
        stage: sorted(
            {
                label
                for case in outcomes
                if case.get("stages", {}).get(key)
                for label in case.get("types", [])
            }
        )
        for stage, key in (
            ("automatic_typing", "typing"),
            ("native_charges", "charges"),
            ("model_definition", "model"),
            ("independent_numerical", "independent_numerical"),
        )
    }
    result["bounded_fixture_coverage_meaning"] = (
        "Types exercised successfully in at least one declared fixture; not universal coverage of every environment or interaction involving that type."
    )
    return result
