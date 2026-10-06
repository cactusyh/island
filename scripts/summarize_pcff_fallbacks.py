"""Assemble separate J2 receipts and all-source gaps without rewriting evidence.

Success means bounded J2 assignment, kernels and workflow gates. The full-source
flag remains false; --require-full-source returns nonzero.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from pcff_fallback_gates import (
    assignment_and_term_gate,
    case_map,
    corrected_case,
    operational_gate,
)

from island.charge_references.records import unpack
from island.workflows.storage import checksum, json_bytes, publish


def load(path):
    return json.loads(path.read_text())


def hashes(root):
    return {
        str(p.relative_to(root)): checksum(p.read_bytes())
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-full-source", action="store_true")
    args = parser.parse_args()
    root = args.root
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    baseline = load(root / "acceptance/report.json")
    checked = load(root / "corrected-checks/report.json")
    workflow = load(root / "workflow-checked/outcome.json")
    historical = load(root / "historical-readonly.json")
    ledger = load(root / "acceptance/coverage.json")
    case_map(baseline["cases"])
    correction = corrected_case(baseline, checked)
    checked_gate = assignment_and_term_gate(
        checked, root / "corrected-checks", args.source
    )
    gaps = []
    family_counts = Counter()
    cases = []
    for result in baseline["cases"]:
        name = result["case"]
        work = root / "acceptance" / name
        typing = unpack((work / "typing.json").read_text())
        case = {
            k: v
            for k, v in result.items()
            if k not in ("fallback_assignments", "diagnostics")
        }
        case["types"] = sorted(set(typing["assignments"].values()))
        if name == "dichlorine":
            if "error" in case:
                case["original_harness_error"] = case.pop("error")
            case.update(
                {
                    k: v
                    for k, v in correction.items()
                    if k not in ("fallback_assignments", "diagnostics")
                }
            )
        cases.append(case)
        if result["status"] == "charge_missing":
            for diag in result["diagnostics"]:
                gaps.append(
                    {
                        "case": name,
                        "stage": "native_charge",
                        "classification": "source_parameter_missing",
                        **diag,
                    }
                )
        if result["status"] == "parameter_missing":
            assignment = unpack((work / "parameters.json").read_text())
            blocked = {r["id"] for r in result["diagnostics"]}
            for a in assignment["assignments"]:
                if a["id"] not in blocked:
                    continue
                family_counts[a["family"]] += 1
                gaps.append(
                    {
                        "case": name,
                        "stage": "model",
                        "classification": "interpretation_unresolved"
                        if a["status"] == "ambiguous"
                        else "source_parameter_missing",
                        "source_search_scope": "direct/ordinary family exact and leading wildcards; automatic positional bases where declared; no auto cross-term family exists",
                        **a,
                    }
                )
    # Keep all 133 labels and avoid promoting a whole type from one passing fixture.
    for row in ledger["type_rows"]:
        row["fixtures"] = [
            {
                "case": c["case"],
                "typing": c.get("typed", False),
                "charges": c.get("charges_complete", False),
                "model": c.get("model_complete", False),
                "independently_verified_in_j2": c["case"] == "dichlorine"
                and checked_gate["passed"],
            }
            for c in cases
            if row["type"] in c["types"]
        ]
        if not row["rule_implemented"]:
            label = row["type"]
            element = row["element"]
            if label in ("h", "h*", "o", "n", "n+", "ho2", "s", "o=", "oe"):
                action = "Establish specificity/aliases against specialized labels with independently typed fixtures; generic label is not a catch-all. Preserve oh/ho and hn2 distinctions."
            elif element in (
                "Ag",
                "Al",
                "Au",
                "Cr",
                "Cu",
                "Fe",
                "K",
                "Li",
                "Mo",
                "Na",
                "Ni",
                "Pb",
                "Pd",
                "Pt",
                "Sn",
                "W",
                "Ca",
            ):
                action = "Establish oxidation/base-charge and coordination interpretation for this exact source type; separately check zeolite versus elemental/ionic records before graph rules."
            elif label == "dw":
                action = "Add explicit isotope/mass semantics and heavy-water typed reference; current standard-mass guard must remain until versioned."
            else:
                action = "Author bond-order/formal-charge/ring/neighbor predicates for the exact source description and resolve competing labels using independently typed reference data; verify native increments and every coupling."
            row["next_action"] = {
                "source_type": label,
                "source_row": row["source_row"],
                "chemical_description": row["description"],
                "action": action,
            }
    term_checks = checked.get("term_checks", [])
    for row in ledger["families"]:
        matching = [
            r for r in term_checks if r["family"] == row["family"] and "source_row" in r
        ]
        if row["namespace"] == "cff91_auto":
            row["issues"] = [
                "J2 supplementation is explicit and versioned; no absent cross-term inference. See exact unresolved interactions."
            ]
        if matching and checked_gate["passed"]:
            row["stages"]["parameter_resolution"] = "implemented_and_verified"
            row["stages"]["energy_force"] = "implemented_and_verified"
            row["stages"]["independent_verification"] = "implemented_and_verified"
            row["verification_scope"] = {
                "meaning": "representative selection/equation verification, not every source row or all chemical applicability",
                "source_lines": [r["source_row"]["line"] for r in matching],
            }
            if row["family"] == "wilson_out_of_plane":
                row["issues"].append(
                    "Nonzero equilibrium remains unsupported; verified subset is chi0=0."
                )
    ledger["j2_fixture_counts"] = {
        "denominator": len(cases),
        "typing": sum(c.get("typed", False) for c in cases),
        "native_charges": sum(c.get("charges_complete", False) for c in cases),
        "complete_model": sum(c.get("model_complete", False) for c in cases),
        "new_independent_whole_system": 1 if checked_gate["passed"] else 0,
    }
    ledger["full_source_complete"] = False
    ledger["remaining_rule_labels"] = [
        r["type"] for r in ledger["type_rows"] if not r["rule_implemented"]
    ]
    evidence = {
        "schema": "island_j2_evidence_v1",
        "base": "758449bb125ed706d2e982b34b4a2c2d8f27d56c",
        "source": baseline["inventory"],
        "declaration": load(root / "declaration.json"),
        "cases": cases,
        "term_checks": term_checks,
        "workflow": workflow,
        "historical": historical,
        "offline_completed": load(root / "offline-completed.json"),
        "environment": load(root / "environment.json"),
        "independent_assignment": load(root / "independent-assignment.json"),
        "msi2lmp_control": load(root / "msi2lmp-control/outcome.json"),
        "original_failed_receipts": {
            "audit": baseline.get("term_failure"),
            "workflow": load(root / "workflow/outcome.json"),
        },
        "assignment_and_term_validation": checked_gate,
        "j2_operational_gate": operational_gate(checked_gate, workflow, historical),
        "full_source_complete": False,
        "artifacts": hashes(root),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(out / "phase_4j2.json", json_bytes(evidence))
    publish(out / "phase_4j2_coverage.json", json_bytes(ledger))
    publish(
        out / "phase_4j2_gaps.json",
        json_bytes(
            {
                "remaining_interactions": gaps,
                "counts_by_family": dict(family_counts),
                "remaining_types": ledger["remaining_rule_labels"],
                "full_source_complete": False,
            }
        ),
    )
    print(
        "J2 operational gate:",
        evidence["j2_operational_gate"],
        "full-source gate: False",
    )
    return 0 if evidence["j2_operational_gate"] and not args.require_full_source else 1


if __name__ == "__main__":
    raise SystemExit(main())
